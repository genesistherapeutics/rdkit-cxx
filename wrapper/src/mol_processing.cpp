// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_processing.h"

#include <map>
#include <mutex>
#include <streambuf>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <GraphMol/FileParsers/GeneralFileReader.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDConfig.h>
#ifdef RDK_USE_BOOST_IOSTREAMS
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#endif

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/mol_processing.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/streams.h"
#include "rdkit-cxx/threads.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace FP = RDKit::v2::FileParsers;
using Options = RDKit::GeneralMolSupplier::SupplierOptions;
using Fingerprints = std::vector<std::unique_ptr<ExplicitBitVect>>;

Options to_options(const SupplierOptionsData &data) {
  // RDKit's SMILES suppliers index each line's fields by the name column
  // whenever it is not -1 and is less than their number
  // (SmilesMolSupplier::processLine and
  // MultithreadedSmilesMolSupplier::processMoleculeRecord), reading before
  // the fields for a column below -1, as smiles_params in suppliers.cpp
  // refuses. A smiles column below 0 is refused by RDKit itself.
  if (data.name_column < -1) {
    throw ValueErrorException("nameColumn must be -1 or more");
  }
  Options options;
  options.numWriterThreads = data.num_threads;
  options.sanitize = data.sanitize;
  options.removeHs = data.remove_hs;
  options.strictParsing = data.strict_parsing;
  options.delimiter = std::string(data.delimiter);
  options.smilesColumn = data.smiles_column;
  options.nameColumn = data.name_column;
  options.titleLine = data.title_line;
  options.nameRecord = std::string(data.name_record);
  options.confId2D = data.conf_id_2d;
  options.confId3D = data.conf_id_3d;
  return options;
}

#ifdef RDK_USE_BOOST_IOSTREAMS
// The decompressed bytes of a gzip file, read as RDKit's gzstream
// (RDStreams/streams.h) reads them, through Boost's gzip_decompressor. Boost
// throws for bytes that are not gzip, or end too soon; RDKit's istreams catch
// that and read it as the end of the input (badbit), and so does this
// streambuf.
class GzipFileBuf : public std::streambuf {
 public:
  explicit GzipFileBuf(const std::string &path) {
    d_chain.push(boost::iostreams::gzip_decompressor());
    d_chain.push(boost::iostreams::file_source(
        path, std::ios_base::in | std::ios_base::binary));
  }

 protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    if (d_done) {
      return traits_type::eof();
    }
    std::streamsize n = 0;
    try {
      n = d_chain.sgetn(d_buffer.data(),
                        static_cast<std::streamsize>(d_buffer.size()));
    } catch (...) {
      n = 0;
    }
    if (n <= 0) {
      d_done = true;
      return traits_type::eof();
    }
    setg(d_buffer.data(), d_buffer.data(), d_buffer.data() + n);
    return traits_type::to_int_type(*gptr());
  }

 private:
  boost::iostreams::filtering_streambuf<boost::iostreams::input> d_chain;
  std::vector<char> d_buffer = std::vector<char>(std::size_t{1} << 16);
  bool d_done = false;
};

// What RDKit reads a gzip file through: the checked stream's bytes, from a
// streambuf that, as gzstream's (a filtering_streambuf over a gzip filter),
// throws for every seek and tell, which RDKit's istreams turn into badbit.
// RDKit's suppliers that seek (SmilesMolSupplier, TDTMolSupplier) then read
// nothing from it, as from gzstream, where over a checked stream, which
// reports its position, a failed seek alone would leave SmilesMolSupplier
// looping forever (skipComments ends only at eof or badbit).
class UnseekableBuf : public std::streambuf {
 public:
  explicit UnseekableBuf(std::streambuf &source) : d_source(source) {}

 protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    auto n = d_source.sgetn(d_buffer.data(),
                            static_cast<std::streamsize>(d_buffer.size()));
    if (n <= 0) {
      return traits_type::eof();
    }
    setg(d_buffer.data(), d_buffer.data(), d_buffer.data() + n);
    return traits_type::to_int_type(*gptr());
  }
  pos_type seekoff(off_type, std::ios_base::seekdir,
                   std::ios_base::openmode) override {
    throw std::ios_base::failure("no random access");
  }
  pos_type seekpos(pos_type, std::ios_base::openmode) override {
    throw std::ios_base::failure("no random access");
  }

 private:
  std::streambuf &d_source;
  std::vector<char> d_buffer = std::vector<char>(std::size_t{1} << 16);
};
#endif

// The stream RDKit reads, and the checked stream under it (U12).
struct Input {
  std::unique_ptr<CheckedIStream> checked;
  std::unique_ptr<std::streambuf> unseekable;
  std::unique_ptr<std::istream> outer;
  std::istream *stream() const { return outer ? outer.get() : checked.get(); }
  // Throws ValueErrorException once RDKit has asked for refused input.
  void check() const { checked->check(); }
};

// Opens the file as GeneralMolSupplier::getSupplier does, with its errors.
Input open_input(const std::string &path, const std::string &compression,
                 InputScanner scanner) {
  Input in;
  if (compression.empty()) {
    in.checked = open_checked_file(path, true, std::move(scanner));
    return in;
  }
#ifdef RDK_USE_BOOST_IOSTREAMS
  in.checked = std::make_unique<CheckedIStream>(
      std::make_unique<GzipFileBuf>(path), std::move(scanner));
  in.unseekable = std::make_unique<UnseekableBuf>(*in.checked->rdbuf());
  in.outer = std::make_unique<std::istream>(in.unseekable.get());
  // gzstream opens the file without reporting a failure, so a file it cannot
  // open reads as empty, which getSupplier reports as invalid.
  if (in.outer->peek() == std::char_traits<char>::eof()) {
    throw RDKit::BadFileException("Invalid input file " + path);
  }
  return in;
#else
  throw RDKit::BadFileException(
      "compressed files are only supported if the RDKit is built with "
      "boost::iostreams support");
#endif
}

// Calls into RDKit, then checks the input, as checked in suppliers.cpp:
// refused input explains whatever RDKit did with it, including an exception,
// and once the input has failed RDKit is not called again.
template <class F>
auto checked(const Input &in, F &&f) -> decltype(f()) {
  in.check();
  try {
    if constexpr (std::is_void_v<decltype(f())>) {
      f();
      in.check();
    } else {
      auto result = f();
      in.check();
      return result;
    }
  } catch (...) {
    in.check();
    throw;
  }
}

// Fingerprints one molecule from the file, after the checks the generator's
// own calls make (check_fingerprint_mol), and with its matrix caches cleared
// around the call (U6).
template <class Gen>
struct Fingerprinter {
  const Gen &gen;
  // The generator to check molecules for; null for the default one, whose
  // fingerprints trust nothing check_fingerprint_mol checks.
  const FpGenerator *checked_gen;

  std::unique_ptr<ExplicitBitVect> operator()(RDKit::ROMol &mol) const {
    if (checked_gen) {
      check_fingerprint_mol(*checked_gen, mol);
    }
    MatrixCacheGuard guard(mol);
    RDKit::FingerprintFuncArguments args;
    return gen.getFingerprint(mol, args);
  }
};

// As worker in MolProcessing.cpp: every record read in order, with null for a
// record RDKit cannot read unless it is the last. An exception from the
// generator reaches the caller.
template <class S, class F>
Fingerprints read_in_order(const Input &in, S &suppl, const F &fp) {
  Fingerprints results;
  while (!suppl.atEnd()) {
    auto mol = checked(in, [&] { return suppl.next(); });
    if (mol) {
      results.push_back(fp(*mol));
    } else if (!suppl.atEnd()) {
      results.push_back(nullptr);
    }
  }
  in.check();
  return results;
}

// The results of mtWorker in MolProcessing.cpp: in input order, up to the
// last record that has a fingerprint, which a record RDKit cannot read lacks.
using ByRecord = std::map<unsigned, std::unique_ptr<ExplicitBitVect>>;

Fingerprints in_record_order(ByRecord &by_record) {
  Fingerprints results;
  if (!by_record.empty()) {
    // Record ids start at 1 (MultithreadedMolSupplier::d_currentRecordId);
    // RDKit's own loop writes at id - 1 without a check.
    results.resize(by_record.rbegin()->first);
  }
  for (auto &[record, result] : by_record) {
    if (record >= 1) {
      results[record - 1] = std::move(result);
    }
  }
  return results;
}

// As mtWorker: RDKit's parsing threads fingerprint each molecule they read
// (its write callback). A molecule whose fingerprint throws has none: RDKit's
// threads catch everything their callback throws
// (MultithreadedMolSupplier::writer), so nothing escapes them to call
// std::terminate (CONVENTIONS.md S6).
template <class S, class F>
Fingerprints read_in_threads(const Input &in, S &suppl, const F &fp,
                             std::mutex &mutex, ByRecord &by_record) {
  suppl.setWriteCallback(
      [&](RDKit::RWMol &mol, const std::string &, unsigned int record) {
        auto result = fp(mol);
        std::lock_guard<std::mutex> lock(mutex);
        by_record[record] = std::move(result);
      });
  while (!suppl.atEnd()) {
    checked(in, [&] { suppl.next(); });
  }
  in.check();
  // Every thread has ended: atEnd is true only once the last one has.
  std::lock_guard<std::mutex> lock(mutex);
  return in_record_order(by_record);
}

// As mtWorker, for a generator that may construct a ControlCHandler, which
// RDKit's threads would construct under no lock (U4): they only parse, and
// this thread, which holds control_c_mutex, fingerprints the molecules as
// they come, skipping any whose fingerprint throws, as RDKit's threads do.
template <class S, class F>
Fingerprints read_in_threads_fingerprint_here(const Input &in, S &suppl,
                                              const F &fp) {
  ByRecord by_record;
  while (!suppl.atEnd()) {
    auto mol = checked(in, [&] { return suppl.next(); });
    if (!mol) {
      continue;
    }
    auto record = suppl.getLastRecordId();
    try {
      by_record[record] = fp(*mol);
    } catch (...) {
      // As MultithreadedMolSupplier::writer, which drops the molecule.
    }
  }
  in.check();
  return in_record_order(by_record);
}

// The parameters getSupplier reads each format with.
FP::MolFileParserParams sdf_params(const Options &options) {
  FP::MolFileParserParams params;
  params.sanitize = options.sanitize;
  params.removeHs = options.removeHs;
  params.strictParsing = options.strictParsing;
  return params;
}

FP::TDTMolSupplierParams tdt_params(const Options &options) {
  FP::TDTMolSupplierParams params;
  params.nameRecord = options.nameRecord;
  params.confId2D = options.confId2D;
  params.confId3D = options.confId3D;
  params.parseParameters.sanitize = options.sanitize;
  return params;
}

FP::SmilesMolSupplierParams smiles_params(const Options &options) {
  FP::SmilesMolSupplierParams params;
  params.delimiter = options.delimiter;
  params.smilesColumn = options.smilesColumn;
  params.nameColumn = options.nameColumn;
  params.titleLine = options.titleLine;
  params.parseParameters.sanitize = options.sanitize;
  return params;
}

// `control_c` holds control_c_mutex if the generator may construct a
// ControlCHandler (lock_control_c_for), for the whole call.
template <class Gen>
Fingerprints run(const std::string &path, const Gen &gen,
                 const FpGenerator *checked_gen,
                 std::unique_lock<ControlCMutex> control_c,
                 const Options &options) {
  std::string format;
  std::string compression;
  RDKit::GeneralMolSupplier::determineFormat(path, format, compression);
  // The records' relevant cycles are checked as they are read (streams.h).
  InputScanner scanner;
  if (format == "sdf") {
    scanner = InputScanner::mol_records(sdf_params(options), true);
  } else if (format == "tdt") {
    scanner = InputScanner::tdt_lines(tdt_params(options));
  } else if (format != "mae") {
    scanner = InputScanner::smiles_lines(smiles_params(options));
  }
  auto in = open_input(path, compression, std::move(scanner));
  if (format == "mae") {
    // getSupplier reads MAE files with RDKit's MaeMolSupplier, which trusts
    // maeparser's parse (see MaeMolSupplier in inventory/manifest.toml).
    throw ValueErrorException(
        "MAE files are not supported: RDKit's MaeMolSupplier trusts what "
        "maeparser parses");
  }
  // RDKit negates the count in int (getNumThreadsToUse), which is undefined
  // for INT_MIN; num_threads_to_use computes in 64 bits.
  if (options.numWriterThreads > static_cast<int>(kMaxSupplierThreads)) {
    throw ValueErrorException("numThreads must be at most " +
                              std::to_string(kMaxSupplierThreads));
  }
  unsigned num_threads = num_threads_to_use(options.numWriterThreads);
  bool threaded = num_threads > 1 && format != "tdt";
  Fingerprinter<Gen> fp{gen, checked_gen};
  auto read_threaded = [&](auto &suppl, std::mutex &mutex,
                           ByRecord &by_record) {
    return control_c.owns_lock()
               ? read_in_threads_fingerprint_here(in, suppl, fp)
               : read_in_threads(in, suppl, fp, mutex, by_record);
  };

  FP::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = num_threads;
  // Declared before the suppliers, whose threads write them until the
  // suppliers are destroyed.
  std::mutex mutex;
  ByRecord by_record;
  // The suppliers do not own the stream, which outlives them (as build in
  // suppliers.cpp).
  if (format == "sdf") {
    auto parse_params = sdf_params(options);
    if (threaded) {
      auto suppl = checked(in, [&] {
        return std::make_unique<FP::MultithreadedSDMolSupplier>(
            in.stream(), false, params, parse_params);
      });
      return read_threaded(*suppl, mutex, by_record);
    }
    auto suppl = checked(in, [&] {
      return std::make_unique<FP::ForwardSDMolSupplier>(in.stream(), false,
                                                        parse_params);
    });
    return read_in_order(in, *suppl, fp);
  }
  if (format == "tdt") {
    auto parse_params = tdt_params(options);
    auto suppl = checked(in, [&] {
      return std::make_unique<FP::TDTMolSupplier>(in.stream(), false,
                                                  parse_params);
    });
    return read_in_order(in, *suppl, fp);
  }
  // smi, csv, txt, and tsv; determineFormat knows no other.
  auto parse_params = smiles_params(options);
  if (threaded) {
    auto suppl = checked(in, [&] {
      return std::make_unique<FP::MultithreadedSmilesMolSupplier>(
          in.stream(), false, params, parse_params);
    });
    return read_threaded(*suppl, mutex, by_record);
  }
  auto suppl = checked(in, [&] {
    return std::make_unique<FP::SmilesMolSupplier>(in.stream(), false,
                                                   parse_params);
  });
  return read_in_order(in, *suppl, fp);
}

std::unique_ptr<FingerprintList> to_list(Fingerprints fps) {
  auto list = std::make_unique<FingerprintList>();
  for (auto &fp : fps) {
    if (fp) {
      list->bit_vects.push_back(wrap_bit_vect(std::move(fp)));
    } else {
      list->bit_vects.push_back(nullptr);
    }
  }
  return list;
}

}  // namespace

SupplierOptionsData supplier_options_default() {
  Options options;
  SupplierOptionsData data;
  data.num_threads = options.numWriterThreads;
  data.sanitize = options.sanitize;
  data.remove_hs = options.removeHs;
  data.strict_parsing = options.strictParsing;
  data.delimiter = rust::String(options.delimiter);
  data.smiles_column = options.smilesColumn;
  data.name_column = options.nameColumn;
  data.title_line = options.titleLine;
  data.name_record = rust::String(options.nameRecord);
  data.conf_id_2d = options.confId2D;
  data.conf_id_3d = options.confId3D;
  return data;
}

// Python's first overload: with no generator, RDKit's default, a 32-bit
// Morgan generator of radius 3, which may construct a ControlCHandler only
// through the CIP labeler (lock_control_c_for_cip).
std::unique_ptr<FingerprintList> get_fingerprints_for_mols_in_file(
    rust::Slice<const std::uint8_t> path, const SupplierOptionsData &options) {
  auto opts = to_options(options);
  std::unique_ptr<RDKit::FingerprintGenerator<std::uint32_t>> morgan(
      RDKit::MorganFingerprint::getMorganGenerator<std::uint32_t>(3));
  return to_list(
      run(to_path(path), *morgan, nullptr, lock_control_c_for_cip(), opts));
}

std::unique_ptr<FingerprintList>
get_fingerprints_for_mols_in_file_with_generator(
    rust::Slice<const std::uint8_t> path, const FpGenerator &generator,
    const SupplierOptionsData &options) {
  auto opts = to_options(options);
  return to_list(run(to_path(path), generator, &generator,
                     lock_control_c_for(generator), opts));
}

}  // namespace rdkit_cxx
