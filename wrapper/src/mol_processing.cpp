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

#include <GraphMol/Chirality.h>
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
Input open_input(const std::string &path, const std::string &compression) {
  Input in;
  if (compression.empty()) {
    in.checked = open_checked_file(path, true);
    return in;
  }
#ifdef RDK_USE_BOOST_IOSTREAMS
  in.checked =
      std::make_unique<CheckedIStream>(std::make_unique<GzipFileBuf>(path));
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
// around the call (U6). `lock` holds control_c_mutex for the call, for a
// generator that may construct a ControlCHandler, when the call runs on one
// of RDKit's threads.
template <class Gen>
struct Fingerprinter {
  const Gen &gen;
  // The generator to check molecules for; null for the default one, whose
  // fingerprints trust nothing check_fingerprint_mol checks.
  const FpGenerator *checked_gen;
  bool lock;

  std::unique_ptr<ExplicitBitVect> operator()(RDKit::ROMol &mol) const {
    if (checked_gen) {
      check_fingerprint_mol(*checked_gen, mol);
    }
    MatrixCacheGuard guard(mol);
    std::unique_lock<std::mutex> control_c;
    if (lock) {
      control_c = std::unique_lock<std::mutex>(control_c_mutex());
    }
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

// As mtWorker in MolProcessing.cpp: RDKit's parsing threads fingerprint each
// molecule they read (its writer callback), and the results come back in
// input order, up to the last record that has one. A record RDKit cannot read
// has none, and nor does one whose fingerprint throws: RDKit's threads catch
// everything their callback throws (MultithreadedMolSupplier::writer), so
// nothing escapes them to call std::terminate.
template <class S, class F>
Fingerprints read_in_threads(
    const Input &in, S &suppl, const F &fp, std::mutex &mutex,
    std::map<unsigned, std::unique_ptr<ExplicitBitVect>> &by_record) {
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

template <class Gen>
Fingerprints run(const std::string &path, const Gen &gen,
                 const FpGenerator *checked_gen, bool needs_lock,
                 const Options &options) {
  std::string format;
  std::string compression;
  RDKit::GeneralMolSupplier::determineFormat(path, format, compression);
  auto in = open_input(path, compression);
  if (format == "mae") {
    // getSupplier reads MAE files with RDKit's MaeMolSupplier, which trusts
    // maeparser's parse (see MaeMolSupplier in inventory/manifest.toml).
    throw ValueErrorException(
        "MAE files are not supported: RDKit's MaeMolSupplier trusts what "
        "maeparser parses");
  }
  // RDKit negates the count in int (getNumThreadsToUse), which is undefined
  // for INT_MIN; num_threads_to_use computes in 64 bits.
  unsigned num_threads = num_threads_to_use(options.numWriterThreads);
  bool threaded = num_threads > 1 && format != "tdt";
  if (needs_lock) {
    // A closure RDKit calls inside a match may not wait for threads that
    // need the lock it holds.
    check_not_in_control_c_callback();
  }
  // RDKit's threads take the lock for each fingerprint; this thread, which
  // the lock can order, takes it for the whole call.
  std::unique_lock<std::mutex> control_c;
  if (needs_lock && !threaded) {
    control_c = std::unique_lock<std::mutex>(control_c_mutex());
  }
  Fingerprinter<Gen> fp{gen, checked_gen, needs_lock && threaded};

  FP::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = num_threads;
  // Declared before the suppliers, whose threads write them until the
  // suppliers are destroyed.
  std::mutex mutex;
  std::map<unsigned, std::unique_ptr<ExplicitBitVect>> by_record;
  // The suppliers do not own the stream, which outlives them (as build in
  // suppliers.cpp).
  if (format == "sdf") {
    FP::MolFileParserParams parse_params;
    parse_params.sanitize = options.sanitize;
    parse_params.removeHs = options.removeHs;
    parse_params.strictParsing = options.strictParsing;
    if (threaded) {
      auto suppl = checked(in, [&] {
        return std::make_unique<FP::MultithreadedSDMolSupplier>(
            in.stream(), false, params, parse_params);
      });
      return read_in_threads(in, *suppl, fp, mutex, by_record);
    }
    auto suppl = checked(in, [&] {
      return std::make_unique<FP::ForwardSDMolSupplier>(in.stream(), false,
                                                        parse_params);
    });
    return read_in_order(in, *suppl, fp);
  }
  if (format == "tdt") {
    FP::TDTMolSupplierParams parse_params;
    parse_params.nameRecord = options.nameRecord;
    parse_params.confId2D = options.confId2D;
    parse_params.confId3D = options.confId3D;
    parse_params.parseParameters.sanitize = options.sanitize;
    auto suppl = checked(in, [&] {
      return std::make_unique<FP::TDTMolSupplier>(in.stream(), false,
                                                  parse_params);
    });
    return read_in_order(in, *suppl, fp);
  }
  // smi, csv, txt, and tsv; determineFormat knows no other.
  FP::SmilesMolSupplierParams parse_params;
  parse_params.delimiter = options.delimiter;
  parse_params.smilesColumn = options.smilesColumn;
  parse_params.nameColumn = options.nameColumn;
  parse_params.titleLine = options.titleLine;
  parse_params.parseParameters.sanitize = options.sanitize;
  if (threaded) {
    auto suppl = checked(in, [&] {
      return std::make_unique<FP::MultithreadedSmilesMolSupplier>(
          in.stream(), false, params, parse_params);
    });
    return read_in_threads(in, *suppl, fp, mutex, by_record);
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
  bool needs_lock = !RDKit::Chirality::getUseLegacyStereoPerception();
  return to_list(run(to_path(path), *morgan, nullptr, needs_lock, opts));
}

// A generator may construct a ControlCHandler when it finds feature
// invariants, by substructure matches (MorganFeatureAtomInvGenerator, which
// names itself in the generator's info string), and when it runs the CIP
// labeler (lock_control_c_for_cip).
std::unique_ptr<FingerprintList>
get_fingerprints_for_mols_in_file_with_generator(
    rust::Slice<const std::uint8_t> path, const FpGenerator &generator,
    const SupplierOptionsData &options) {
  auto opts = to_options(options);
  bool needs_lock =
      generator.infoString().find("MorganFeatureInvariantGenerator") !=
          std::string::npos ||
      !RDKit::Chirality::getUseLegacyStereoPerception();
  return to_list(run(to_path(path), generator, &generator, needs_lock, opts));
}

}  // namespace rdkit_cxx
