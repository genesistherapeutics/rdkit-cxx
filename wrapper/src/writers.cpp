// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/writers.h"

#include <algorithm>
#include <climits>
#include <iostream>
#include <mutex>
#include <utility>

#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/types.h>

#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/writers.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

// A streambuf that writes to a Rust `Write` value (Python's file objects)
// through a buffer it hands over when full and on sync. A failed write is
// kept, and fails every later one; WriterOutput::check reports it.
class RustWriteBuf : public std::streambuf {
 public:
  explicit RustWriteBuf(rust::Box<RustWriter> writer)
      : d_writer(std::move(writer)), d_buffer(new char[kSize]) {
    setp(d_buffer.get(), d_buffer.get() + kSize);
  }
  bool failed() const { return !d_error.empty(); }
  const std::string &error() const { return d_error; }

 protected:
  int_type overflow(int_type c) override {
    if (!drain()) {
      return traits_type::eof();
    }
    if (!traits_type::eq_int_type(c, traits_type::eof())) {
      *pptr() = traits_type::to_char_type(c);
      pbump(1);
    }
    return traits_type::not_eof(c);
  }

  int sync() override {
    if (!drain()) {
      return -1;
    }
    try {
      rust_writer_flush(*d_writer);
    } catch (const std::exception &e) {
      d_error = e.what();
      return -1;
    }
    return 0;
  }

 private:
  // Hands the buffered bytes to the writer.
  bool drain() {
    if (failed()) {
      return false;
    }
    auto n = static_cast<std::size_t>(pptr() - pbase());
    setp(d_buffer.get(), d_buffer.get() + kSize);
    if (n == 0) {
      return true;
    }
    try {
      rust_writer_write(
          *d_writer,
          rust::Slice<const std::uint8_t>(
              reinterpret_cast<const std::uint8_t *>(d_buffer.get()), n));
    } catch (const std::exception &e) {
      d_error = e.what();
      return false;
    }
    return true;
  }

  static constexpr std::size_t kSize = std::size_t{1} << 16;
  rust::Box<RustWriter> d_writer;
  std::unique_ptr<char[]> d_buffer;
  std::string d_error;
};

void WriterOutput::check() const {
  if (rust && rust->failed()) {
    throw RDKit::BadFileException("could not write the output: " +
                                  rust->error());
  }
  if (stream && stream->fail()) {
    throw RDKit::BadFileException("could not write the output");
  }
}

namespace {

// An ostream that owns its RustWriteBuf.
class RustOStream : public std::ostream {
 public:
  explicit RustOStream(rust::Box<RustWriter> writer)
      : std::ostream(nullptr), d_buf(std::move(writer)) {
    rdbuf(&d_buf);
  }
  const RustWriteBuf &buf() const { return d_buf; }

 private:
  RustWriteBuf d_buf;
};

// Opens `name` as RDKit's writers do: "Bad output file <name>" if it cannot
// be opened, and, when `dash_is_stdout`, standard output for "-".
WriterOutput open_output(rust::Slice<const std::uint8_t> path,
                         bool dash_is_stdout) {
  auto name = to_path(path);
  WriterOutput out;
  if (dash_is_stdout && name == "-") {
    out.stream =
        std::shared_ptr<std::ostream>(&std::cout, [](std::ostream *) {});
    return out;
  }
  auto file = std::make_shared<std::ofstream>(name);
  if (!*file || file->bad()) {
    throw RDKit::BadFileException("Bad output file " + name);
  }
  out.file = file.get();
  out.stream = std::move(file);
  return out;
}

WriterOutput rust_output(rust::Box<RustWriter> writer) {
  auto stream = std::make_shared<RustOStream>(std::move(writer));
  WriterOutput out;
  out.rust = &stream->buf();
  out.stream = std::move(stream);
  return out;
}

// Calls into RDKit with the output checked before and after; a failed output
// explains whatever RDKit did, including an exception.
template <class W, class F>
void checked(W &w, F &&f) {
  w.out.check();
  try {
    f();
  } catch (...) {
    w.out.check();
    throw;
  }
  w.out.check();
}

template <class W>
void flush_writer(W &w) {
  checked(w, [&] { w.writer->flush(); });
}

// As RDKit's close, which flushes and lets go of the stream; the wrapper then
// closes the file it opened, which reports a failed write. Closing again does
// nothing, as for RDKit.
template <class W>
void close_writer(W &w) {
  checked(w, [&] {
    w.writer->close();
    if (w.out.file && w.out.file->is_open()) {
      w.out.file->close();
    }
  });
}

RDKit::STR_VECT to_strings(const rust::Vec<rust::String> &props) {
  RDKit::STR_VECT result;
  result.reserve(props.size());
  for (const auto &prop : props) {
    result.emplace_back(prop);
  }
  return result;
}

// The SD and Maestro writers compute 2D coordinates for a molecule without a
// conformer (prepareMol in MolFileWriter.cpp, _MolToMaeCtBlock in
// MaeWriter.cpp), so they hold the depiction lock then (U4). They use
// compute2DCoords, which trusts the molecule (depiction_checks.h).
std::unique_lock<std::mutex> lock_if_depicting(const RWMol &mol) {
  if (mol.getNumConformers() == 0) {
    auto lock = lock_depiction();
    check_compute_2d_coords_mol(mol, false);
    return lock;
  }
  return {};
}

}  // namespace

// SDWriter

namespace {

// RDKit's SD writer numbers each property it writes as molid + 1 in int
// (_writePropToStream in SDWriter.cpp), which overflows for INT_MAX.
void check_molid(std::int64_t molid) {
  if (molid >= INT_MAX) {
    throw ValueErrorException(
        "the molecule id must be below INT_MAX: RDKit writes it plus 1 in int");
  }
}

std::unique_ptr<SDWriter> make_sd_writer(WriterOutput out) {
  auto w = std::make_unique<SDWriter>();
  w->writer = std::make_unique<RDKit::SDWriter>(out.stream.get(), false);
  w->out = std::move(out);
  return w;
}

}  // namespace

std::unique_ptr<SDWriter> sd_writer_new(rust::Slice<const std::uint8_t> path) {
  return make_sd_writer(open_output(path, true));
}

std::unique_ptr<SDWriter> sd_writer_from_writer(rust::Box<RustWriter> writer) {
  return make_sd_writer(rust_output(std::move(writer)));
}

void sd_writer_write(SDWriter &w, const RWMol &mol, std::int32_t conf_id) {
  check_molid(w.writer->numMols());
  checked(w, [&] {
    auto lock = lock_if_depicting(mol);
    w.writer->write(mol, conf_id);
  });
}

void sd_writer_flush(SDWriter &w) { flush_writer(w); }

void sd_writer_close(SDWriter &w) { close_writer(w); }

std::uint32_t sd_writer_num_mols(const SDWriter &w) {
  return w.writer->numMols();
}

void sd_writer_set_props(SDWriter &w, const rust::Vec<rust::String> &props) {
  w.writer->setProps(to_strings(props));
}

bool sd_writer_get_force_v3000(const SDWriter &w) {
  return w.writer->getForceV3000();
}

void sd_writer_set_force_v3000(SDWriter &w, bool val) {
  w.writer->setForceV3000(val);
}

bool sd_writer_get_kekulize(const SDWriter &w) {
  return w.writer->getKekulize();
}

void sd_writer_set_kekulize(SDWriter &w, bool val) {
  w.writer->setKekulize(val);
}

std::unique_ptr<std::string> sd_writer_get_text(const RWMol &mol,
                                                std::int32_t conf_id,
                                                bool kekulize, bool force_v3000,
                                                std::int32_t molid) {
  check_molid(molid);
  auto lock = lock_if_depicting(mol);
  return std::make_unique<std::string>(
      RDKit::SDWriter::getText(mol, conf_id, kekulize, force_v3000, molid));
}

// SmilesWriter

namespace {

std::unique_ptr<SmilesWriter> make_smiles_writer(
    WriterOutput out, const SmilesWriterParamsData &params) {
  auto w = std::make_unique<SmilesWriter>();
  w->writer = std::make_unique<RDKit::SmilesWriter>(
      out.stream.get(), std::string(params.delimiter),
      std::string(params.name_header), params.include_header, false,
      params.isomeric_smiles, params.kekule_smiles);
  w->out = std::move(out);
  return w;
}

}  // namespace

std::unique_ptr<SmilesWriter> smiles_writer_new(
    rust::Slice<const std::uint8_t> path,
    const SmilesWriterParamsData &params) {
  return make_smiles_writer(open_output(path, true), params);
}

std::unique_ptr<SmilesWriter> smiles_writer_from_writer(
    rust::Box<RustWriter> writer, const SmilesWriterParamsData &params) {
  return make_smiles_writer(rust_output(std::move(writer)), params);
}

// The SMILES and TDT writers write SMILES (mol_to_smiles).
void smiles_writer_write(SmilesWriter &w, const RWMol &mol,
                         std::int32_t conf_id) {
  check_removals(mol, Removes::Fragments);
  checked(w, [&] { w.writer->write(mol, conf_id); });
}

void smiles_writer_flush(SmilesWriter &w) { flush_writer(w); }

void smiles_writer_close(SmilesWriter &w) { close_writer(w); }

std::uint32_t smiles_writer_num_mols(const SmilesWriter &w) {
  return w.writer->numMols();
}

void smiles_writer_set_props(SmilesWriter &w,
                             const rust::Vec<rust::String> &props) {
  w.writer->setProps(to_strings(props));
}

// TDTWriter

namespace {

std::unique_ptr<TDTWriter> make_tdt_writer(WriterOutput out) {
  auto w = std::make_unique<TDTWriter>();
  w->writer = std::make_unique<RDKit::TDTWriter>(out.stream.get(), false);
  w->out = std::move(out);
  return w;
}

}  // namespace

std::unique_ptr<TDTWriter> tdt_writer_new(
    rust::Slice<const std::uint8_t> path) {
  return make_tdt_writer(open_output(path, true));
}

std::unique_ptr<TDTWriter> tdt_writer_from_writer(
    rust::Box<RustWriter> writer) {
  return make_tdt_writer(rust_output(std::move(writer)));
}

void tdt_writer_write(TDTWriter &w, const RWMol &mol, std::int32_t conf_id) {
  check_removals(mol, Removes::Fragments);
  checked(w, [&] { w.writer->write(mol, conf_id); });
}

void tdt_writer_flush(TDTWriter &w) { flush_writer(w); }

// RDKit's close also ends the last record with its "|" line.
void tdt_writer_close(TDTWriter &w) { close_writer(w); }

std::uint32_t tdt_writer_num_mols(const TDTWriter &w) {
  return w.writer->numMols();
}

void tdt_writer_set_props(TDTWriter &w, const rust::Vec<rust::String> &props) {
  w.writer->setProps(to_strings(props));
}

std::uint32_t tdt_writer_get_num_digits(const TDTWriter &w) {
  return w.writer->getNumDigits();
}

// RDKit formats each coordinate with std::setprecision(numDigits), which
// libstdc++ formats in a stack buffer as long as its output.
void tdt_writer_set_num_digits(TDTWriter &w, std::uint32_t num_digits) {
  check_precision(num_digits);
  w.writer->setNumDigits(num_digits);
}

bool tdt_writer_get_write_2d(const TDTWriter &w) {
  return w.writer->getWrite2D();
}

void tdt_writer_set_write_2d(TDTWriter &w, bool state) {
  w.writer->setWrite2D(state);
}

bool tdt_writer_get_write_names(const TDTWriter &w) {
  return w.writer->getWriteNames();
}

void tdt_writer_set_write_names(TDTWriter &w, bool state) {
  w.writer->setWriteNames(state);
}

// PDBWriter

namespace {

std::unique_ptr<PDBWriter> make_pdb_writer(WriterOutput out,
                                           std::uint32_t flavor) {
  auto w = std::make_unique<PDBWriter>();
  w->writer =
      std::make_unique<RDKit::PDBWriter>(out.stream.get(), false, flavor);
  w->out = std::move(out);
  return w;
}

}  // namespace

std::unique_ptr<PDBWriter> pdb_writer_new(rust::Slice<const std::uint8_t> path,
                                          std::uint32_t flavor) {
  return make_pdb_writer(open_output(path, true), flavor);
}

std::unique_ptr<PDBWriter> pdb_writer_from_writer(rust::Box<RustWriter> writer,
                                                  std::uint32_t flavor) {
  return make_pdb_writer(rust_output(std::move(writer)), flavor);
}

void pdb_writer_write(PDBWriter &w, const RWMol &mol, std::int32_t conf_id) {
  checked(w, [&] { w.writer->write(mol, conf_id); });
}

void pdb_writer_flush(PDBWriter &w) { flush_writer(w); }

void pdb_writer_close(PDBWriter &w) { close_writer(w); }

std::uint32_t pdb_writer_num_mols(const PDBWriter &w) {
  return w.writer->numMols();
}

// MaeWriter

namespace {

// The longest property name the Maestro writer takes. It matches the name of
// every property it exports with std::regex ("[birs]_[^_ ]+_.+",
// copyProperties in MaeWriter.cpp), whose matcher recurses once per
// character, so a long name overflows the stack (Python crashes at 40,000).
constexpr std::size_t kMaxMaeProp = 1000;

void check_mae_props(const RDKit::RDProps &obj, const RDKit::STR_VECT &props) {
  RDKit::STR_VECT computed;
  obj.getPropIfPresent(RDKit::detail::computedPropName, computed);
  for (const auto &item : obj.getDict().getData()) {
    if (item.key.size() <= kMaxMaeProp ||
        item.key == RDKit::detail::computedPropName ||
        std::find(computed.begin(), computed.end(), item.key) !=
            computed.end() ||
        (!props.empty() &&
         std::find(props.begin(), props.end(), item.key) == props.end())) {
      continue;
    }
    throw ValueErrorException(
        "a property name is longer than " + std::to_string(kMaxMaeProp) +
        " characters; RDKit's Maestro writer matches it with a regex that "
        "recurses once per character");
  }
}

// The Maestro writer casts every atom's monomer info to AtomPDBResidueInfo
// (mapAtom in MaeWriter.cpp), whatever its type, which reads past a plain
// AtomMonomerInfo; and it matches property names as above.
void check_mae(const RWMol &mol, const RDKit::STR_VECT &props) {
  for (const auto atom : mol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    if (info && !dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info)) {
      throw ValueErrorException(
          "atom " + std::to_string(atom->getIdx()) +
          " has monomer info that is not an AtomPDBResidueInfo, which RDKit's "
          "Maestro writer would read as one");
    }
  }
  check_mae_props(mol, props);
  for (const auto atom : mol.atoms()) {
    check_mae_props(*atom, props);
  }
  for (const auto bond : mol.bonds()) {
    check_mae_props(*bond, props);
  }
}

// RDKit 2024.09's MaeWriter(std::shared_ptr<std::ostream>) checks the
// pointer it has just moved from, and always fails; MaeWriter(std::ostream *)
// takes ownership of its stream, and deletes it when it closes. So it gets a
// stream of its own over the output's buffer, which `source` keeps alive
// past the writer, and the output's checks read that stream until the
// writer closes (mae_writer_close), then the source.
std::unique_ptr<MaeWriter> make_mae_writer(WriterOutput out) {
  auto w = std::make_unique<MaeWriter>();
  auto *stream = new std::ostream(out.stream->rdbuf());
  w->source = out.stream;
  // The writer deletes the stream, even when its constructor throws.
  w->writer = std::make_unique<RDKit::MaeWriter>(stream);
  out.stream = std::shared_ptr<std::ostream>(stream, [](std::ostream *) {});
  w->out = std::move(out);
  return w;
}

}  // namespace

std::unique_ptr<MaeWriter> mae_writer_new(
    rust::Slice<const std::uint8_t> path) {
  return make_mae_writer(open_output(path, false));
}

std::unique_ptr<MaeWriter> mae_writer_from_writer(
    rust::Box<RustWriter> writer) {
  return make_mae_writer(rust_output(std::move(writer)));
}

void mae_writer_write(MaeWriter &w, const RWMol &mol, std::int32_t conf_id) {
  check_mae(mol, w.props);
  checked(w, [&] {
    auto lock = lock_if_depicting(mol);
    w.writer->write(mol, conf_id);
  });
}

void mae_writer_flush(MaeWriter &w) { flush_writer(w); }

// As close_writer, but the writer deletes its stream as it closes, so the
// output's checks read the source from then on.
void mae_writer_close(MaeWriter &w) {
  w.out.check();
  try {
    w.writer->close();
  } catch (...) {
    w.out.stream = w.source;
    w.out.check();
    throw;
  }
  w.out.stream = w.source;
  checked(w, [&] {
    if (w.out.file && w.out.file->is_open()) {
      w.out.file->close();
    }
  });
}

std::uint32_t mae_writer_num_mols(const MaeWriter &w) {
  return w.writer->numMols();
}

void mae_writer_set_props(MaeWriter &w, const rust::Vec<rust::String> &props) {
  w.props = to_strings(props);
  w.writer->setProps(w.props);
}

std::unique_ptr<std::string> mae_writer_get_text(
    const RWMol &mol, std::int32_t conf_id,
    const rust::Vec<rust::String> &props) {
  auto names = to_strings(props);
  check_mae(mol, names);
  auto lock = lock_if_depicting(mol);
  return std::make_unique<std::string>(
      RDKit::MaeWriter::getText(mol, conf_id, names));
}

}  // namespace rdkit_cxx
