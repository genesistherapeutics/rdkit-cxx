/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolfiles's molecule writers. Each wrapper owns its output (a file it
// opens, standard output for "-", or a Rust writer) and an RDKit writer that
// writes to it without owning it, and checks the output after every call into
// RDKit: RDKit's writers ignore a stream that fails.

#pragma once

#include <cstdint>
#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

// MolWriters.h declares MaeWriter when RDConfig.h, included first, defines
// RDK_BUILD_MAEPARSER_SUPPORT.
// clang-format off
#include <RDGeneral/RDConfig.h>
#include <GraphMol/FileParsers/MolWriters.h>
// clang-format on

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct RustWriter;
class RustWriteBuf;

// The output of a writer.
struct WriterOutput {
  std::shared_ptr<std::ostream> stream;
  // The file, when the wrapper opened one; the Rust writer, if any.
  std::ofstream *file = nullptr;
  const RustWriteBuf *rust = nullptr;
  // A stream of the RDKit writer's own over `stream`'s buffer, if it has one
  // (MaeWriter, in writers.cpp), while the writer holds it.
  const std::ostream *writer_stream = nullptr;
  // Throws BadFileException if the output has failed.
  void check() const;
};

// Each wrapper declares its output before the writer, which is then
// destroyed (and closed) first.
class SDWriter {
 public:
  WriterOutput out;
  std::unique_ptr<RDKit::SDWriter> writer;
};

class SmilesWriter {
 public:
  WriterOutput out;
  std::unique_ptr<RDKit::SmilesWriter> writer;
};

class TDTWriter {
 public:
  WriterOutput out;
  std::unique_ptr<RDKit::TDTWriter> writer;
};

class PDBWriter {
 public:
  WriterOutput out;
  std::unique_ptr<RDKit::PDBWriter> writer;
};

class MaeWriter {
 public:
  WriterOutput out;
  std::unique_ptr<RDKit::MaeWriter> writer;
  // The properties setProps gave the writer, which its checks need.
  std::vector<std::string> props;
};

// Shared struct, defined by the bridge (src/ffi/writers.rs).
struct SmilesWriterParamsData;

// SDWriter; "-" is standard output.
std::unique_ptr<SDWriter> sd_writer_new(rust::Slice<const std::uint8_t> path);
std::unique_ptr<SDWriter> sd_writer_from_writer(rust::Box<RustWriter> writer);
void sd_writer_write(SDWriter &w, const RWMol &mol, std::int32_t conf_id);
void sd_writer_flush(SDWriter &w);
void sd_writer_close(SDWriter &w);
std::uint32_t sd_writer_num_mols(const SDWriter &w);
void sd_writer_set_props(SDWriter &w, const rust::Vec<rust::String> &props);
bool sd_writer_get_force_v3000(const SDWriter &w);
void sd_writer_set_force_v3000(SDWriter &w, bool val);
bool sd_writer_get_kekulize(const SDWriter &w);
void sd_writer_set_kekulize(SDWriter &w, bool val);
std::unique_ptr<std::string> sd_writer_get_text(const RWMol &mol,
                                                std::int32_t conf_id,
                                                bool kekulize, bool force_v3000,
                                                std::int32_t molid);

// SmilesWriter; "-" is standard output.
std::unique_ptr<SmilesWriter> smiles_writer_new(
    rust::Slice<const std::uint8_t> path, const SmilesWriterParamsData &params);
std::unique_ptr<SmilesWriter> smiles_writer_from_writer(
    rust::Box<RustWriter> writer, const SmilesWriterParamsData &params);
void smiles_writer_write(SmilesWriter &w, const RWMol &mol,
                         std::int32_t conf_id);
void smiles_writer_flush(SmilesWriter &w);
void smiles_writer_close(SmilesWriter &w);
std::uint32_t smiles_writer_num_mols(const SmilesWriter &w);
void smiles_writer_set_props(SmilesWriter &w,
                             const rust::Vec<rust::String> &props);

// TDTWriter; "-" is standard output.
std::unique_ptr<TDTWriter> tdt_writer_new(rust::Slice<const std::uint8_t> path);
std::unique_ptr<TDTWriter> tdt_writer_from_writer(rust::Box<RustWriter> writer);
void tdt_writer_write(TDTWriter &w, const RWMol &mol, std::int32_t conf_id);
void tdt_writer_flush(TDTWriter &w);
void tdt_writer_close(TDTWriter &w);
std::uint32_t tdt_writer_num_mols(const TDTWriter &w);
void tdt_writer_set_props(TDTWriter &w, const rust::Vec<rust::String> &props);
std::uint32_t tdt_writer_get_num_digits(const TDTWriter &w);
// A number of digits past kMaxPrecision is refused (file_checks.h).
void tdt_writer_set_num_digits(TDTWriter &w, std::uint32_t num_digits);
bool tdt_writer_get_write_2d(const TDTWriter &w);
void tdt_writer_set_write_2d(TDTWriter &w, bool state);
bool tdt_writer_get_write_names(const TDTWriter &w);
void tdt_writer_set_write_names(TDTWriter &w, bool state);

// PDBWriter; "-" is standard output.
std::unique_ptr<PDBWriter> pdb_writer_new(rust::Slice<const std::uint8_t> path,
                                          std::uint32_t flavor);
std::unique_ptr<PDBWriter> pdb_writer_from_writer(rust::Box<RustWriter> writer,
                                                  std::uint32_t flavor);
void pdb_writer_write(PDBWriter &w, const RWMol &mol, std::int32_t conf_id);
void pdb_writer_flush(PDBWriter &w);
void pdb_writer_close(PDBWriter &w);
std::uint32_t pdb_writer_num_mols(const PDBWriter &w);

// MaeWriter; "-" is a file name, as for RDKit's.
std::unique_ptr<MaeWriter> mae_writer_new(rust::Slice<const std::uint8_t> path);
std::unique_ptr<MaeWriter> mae_writer_from_writer(rust::Box<RustWriter> writer);
void mae_writer_write(MaeWriter &w, const RWMol &mol, std::int32_t conf_id);
void mae_writer_flush(MaeWriter &w);
void mae_writer_close(MaeWriter &w);
std::uint32_t mae_writer_num_mols(const MaeWriter &w);
void mae_writer_set_props(MaeWriter &w, const rust::Vec<rust::String> &props);
std::unique_ptr<std::string> mae_writer_get_text(
    const RWMol &mol, std::int32_t conf_id,
    const rust::Vec<rust::String> &props);

}  // namespace rdkit_cxx
