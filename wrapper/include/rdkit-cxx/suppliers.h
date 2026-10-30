/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolfiles's molecule suppliers. Each wrapper owns an RDKit supplier that
// reads from a checked stream (CONVENTIONS.md U12), and checks the stream after
// every call into RDKit. RDKit's MaeMolSupplier is not bound (see
// inventory/manifest.toml).

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MultithreadedSDMolSupplier.h>
#include <GraphMol/FileParsers/MultithreadedSmilesMolSupplier.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/streams.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct RustReader;
class RustReadBuf;

// The input of a supplier: the stream RDKit reads, and the Rust reader under
// it, if any. The wrapper owns the stream, not the RDKit supplier (see build
// in suppliers.cpp); each wrapper declares it before the supplier, which is
// then destroyed first.
struct SupplierInput {
  std::unique_ptr<CheckedIStream> stream;
  const RustReadBuf *reader = nullptr;
  // Throws BadFileException for a failed Rust reader, then
  // ValueErrorException for refused input.
  void check() const;
};

class SDMolSupplier {
 public:
  SupplierInput input;
  std::unique_ptr<RDKit::v2::FileParsers::SDMolSupplier> suppl;
};

class ForwardSDMolSupplier {
 public:
  SupplierInput input;
  std::unique_ptr<RDKit::v2::FileParsers::ForwardSDMolSupplier> suppl;
};

class SmilesMolSupplier {
 public:
  SupplierInput input;
  std::unique_ptr<RDKit::v2::FileParsers::SmilesMolSupplier> suppl;
};

// RDKit's TDTMolSupplier is not used across calls (see suppliers.cpp): the
// wrapper keeps the text and where its records start, and reads each record
// with a fresh RDKit supplier.
class TDTMolSupplier {
 public:
  std::string text;
  // The offsets of the lines that start with "$SMI<".
  std::vector<std::size_t> starts;
  RDKit::v2::FileParsers::TDTMolSupplierParams params;
  // The record next() reads, and whether the supplier is at its end.
  std::size_t next = 0;
  bool at_end = true;
};

class MultithreadedSDMolSupplier {
 public:
  ~MultithreadedSDMolSupplier();
  SupplierInput input;
  std::unique_ptr<RDKit::v2::FileParsers::MultithreadedSDMolSupplier> suppl;
  bool started = false;
  // Whether starting its threads failed (suppliers.cpp, mt_next).
  bool broken = false;
};

class MultithreadedSmilesMolSupplier {
 public:
  ~MultithreadedSmilesMolSupplier();
  SupplierInput input;
  std::unique_ptr<RDKit::v2::FileParsers::MultithreadedSmilesMolSupplier> suppl;
  bool broken = false;
};

// Shared structs, defined by the bridge (src/ffi/suppliers.rs).
struct SmilesMolSupplierParamsData;
struct TDTMolSupplierParamsData;
struct MultithreadedParamsData;

// SDMolSupplier. Text is copied; files are read through checked streams.
std::unique_ptr<SDMolSupplier> sd_mol_supplier_new();
std::unique_ptr<SDMolSupplier> sd_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing);
void sd_mol_supplier_set_data(SDMolSupplier &s,
                              rust::Slice<const std::uint8_t> data,
                              bool sanitize, bool remove_hs,
                              bool strict_parsing);
// Null for a record RDKit cannot read.
std::unique_ptr<RWMol> sd_mol_supplier_next(SDMolSupplier &s);
// Throws IndexError past the last record.
std::unique_ptr<RWMol> sd_mol_supplier_get(SDMolSupplier &s, std::uint32_t idx);
std::unique_ptr<std::string> sd_mol_supplier_get_item_text(SDMolSupplier &s,
                                                           std::uint32_t idx);
std::uint32_t sd_mol_supplier_length(SDMolSupplier &s);
bool sd_mol_supplier_at_end(SDMolSupplier &s);
void sd_mol_supplier_reset(SDMolSupplier &s);
bool sd_mol_supplier_get_process_property_lists(const SDMolSupplier &s);
void sd_mol_supplier_set_process_property_lists(SDMolSupplier &s, bool val);

// ForwardSDMolSupplier, from a file or a Rust reader.
std::unique_ptr<ForwardSDMolSupplier> forward_sd_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing);
std::unique_ptr<ForwardSDMolSupplier> forward_sd_mol_supplier_from_reader(
    rust::Box<RustReader> reader, bool sanitize, bool remove_hs,
    bool strict_parsing);
std::unique_ptr<RWMol> forward_sd_mol_supplier_next(ForwardSDMolSupplier &s);
bool forward_sd_mol_supplier_at_end(ForwardSDMolSupplier &s);
bool forward_sd_mol_supplier_get_eof_hit_on_read(const ForwardSDMolSupplier &s);
bool forward_sd_mol_supplier_get_process_property_lists(
    const ForwardSDMolSupplier &s);
void forward_sd_mol_supplier_set_process_property_lists(ForwardSDMolSupplier &s,
                                                        bool val);

// SmilesMolSupplier; a name column below -1 is refused.
std::unique_ptr<SmilesMolSupplier> smiles_mol_supplier_new();
std::unique_ptr<SmilesMolSupplier> smiles_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path,
    const SmilesMolSupplierParamsData &params);
void smiles_mol_supplier_set_data(SmilesMolSupplier &s, rust::Str data,
                                  const SmilesMolSupplierParamsData &params);
std::unique_ptr<RWMol> smiles_mol_supplier_next(SmilesMolSupplier &s);
std::unique_ptr<RWMol> smiles_mol_supplier_get(SmilesMolSupplier &s,
                                               std::uint32_t idx);
std::unique_ptr<std::string> smiles_mol_supplier_get_item_text(
    SmilesMolSupplier &s, std::uint32_t idx);
std::uint32_t smiles_mol_supplier_length(SmilesMolSupplier &s);
bool smiles_mol_supplier_at_end(SmilesMolSupplier &s);
void smiles_mol_supplier_reset(SmilesMolSupplier &s);

// TDTMolSupplier, which reads its input into memory.
std::unique_ptr<TDTMolSupplier> tdt_mol_supplier_new();
std::unique_ptr<TDTMolSupplier> tdt_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path,
    const TDTMolSupplierParamsData &params);
void tdt_mol_supplier_set_data(TDTMolSupplier &s, rust::Str data,
                               const TDTMolSupplierParamsData &params);
std::unique_ptr<RWMol> tdt_mol_supplier_next(TDTMolSupplier &s);
std::unique_ptr<RWMol> tdt_mol_supplier_get(TDTMolSupplier &s,
                                            std::uint32_t idx);
std::unique_ptr<std::string> tdt_mol_supplier_get_item_text(TDTMolSupplier &s,
                                                            std::uint32_t idx);
std::uint32_t tdt_mol_supplier_length(TDTMolSupplier &s);
bool tdt_mol_supplier_at_end(TDTMolSupplier &s);
void tdt_mol_supplier_reset(TDTMolSupplier &s);

// The multithreaded suppliers. next returns null at the end, and sets
// `record` to the record id, or 0 at the end (a null molecule with a nonzero
// id is a record RDKit could not read).
std::unique_ptr<MultithreadedSDMolSupplier> mt_sd_mol_supplier_new();
std::unique_ptr<MultithreadedSDMolSupplier> mt_sd_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing, const MultithreadedParamsData &params);
std::unique_ptr<RWMol> mt_sd_mol_supplier_next(MultithreadedSDMolSupplier &s,
                                               std::uint32_t &record);
bool mt_sd_mol_supplier_at_end(MultithreadedSDMolSupplier &s);
std::uint32_t mt_sd_mol_supplier_get_last_record_id(
    const MultithreadedSDMolSupplier &s);
std::unique_ptr<std::string> mt_sd_mol_supplier_get_last_item_text(
    const MultithreadedSDMolSupplier &s);
bool mt_sd_mol_supplier_get_process_property_lists(
    const MultithreadedSDMolSupplier &s);
void mt_sd_mol_supplier_set_process_property_lists(
    MultithreadedSDMolSupplier &s, bool val);

std::unique_ptr<MultithreadedSmilesMolSupplier> mt_smiles_mol_supplier_new();
std::unique_ptr<MultithreadedSmilesMolSupplier>
mt_smiles_mol_supplier_from_file(rust::Slice<const std::uint8_t> path,
                                 const SmilesMolSupplierParamsData &smiles,
                                 const MultithreadedParamsData &params);
std::unique_ptr<RWMol> mt_smiles_mol_supplier_next(
    MultithreadedSmilesMolSupplier &s, std::uint32_t &record);
bool mt_smiles_mol_supplier_at_end(MultithreadedSmilesMolSupplier &s);
std::uint32_t mt_smiles_mol_supplier_get_last_record_id(
    const MultithreadedSmilesMolSupplier &s);
std::unique_ptr<std::string> mt_smiles_mol_supplier_get_last_item_text(
    const MultithreadedSmilesMolSupplier &s);

}  // namespace rdkit_cxx
