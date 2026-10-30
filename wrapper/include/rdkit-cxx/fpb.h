/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// FPB fingerprint files. RDKit trusts every size, count, and offset in an
// FPB file, so the wrapper reads the file itself, checks its chunks, and
// gives RDKit only the checked bytes (CONVENTIONS.md U2). Fingerprints cross
// as bytes, checked against the file's fingerprint size and padded to its
// storage size, which is what RDKit reads.

#pragma once

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <DataStructs/FPBReader.h>

#include "rdkit-cxx/data_structs.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// What the checks learn about a file, for the reads RDKit leaves unchecked.
struct FpbLayout {
  std::uint32_t num_bytes = 0;     // bytes in a fingerprint
  std::uint32_t storage_size = 0;  // bytes per fingerprint in the file
  std::uint64_t num_pop_counts = 0;
};

// An RDKit::FPBReader over a checked copy of its file.
class FPBReader {
 public:
  FPBReader(std::string data, bool lazy);

  void init();
  const RDKit::FPBReader &rdkit() const { return *d_reader; }
  const FpbLayout &layout() const { return d_layout; }

  // The query padded to the storage size; throws unless it has num_bytes.
  std::vector<std::uint8_t> query(rust::Slice<const std::uint8_t> bytes) const;
  // Whether RDKit's popcount screen indexes its table in range for these
  // bounds, which it computes, and casts, unchecked.
  bool screen_in_range(double min_count, double max_count) const;

 private:
  FpbLayout d_layout;
  bool d_lazy;
  std::istringstream *d_stream;  // owned by d_reader
  std::unique_ptr<RDKit::FPBReader> d_reader;
};

// RDKit's MultiFPBReader searches its readers with the popcount screen
// always on, so the wrapper searches them itself, as RDKit does but through
// the checks above. The readers stay owned by Rust
// (src/data_structs/c_data_structs/fpb.rs), which keeps them alive while
// this holds their addresses.
class MultiFPBReader {
 public:
  explicit MultiFPBReader(bool init_on_search)
      : d_init_on_search(init_on_search) {}

  std::uint32_t add_reader(FPBReader &reader);
  void init();
  std::uint32_t num_bits() const;
  void check_searchable() const;
  bool init_on_search() const { return d_init_on_search; }
  const std::vector<FPBReader *> &readers() const { return d_readers; }

 private:
  std::vector<FPBReader *> d_readers;
  bool d_init_on_search;
  bool d_init = false;
};

// Shared structs, defined by the bridge (src/ffi/fpb.rs).
struct FPBNeighbor;
struct MultiFPBNeighbor;
struct MultiFPBContaining;

std::unique_ptr<FPBReader> fpb_reader_new(rust::Slice<const std::uint8_t> path,
                                          bool lazy);
void fpb_reader_init(FPBReader &reader);
std::uint32_t fpb_reader_length(const FPBReader &reader);
std::uint32_t fpb_reader_get_num_bits(const FPBReader &reader);
std::unique_ptr<BitVect> fpb_reader_get_fp(const FPBReader &reader,
                                           std::uint32_t idx);
rust::Vec<std::uint8_t> fpb_reader_get_bytes(const FPBReader &reader,
                                             std::uint32_t which);
std::unique_ptr<std::string> fpb_reader_get_id(const FPBReader &reader,
                                               std::uint32_t idx);
double fpb_reader_get_tanimoto(const FPBReader &reader, std::uint32_t which,
                               rust::Slice<const std::uint8_t> bytes);
double fpb_reader_get_tversky(const FPBReader &reader, std::uint32_t which,
                              rust::Slice<const std::uint8_t> bytes, double ca,
                              double cb);
rust::Vec<FPBNeighbor> fpb_reader_get_tanimoto_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv,
    double threshold);
rust::Vec<FPBNeighbor> fpb_reader_get_tversky_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold);
rust::Vec<std::uint32_t> fpb_reader_get_containing_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv);

std::unique_ptr<MultiFPBReader> multi_fpb_reader_new(bool init_on_search);
std::uint32_t multi_fpb_reader_add_reader(MultiFPBReader &multi,
                                          FPBReader &reader);
void multi_fpb_reader_init(MultiFPBReader &multi);
std::uint32_t multi_fpb_reader_get_num_bits(const MultiFPBReader &multi);
rust::Vec<MultiFPBNeighbor> multi_fpb_reader_get_tanimoto_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv,
    double threshold, std::uint32_t num_threads);
rust::Vec<MultiFPBNeighbor> multi_fpb_reader_get_tversky_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold, std::uint32_t num_threads);
rust::Vec<MultiFPBContaining> multi_fpb_reader_get_containing_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv,
    std::uint32_t num_threads);

}  // namespace rdkit_cxx
