/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// FPB fingerprint files. Fingerprints cross as bytes, which are checked
// against the file's size; RDKit reads them without a length.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <DataStructs/FPBReader.h>
#include <DataStructs/MultiFPBReader.h>

#include "rdkit-cxx/data_structs.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using FPBReader = RDKit::FPBReader;
using MultiFPBReader = RDKit::MultiFPBReader;

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

// The readers stay owned by Rust (src/data_structs/c_data_structs/fpb.rs),
// which keeps them alive while the MultiFPBReader holds their addresses.
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
