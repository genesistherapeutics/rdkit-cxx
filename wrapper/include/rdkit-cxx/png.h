/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolfiles's PNG metadata functions.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridge (src/ffi/png.rs).
struct PNGMetadataParamsData;
struct SmilesParserParamsData;
struct SmilesWriteParamsData;

// Metadata as parallel key and value lists, in file order.
class PNGMetadata {
 public:
  std::vector<std::pair<std::string, std::string>> items;
};

// Checks the chunk structure RDKit walks without checking it (check_png in
// png.cpp); the file forms read the file once and check those bytes (U2).
std::unique_ptr<RWMol> mol_from_png_string(
    rust::Slice<const std::uint8_t> png, const SmilesParserParamsData &params);
std::unique_ptr<RWMol> mol_from_png_file(rust::Slice<const std::uint8_t> path,
                                         const SmilesParserParamsData &params);
// Null entries are records RDKit could not read.
std::unique_ptr<MolList> mols_from_png_string(
    rust::Slice<const std::uint8_t> png, rust::Str tag,
    const SmilesParserParamsData &params);
std::unique_ptr<MolList> mols_from_png_file(
    rust::Slice<const std::uint8_t> path, rust::Str tag,
    const SmilesParserParamsData &params);
rust::Vec<std::uint8_t> mol_metadata_to_png_string(
    const RWMol &mol, rust::Slice<const std::uint8_t> png,
    const PNGMetadataParamsData &params);
rust::Vec<std::uint8_t> mol_metadata_to_png_file(
    const RWMol &mol, rust::Slice<const std::uint8_t> path,
    const PNGMetadataParamsData &params);
rust::Vec<std::uint8_t> add_metadata_to_png_string(
    const rust::Vec<rust::String> &keys, const rust::Vec<rust::String> &values,
    rust::Slice<const std::uint8_t> png);
rust::Vec<std::uint8_t> add_metadata_to_png_file(
    const rust::Vec<rust::String> &keys, const rust::Vec<rust::String> &values,
    rust::Slice<const std::uint8_t> path);
std::unique_ptr<PNGMetadata> metadata_from_png_string(
    rust::Slice<const std::uint8_t> png);
std::unique_ptr<PNGMetadata> metadata_from_png_file(
    rust::Slice<const std::uint8_t> path);
std::size_t png_metadata_size(const PNGMetadata &metadata);
const std::string &png_metadata_key(const PNGMetadata &metadata, std::size_t i);
rust::Vec<std::uint8_t> png_metadata_value(const PNGMetadata &metadata,
                                           std::size_t i);

}  // namespace rdkit_cxx
