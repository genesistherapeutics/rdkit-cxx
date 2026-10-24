/*
 * Copyright (c) 2026 Steven Kearnes. All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridge (src/ffi/molfiles.rs).
struct SmilesParserParamsData;
struct SmartsParserParamsData;
struct SmilesWriteParamsData;
struct MolWriterParamsData;

// Each returns the C++ class's defaults (CONVENTIONS.md A3).
SmilesParserParamsData smiles_parser_params_default();
SmartsParserParamsData smarts_parser_params_default();
SmilesWriteParamsData smiles_write_params_default();
MolWriterParamsData mol_writer_params_default();

std::unique_ptr<RWMol> mol_from_smiles(rust::Str smiles,
                                       const SmilesParserParamsData &params);
std::unique_ptr<std::string> mol_to_smiles(const RWMol &mol,
                                           const SmilesWriteParamsData &params);
std::unique_ptr<RWMol> mol_from_smarts(rust::Str smarts,
                                       const SmartsParserParamsData &params);
std::unique_ptr<std::string> mol_to_smarts(const RWMol &mol,
                                           const SmilesWriteParamsData &params);

std::unique_ptr<RWMol> mol_from_mol_block(rust::Str mol_block, bool sanitize,
                                          bool remove_hs, bool strict_parsing);
std::unique_ptr<RWMol> mol_from_mol_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs,
                                         bool strict_parsing);
std::unique_ptr<std::string> mol_to_mol_block(
    const RWMol &mol, const MolWriterParamsData &params);
void mol_to_mol_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     const MolWriterParamsData &params);

}  // namespace rdkit_cxx
