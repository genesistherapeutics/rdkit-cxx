/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolInterchange: molecules to and from RDKit's JSON (commonchem and
// rdkitjson).
//
// RDKit's JSON parser (MolInterchange/Parser.cpp) trusts the indices, counts,
// and query trees it reads, and no fuzzer covers it, so the wrapper parses
// the JSON itself first, with the same boost::json parser, checks what RDKit
// would trust (check_json_document in mol_interchange.cpp), and hands RDKit
// the checked document (CONVENTIONS.md U2).

#pragma once

#include <memory>
#include <string>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridge (src/ffi/mol_interchange.rs).
struct JSONParseParametersData;
struct JSONWriteParametersData;

JSONParseParametersData json_parse_parameters_default();
JSONWriteParametersData json_write_parameters_default();

// As JSONDataToMols, after the checks. Malformed JSON (boost::json's
// errors and RDKit's FileParseException) is a FileParseException.
std::unique_ptr<MolList> json_to_mols(rust::Str json_block,
                                      const JSONParseParametersData &params);
std::unique_ptr<std::string> mol_to_json(const RWMol &mol,
                                         const JSONWriteParametersData &params);
std::unique_ptr<std::string> mols_to_json(
    BorrowedMols &mols, const JSONWriteParametersData &params);

}  // namespace rdkit_cxx
