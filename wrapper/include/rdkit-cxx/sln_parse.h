/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdSLNParse: Sybyl line notation.

#pragma once

#include <memory>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Null when the SLN does not parse.
std::unique_ptr<RWMol> mol_from_sln(rust::Str sln, bool sanitize,
                                    bool debug_parser);
std::unique_ptr<RWMol> mol_from_query_sln(rust::Str sln, bool merge_hs,
                                          bool debug_parser);

}  // namespace rdkit_cxx
