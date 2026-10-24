/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolops's query adjustment: AdjustQueryProperties and the generic-group
// queries.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolOps.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A C enum, which cxx checks against the bridge's values.
using AdjustQueryWhichFlags = RDKit::MolOps::AdjustQueryWhichFlags;

// A shared struct, defined by the bridge (src/ffi/adjust_query.rs).
struct AdjustQueryParametersData;

AdjustQueryParametersData adjust_query_parameters_default();
AdjustQueryParametersData adjust_query_parameters_no_adjustments();
std::unique_ptr<RWMol> adjust_query_properties(
    const RWMol &mol, const AdjustQueryParametersData &params);
std::unique_ptr<RWMol> adjust_query_properties_with_generic_groups(
    const RWMol &mol, const AdjustQueryParametersData &params);
void set_generic_queries_from_properties(RWMol &mol, bool use_atom_labels,
                                         bool use_s_groups);
void convert_generic_queries_to_substance_groups(RWMol &mol);

}  // namespace rdkit_cxx
