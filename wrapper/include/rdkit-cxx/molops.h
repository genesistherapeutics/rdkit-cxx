/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolOps.h>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A shared struct, defined by the bridge (src/ffi/molops.rs).
struct AddHsParametersData;

// RDKit's BETTER_ENUM is a plain enum with a fixed underlying type unless
// USE_BETTER_ENUMS is defined, so cxx can check its values.
using SanitizeFlags = RDKit::MolOps::SanitizeFlags;

std::uint32_t sanitize_mol(RWMol &mol, std::uint32_t sanitize_ops,
                           bool catch_errors);
void kekulize(RWMol &mol, bool clear_aromatic_flags, bool canonical);
AddHsParametersData add_hs_parameters_default();
std::unique_ptr<RWMol> add_hs(const RWMol &mol,
                              const AddHsParametersData &params);
rust::Vec<IndexList> get_mol_frags(const RWMol &mol);
std::unique_ptr<MolList> get_mol_frags_as_mols(const RWMol &mol,
                                               bool sanitize_frags);

std::size_t mol_list_size(const MolList &list);
std::unique_ptr<RWMol> mol_list_take(MolList &list, std::size_t i);

}  // namespace rdkit_cxx
