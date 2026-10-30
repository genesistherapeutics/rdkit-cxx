/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMMPA: fragmentation for matched molecular pair analysis.

#pragma once

#include <cstdint>
#include <memory>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Each result is two molecules in the list: the core (null for a single
// cut) and the side chains. Empty when the pattern matches no bond or more
// than maxCutBonds.
std::unique_ptr<MolList> fragment_mol(const RWMol &mol, std::uint32_t min_cuts,
                                      std::uint32_t max_cuts,
                                      std::uint32_t max_cut_bonds,
                                      rust::Str pattern);
std::unique_ptr<MolList> fragment_mol_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bonds_to_cut,
    std::uint32_t min_cuts, std::uint32_t max_cuts);

}  // namespace rdkit_cxx
