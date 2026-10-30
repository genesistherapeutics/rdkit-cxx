/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for GraphMol/DetermineBonds (rdkit.Chem.rdDetermineBonds).

#pragma once

#include <cstddef>
#include <cstdint>

#include "rdkit-cxx/mol.h"

namespace rdkit_cxx {

void determine_connectivity(RWMol &mol, bool use_hueckel, std::int32_t charge,
                            double cov_factor, bool use_vdw);
void determine_bond_orders(RWMol &mol, std::int32_t charge,
                           bool allow_charged_fragments, bool embed_chiral,
                           bool use_atom_map, std::size_t max_iterations);
void determine_bonds(RWMol &mol, bool use_hueckel, std::int32_t charge,
                     double cov_factor, bool allow_charged_fragments,
                     bool embed_chiral, bool use_atom_map, bool use_vdw,
                     std::size_t max_iterations);
bool hueckel_enabled();

}  // namespace rdkit_cxx
