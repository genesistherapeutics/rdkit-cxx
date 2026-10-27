/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolops's ring finding, aromaticity, cleanup, and other perception steps,
// most of them the steps sanitization runs.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolOps.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A C enum, which cxx checks against the bridge's values.
using AromaticityModel = RDKit::MolOps::AromaticityModel;

rust::Vec<IndexList> get_sssr(RWMol &mol, bool include_dative_bonds);
rust::Vec<IndexList> get_symm_sssr(RWMol &mol, bool include_dative_bonds);
void fast_find_rings(RWMol &mol);
void find_ring_families(RWMol &mol);

void cleanup(RWMol &mol);
void cleanup_chirality(RWMol &mol);
void cleanup_atropisomers(RWMol &mol);
void cleanup_organometallics(RWMol &mol);
void set_aromaticity(RWMol &mol, AromaticityModel model);
void set_conjugation(RWMol &mol);
void set_hybridization(RWMol &mol);
void assign_radicals(RWMol &mol);
void kekulize_if_possible(RWMol &mol, bool clear_aromatic_flags);
std::unique_ptr<RWMol> haptic_bonds_to_dative(const RWMol &mol);
std::unique_ptr<RWMol> dative_bonds_to_haptic(const RWMol &mol);
bool needs_hs(const RWMol &mol);
std::int32_t count_atom_elec(const Atom &atom);
bool atom_has_conjugated_bond(const Atom &atom);
std::int32_t get_formal_charge(const RWMol &mol);
rust::Vec<std::uint32_t> get_shortest_path(const RWMol &mol, std::uint32_t aid1,
                                           std::uint32_t aid2);
void set_terminal_atom_coords(RWMol &mol, std::uint32_t idx,
                              std::uint32_t other_idx);

}  // namespace rdkit_cxx
