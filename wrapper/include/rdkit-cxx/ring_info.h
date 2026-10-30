/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RingInfo's methods. Each takes the molecule, which owns the RingInfo and
// bounds the indices (CONVENTIONS.md E3).

#pragma once

#include <cstdint>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

unsigned ring_info_num_rings(const RWMol &mol);
rust::Vec<IndexList> ring_info_atom_rings(const RWMol &mol);
rust::Vec<IndexList> ring_info_bond_rings(const RWMol &mol);
unsigned ring_info_num_atom_rings(const RWMol &mol, unsigned idx);
unsigned ring_info_num_bond_rings(const RWMol &mol, unsigned idx);
bool ring_info_is_atom_in_ring_of_size(const RWMol &mol, unsigned idx,
                                       unsigned size);
bool ring_info_is_bond_in_ring_of_size(const RWMol &mol, unsigned idx,
                                       unsigned size);
rust::Vec<std::int32_t> ring_info_atom_members(const RWMol &mol, unsigned idx);
rust::Vec<std::int32_t> ring_info_bond_members(const RWMol &mol, unsigned idx);
rust::Vec<std::int32_t> ring_info_atom_ring_sizes(const RWMol &mol,
                                                  unsigned idx);
rust::Vec<std::int32_t> ring_info_bond_ring_sizes(const RWMol &mol,
                                                  unsigned idx);
unsigned ring_info_min_atom_ring_size(const RWMol &mol, unsigned idx);
unsigned ring_info_min_bond_ring_size(const RWMol &mol, unsigned idx);
bool ring_info_are_atoms_in_same_ring(const RWMol &mol, unsigned idx1,
                                      unsigned idx2);
bool ring_info_are_atoms_in_same_ring_of_size(const RWMol &mol, unsigned idx1,
                                              unsigned idx2, unsigned size);
bool ring_info_are_bonds_in_same_ring(const RWMol &mol, unsigned idx1,
                                      unsigned idx2);
bool ring_info_are_bonds_in_same_ring_of_size(const RWMol &mol, unsigned idx1,
                                              unsigned idx2, unsigned size);
// The fusion queries compute and cache fused-ring data on first use.
bool ring_info_is_ring_fused(const RWMol &mol, unsigned ring_idx);
bool ring_info_are_rings_fused(const RWMol &mol, unsigned ring1_idx,
                               unsigned ring2_idx);
unsigned ring_info_num_fused_bonds(const RWMol &mol, unsigned ring_idx);
bool ring_info_are_ring_families_initialized(const RWMol &mol);
unsigned ring_info_num_ring_families(const RWMol &mol);
unsigned ring_info_num_relevant_cycles(const RWMol &mol);
rust::Vec<IndexList> ring_info_atom_ring_families(const RWMol &mol);
rust::Vec<IndexList> ring_info_bond_ring_families(const RWMol &mol);
// As addRing in Code/GraphMol/Wrap/RingInfo.cpp.
void ring_info_add_ring(RWMol &mol, rust::Slice<const std::uint32_t> atom_ids,
                        rust::Slice<const std::uint32_t> bond_ids);

}  // namespace rdkit_cxx
