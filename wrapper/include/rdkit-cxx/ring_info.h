// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause
//
// RingInfo's methods. Each takes the molecule, which owns the RingInfo and
// bounds the indices (CONVENTIONS.md E3).

#pragma once

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

}  // namespace rdkit_cxx
