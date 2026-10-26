/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

using Bond = RDKit::Bond;
using BondType = RDKit::Bond::BondType;
using BondDir = RDKit::Bond::BondDir;
using BondStereo = RDKit::Bond::BondStereo;

unsigned bond_get_idx(const Bond &bond);
BondType bond_get_bond_type(const Bond &bond);
double bond_get_bond_type_as_double(const Bond &bond);
unsigned bond_get_begin_atom_idx(const Bond &bond);
unsigned bond_get_end_atom_idx(const Bond &bond);
unsigned bond_get_other_atom_idx(const Bond &bond, unsigned this_idx);
bool bond_get_is_aromatic(const Bond &bond);
bool bond_get_is_conjugated(const Bond &bond);
BondDir bond_get_bond_dir(const Bond &bond);
BondStereo bond_get_stereo(const Bond &bond);
rust::Vec<std::int32_t> bond_get_stereo_atoms(const Bond &bond);
bool bond_has_owning_mol(const Bond &bond);
bool bond_has_query(const Bond &bond);
bool bond_is_in_ring(const Bond &bond);
bool bond_is_in_ring_size(const Bond &bond, int size);
double bond_get_valence_contrib(const Bond &bond, const RDKit::Atom &at);
// As BondGetSmarts in Code/GraphMol/Wrap/Bond.cpp.
std::unique_ptr<std::string> bond_get_smarts(const Bond &bond,
                                             bool all_bonds_explicit);

void bond_set_bond_type(Bond &bond, BondType b_t);
void bond_set_is_aromatic(Bond &bond, bool what);
void bond_set_is_conjugated(Bond &bond, bool what);
void bond_set_bond_dir(Bond &bond, BondDir what);
void bond_set_stereo(Bond &bond, BondStereo what);
void bond_set_stereo_atoms(Bond &bond, unsigned bgn_idx, unsigned end_idx);
bool bond_invert_chirality(Bond &bond);
std::unique_ptr<Bond> bond_clone(const Bond &bond);

}  // namespace rdkit_cxx
