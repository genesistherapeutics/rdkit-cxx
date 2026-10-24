/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <GraphMol/Bond.h>

namespace rdkit_cxx {

using Bond = RDKit::Bond;
using BondType = RDKit::Bond::BondType;

unsigned bond_get_idx(const Bond &bond);
BondType bond_get_bond_type(const Bond &bond);
double bond_get_bond_type_as_double(const Bond &bond);
unsigned bond_get_begin_atom_idx(const Bond &bond);
unsigned bond_get_end_atom_idx(const Bond &bond);
unsigned bond_get_other_atom_idx(const Bond &bond, unsigned this_idx);
bool bond_get_is_aromatic(const Bond &bond);
bool bond_get_is_conjugated(const Bond &bond);

void bond_set_bond_type(Bond &bond, BondType b_t);
void bond_set_is_aromatic(Bond &bond, bool what);
void bond_set_is_conjugated(Bond &bond, bool what);

}  // namespace rdkit_cxx
