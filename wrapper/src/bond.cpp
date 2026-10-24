// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#include "rdkit-cxx/bond.h"

namespace rdkit_cxx {

unsigned bond_get_idx(const Bond &bond) { return bond.getIdx(); }

BondType bond_get_bond_type(const Bond &bond) { return bond.getBondType(); }

double bond_get_bond_type_as_double(const Bond &bond) {
  return bond.getBondTypeAsDouble();
}

unsigned bond_get_begin_atom_idx(const Bond &bond) {
  return bond.getBeginAtomIdx();
}

unsigned bond_get_end_atom_idx(const Bond &bond) {
  return bond.getEndAtomIdx();
}

unsigned bond_get_other_atom_idx(const Bond &bond, unsigned this_idx) {
  return bond.getOtherAtomIdx(this_idx);
}

bool bond_get_is_aromatic(const Bond &bond) { return bond.getIsAromatic(); }

bool bond_get_is_conjugated(const Bond &bond) {
  return bond.getIsConjugated();
}

}  // namespace rdkit_cxx
