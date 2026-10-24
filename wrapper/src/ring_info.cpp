// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/ring_info.h"

#include <GraphMol/RingInfo.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

rust::Vec<IndexList> index_lists(const RDKit::RingInfo::VECT_INT_VECT &rings) {
  rust::Vec<IndexList> result;
  for (const auto &ring : rings) {
    IndexList list;
    for (auto idx : ring) {
      list.indices.push_back(static_cast<std::uint32_t>(idx));
    }
    result.push_back(std::move(list));
  }
  return result;
}

}  // namespace

unsigned ring_info_num_rings(const RWMol &mol) {
  return mol.getRingInfo()->numRings();
}

rust::Vec<IndexList> ring_info_atom_rings(const RWMol &mol) {
  return index_lists(mol.getRingInfo()->atomRings());
}

rust::Vec<IndexList> ring_info_bond_rings(const RWMol &mol) {
  return index_lists(mol.getRingInfo()->bondRings());
}

unsigned ring_info_num_atom_rings(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  return mol.getRingInfo()->numAtomRings(idx);
}

unsigned ring_info_num_bond_rings(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumBonds());
  return mol.getRingInfo()->numBondRings(idx);
}

bool ring_info_is_atom_in_ring_of_size(const RWMol &mol, unsigned idx,
                                       unsigned size) {
  check_index(idx, mol.getNumAtoms());
  return mol.getRingInfo()->isAtomInRingOfSize(idx, size);
}

bool ring_info_is_bond_in_ring_of_size(const RWMol &mol, unsigned idx,
                                       unsigned size) {
  check_index(idx, mol.getNumBonds());
  return mol.getRingInfo()->isBondInRingOfSize(idx, size);
}

}  // namespace rdkit_cxx
