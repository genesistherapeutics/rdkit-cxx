// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/ring_info.h"

#include <GraphMol/RingInfo.h>
#include <RDGeneral/Exceptions.h>

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

rust::Vec<std::int32_t> ints(const std::vector<int> &values) {
  rust::Vec<std::int32_t> result;
  for (auto v : values) {
    result.push_back(v);
  }
  return result;
}

RDKit::RingInfo &ring_info(const RWMol &mol) {
  return *const_cast<RWMol &>(mol).getRingInfo();
}

void check_ring(const RWMol &mol, unsigned ring_idx) {
  check_index(ring_idx, mol.getRingInfo()->numRings());
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

rust::Vec<std::int32_t> ring_info_atom_members(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  return ints(mol.getRingInfo()->atomMembers(idx));
}

rust::Vec<std::int32_t> ring_info_bond_members(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumBonds());
  return ints(mol.getRingInfo()->bondMembers(idx));
}

rust::Vec<std::int32_t> ring_info_atom_ring_sizes(const RWMol &mol,
                                                  unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  return ints(mol.getRingInfo()->atomRingSizes(idx));
}

rust::Vec<std::int32_t> ring_info_bond_ring_sizes(const RWMol &mol,
                                                  unsigned idx) {
  check_index(idx, mol.getNumBonds());
  return ints(mol.getRingInfo()->bondRingSizes(idx));
}

unsigned ring_info_min_atom_ring_size(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  return mol.getRingInfo()->minAtomRingSize(idx);
}

unsigned ring_info_min_bond_ring_size(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumBonds());
  return mol.getRingInfo()->minBondRingSize(idx);
}

bool ring_info_are_atoms_in_same_ring(const RWMol &mol, unsigned idx1,
                                      unsigned idx2) {
  check_index(idx1, mol.getNumAtoms());
  check_index(idx2, mol.getNumAtoms());
  return mol.getRingInfo()->areAtomsInSameRing(idx1, idx2);
}

bool ring_info_are_atoms_in_same_ring_of_size(const RWMol &mol, unsigned idx1,
                                              unsigned idx2, unsigned size) {
  check_index(idx1, mol.getNumAtoms());
  check_index(idx2, mol.getNumAtoms());
  return mol.getRingInfo()->areAtomsInSameRingOfSize(idx1, idx2, size);
}

bool ring_info_are_bonds_in_same_ring(const RWMol &mol, unsigned idx1,
                                      unsigned idx2) {
  check_index(idx1, mol.getNumBonds());
  check_index(idx2, mol.getNumBonds());
  return mol.getRingInfo()->areBondsInSameRing(idx1, idx2);
}

bool ring_info_are_bonds_in_same_ring_of_size(const RWMol &mol, unsigned idx1,
                                              unsigned idx2, unsigned size) {
  check_index(idx1, mol.getNumBonds());
  check_index(idx2, mol.getNumBonds());
  return mol.getRingInfo()->areBondsInSameRingOfSize(idx1, idx2, size);
}

bool ring_info_is_ring_fused(const RWMol &mol, unsigned ring_idx) {
  check_ring(mol, ring_idx);
  return ring_info(mol).isRingFused(ring_idx);
}

bool ring_info_are_rings_fused(const RWMol &mol, unsigned ring1_idx,
                               unsigned ring2_idx) {
  check_ring(mol, ring1_idx);
  check_ring(mol, ring2_idx);
  return ring_info(mol).areRingsFused(ring1_idx, ring2_idx);
}

unsigned ring_info_num_fused_bonds(const RWMol &mol, unsigned ring_idx) {
  check_ring(mol, ring_idx);
  return ring_info(mol).numFusedBonds(ring_idx);
}

bool ring_info_are_ring_families_initialized(const RWMol &mol) {
  return mol.getRingInfo()->areRingFamiliesInitialized();
}

unsigned ring_info_num_ring_families(const RWMol &mol) {
  return mol.getRingInfo()->numRingFamilies();
}

unsigned ring_info_num_relevant_cycles(const RWMol &mol) {
  return mol.getRingInfo()->numRelevantCycles();
}

rust::Vec<IndexList> ring_info_atom_ring_families(const RWMol &mol) {
  return index_lists(mol.getRingInfo()->atomRingFamilies());
}

rust::Vec<IndexList> ring_info_bond_ring_families(const RWMol &mol) {
  return index_lists(mol.getRingInfo()->bondRingFamilies());
}

void ring_info_add_ring(RWMol &mol, rust::Slice<const std::uint32_t> atom_ids,
                        rust::Slice<const std::uint32_t> bond_ids) {
  if (atom_ids.size() != bond_ids.size()) {
    throw ValueErrorException("list sizes must match");
  }
  RDKit::INT_VECT atoms;
  for (auto idx : atom_ids) {
    check_index(idx, mol.getNumAtoms());
    atoms.push_back(idx);
  }
  RDKit::INT_VECT bonds;
  for (auto idx : bond_ids) {
    check_index(idx, mol.getNumBonds());
    bonds.push_back(idx);
  }
  auto *info = mol.getRingInfo();
  if (!info->isInitialized()) {
    info->initialize();
  }
  info->addRing(atoms, bonds);
}

}  // namespace rdkit_cxx
