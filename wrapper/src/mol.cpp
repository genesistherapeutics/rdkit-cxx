// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#include "rdkit-cxx/mol.h"

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

std::unique_ptr<RWMol> mol_new() { return std::make_unique<RWMol>(); }

std::unique_ptr<RWMol> mol_clone(const RWMol &mol) {
  return std::make_unique<RWMol>(mol);
}

unsigned mol_get_num_atoms(const RWMol &mol) { return mol.getNumAtoms(); }

unsigned mol_get_num_atoms_with_only_explicit(const RWMol &mol,
                                              bool only_explicit) {
  return mol.getNumAtoms(only_explicit);
}

unsigned mol_get_num_heavy_atoms(const RWMol &mol) {
  return mol.getNumHeavyAtoms();
}

unsigned mol_get_num_bonds(const RWMol &mol, bool only_heavy) {
  return mol.getNumBonds(only_heavy);
}

const Atom &mol_get_atom_with_idx(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  return *mol.getAtomWithIdx(idx);
}

const Bond &mol_get_bond_with_idx(const RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumBonds());
  return *mol.getBondWithIdx(idx);
}

std::int64_t mol_get_bond_between_atoms(const RWMol &mol, unsigned idx1,
                                        unsigned idx2) {
  check_index(idx1, mol.getNumAtoms());
  check_index(idx2, mol.getNumAtoms());
  const auto *bond = mol.getBondBetweenAtoms(idx1, idx2);
  return bond ? static_cast<std::int64_t>(bond->getIdx()) : -1;
}

}  // namespace rdkit_cxx
