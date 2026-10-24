// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#include "rdkit-cxx/mol.h"

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

}  // namespace rdkit_cxx
