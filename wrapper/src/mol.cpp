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

Atom &mol_get_atom_with_idx_mut(RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  return *mol.getAtomWithIdx(idx);
}

Bond &mol_get_bond_with_idx_mut(RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumBonds());
  return *mol.getBondWithIdx(idx);
}

void mol_update_property_cache(RWMol &mol, bool strict) {
  mol.updatePropertyCache(strict);
}

// As ReadWriteMol in Code/GraphMol/nbWrap/Mol.cpp: the molecule copies the
// atom, and the bond functions return the new number of bonds.
unsigned rwmol_add_atom(RWMol &mol, const Atom &atom) {
  return mol.addAtom(atom.copy(), true, true);
}

unsigned rwmol_add_bond(RWMol &mol, unsigned begin_atom_idx,
                        unsigned end_atom_idx, BondType order) {
  check_index(begin_atom_idx, mol.getNumAtoms());
  check_index(end_atom_idx, mol.getNumAtoms());
  return mol.addBond(begin_atom_idx, end_atom_idx, order);
}

void rwmol_remove_atom(RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  mol.removeAtom(idx);
}

void rwmol_remove_bond(RWMol &mol, unsigned idx1, unsigned idx2) {
  check_index(idx1, mol.getNumAtoms());
  check_index(idx2, mol.getNumAtoms());
  mol.removeBond(idx1, idx2);
}

void rwmol_replace_atom(RWMol &mol, unsigned index, const Atom &new_atom,
                        bool update_label, bool preserve_props) {
  check_index(index, mol.getNumAtoms());
  // replaceAtom copies the atom.
  mol.replaceAtom(index, const_cast<Atom *>(&new_atom), update_label,
                  preserve_props);
}

}  // namespace rdkit_cxx
