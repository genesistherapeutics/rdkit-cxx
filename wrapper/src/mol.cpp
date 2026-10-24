// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol.h"

#include <iostream>
#include <string>

#include <GraphMol/MolPickler.h>

// set_default_pickle_properties is safe Rust only because RDKit guards the
// default with a mutex in thread-safe builds (CONVENTIONS.md S2).
#ifndef RDK_BUILD_THREADSAFE_SSS
#error "rdkit-cxx needs an RDKit built with RDK_BUILD_THREADSAFE_SSS"
#endif

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

std::unique_ptr<RWMol> mol_new() { return std::make_unique<RWMol>(); }

std::unique_ptr<RWMol> mol_from_pkl_string(rust::Slice<const std::uint8_t> pkl,
                                           bool has_flags, unsigned flags) {
  std::string text(reinterpret_cast<const char *>(pkl.data()), pkl.size());
  auto mol = std::make_unique<RWMol>();
  if (has_flags) {
    RDKit::MolPickler::molFromPickle(text, mol.get(), flags);
  } else {
    RDKit::MolPickler::molFromPickle(text, mol.get());
  }
  return mol;
}

rust::Vec<std::uint8_t> mol_to_binary(const RWMol &mol, bool has_flags,
                                      unsigned flags) {
  std::string res;
  if (has_flags) {
    RDKit::MolPickler::pickleMol(mol, res, flags);
  } else {
    RDKit::MolPickler::pickleMol(mol, res);
  }
  rust::Vec<std::uint8_t> bytes;
  bytes.reserve(res.size());
  for (char c : res) {
    bytes.push_back(static_cast<std::uint8_t>(c));
  }
  return bytes;
}

// RDKit guards the default with a mutex in thread-safe builds
// (RDK_BUILD_THREADSAFE_SSS), as conda-forge's is.
unsigned get_default_pickle_properties() {
  return RDKit::MolPickler::getDefaultPickleProperties();
}

void set_default_pickle_properties(unsigned props) {
  RDKit::MolPickler::setDefaultPickleProperties(props);
}

rust::Vec<std::uint32_t> property_pickle_option_values() {
  using namespace RDKit::PicklerOps;
  rust::Vec<std::uint32_t> values;
  for (unsigned v :
       {NoProps, MolProps, AtomProps, BondProps, QueryAtomData, PrivateProps,
        ComputedProps, AllProps, CoordsAsDouble, NoConformers}) {
    values.push_back(v);
  }
  return values;
}

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

void mol_clear_property_cache(RWMol &mol) { mol.clearPropertyCache(); }

bool mol_needs_update_property_cache(const RWMol &mol) {
  return mol.needsUpdatePropertyCache();
}

bool mol_has_query(const RWMol &mol) { return mol.hasQuery(); }

std::unique_ptr<std::string> mol_get_name(const RWMol &mol) {
  return std::make_unique<std::string>(mol.getName());
}

void mol_set_name(RWMol &mol, rust::Str name) {
  mol.setName(std::string(name));
}

void mol_debug(const RWMol &mol, bool use_stdout) {
  mol.debugMol(use_stdout ? std::cout : std::cerr);
}

// As GetAromaticAtoms in Code/GraphMol/nbWrap/Mol.cpp, which matches an
// aromaticity query.
rust::Vec<std::uint32_t> mol_get_aromatic_atoms(const RWMol &mol) {
  rust::Vec<std::uint32_t> indices;
  for (const auto *atom : mol.atoms()) {
    if (atom->getIsAromatic()) {
      indices.push_back(atom->getIdx());
    }
  }
  return indices;
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

void rwmol_replace_bond(RWMol &mol, unsigned index, const Bond &new_bond,
                        bool preserve_props, bool keep_s_groups) {
  check_index(index, mol.getNumBonds());
  // replaceBond copies the bond.
  mol.replaceBond(index, const_cast<Bond *>(&new_bond), preserve_props,
                  keep_s_groups);
}

void rwmol_insert_mol(RWMol &mol, const RWMol &other) { mol.insertMol(other); }

void rwmol_begin_batch_edit(RWMol &mol) { mol.beginBatchEdit(); }

void rwmol_commit_batch_edit(RWMol &mol) { mol.commitBatchEdit(); }

void rwmol_rollback_batch_edit(RWMol &mol) { mol.rollbackBatchEdit(); }

}  // namespace rdkit_cxx
