// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/MolPickler.h>

// set_default_pickle_properties is safe Rust only because RDKit guards the
// default with a mutex in thread-safe builds (CONVENTIONS.md S2).
#ifndef RDK_BUILD_THREADSAFE_SSS
#error "rdkit-cxx needs an RDKit built with RDK_BUILD_THREADSAFE_SSS"
#endif

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"

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
  check_pickled_mol(*mol);
  return mol;
}

rust::Vec<std::uint8_t> mol_to_binary(const RWMol &mol, bool has_flags,
                                      unsigned flags) {
  check_picklable(mol);
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

bool mol_needs_update_property_cache(const RWMol &mol) {
  return mol.needsUpdatePropertyCache();
}

bool mol_has_query(const RWMol &mol) { return mol.hasQuery(); }

void mol_debug(const RWMol &mol, bool use_stdout) {
  mol.debugMol(use_stdout ? std::cout : std::cerr);
}

// As MolGetAromaticAtoms in Code/GraphMol/Wrap/Mol.cpp, which matches an
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

// As ReadWriteMol in Code/GraphMol/Wrap/Mol.cpp: the molecule copies the
// atom, and the bond functions return the new number of bonds. RDKit keeps
// cached matrices through these edits, so they clear them (U6).
unsigned rwmol_add_atom(RWMol &mol, const Atom &atom) {
  auto idx = mol.addAtom(atom.copy(), true, true);
  clear_matrix_caches(mol);
  return idx;
}

unsigned rwmol_add_bond(RWMol &mol, unsigned begin_atom_idx,
                        unsigned end_atom_idx, BondType order) {
  check_index(begin_atom_idx, mol.getNumAtoms());
  check_index(end_atom_idx, mol.getNumAtoms());
  auto num_bonds = mol.addBond(begin_atom_idx, end_atom_idx, order);
  clear_matrix_caches(mol);
  return num_bonds;
}

namespace {

// RWMol keeps a batch edit's removals in ROMol's private dp_delAtoms and
// dp_delBonds, which nothing public reads. An explicit instantiation may
// name private members ([temp.spec.general]), so Steal<Tag, &C::member>
// defines a friend that returns the pointer to the member, found through its
// tag (as in filter_catalog.cpp).
using DeletionsPtr = std::unique_ptr<boost::dynamic_bitset<>>;

template <class Tag, typename Tag::type Member>
struct Steal {
  friend typename Tag::type member(Tag) { return Member; }
};

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif
struct DelAtoms {
  using type = DeletionsPtr RDKit::ROMol::*;
  friend type member(DelAtoms);
};
template struct Steal<DelAtoms, &RDKit::ROMol::dp_delAtoms>;
struct DelBonds {
  using type = DeletionsPtr RDKit::ROMol::*;
  friend type member(DelBonds);
};
template struct Steal<DelBonds, &RDKit::ROMol::dp_delBonds>;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

struct BatchEdit {
  static const boost::dynamic_bitset<> *atoms(const RDKit::ROMol &mol) {
    return (mol.*member(DelAtoms())).get();
  }
  static const boost::dynamic_bitset<> *bonds(const RDKit::ROMol &mol) {
    return (mol.*member(DelBonds())).get();
  }
};

bool in_batch_edit(const RWMol &mol) {
  return BatchEdit::atoms(mol) || BatchEdit::bonds(mol);
}

// The checks of check_removals but check_stereo_group_removals, whose
// pointers the removals here take out of the groups themselves.
void check_edit_removals(const RWMol &mol) {
  check_bond_end_points(mol);
  check_substance_group_parents(mol);
}

// RDKit 2024.03's removeBond and batchRemoveBonds (RWMol.cpp) delete bonds
// and leave the stereo groups that hold them pointing at them, and
// removeAtomFromGroups erases only an atom's first copy from a group
// (check_stereo_group_removals), so this takes every copy of the atoms and
// bonds about to go out of the groups first, as 2026.03's
// removeBondFromGroups and removeAtomFromGroups do. removeAtomFromGroups
// then also drops the groups left with no atoms, which 2026.03 keeps while
// they hold bonds; their bonds stay, so they are put back afterwards.
template <typename Remove>
void remove_scrubbing_groups(RWMol &mol,
                             const std::vector<const RDKit::Atom *> &atoms,
                             const std::vector<const RDKit::Bond *> &bonds,
                             Remove remove) {
  remove_from_stereo_groups(mol, atoms, bonds);
  const auto groups = mol.getStereoGroups();
  const bool restore =
      std::any_of(groups.begin(), groups.end(),
                  [](const auto &group) { return group.getAtoms().empty(); });
  remove();
  if (restore) {
    mol.setStereoGroups(groups);
  }
}

// The atom and its bonds.
void add_atom_removal(const RWMol &mol, const RDKit::Atom *atom,
                      std::vector<const RDKit::Atom *> &atoms,
                      std::vector<const RDKit::Bond *> &bonds) {
  atoms.push_back(atom);
  for (const auto *bond : mol.atomBonds(atom)) {
    bonds.push_back(bond);
  }
}

}  // namespace

// RWMol's removals parse the bonds' _MolFileBondEndPts and follow the
// substance groups' PARENT chains unchecked (check_removals), so they are
// checked before anything is removed, and leave the stereo groups pointing
// at what they delete (check_stereo_group_removals), so what goes is taken
// out of the groups first. In a batch the removals happen at the commit,
// which checks again and takes them out of the groups then.
void rwmol_remove_atom(RWMol &mol, unsigned idx) {
  check_index(idx, mol.getNumAtoms());
  check_edit_removals(mol);
  if (in_batch_edit(mol)) {
    mol.removeAtom(idx);
    return;
  }
  std::vector<const RDKit::Atom *> atoms;
  std::vector<const RDKit::Bond *> bonds;
  add_atom_removal(mol, mol.getAtomWithIdx(idx), atoms, bonds);
  remove_scrubbing_groups(mol, atoms, bonds, [&] { mol.removeAtom(idx); });
}

void rwmol_remove_bond(RWMol &mol, unsigned idx1, unsigned idx2) {
  check_index(idx1, mol.getNumAtoms());
  check_index(idx2, mol.getNumAtoms());
  check_edit_removals(mol);
  const auto *bond = mol.getBondBetweenAtoms(idx1, idx2);
  if (bond && !in_batch_edit(mol)) {
    remove_scrubbing_groups(mol, {}, {bond},
                            [&] { mol.removeBond(idx1, idx2); });
  } else {
    mol.removeBond(idx1, idx2);
  }
  clear_matrix_caches(mol);
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
  if (!keep_s_groups) {
    // replaceBond then removes the substance groups that hold the bond as
    // a bond removal does, following their PARENT chains unchecked.
    check_edit_removals(mol);
  }
  // replaceBond copies the bond. RDKit's would leave a stereo group that
  // held the old bond pointing at it, deleted (replace_bond).
  replace_bond(mol, index, const_cast<Bond *>(&new_bond), preserve_props,
               keep_s_groups);
}

void rwmol_insert_mol(RWMol &mol, const RWMol &other) {
  mol.insertMol(other);
  clear_matrix_caches(mol);
}

void rwmol_begin_batch_edit(RWMol &mol) { mol.beginBatchEdit(); }

// A failed check leaves the batch open and the molecule as it was. The
// batch's bitsets can be shorter than the molecule, which grows during it.
void rwmol_commit_batch_edit(RWMol &mol) {
  check_edit_removals(mol);
  std::vector<const RDKit::Atom *> atoms;
  std::vector<const RDKit::Bond *> bonds;
  if (const auto *del_atoms = BatchEdit::atoms(mol)) {
    for (unsigned i = 0; i < del_atoms->size() && i < mol.getNumAtoms(); ++i) {
      if ((*del_atoms)[i]) {
        add_atom_removal(mol, mol.getAtomWithIdx(i), atoms, bonds);
      }
    }
  }
  if (const auto *del_bonds = BatchEdit::bonds(mol)) {
    for (unsigned i = 0; i < del_bonds->size() && i < mol.getNumBonds(); ++i) {
      if ((*del_bonds)[i]) {
        bonds.push_back(mol.getBondWithIdx(i));
      }
    }
  }
  remove_scrubbing_groups(mol, atoms, bonds, [&] { mol.commitBatchEdit(); });
}

void rwmol_rollback_batch_edit(RWMol &mol) { mol.rollbackBatchEdit(); }

}  // namespace rdkit_cxx
