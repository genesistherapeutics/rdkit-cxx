/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/RWMol.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/bond.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Every Rust molecule owns an RWMol (CONVENTIONS.md W1).
using RWMol = RDKit::RWMol;

std::unique_ptr<RWMol> mol_new();
// Pickles, as MolToBinary and the pickle constructors in Wrap/Mol.cpp.
std::unique_ptr<RWMol> mol_from_pkl_string(rust::Slice<const std::uint8_t> pkl,
                                           bool has_flags, unsigned flags);
rust::Vec<std::uint8_t> mol_to_binary(const RWMol &mol, bool has_flags,
                                      unsigned flags);
unsigned get_default_pickle_properties();
void set_default_pickle_properties(unsigned props);
// PropertyPickleOptions' values in a fixed order, for a test that pins them.
rust::Vec<std::uint32_t> property_pickle_option_values();
std::unique_ptr<RWMol> mol_clone(const RWMol &mol);
unsigned mol_get_num_atoms(const RWMol &mol);
unsigned mol_get_num_atoms_with_only_explicit(const RWMol &mol,
                                              bool only_explicit);
unsigned mol_get_num_heavy_atoms(const RWMol &mol);
unsigned mol_get_num_bonds(const RWMol &mol, bool only_heavy);
// The accessors check their indices (E3).
const Atom &mol_get_atom_with_idx(const RWMol &mol, unsigned idx);
const Bond &mol_get_bond_with_idx(const RWMol &mol, unsigned idx);
// The bond's index, or -1 when the atoms are not bonded.
std::int64_t mol_get_bond_between_atoms(const RWMol &mol, unsigned idx1,
                                        unsigned idx2);
Atom &mol_get_atom_with_idx_mut(RWMol &mol, unsigned idx);
Bond &mol_get_bond_with_idx_mut(RWMol &mol, unsigned idx);
void mol_update_property_cache(RWMol &mol, bool strict);
void mol_clear_property_cache(RWMol &mol);
bool mol_needs_update_property_cache(const RWMol &mol);
bool mol_has_query(const RWMol &mol);
void mol_debug(const RWMol &mol, bool use_stdout);
rust::Vec<std::uint32_t> mol_get_aromatic_atoms(const RWMol &mol);

// RWMol's editing methods.
unsigned rwmol_add_atom(RWMol &mol, const Atom &atom);
unsigned rwmol_add_bond(RWMol &mol, unsigned begin_atom_idx,
                        unsigned end_atom_idx, BondType order);
void rwmol_remove_atom(RWMol &mol, unsigned idx);
void rwmol_remove_bond(RWMol &mol, unsigned idx1, unsigned idx2);
void rwmol_replace_atom(RWMol &mol, unsigned index, const Atom &new_atom,
                        bool update_label, bool preserve_props);
void rwmol_replace_bond(RWMol &mol, unsigned index, const Bond &new_bond,
                        bool preserve_props, bool keep_s_groups);
void rwmol_insert_mol(RWMol &mol, const RWMol &other);
void rwmol_begin_batch_edit(RWMol &mol);
void rwmol_commit_batch_edit(RWMol &mol);
void rwmol_rollback_batch_edit(RWMol &mol);

}  // namespace rdkit_cxx
