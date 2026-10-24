// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/RWMol.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/bond.h"

namespace rdkit_cxx {

// Every Rust molecule owns an RWMol (CONVENTIONS.md W1).
using RWMol = RDKit::RWMol;

std::unique_ptr<RWMol> mol_new();
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

// RWMol's editing methods.
unsigned rwmol_add_atom(RWMol &mol, const Atom &atom);
unsigned rwmol_add_bond(RWMol &mol, unsigned begin_atom_idx,
                        unsigned end_atom_idx, BondType order);
void rwmol_remove_atom(RWMol &mol, unsigned idx);
void rwmol_remove_bond(RWMol &mol, unsigned idx1, unsigned idx2);
void rwmol_replace_atom(RWMol &mol, unsigned index, const Atom &new_atom,
                        bool update_label, bool preserve_props);

}  // namespace rdkit_cxx
