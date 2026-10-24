/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/Atom.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

using Atom = RDKit::Atom;
using HybridizationType = RDKit::Atom::HybridizationType;
using ChiralType = RDKit::Atom::ChiralType;
using ValenceType = RDKit::Atom::ValenceType;

unsigned atom_noatom();

std::unique_ptr<Atom> atom_new(rust::Str what);
std::unique_ptr<Atom> atom_from_num(unsigned num);
std::unique_ptr<Atom> atom_clone(const Atom &atom);

// Getters that need no owning molecule.
unsigned atom_get_idx(const Atom &atom);
int atom_get_atomic_num(const Atom &atom);
std::unique_ptr<std::string> atom_get_symbol(const Atom &atom);
int atom_get_formal_charge(const Atom &atom);
bool atom_get_is_aromatic(const Atom &atom);
unsigned atom_get_isotope(const Atom &atom);
unsigned atom_get_num_explicit_hs(const Atom &atom);
HybridizationType atom_get_hybridization(const Atom &atom);
ChiralType atom_get_chiral_tag(const Atom &atom);
int atom_get_atom_map_num(const Atom &atom);
double atom_get_mass(const Atom &atom);
bool atom_get_no_implicit(const Atom &atom);
unsigned atom_get_num_radical_electrons(const Atom &atom);
bool atom_has_owning_mol(const Atom &atom);
bool atom_has_query(const Atom &atom);
bool atom_needs_update_property_cache(const Atom &atom);
unsigned atom_get_valence(const Atom &atom, ValenceType which);
// As AtomGetSmarts in Code/GraphMol/nbWrap/Atom.cpp.
std::unique_ptr<std::string> atom_get_smarts(const Atom &atom, bool do_kekule,
                                             bool all_hs_explicit,
                                             bool isomeric_smiles);

// Setters.
void atom_set_atomic_num(Atom &atom, int new_num);
void atom_set_formal_charge(Atom &atom, int what);
void atom_set_is_aromatic(Atom &atom, bool what);
void atom_set_isotope(Atom &atom, unsigned what);
void atom_set_num_explicit_hs(Atom &atom, unsigned what);
void atom_set_no_implicit(Atom &atom, bool what);
void atom_set_hybridization(Atom &atom, HybridizationType what);
void atom_set_chiral_tag(Atom &atom, ChiralType what);
void atom_set_atom_map_num(Atom &atom, int mapno, bool strict);
void atom_set_num_radical_electrons(Atom &atom, unsigned num);
void atom_clear_property_cache(Atom &atom);
bool atom_invert_chirality(Atom &atom);

// Getters that need the owning molecule (CONVENTIONS.md W4).
unsigned atom_get_degree(const Atom &atom);
unsigned atom_get_total_degree(const Atom &atom);
unsigned atom_get_total_num_hs(const Atom &atom, bool include_neighbors);
unsigned atom_get_num_implicit_hs(const Atom &atom);
unsigned atom_get_total_valence(const Atom &atom);
bool atom_is_in_ring(const Atom &atom);
bool atom_is_in_ring_size(const Atom &atom, int size);
bool atom_has_valence_violation(const Atom &atom);

// The module-level atom label functions in Code/GraphMol/nbWrap/Atom.cpp.
std::unique_ptr<std::string> atom_get_alias(const Atom &atom);
void atom_set_alias(Atom &atom, rust::Str rlabel);
int atom_get_r_label(const Atom &atom);
void atom_set_r_label(Atom &atom, int rlabel);
std::unique_ptr<std::string> atom_get_value(const Atom &atom);
void atom_set_value(Atom &atom, rust::Str rlabel);
std::unique_ptr<std::string> atom_get_supplemental_smiles_label(
    const Atom &atom);
void atom_set_supplemental_smiles_label(Atom &atom, rust::Str label);
unsigned atom_get_num_pi_electrons(const Atom &atom);
rust::Vec<std::uint32_t> atom_get_neighbors(const Atom &atom);
rust::Vec<std::uint32_t> atom_get_bonds(const Atom &atom);
void atom_update_property_cache(Atom &atom, bool strict);

}  // namespace rdkit_cxx
