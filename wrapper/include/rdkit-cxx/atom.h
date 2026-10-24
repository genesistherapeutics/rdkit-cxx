// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

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

// Setters.
void atom_set_atomic_num(Atom &atom, int new_num);
void atom_set_formal_charge(Atom &atom, int what);
void atom_set_is_aromatic(Atom &atom, bool what);
void atom_set_isotope(Atom &atom, unsigned what);
void atom_set_num_explicit_hs(Atom &atom, unsigned what);
void atom_set_no_implicit(Atom &atom, bool what);
void atom_set_hybridization(Atom &atom, HybridizationType what);
void atom_set_chiral_tag(Atom &atom, ChiralType what);

// Getters that need the owning molecule (CONVENTIONS.md W4).
unsigned atom_get_degree(const Atom &atom);
unsigned atom_get_total_degree(const Atom &atom);
unsigned atom_get_total_num_hs(const Atom &atom, bool include_neighbors);
unsigned atom_get_num_implicit_hs(const Atom &atom);
unsigned atom_get_total_valence(const Atom &atom);
bool atom_is_in_ring(const Atom &atom);
rust::Vec<std::uint32_t> atom_get_neighbors(const Atom &atom);
rust::Vec<std::uint32_t> atom_get_bonds(const Atom &atom);
void atom_update_property_cache(Atom &atom, bool strict);

}  // namespace rdkit_cxx
