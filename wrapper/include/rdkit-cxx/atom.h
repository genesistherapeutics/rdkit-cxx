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

// Getters that need the owning molecule (CONVENTIONS.md W4).
unsigned atom_get_degree(const Atom &atom);
unsigned atom_get_total_degree(const Atom &atom);
unsigned atom_get_total_num_hs(const Atom &atom, bool include_neighbors);
unsigned atom_get_num_implicit_hs(const Atom &atom);
unsigned atom_get_total_valence(const Atom &atom);
bool atom_is_in_ring(const Atom &atom);
rust::Vec<std::uint32_t> atom_get_neighbors(const Atom &atom);
rust::Vec<std::uint32_t> atom_get_bonds(const Atom &atom);

}  // namespace rdkit_cxx
