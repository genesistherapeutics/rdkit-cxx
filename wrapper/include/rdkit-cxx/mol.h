// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include <GraphMol/RWMol.h>

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

}  // namespace rdkit_cxx
