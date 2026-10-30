/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Query atoms and bonds. Rust has no separate query types: these check the
// C++ type and throw ValueErrorException for a plain atom or bond where
// Python's signature asks for a QueryAtom or QueryBond.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/bond.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

std::unique_ptr<std::string> atom_describe_query(const Atom &atom);
std::unique_ptr<std::string> atom_get_query_type(const Atom &atom);
bool atom_match(const Atom &atom, const Atom &other);
// As expandQuery and setQuery in Code/GraphMol/Wrap/Atom.cpp. `how` is a
// Queries::CompositeQueryType.
void atom_expand_query(Atom &atom, const Atom &other, int how,
                       bool maintain_order);
void atom_set_query(Atom &atom, const Atom &other);

std::unique_ptr<std::string> bond_describe_query(const Bond &bond);
bool bond_match(const Bond &bond, const Bond &other);
// As expandQuery and setQuery in Code/GraphMol/Wrap/Bond.cpp.
void bond_expand_query(Bond &bond, const Bond &other, int how,
                       bool maintain_order);
void bond_set_query(Bond &bond, const Bond &other);

// As MolGetQueryAtoms in Code/GraphMol/Wrap/Mol.cpp, which finds the
// matches of recursive queries first; sorted by index.
rust::Vec<std::uint32_t> mol_get_atoms_matching_query(const RWMol &mol,
                                                      const Atom &query);

}  // namespace rdkit_cxx
