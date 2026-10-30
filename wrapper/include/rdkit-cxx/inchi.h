/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdinchi: InChI strings and keys, through the InChI library RDKit links.

#pragma once

#include <memory>
#include <string>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct InchiToMolData;
struct MolToInchiData;

// As InchiToMol. The molecule is null when InChI cannot read the string.
InchiToMolData inchi_to_mol(rust::Str inchi, bool sanitize, bool remove_hs);
// As MolToInchi, after the checks in inchi.cpp; an empty options string is
// none, as Python passes it.
MolToInchiData mol_to_inchi(const RWMol &mol, rust::Str options);
// As MolToInchiKey, with mol_to_inchi's checks.
std::unique_ptr<std::string> mol_to_inchi_key(const RWMol &mol,
                                              rust::Str options);
std::unique_ptr<std::string> inchi_to_inchi_key(rust::Str inchi);
std::unique_ptr<std::string> get_inchi_version();

}  // namespace rdkit_cxx
