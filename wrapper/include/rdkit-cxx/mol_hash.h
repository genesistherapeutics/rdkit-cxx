/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolHash: NextMove's molecular hashes.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/MolHash/MolHash.h>

#include "rdkit-cxx/mol.h"

namespace rdkit_cxx {

using HashFunction = RDKit::MolHash::HashFunction;

std::unique_ptr<std::string> mol_hash(const RWMol &mol, HashFunction func,
                                      bool use_cx_smiles,
                                      std::uint32_t cx_flags_to_skip);

}  // namespace rdkit_cxx
