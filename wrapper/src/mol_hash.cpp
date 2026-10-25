// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_hash.h"

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

std::unique_ptr<std::string> mol_hash(const RWMol &mol, HashFunction func,
                                      bool use_cx_smiles,
                                      std::uint32_t cx_flags_to_skip) {
  // As Code/GraphMol/MolHash/nbWrap/rdMolHash.cpp, which hashes a copy:
  // MolHash edits the molecule it is given.
  RWMol copy(mol);
  return std::make_unique<std::string>(
      RDKit::MolHash::MolHash(&copy, func, use_cx_smiles, cx_flags_to_skip));
}

}  // namespace rdkit_cxx
