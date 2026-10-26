// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_hash.h"

#include <mutex>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

std::unique_ptr<std::string> mol_hash(const RWMol &mol, HashFunction func,
                                      bool use_cx_smiles,
                                      std::uint32_t cx_flags_to_skip) {
  // These hashes sum the formal charges in int (hashfunctions.cpp), and the
  // formula hash negates a negative sum. The tautomer and protomer hashes
  // also sum hydrogen counts in int and subtract the charge from them;
  // MolHash first moves each atom's hydrogens into its 8-bit explicit count
  // (NMRDKitSanitizeHydrogens), so each atom adds at most 255 (U9).
  switch (func) {
    case HashFunction::MolFormula:
    case HashFunction::Mesomer:
    case HashFunction::RedoxPair:
    case HashFunction::NetCharge:
      check_formal_charge_sum(mol);
      break;
    case HashFunction::HetAtomTautomer:
    case HashFunction::HetAtomTautomerv2:
    case HashFunction::HetAtomProtomer:
    case HashFunction::HetAtomProtomerv2:
      check_formal_charge_sum(mol, std::int64_t{255} * mol.getNumAtoms());
      break;
    default:
      break;
  }
  // As Code/GraphMol/MolHash/nbWrap/rdMolHash.cpp, which hashes a copy:
  // MolHash edits the molecule it is given.
  RWMol copy(mol);
  // Some hashes match substructures (control_c_lock.h).
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return std::make_unique<std::string>(
      RDKit::MolHash::MolHash(&copy, func, use_cx_smiles, cx_flags_to_skip));
}

}  // namespace rdkit_cxx
