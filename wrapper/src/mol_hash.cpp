// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_hash.h"

#include <mutex>
#include <string>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

// The most atoms the ExtendedMurcko hash takes, as on rdkit-2024.03: each
// alloca stays within a 4 KiB page, with room for the frame's other locals.
constexpr unsigned kMaxExtendedMurckoAtoms = 4000;

}  // namespace

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
    case HashFunction::ExtendedMurcko:
      // RDKit 2025.09's ExtendedMurckoScaffold and DepthFirstSearchForRing
      // each take a byte per atom from the stack with alloca
      // (hashfunctions.cpp), and conda-forge's build has no stack-clash
      // protection, so a frame that grows by a page or more can step past a
      // thread stack's one-page guard into other memory (RDKit 2026.03 uses
      // vectors).
      if (mol.getNumAtoms() > kMaxExtendedMurckoAtoms) {
        throw ValueErrorException(
            "the ExtendedMurcko hash takes at most " +
            std::to_string(kMaxExtendedMurckoAtoms) +
            " atoms: RDKit 2025.09 allocates a byte per atom on the stack");
      }
      break;
    default:
      break;
  }
  // The hashes write SMILES, and the scaffold and regioisomer hashes remove
  // atoms or bonds, both of which trust what check_removals checks.
  check_removals(mol);
  // As MolHashHelper in Code/GraphMol/MolHash/Wrap/rdMolHash.cpp, which
  // hashes a copy:
  // MolHash edits the molecule it is given.
  RWMol copy(mol);
  // Some hashes match substructures (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return std::make_unique<std::string>(
      RDKit::MolHash::MolHash(&copy, func, use_cx_smiles, cx_flags_to_skip));
}

}  // namespace rdkit_cxx
