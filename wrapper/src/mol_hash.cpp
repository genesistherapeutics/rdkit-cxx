// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_hash.h"

#include <algorithm>
#include <mutex>
#include <string>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"

namespace rdkit_cxx {

namespace {

// The hashes that write SMILES first assign stereochemistry, with cleanIt
// and force, which symmetrizes what the hash leaves (hashfunctions.cpp). The
// scaffold and regioisomer hashes remove only atoms of at most one bond and
// bonds in no ring, which leaves the biconnected blocks, and so the relevant
// cycles, as they were. The graph and mesomer hashes make every bond single.
// The tautomer hashes make some bonds single or aromatic: a molecule above
// the bound with bonds RDKit leaves out of its ring graph (zero-order,
// dative, hydrogen), which they could bring in, is refused for them.
void check_hash_relevant_cycles(const RWMol &mol, HashFunction func) {
  switch (func) {
    case HashFunction::CanonicalSmiles:
    case HashFunction::MurckoScaffold:
    case HashFunction::ExtendedMurcko:
    case HashFunction::Regioisomer:
      check_graph_relevant_cycles(mol);
      return;
    case HashFunction::AnonymousGraph:
    case HashFunction::ElementGraph:
    case HashFunction::Mesomer:
    case HashFunction::RedoxPair: {
      if (independent_cycles(mol) <= MAX_INDEPENDENT_CYCLES) {
        return;
      }
      RWMol single(mol);
      for (auto *bond : single.bonds()) {
        bond->setBondType(RDKit::Bond::SINGLE);
      }
      check_graph_relevant_cycles(single);
      return;
    }
    case HashFunction::HetAtomTautomer:
    case HashFunction::HetAtomTautomerv2:
    case HashFunction::HetAtomProtomer:
    case HashFunction::HetAtomProtomerv2: {
      bool outside = std::any_of(
          mol.bonds().begin(), mol.bonds().end(), [](const auto *bond) {
            auto type = bond->getBondType();
            return type == RDKit::Bond::ZERO || RDKit::isDative(type) ||
                   type == RDKit::Bond::HYDROGEN;
          });
      if (!outside) {
        check_graph_relevant_cycles(mol);
      } else if (independent_cycles(mol) > MAX_INDEPENDENT_CYCLES) {
        throw ValueErrorException(
            "the tautomer hashes of a molecule with more than " +
            std::to_string(MAX_INDEPENDENT_CYCLES) +
            " independent cycles and zero-order, dative or hydrogen bonds "
            "cannot be checked for their relevant cycles");
      }
      return;
    }
    default:
      return;
  }
}

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
    default:
      break;
  }
  // The hashes write SMILES, and the scaffold and regioisomer hashes remove
  // atoms or bonds, both of which trust what check_removals checks.
  check_removals(mol);
  check_hash_relevant_cycles(mol, func);
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
