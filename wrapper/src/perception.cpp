/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/perception.h"

#include <limits>
#include <vector>

#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/perception.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

rust::Vec<IndexList> get_sssr(RWMol &mol, bool include_dative_bonds,
                              bool include_hydrogen_bonds) {
  std::vector<std::vector<int>> rings;
  RDKit::MolOps::findSSSR(mol, rings, include_dative_bonds,
                          include_hydrogen_bonds);
  return to_index_lists(rings);
}

rust::Vec<IndexList> get_symm_sssr(RWMol &mol, bool include_dative_bonds,
                                   bool include_hydrogen_bonds) {
  std::vector<std::vector<int>> rings;
  RDKit::MolOps::symmetrizeSSSR(mol, rings, include_dative_bonds,
                                include_hydrogen_bonds);
  return to_index_lists(rings);
}

void fast_find_rings(RWMol &mol) { RDKit::MolOps::fastFindRings(mol); }

void find_ring_families(RWMol &mol, bool include_dative_bonds,
                        bool include_hydrogen_bonds) {
  RDKit::MolOps::findRingFamilies(mol, include_dative_bonds,
                                  include_hydrogen_bonds);
}

void cleanup(RWMol &mol) { RDKit::MolOps::cleanUp(mol); }

void cleanup_chirality(RWMol &mol) { RDKit::MolOps::cleanupChirality(mol); }

void cleanup_atropisomers(RWMol &mol) {
  RDKit::MolOps::cleanupAtropisomers(mol);
}

void cleanup_organometallics(RWMol &mol) {
  RDKit::MolOps::cleanUpOrganometallics(mol);
}

// RDKit's own check rejects AROMATICITY_CUSTOM, which needs a C++ function.
void set_aromaticity(RWMol &mol, AromaticityModel model) {
  RDKit::MolOps::setAromaticity(mol, model);
}

void set_conjugation(RWMol &mol) { RDKit::MolOps::setConjugation(mol); }

void set_hybridization(RWMol &mol) { RDKit::MolOps::setHybridization(mol); }

void assign_radicals(RWMol &mol) { RDKit::MolOps::assignRadicals(mol); }

void kekulize_if_possible(RWMol &mol, bool clear_aromatic_flags,
                          bool canonical) {
  RDKit::MolOps::KekulizeIfPossible(mol, clear_aromatic_flags, canonical);
}

// Both remove atoms or bonds from the copy, which trusts what
// check_removals checks.
std::unique_ptr<RWMol> haptic_bonds_to_dative(const RWMol &mol) {
  check_removals(mol);
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::hapticBondsToDative(*res);
  return res;
}

std::unique_ptr<RWMol> dative_bonds_to_haptic(const RWMol &mol) {
  check_removals(mol);
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::dativeBondsToHaptic(*res);
  return res;
}

bool needs_hs(const RWMol &mol) { return RDKit::MolOps::needsHs(mol); }

// RDKit's own checks reject atoms with no molecule.
std::int32_t count_atom_elec(const Atom &atom) {
  return RDKit::MolOps::countAtomElec(&atom);
}

bool atom_has_conjugated_bond(const Atom &atom) {
  return RDKit::MolOps::atomHasConjugatedBond(&atom);
}

// getFormalCharge sums the atoms' charges in an int accumulator, in atom
// order (std::accumulate in MolOps.cpp). Each charge is 8 bits, so the sum
// leaves int's range only past 2^24 atoms, and overflows there (U9).
std::int32_t get_formal_charge(const RWMol &mol) {
  std::int64_t sum = 0;
  for (const auto *atom : mol.atoms()) {
    sum += atom->getFormalCharge();
    if (sum < std::numeric_limits<int>::min() ||
        sum > std::numeric_limits<int>::max()) {
      throw ValueErrorException(
          "RDKit sums formal charges in int, which overflows for this "
          "molecule");
    }
  }
  return RDKit::MolOps::getFormalCharge(mol);
}

rust::Vec<std::uint32_t> get_shortest_path(const RWMol &mol, std::uint32_t aid1,
                                           std::uint32_t aid2) {
  check_index(aid1, mol.getNumAtoms());
  check_index(aid2, mol.getNumAtoms());
  rust::Vec<std::uint32_t> result;
  for (auto idx : RDKit::MolOps::getShortestPath(mol, static_cast<int>(aid1),
                                                 static_cast<int>(aid2))) {
    result.push_back(static_cast<std::uint32_t>(idx));
  }
  return result;
}

void set_terminal_atom_coords(RWMol &mol, std::uint32_t idx,
                              std::uint32_t other_idx) {
  check_index(idx, mol.getNumAtoms());
  check_index(other_idx, mol.getNumAtoms());
  RDKit::MolOps::setTerminalAtomCoords(mol, idx, other_idx);
}

}  // namespace rdkit_cxx
