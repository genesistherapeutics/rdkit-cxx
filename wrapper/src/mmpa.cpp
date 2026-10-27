// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mmpa.h"

#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/MMPA/MMPA.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

using Results = std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>>;

std::unique_ptr<RWMol> copy(const RDKit::ROMOL_SPTR &mol) {
  if (!mol) {
    return nullptr;
  }
  auto result = std::make_unique<RWMol>(*mol);
  clear_matrix_caches(*result);
  return result;
}

std::unique_ptr<MolList> to_list(const Results &results) {
  auto list = std::make_unique<MolList>();
  for (const auto &[core, side_chains] : results) {
    list->push(copy(core));
    list->push(copy(side_chains));
  }
  return list;
}

}  // namespace

std::unique_ptr<MolList> fragment_mol(const RWMol &mol, std::uint32_t min_cuts,
                                      std::uint32_t max_cuts,
                                      std::uint32_t max_cut_bonds,
                                      rust::Str pattern) {
  // fragmentMol matches the pattern unchecked: SmartsToMol's null for bad
  // SMARTS is dereferenced, and each match's first two atoms are read as
  // the bond to cut.
  std::string smarts(pattern);
  std::unique_ptr<RDKit::RWMol> query(RDKit::SmartsToMol(smarts));
  if (!query) {
    throw RDKit::SmilesParse::SmilesParseException(
        "could not parse the pattern SMARTS '" + smarts + "'");
  }
  if (query->getNumAtoms() < 2) {
    throw ValueErrorException("the pattern must match a bond: two atoms");
  }
  // fragmentMol cuts bonds of a copy (RWMol::removeBond in addResult,
  // MMPA.cpp) and writes the pieces as SMILES, which trusts what
  // check_removals checks.
  //
  // RDKit 2024.03's computeMorganCodeHash (MMPA.cpp), which addResult runs
  // to deduplicate results, left-shifts each atom's formal charge by 16
  // (getFormalCharge() << 16). For an anion that shifts a negative int,
  // which C++17 (RDKit 2024.03's standard) leaves undefined; GCC and Clang
  // define it as the two's-complement shift (C++20 does too, and 2026.03
  // rewrote it with std::abs), so anions are not refused.
  check_removals(mol);
  Results results;
  // fragmentMol matches the pattern (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  if (!RDKit::MMPA::fragmentMol(mol, results, min_cuts, max_cuts, max_cut_bonds,
                                smarts)) {
    results.clear();
  }
  return to_list(results);
}

std::unique_ptr<MolList> fragment_mol_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bonds_to_cut,
    std::uint32_t min_cuts, std::uint32_t max_cuts) {
  // As fragmentMolHelper3 in Code/GraphMol/MMPA/Wrap/rdMMPA.cpp.
  if (bonds_to_cut.empty()) {
    throw ValueErrorException("bondsToCut must be non-empty");
  }
  for (auto idx : bonds_to_cut) {
    check_index(idx, mol.getNumBonds());
  }
  // As fragment_mol: removals, and computeMorganCodeHash's shift of a
  // negative charge.
  check_removals(mol);
  Results results;
  RDKit::MMPA::fragmentMol(
      mol, results,
      std::vector<unsigned int>(bonds_to_cut.begin(), bonds_to_cut.end()),
      min_cuts, max_cuts);
  return to_list(results);
}

}  // namespace rdkit_cxx
