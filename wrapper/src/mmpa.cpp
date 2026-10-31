// Copyright (c) 2026 Genesis Molecular AI, Inc.
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
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"

namespace rdkit_cxx {

namespace {

using Results = std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>>;

// Whether the atoms are joined by a path but for their bond.
bool joined_without(const RDKit::ROMol &mol, const RDKit::Bond &bond) {
  std::vector<char> seen(mol.getNumAtoms(), 0);
  std::vector<unsigned> stack{bond.getBeginAtomIdx()};
  seen[stack.back()] = 1;
  while (!stack.empty()) {
    const auto *atom = mol.getAtomWithIdx(stack.back());
    stack.pop_back();
    for (const auto *next : mol.atomBonds(atom)) {
      if (next == &bond) {
        continue;
      }
      auto other = next->getOtherAtomIdx(atom->getIdx());
      if (!seen[other]) {
        seen[other] = 1;
        stack.push_back(other);
      }
    }
  }
  return seen[bond.getEndAtomIdx()];
}

// fragmentMol writes each cut's side chains with MolToSmiles (addResult in
// MMPA.cpp), which perceives their stereochemistry. With one bond cut, the
// side chains are the molecule less that bond. With more, they are the
// pieces only one cut bond leaves attached: each is cut off by that bond
// alone, so its biconnected blocks are some of the molecule's, and its
// relevant cycles, which add up over the blocks, at most the molecule's. So
// is the molecule less a bond in no ring.
void check_cuts_relevant_cycles(
    const RDKit::ROMol &mol,
    const std::vector<std::pair<unsigned, unsigned>> &cuts,
    std::uint32_t min_cuts) {
  if (independent_cycles(mol) <= MAX_INDEPENDENT_CYCLES) {
    return;
  }
  check_stereo_relevant_cycles(mol);
  if (min_cuts != 1) {
    return;
  }
  for (const auto &[a, b] : cuts) {
    const auto *bond = a < mol.getNumAtoms() && b < mol.getNumAtoms()
                           ? mol.getBondBetweenAtoms(a, b)
                           : nullptr;
    if (bond && joined_without(mol, *bond)) {
      RDKit::RWMol cut(mol);
      cut.removeBond(a, b);
      check_stereo_relevant_cycles(cut);
    }
  }
}

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
  // fragmentMol cuts bonds of a copy (fragmentOnBonds) and writes the
  // pieces as SMILES, which trusts what check_removals checks.
  check_removals(mol);
  Results results;
  // fragmentMol matches the pattern (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  if (independent_cycles(mol) > MAX_INDEPENDENT_CYCLES) {
    // The bonds fragmentMol cuts, matched as it matches them.
    std::vector<RDKit::MatchVectType> matches;
    RDKit::SubstructMatch(mol, *query, matches);
    if (matches.size() <= max_cut_bonds) {
      std::vector<std::pair<unsigned, unsigned>> cuts;
      for (const auto &match : matches) {
        cuts.emplace_back(match[0].second, match[1].second);
      }
      check_cuts_relevant_cycles(mol, cuts, min_cuts);
    }
  }
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
  check_removals(mol);
  std::vector<std::pair<unsigned, unsigned>> cuts;
  for (auto idx : bonds_to_cut) {
    const auto *bond = mol.getBondWithIdx(idx);
    cuts.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
  }
  check_cuts_relevant_cycles(mol, cuts, min_cuts);
  Results results;
  RDKit::MMPA::fragmentMol(
      mol, results,
      std::vector<unsigned int>(bonds_to_cut.begin(), bonds_to_cut.end()),
      min_cuts, max_cuts);
  return to_list(results);
}

}  // namespace rdkit_cxx
