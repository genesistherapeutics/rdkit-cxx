/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdRascalMCES: maximum common edge substructures (RASCAL) and the
// clusterings built on them.

#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using RascalResult = RDKit::RascalMCES::RascalResult;

// Results a call returns, moved to Rust one at a time.
class RascalResultList {
 public:
  std::vector<RascalResult> results;
};

// Shared structs, defined by the bridges (src/ffi/rascal_mces.rs,
// src/ffi/chirality.rs, and src/ffi/geometry.rs).
struct RascalOptionsData;
struct RascalClusterOptionsData;
struct IndexPair;
struct IndexList;

// The options a RascalOptionsData holds, checked for what RDKit trusts about
// them (rascal_mces.cpp says what); throws ValueErrorException or
// SmilesParseException.
RDKit::RascalMCES::RascalOptions rascal_options_from(
    const RascalOptionsData &data);

// tier1Sim and tier2Sim (RascalMCES.cpp) square, in int, a count of matched
// atoms and bonds, which the atoms plus bonds of the smaller molecule bound;
// RascalResult reads each molecule's distance matrix at indices it computes
// in int (applyMaxFragSep, calcMaxDeltaAtomAtomDistScore).
inline constexpr std::int64_t kRascalMaxSquared = 46340;

// Checks molecules whose MCESs RDKit will compute, pairwise: each one's
// matrices (check_matrix_size) and what its results' removals trust
// (check_removals), and at most one with more than kRascalMaxSquared atoms
// and bonds. Throws ValueErrorException.
void check_rascal_mols(const std::vector<const RDKit::ROMol *> &mols);

RascalOptionsData rascal_options_default();
RascalClusterOptionsData rascal_cluster_options_default();

std::unique_ptr<RascalResultList> find_mces(const RWMol &mol1,
                                            const RWMol &mol2,
                                            const RascalOptionsData &opts);
std::size_t rascal_result_list_size(const RascalResultList &list);
std::unique_ptr<RascalResult> rascal_result_list_take(RascalResultList &list,
                                                      std::size_t i);

std::unique_ptr<RascalResult> rascal_result_clone(const RascalResult &result);
std::unique_ptr<std::string> rascal_result_smarts_string(
    const RascalResult &result);
rust::Vec<IndexPair> rascal_result_bond_matches(const RascalResult &result);
rust::Vec<IndexPair> rascal_result_atom_matches(const RascalResult &result);
void rascal_result_largest_fragment_only(RascalResult &result);
double rascal_result_similarity(const RascalResult &result);
int rascal_result_num_fragments(const RascalResult &result);
std::uint32_t rascal_result_largest_fragment_size(const RascalResult &result);
double rascal_result_tier1_sim(const RascalResult &result);
double rascal_result_tier2_sim(const RascalResult &result);
bool rascal_result_timed_out(const RascalResult &result);

rust::Vec<IndexList> rascal_cluster(BorrowedMols &mols,
                                    const RascalClusterOptionsData &opts);
rust::Vec<IndexList> rascal_butina_cluster(
    BorrowedMols &mols, const RascalClusterOptionsData &opts);

}  // namespace rdkit_cxx
