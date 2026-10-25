/*
 * Copyright (c) 2026 Steven Kearnes
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
