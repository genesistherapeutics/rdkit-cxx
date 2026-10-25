/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdRGroupDecomposition: decomposing molecules into a core and R groups.
//
// RDKit's decomposition keeps state of its own and process-wide state:
// the GA draws from GarethUtil's RandomUtil singleton, which it also
// seeds, and fingerprint-variance scoring rewrites a function-local static
// RGroupData (dummyHydrogenFingerprint in RGroupFingerprintScore.cpp), both
// without a lock. Every call that can score or run the GA holds
// rgroup_mutex(), a module lock (U4), and RDKit's parallel GA runs, whose
// threads share both, run one after another on the calling thread. Making
// cores, adding molecules, and processing match substructures (and the
// cores' MCS and enumeration do), so those calls hold control_c_mutex too,
// taken after rgroup_mutex.

#pragma once

#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using RGroupDecomposition = RDKit::RGroupDecomposition;
using RGroupLabels = RDKit::RGroupLabels;
using RGroupMatching = RDKit::RGroupMatching;
using RGroupLabelling = RDKit::RGroupLabelling;
using RGroupCoreAlignment = RDKit::RGroupCoreAlignment;
using RGroupScore = RDKit::RGroupScore;

// Shared structs, defined by the bridge (src/ffi/rgroup.rs).
struct RGroupDecompositionParametersData;
struct RGroupTable;
struct RGroupDecomposeResult;
struct ProcessAndScoreResult;
struct MatchingCoreResult;

std::mutex &rgroup_mutex();

RGroupDecompositionParametersData r_group_decomposition_parameters_default();

// The cores are copied.
std::unique_ptr<RGroupDecomposition> r_group_decomposition_new(
    BorrowedMols &cores, const RGroupDecompositionParametersData &params);
std::int32_t r_group_decomposition_add(RGroupDecomposition &decomp,
                                       const RWMol &mol);
MatchingCoreResult r_group_decomposition_get_matching_core_idx(
    RGroupDecomposition &decomp, const RWMol &mol);
bool r_group_decomposition_process(RGroupDecomposition &decomp);
ProcessAndScoreResult r_group_decomposition_process_and_score(
    RGroupDecomposition &decomp);
std::unique_ptr<std::vector<std::string>>
r_group_decomposition_get_r_group_labels(RGroupDecomposition &decomp);
// Copies of the molecules, named by their labels, with the number in each
// row or column.
RGroupTable r_group_decomposition_get_r_groups_as_rows(
    RGroupDecomposition &decomp);
RGroupTable r_group_decomposition_get_r_groups_as_columns(
    RGroupDecomposition &decomp);

// Python's RGroupDecompose (RGroupDecomp in rdRGroupComposition.cpp).
RGroupDecomposeResult r_group_decompose(
    BorrowedMols &cores, BorrowedMols &mols,
    const RGroupDecompositionParametersData &params, bool as_rows);

void relabel_mapped_dummies(RWMol &mol, std::uint32_t input_labels,
                            std::uint32_t output_labels);

}  // namespace rdkit_cxx
