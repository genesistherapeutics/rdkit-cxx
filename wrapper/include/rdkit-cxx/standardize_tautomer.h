/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's tautomer enumeration and scoring (Tautomer.h), and
// callbacks and scoring functions that call Rust closures.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/MolStandardize/Tautomer.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Defined by the bridges (src/ffi/standardize.rs and
// src/ffi/standardize_tautomer.rs).
struct CleanupParametersData;
struct SubstructTermData;
struct TautomerCallback;
struct ScoreFunction;

using TautomerEnumerator = RDKit::MolStandardize::TautomerEnumerator;
using TautomerEnumeratorResult =
    RDKit::MolStandardize::TautomerEnumeratorResult;
using TautomerEnumeratorStatus =
    RDKit::MolStandardize::TautomerEnumeratorStatus;

// A tautomer enumerator callback, shared with the enumerators it is set on.
struct TautomerCallbackHandle {
  std::shared_ptr<RDKit::MolStandardize::TautomerEnumeratorCallback> ptr;
};

std::unique_ptr<TautomerEnumerator> tautomer_enumerator_default();
// Checks the transforms (standardize.h's checked_enumerator).
std::unique_ptr<TautomerEnumerator> tautomer_enumerator_new(
    const CleanupParametersData &params);
std::unique_ptr<TautomerEnumerator> tautomer_enumerator_v1();
// Copies share the catalog and the callback.
std::unique_ptr<TautomerEnumerator> tautomer_enumerator_clone(
    const TautomerEnumerator &enumerator);

std::uint32_t tautomer_enumerator_get_max_tautomers(
    const TautomerEnumerator &enumerator);
void tautomer_enumerator_set_max_tautomers(TautomerEnumerator &enumerator,
                                           std::uint32_t max_tautomers);
std::uint32_t tautomer_enumerator_get_max_transforms(
    const TautomerEnumerator &enumerator);
void tautomer_enumerator_set_max_transforms(TautomerEnumerator &enumerator,
                                            std::uint32_t max_transforms);
bool tautomer_enumerator_get_remove_sp3_stereo(
    const TautomerEnumerator &enumerator);
void tautomer_enumerator_set_remove_sp3_stereo(TautomerEnumerator &enumerator,
                                               bool value);
bool tautomer_enumerator_get_remove_bond_stereo(
    const TautomerEnumerator &enumerator);
void tautomer_enumerator_set_remove_bond_stereo(TautomerEnumerator &enumerator,
                                                bool value);
bool tautomer_enumerator_get_reassign_stereo(
    const TautomerEnumerator &enumerator);
void tautomer_enumerator_set_reassign_stereo(TautomerEnumerator &enumerator,
                                             bool value);
void tautomer_enumerator_set_callback(TautomerEnumerator &enumerator,
                                      const TautomerCallbackHandle &callback);
void tautomer_enumerator_clear_callback(TautomerEnumerator &enumerator);

std::unique_ptr<TautomerEnumeratorResult> tautomer_enumerator_enumerate(
    const TautomerEnumerator &enumerator, const RWMol &mol);
// The default score is checked first (standardize.h's
// check_tautomer_score); a closure's score of INT_MIN counts as INT_MIN + 1:
// RDKit starts its search for the best tautomer below every other score,
// and dereferences null when no tautomer beats it.
std::unique_ptr<RWMol> tautomer_enumerator_canonicalize(
    const TautomerEnumerator &enumerator, const RWMol &mol);
std::unique_ptr<RWMol> tautomer_enumerator_canonicalize_with(
    const TautomerEnumerator &enumerator, const RWMol &mol,
    const ScoreFunction &score);
// Both pickCanonical overloads dereference the best tautomer, null for no
// tautomers, so an empty list or result is ValueErrorException.
std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical(
    const TautomerEnumerator &enumerator, BorrowedMols &mols);
std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical_with(
    const TautomerEnumerator &enumerator, BorrowedMols &mols,
    const ScoreFunction &score);
std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical_result(
    const TautomerEnumerator &enumerator,
    const TautomerEnumeratorResult &result);
std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical_result_with(
    const TautomerEnumerator &enumerator,
    const TautomerEnumeratorResult &result, const ScoreFunction &score);

std::unique_ptr<TautomerCallbackHandle> tautomer_callback_new(
    rust::Box<TautomerCallback> callback);
std::unique_ptr<TautomerCallbackHandle> tautomer_callback_clone(
    const TautomerCallbackHandle &handle);
bool tautomer_callback_call(const TautomerCallbackHandle &handle,
                            const RWMol &mol,
                            const TautomerEnumeratorResult &result);

std::size_t tautomer_result_size(const TautomerEnumeratorResult &result);
// Checks the position (E3).
std::unique_ptr<RWMol> tautomer_result_tautomer(
    const TautomerEnumeratorResult &result, std::size_t i);
// The kekulized form RDKit made while enumerating, or else one made now,
// as Tautomer::getKekulized makes it, without caching it.
std::unique_ptr<RWMol> tautomer_result_kekulized(
    const TautomerEnumeratorResult &result, std::size_t i);
std::unique_ptr<std::vector<std::string>> tautomer_result_smiles(
    const TautomerEnumeratorResult &result);
TautomerEnumeratorStatus tautomer_result_status(
    const TautomerEnumeratorResult &result);
rust::Vec<std::uint32_t> tautomer_result_modified_atoms(
    const TautomerEnumeratorResult &result);
rust::Vec<std::uint32_t> tautomer_result_modified_bonds(
    const TautomerEnumeratorResult &result);

// The scores. scoreRings and the default tautomer score are checked first
// (check_tautomer_score), and so are custom terms: their scores times their
// match counts (at most 1000) are summed in int.
std::int32_t score_rings(const RWMol &mol);
std::int32_t score_hetero_hs(const RWMol &mol);
std::int32_t score_substructs(const RWMol &mol);
std::int32_t score_substructs_with_terms(
    const RWMol &mol, rust::Slice<const SubstructTermData> terms);
std::int32_t score_tautomer(const RWMol &mol);
rust::Vec<SubstructTermData> default_tautomer_score_substructs();
std::unique_ptr<std::string> tautomer_score_version();

}  // namespace rdkit_cxx
