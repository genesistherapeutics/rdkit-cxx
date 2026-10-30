/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdChemReactions' EnumerateLibrary and its enumeration strategies.
//
// A library owns copies of its building blocks: Python's binding makes the
// library share Python's molecules without keeping them alive, so a library
// whose molecules Python has collected reads freed memory. The strategies'
// states are checked wherever they come from outside RDKit's own
// constructors (check_state in enumerate.cpp).

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/ChemReactions/Enumerate/Enumerate.h>
#include <GraphMol/ChemReactions/Enumerate/EnumerationStrategyBase.h>

#include "rdkit-cxx/chem_reactions.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using EnumerateLibrary = RDKit::EnumerateLibrary;
using EnumerationStrategyBase = RDKit::EnumerationStrategyBase;

// Shared struct, defined by the bridge (src/ffi/enumerate.rs).
struct EnumerationParamsData;

// The strategies by kind: 0 CartesianProduct, 1 RandomSample,
// 2 RandomSampleAllBBs, 3 EvenSamplePairs.
std::unique_ptr<EnumerationStrategyBase> enumeration_strategy_new(
    std::uint8_t kind);
std::unique_ptr<EnumerationStrategyBase> enumeration_strategy_clone(
    const EnumerationStrategyBase &strategy);
// Python's Initialize. The molecules are lent for the call, in groups of
// `group_sizes`; the strategies read only how many there are.
void enumeration_strategy_initialize(
    EnumerationStrategyBase &strategy, const ChemicalReaction &rxn,
    BorrowedMols &mols, rust::Slice<const std::uint32_t> group_sizes);
rust::Vec<std::uint64_t> enumeration_strategy_get_position(
    const EnumerationStrategyBase &strategy);
std::uint64_t enumeration_strategy_get_num_permutations(
    const EnumerationStrategyBase &strategy);
std::unique_ptr<std::string> enumeration_strategy_type(
    const EnumerationStrategyBase &strategy);
bool enumeration_strategy_has_next(const EnumerationStrategyBase &strategy);
rust::Vec<std::uint64_t> enumeration_strategy_next(
    EnumerationStrategyBase &strategy);
bool enumeration_strategy_skip(EnumerationStrategyBase &strategy,
                               std::uint64_t skip_count);
std::unique_ptr<std::string> even_sample_pairs_strategy_stats(
    const EnumerationStrategyBase &strategy);

EnumerationParamsData enumeration_params_default();
bool enumerate_library_can_serialize();

std::unique_ptr<EnumerateLibrary> enumerate_library_new();
// Copy the building blocks.
std::unique_ptr<EnumerateLibrary> enumerate_library_new_with_params(
    const ChemicalReaction &rxn, BorrowedMols &mols,
    rust::Slice<const std::uint32_t> group_sizes,
    const EnumerationParamsData &params);
std::unique_ptr<EnumerateLibrary> enumerate_library_new_with_enumerator(
    const ChemicalReaction &rxn, BorrowedMols &mols,
    rust::Slice<const std::uint32_t> group_sizes,
    const EnumerationStrategyBase &strategy,
    const EnumerationParamsData &params);
// The building blocks the reaction's templates match, copied.
std::unique_ptr<ProductSets> enumerate_library_get_reagents(
    const EnumerateLibrary &lib);
std::unique_ptr<ChemicalReaction> enumerate_library_get_reaction(
    const EnumerateLibrary &lib);
std::unique_ptr<EnumerationStrategyBase> enumerate_library_get_enumerator(
    const EnumerateLibrary &lib);
rust::Vec<std::uint64_t> enumerate_library_get_position(
    const EnumerateLibrary &lib);
bool enumerate_library_has_next(const EnumerateLibrary &lib);
// EnumerateLibrary::next, with the strategy's indices checked before they
// index the building blocks.
std::unique_ptr<ProductSets> enumerate_library_next(EnumerateLibrary &lib);
std::unique_ptr<std::string> enumerate_library_get_state(
    const EnumerateLibrary &lib);
void enumerate_library_set_state(EnumerateLibrary &lib, rust::Str state);
void enumerate_library_reset_state(EnumerateLibrary &lib);
rust::Vec<std::uint8_t> enumerate_library_serialize(
    const EnumerateLibrary &lib);
bool enumerate_library_is_empty(const EnumerateLibrary &lib);
// A library from Serialize's bytes, checked.
std::unique_ptr<EnumerateLibrary> enumerate_library_from_string(
    rust::Slice<const std::uint8_t> data);

}  // namespace rdkit_cxx
