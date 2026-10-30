/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdSynthonSpaceSearch: synthon spaces (combinatorial libraries such as
// Enamine REAL, held as reactions and their synthons) and their substructure,
// fingerprint, and RASCAL searches.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>

#include "rdkit-cxx/fingerprint_generator.h"
#include "rdkit-cxx/generalized_substruct.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// The space a Rust SynthonSpace owns. RDKit's SynthonSpace can be neither
// copied nor moved, and its readers add to a space in place, so each read
// builds a new space and replaces this one only when the read succeeds.
class SynthonSpace {
 public:
  SynthonSpace();
  std::unique_ptr<RDKit::SynthonSpaceSearch::SynthonSpace> space;
};

// Python's SubstructureResult.
using SearchResults = RDKit::SynthonSpaceSearch::SearchResults;

// Shared structs and the Rust callback, defined by the bridges
// (src/ffi/synthon_space_search.rs and the bridges they name).
struct SynthonSpaceSearchParamsData;
struct SubstructMatchParametersData;
struct RascalOptionsData;
struct SearchCallback;

SynthonSpaceSearchParamsData synthon_space_search_params_default();

std::unique_ptr<SynthonSpace> synthon_space_new();
void synthon_space_read_text_file(SynthonSpace &space,
                                  rust::Slice<const std::uint8_t> path);
void synthon_space_read_db_file(SynthonSpace &space,
                                rust::Slice<const std::uint8_t> path);
void synthon_space_write_db_file(const SynthonSpace &space,
                                 rust::Slice<const std::uint8_t> path);
void synthon_space_write_enumerated_file(const SynthonSpace &space,
                                         rust::Slice<const std::uint8_t> path);
std::uint64_t synthon_space_get_num_reactions(const SynthonSpace &space);
std::uint64_t synthon_space_get_num_products(const SynthonSpace &space);
void synthon_space_summarise(const SynthonSpace &space);
std::unique_ptr<std::string> synthon_space_get_synthon_fingerprint_type(
    const SynthonSpace &space);
void synthon_space_build_synthon_fingerprints(SynthonSpace &space,
                                              const FpGenerator &gen);

std::unique_ptr<SearchResults> synthon_space_substructure_search(
    SynthonSpace &space, const RWMol &query,
    const SubstructMatchParametersData &match_params,
    const SynthonSpaceSearchParamsData &params);
std::unique_ptr<SearchResults> synthon_space_substructure_search_extended(
    SynthonSpace &space, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &match_params,
    const SynthonSpaceSearchParamsData &params);
void synthon_space_substructure_search_incremental(
    SynthonSpace &space, const RWMol &query, SearchCallback &callback,
    const SubstructMatchParametersData &match_params,
    const SynthonSpaceSearchParamsData &params);

std::unique_ptr<SearchResults> synthon_space_fingerprint_search(
    SynthonSpace &space, const RWMol &query, const FpGenerator &gen,
    const SynthonSpaceSearchParamsData &params);
void synthon_space_fingerprint_search_incremental(
    SynthonSpace &space, const RWMol &query, const FpGenerator &gen,
    SearchCallback &callback, const SynthonSpaceSearchParamsData &params);

std::unique_ptr<SearchResults> synthon_space_rascal_search(
    SynthonSpace &space, const RWMol &query,
    const RascalOptionsData &rascal_options,
    const SynthonSpaceSearchParamsData &params);
void synthon_space_rascal_search_incremental(
    SynthonSpace &space, const RWMol &query,
    const RascalOptionsData &rascal_options, SearchCallback &callback,
    const SynthonSpaceSearchParamsData &params);

std::unique_ptr<MolList> substructure_result_get_hit_molecules(
    const SearchResults &results);
std::uint64_t substructure_result_get_max_num_results(
    const SearchResults &results);
bool substructure_result_get_timed_out(const SearchResults &results);
bool substructure_result_get_cancelled(const SearchResults &results);

std::unique_ptr<std::string> formatted_integer_string(std::int64_t value);

}  // namespace rdkit_cxx
