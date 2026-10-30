/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdTautomerQuery: substructure queries that match any tautomer of a molecule.
// RDKit's TautomerQuery holds the template and tautomers; the matching and the
// construction from a molecule are the wrapper's own (tautomer_query.cpp says
// why), so generalized_substruct.cpp shares them.

#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using TautomerQuery = RDKit::TautomerQuery;

// Shared structs, defined by other bridges.
struct IndexList;
struct SubstructMatchParametersData;

// For the wrappers that hold tautomer queries. A query read from a pickle or
// JSON is checked with check_tautomer_query, which throws
// ValueErrorException for one whose indices the matcher would read past;
// tautomer_query_substruct_of is TautomerQuery::substructOf, whose caller
// holds control_c_mutex (control_c_lock.h).
void check_tautomer_query(const TautomerQuery &query);
// check_picklable (mol_checks.h) for the template and each tautomer, which
// TautomerQuery::serialize pickles as molecules.
void check_picklable_tautomer_query(const TautomerQuery &query);
std::vector<RDKit::MatchVectType> tautomer_query_substruct_of(
    const TautomerQuery &query, const RDKit::ROMol &mol,
    const RDKit::SubstructMatchParameters &params,
    std::vector<RDKit::ROMOL_SPTR> *matching_tautomers);
// Fingerprints are ExplicitBitVects, which hold at most INT_MAX bits (U3).
void check_fingerprint_size(std::uint32_t fingerprint_size);
// As TautomerQuery::fromMol(query), with RDKit's default transforms; it
// replaces the template's bonds with replace_bond (mol_checks.h). The
// enumerator matches its transforms, so the caller holds control_c_mutex
// (control_c_lock.h).
std::unique_ptr<TautomerQuery> make_tautomer_query(const RDKit::ROMol &query);

std::unique_ptr<TautomerQuery> tautomer_query_clone(const TautomerQuery &query);
std::unique_ptr<TautomerQuery> tautomer_query_from_mol(
    const RWMol &mol, rust::Slice<const std::uint8_t> tautomer_transform_file);
std::unique_ptr<TautomerQuery> tautomer_query_from_pickle(
    rust::Slice<const std::uint8_t> pickle);
rust::Vec<std::uint8_t> tautomer_query_to_binary(const TautomerQuery &query);

bool tautomer_query_is_substruct_of(const TautomerQuery &query,
                                    const RWMol &target,
                                    const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> tautomer_query_get_substruct_match(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> tautomer_query_get_substruct_matches(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params);
// Each match's tautomer, in the order of the matches, which fill `matches`.
std::unique_ptr<MolList> tautomer_query_get_substruct_matches_with_tautomers(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params, rust::Vec<IndexList> &matches);

std::unique_ptr<BitVect> tautomer_query_pattern_fingerprint_template(
    const TautomerQuery &query, std::uint32_t fingerprint_size);
std::unique_ptr<RWMol> tautomer_query_get_template_molecule(
    const TautomerQuery &query);
rust::Vec<std::uint32_t> tautomer_query_get_modified_atoms(
    const TautomerQuery &query);
rust::Vec<std::uint32_t> tautomer_query_get_modified_bonds(
    const TautomerQuery &query);
std::unique_ptr<MolList> tautomer_query_get_tautomers(
    const TautomerQuery &query);

std::unique_ptr<BitVect> pattern_fingerprint_tautomer_target(
    const RWMol &target, std::uint32_t fingerprint_size);
bool tautomer_query_can_serialize();

}  // namespace rdkit_cxx
