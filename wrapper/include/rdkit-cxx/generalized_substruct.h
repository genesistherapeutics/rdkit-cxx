/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdGeneralizedSubstruct: extended query molecules, which hold a molecule, a
// bundle of them, a tautomer query, or a bundle of tautomer queries.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/GeneralizedSubstruct/XQMol.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using ExtendedQueryMol = RDKit::GeneralizedSubstruct::ExtendedQueryMol;

// For the wrappers that enumerate molecules with MolEnumerator: throws
// ValueErrorException for link nodes the enumerator would read past an
// empty vector on (generalized_substruct.cpp says why).
void check_link_nodes(const RDKit::ROMol &mol);

// Shared structs, defined by other bridges.
struct AdjustQueryParametersData;
struct IndexList;
struct SubstructMatchParametersData;

std::unique_ptr<ExtendedQueryMol> extended_query_mol_from_binary(
    rust::Slice<const std::uint8_t> pickle);
std::unique_ptr<ExtendedQueryMol> extended_query_mol_from_json(rust::Str text);
std::unique_ptr<ExtendedQueryMol> extended_query_mol_clone(
    const ExtendedQueryMol &query);
rust::Vec<std::uint8_t> extended_query_mol_to_binary(
    const ExtendedQueryMol &query);
std::unique_ptr<std::string> extended_query_mol_to_json(
    const ExtendedQueryMol &query);
std::unique_ptr<BitVect> extended_query_mol_pattern_fingerprint_query(
    const ExtendedQueryMol &query, std::uint32_t fingerprint_size);

std::unique_ptr<ExtendedQueryMol> create_extended_query_mol(
    const RWMol &mol, bool do_enumeration, bool do_tautomers,
    bool adjust_query_properties,
    const AdjustQueryParametersData &adjust_query_parameters);

bool mol_has_substruct_match_extended(
    const RWMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> mol_get_substruct_match_extended(
    const RWMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> mol_get_substruct_matches_extended(
    const RWMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &params);
std::unique_ptr<BitVect> pattern_fingerprint_target(
    const RWMol &target, std::uint32_t fingerprint_size);

}  // namespace rdkit_cxx
