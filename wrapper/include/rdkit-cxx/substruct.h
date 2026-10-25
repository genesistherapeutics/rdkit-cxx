/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <GraphMol/Substruct/SubstructMatch.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_bundle.h"
#include "rdkit-cxx/resonance.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridge (src/ffi/substruct.rs).
struct SubstructMatchParametersData;
struct AtomCoordsMatcherData;

// Conversions for other areas' wrappers that hold the parameters.
SubstructMatchParametersData substruct_params_to_data(
    const RDKit::SubstructMatchParameters &params);
RDKit::SubstructMatchParameters substruct_params_from_data(
    const SubstructMatchParametersData &data);

SubstructMatchParametersData substruct_match_parameters_default();
AtomCoordsMatcherData atom_coords_matcher_default();
AtomCoordsMatcherData atom_coords_matcher_new(int ref_conf_id,
                                              int query_conf_id, double tol);
bool atom_coords_matcher_call(const AtomCoordsMatcherData &matcher,
                              const Atom &query_atom, const Atom &target_atom);

// As helpHasSubstructMatch, helpGetSubstructMatch, and
// helpGetSubstructMatches in Code/GraphMol/nbWrap/substructmethods.h, for each
// pairing of molecule and bundle.
bool mol_has_substruct_match(const RWMol &target, const RWMol &query,
                             const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> mol_get_substruct_match(
    const RWMol &target, const RWMol &query,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> mol_get_substruct_matches(
    const RWMol &target, const RWMol &query,
    const SubstructMatchParametersData &params);
bool mol_has_substruct_match_mol_bundle(
    const RWMol &target, const MolBundle &query,
    const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> mol_get_substruct_match_mol_bundle(
    const RWMol &target, const MolBundle &query,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> mol_get_substruct_matches_mol_bundle(
    const RWMol &target, const MolBundle &query,
    const SubstructMatchParametersData &params);
bool mol_bundle_has_substruct_match(const MolBundle &target, const RWMol &query,
                                    const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> mol_bundle_get_substruct_match(
    const MolBundle &target, const RWMol &query,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> mol_bundle_get_substruct_matches(
    const MolBundle &target, const RWMol &query,
    const SubstructMatchParametersData &params);
bool mol_bundle_has_substruct_match_mol_bundle(
    const MolBundle &target, const MolBundle &query,
    const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> mol_bundle_get_substruct_match_mol_bundle(
    const MolBundle &target, const MolBundle &query,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> mol_bundle_get_substruct_matches_mol_bundle(
    const MolBundle &target, const MolBundle &query,
    const SubstructMatchParametersData &params);

// As GetSubstructMatch and GetResonanceSubstructMatches in
// Code/GraphMol/nbWrap/ResonanceMolSupplier.cpp.
rust::Vec<std::uint32_t> resonance_mol_supplier_get_substruct_match(
    ResonanceMolSupplier &suppl, const RWMol &query, bool use_chirality,
    bool use_query_query_matches);
rust::Vec<IndexList> resonance_mol_supplier_get_substruct_matches(
    ResonanceMolSupplier &suppl, const RWMol &query, bool uniquify,
    bool use_chirality, bool use_query_query_matches, std::uint32_t max_matches,
    int num_threads);

}  // namespace rdkit_cxx
