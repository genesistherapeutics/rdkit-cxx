/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_bundle.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A shared struct, defined by the bridge (src/ffi/substruct.rs).
struct SubstructMatchParametersData;

SubstructMatchParametersData substruct_match_parameters_default();

// As helpHasSubstructMatch, helpGetSubstructMatch, and
// helpGetSubstructMatches in Code/GraphMol/nbWrap/substructmethods.h, for each
// pairing of molecule and bundle.
bool mol_has_substruct_match(
    const RWMol &target, const RWMol &query,
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
bool mol_bundle_has_substruct_match(
    const MolBundle &target, const RWMol &query,
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

}  // namespace rdkit_cxx
