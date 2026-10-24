// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A shared struct, defined by the bridge (src/ffi/substruct.rs).
struct SubstructMatchParametersData;

SubstructMatchParametersData substruct_match_parameters_default();

// As helpHasSubstructMatch, helpGetSubstructMatch, and
// helpGetSubstructMatches in Code/GraphMol/nbWrap/substructmethods.h.
bool mol_has_substruct_match(const RWMol &mol, const RWMol &query,
                             const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> mol_get_substruct_match(
    const RWMol &mol, const RWMol &query,
    const SubstructMatchParametersData &params);
rust::Vec<IndexList> mol_get_substruct_matches(
    const RWMol &mol, const RWMol &query,
    const SubstructMatchParametersData &params);

}  // namespace rdkit_cxx
