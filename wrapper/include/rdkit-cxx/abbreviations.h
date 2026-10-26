/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdAbbreviations: abbreviations (and linkers) found in molecules, condensed
// to labeled dummy atoms or labeled with SUP substance groups.

#pragma once

#include <memory>
#include <vector>

#include <GraphMol/Abbreviations/Abbreviations.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared type, defined by the bridge (src/ffi/abbreviations.rs).
struct AbbreviationDefinitionData;

// Definitions Rust lends to one call. Their molecules are shared_ptrs with
// no-op deleters into Rust's molecules (CONVENTIONS.md W8): Rust builds the
// list, makes the call, and drops the list before its borrows end.
class AbbreviationList {
 public:
  std::vector<RDKit::Abbreviations::AbbreviationDefinition> defs;
};

std::unique_ptr<AbbreviationList> abbreviation_list_new();
void abbreviation_list_push(AbbreviationList &list, rust::Str label,
                            rust::Str display_label, rust::Str display_label_w,
                            bool includes_x_bonds);
void abbreviation_list_push_with_mol(AbbreviationList &list, rust::Str label,
                                     rust::Str display_label,
                                     rust::Str display_label_w,
                                     bool includes_x_bonds, const RWMol &mol);

// A default-constructed definition, without a molecule.
AbbreviationDefinitionData abbreviation_definition_default();
// RDKit's definitions, with copies of their molecules.
rust::Vec<AbbreviationDefinitionData> get_default_abbreviations();
rust::Vec<AbbreviationDefinitionData> get_default_linkers();
// As Utils::parseAbbreviations, skipping lines without a SMARTS.
rust::Vec<AbbreviationDefinitionData> parse_abbreviations(
    rust::Str text, bool remove_extra_dummies,
    bool allow_connection_to_dummies);

// As the lambdas in nbWrap/rdAbbreviations.cpp, each of which edits a copy.
std::unique_ptr<RWMol> condense_mol_abbreviations(
    const RWMol &mol, const AbbreviationList &abbrevs, double max_coverage,
    bool sanitize);
std::unique_ptr<RWMol> label_mol_abbreviations(const RWMol &mol,
                                               const AbbreviationList &abbrevs,
                                               double max_coverage);
std::unique_ptr<RWMol> condense_abbreviation_substance_groups(const RWMol &mol);

}  // namespace rdkit_cxx
