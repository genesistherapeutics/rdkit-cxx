/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolops's chemical transformations: substructure replacement, core and
// sidechain replacement, fragmentation, molzip, and recursive queries.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/ChemTransforms/MolFragmenter.h>

#include "rdkit-cxx/bond.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A plain enum class unless USE_BETTER_ENUMS is defined, which the wrappers
// never do (RDGeneral/BetterEnums.h).
using MolzipLabel = RDKit::MolzipLabel;

// Shared structs, defined by the bridges (src/ffi/transforms.rs,
// geometry.rs, chirality.rs).
struct IndexList;
struct IndexPair;
struct Point3D;
struct MolzipParamsData;

std::unique_ptr<RWMol> delete_substructs(const RWMol &mol, const RWMol &query,
                                         bool only_frags, bool use_chirality);
std::unique_ptr<MolList> replace_substructs(
    const RWMol &mol, const RWMol &query, const RWMol &replacement,
    bool replace_all, std::uint32_t replacement_connection_point,
    bool use_chirality);
std::unique_ptr<RWMol> replace_sidechains(const RWMol &mol,
                                          const RWMol &core_query,
                                          bool use_chirality);
std::unique_ptr<RWMol> replace_core_with_matches(
    const RWMol &mol, const RWMol &core, rust::Slice<const IndexPair> matches,
    bool replace_dummies, bool label_by_index, bool require_dummy_match);
std::unique_ptr<RWMol> replace_core(const RWMol &mol, const RWMol &core_query,
                                    bool replace_dummies, bool label_by_index,
                                    bool require_dummy_match,
                                    bool use_chirality);
std::unique_ptr<RWMol> murcko_decompose(const RWMol &mol);
std::unique_ptr<RWMol> combine_mols(const RWMol &mol1, const RWMol &mol2,
                                    const Point3D &offset);
rust::Vec<std::uint32_t> get_most_substituted_core_match(
    const RWMol &mol, const RWMol &core, rust::Slice<const IndexList> matches);
rust::Vec<IndexList> sort_matches_by_degree_of_core_substitution(
    const RWMol &mol, const RWMol &core, rust::Slice<const IndexList> matches);

void add_recursive_query(RWMol &mol, const RWMol &query, std::uint32_t atom_idx,
                         bool preserve_existing_query);
void add_recursive_queries(RWMol &mol, BorrowedMols &queries,
                           rust::Str prop_name);
std::unique_ptr<MolList> parse_mol_query_def_file(
    rust::Slice<const std::uint8_t> path, bool standardize, rust::Str delimiter,
    rust::Str comment, std::uint32_t name_column, std::uint32_t smarts_column);
std::unique_ptr<MolList> parse_mol_query_def_text(
    rust::Str text, bool standardize, rust::Str delimiter, rust::Str comment,
    std::uint32_t name_column, std::uint32_t smarts_column);

std::unique_ptr<RWMol> fragment_on_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    bool add_dummies, bool has_dummy_labels,
    rust::Slice<const IndexPair> dummy_labels, bool has_bond_types,
    rust::Slice<const BondType> bond_types,
    rust::Vec<std::uint32_t> &cuts_per_atom);
std::unique_ptr<MolList> fragment_on_some_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    std::uint32_t num_to_break, bool add_dummies, bool has_dummy_labels,
    rust::Slice<const IndexPair> dummy_labels, bool has_bond_types,
    rust::Slice<const BondType> bond_types, bool return_cuts_per_atom,
    rust::Vec<IndexList> &cuts_per_atom);
std::unique_ptr<RWMol> fragment_on_brics_bonds(const RWMol &mol);

MolzipParamsData molzip_params_default();
std::unique_ptr<RWMol> molzip(const RWMol &a, const RWMol &b,
                              const MolzipParamsData &params);
std::unique_ptr<RWMol> molzip_single(const RWMol &a,
                                     const MolzipParamsData &params);
std::unique_ptr<RWMol> molzip_fragments(BorrowedMols &mols,
                                        const MolzipParamsData &params);
std::unique_ptr<RWMol> molzip_row(BorrowedMols &row,
                                  const MolzipParamsData &params);

}  // namespace rdkit_cxx
