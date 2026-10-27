/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolops's graph matrices, paths, and subgraphs.

#pragma once

#include <cstdint>
#include <map>
#include <memory>

#include <GraphMol/ROMol.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/graph.rs, geometry.rs,
// chirality.rs).
struct IndexList;
struct IndexPair;

rust::Vec<double> get_distance_matrix(const RWMol &mol, bool use_bo,
                                      bool use_atom_wts);
rust::Vec<double> get_3d_distance_matrix(const RWMol &mol, std::int32_t conf_id,
                                         bool use_atom_wts);
rust::Vec<double> get_adjacency_matrix(const RWMol &mol, bool use_bo,
                                       std::int32_t empty_val);

rust::Vec<IndexList> find_all_paths_of_length_n(const RWMol &mol,
                                                std::uint32_t length,
                                                bool use_bonds, bool use_hs,
                                                std::int32_t rooted_at_atom,
                                                bool only_shortest_paths);
rust::Vec<IndexList> find_all_subgraphs_of_length_n(
    const RWMol &mol, std::uint32_t length, bool use_hs,
    std::int32_t rooted_at_atom);
rust::Vec<IndexList> find_all_subgraphs_of_length_m_to_n(
    const RWMol &mol, std::uint32_t min, std::uint32_t max, bool use_hs,
    std::int32_t rooted_at_atom, rust::Vec<std::uint32_t> &counts);
rust::Vec<IndexList> find_unique_subgraphs_of_length_n(
    const RWMol &mol, std::uint32_t length, bool use_hs, bool use_bo,
    std::int32_t rooted_at_atom);
rust::Vec<std::uint32_t> find_atom_environment_of_radius_n(
    const RWMol &mol, std::uint32_t radius, std::uint32_t rooted_at_atom,
    bool use_hs, bool enforce_size, rust::Vec<IndexPair> &atom_map);
std::unique_ptr<RWMol> path_to_submol(const RWMol &mol,
                                      rust::Slice<const std::uint32_t> path,
                                      bool use_query,
                                      rust::Vec<IndexPair> &atom_map);

std::unique_ptr<MolList> split_mol_by_pdb_residues(
    const RWMol &mol, bool has_white_list,
    rust::Slice<const rust::String> white_list, bool negate_list);
std::unique_ptr<MolList> split_mol_by_pdb_chain_id(
    const RWMol &mol, bool has_white_list,
    rust::Slice<const rust::String> white_list, bool negate_list);

std::unique_ptr<RWMol> renumber_atoms(
    const RWMol &mol, rust::Slice<const std::uint32_t> new_order);

}  // namespace rdkit_cxx
