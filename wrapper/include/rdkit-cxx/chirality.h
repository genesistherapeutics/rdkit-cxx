/*
 * Copyright (c) 2026 Steven Kearnes. All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolops's stereochemistry perception, wedging, and enhanced stereo.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>

#include "rdkit-cxx/bond.h"
#include "rdkit-cxx/conformer.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/stereo.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using StereoGroupAbsOptions = RDKit::StereoGroupAbsOptions;

// Shared structs, defined by the bridge (src/ffi/chirality.rs).
struct IndexPair;
struct BondWedgingParametersData;

void assign_stereochemistry(RWMol &mol, bool clean_it, bool force,
                            bool flag_possible_stereo_centers);
void assign_stereochemistry_from_3d(RWMol &mol, std::int32_t conf_id,
                                    bool replace_existing_tags);
void assign_chiral_types_from_bond_dirs(RWMol &mol, std::int32_t conf_id,
                                        bool replace_existing_tags);
void assign_atom_chiral_tags_from_structure(RWMol &mol, std::int32_t conf_id,
                                            bool replace_existing_tags);
void assign_atom_chiral_tags_from_mol_parity(RWMol &mol,
                                             bool replace_existing_tags);
void find_potential_stereo_bonds(RWMol &mol, bool clean_it);
void remove_stereochemistry(RWMol &mol);
void remove_non_explicit_3d_chirality(RWMol &mol);
void detect_bond_stereochemistry(RWMol &mol, std::int32_t conf_id);
void set_double_bond_neighbor_directions(RWMol &mol);
void set_double_bond_neighbor_directions_with_conf(RWMol &mol,
                                                   const Conformer &conf);
void set_bond_stereo_from_directions(RWMol &mol);
rust::Vec<std::uint32_t> compute_atom_cip_ranks(RWMol &mol);
rust::Vec<IndexPair> find_meso_centers(const RWMol &mol, bool include_isotopes,
                                       bool include_atom_maps);

std::unique_ptr<RWMol> canonicalize_stereo_groups(
    const RWMol &mol, StereoGroupAbsOptions output_absolute_groups,
    std::uint32_t max_stereo_groups);
void cleanup_stereo_groups(RWMol &mol);
void simplify_enhanced_stereo(RWMol &mol, bool remove_affected_stereo_groups);
void translate_chiral_flag_to_stereo_groups(
    RWMol &mol, StereoGroupType zero_flag_group_type);
void add_stereo_annotations(RWMol &mol, rust::Str abs_label, rust::Str or_label,
                            rust::Str and_label, rust::Str cip_label,
                            rust::Str bond_label);

void add_wavy_bonds_for_stereo_any(RWMol &mol, bool clear_double_bond_flags,
                                   std::uint32_t add_when_impossible);
// The conformer is the bond's molecule's, as RDKit requires.
void wedge_bond(Bond &bond, std::uint32_t from_atom_idx, std::int32_t conf_id);
BondWedgingParametersData bond_wedging_parameters_default();
void wedge_mol_bonds(RWMol &mol, const Conformer &conformer,
                     const BondWedgingParametersData &params);
void reapply_mol_block_wedging(RWMol &mol, bool all_bond_types, bool verify);

// These settings live in environment variables, which RDKit sets with
// setenv; the Rust setters are unsafe for that reason (CONVENTIONS.md U5).
void set_allow_nontetrahedral_chirality(bool val);
bool get_allow_nontetrahedral_chirality();
void set_use_legacy_stereo_perception(bool val);
bool get_use_legacy_stereo_perception();

void expand_attachment_points(RWMol &mol, bool add_as_queries, bool add_coords);
void collapse_attachment_points(RWMol &mol, bool marked_only);
// Each problem, encoded as the trycatch encodes exceptions (errors.h).
std::unique_ptr<std::vector<std::string>> detect_chemistry_problems(
    const RWMol &mol, std::uint32_t sanitize_ops);

}  // namespace rdkit_cxx
