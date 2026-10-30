/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdDepictor: 2D coordinates. Every function that computes coordinates or
// reads or changes the depiction settings holds lock_depiction()
// (CONVENTIONS.md U4).

#pragma once

#include <cstdint>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct Compute2DCoordsData;
struct Compute2DCoordsMimicDistmatData;
struct ConstrainedDepictionParamsData;
struct IndexPair;

bool is_coord_gen_support_available();
void set_prefer_coord_gen(bool val);
bool get_prefer_coord_gen();
// Sets the preference and returns the one it replaces, in one locked step.
bool replace_prefer_coord_gen(bool val);

void set_ring_system_templates(rust::Slice<const std::uint8_t> path);
void add_ring_system_templates(rust::Slice<const std::uint8_t> path);
void load_default_ring_system_templates();

ConstrainedDepictionParamsData constrained_depiction_params_default();

unsigned compute_2d_coords(RWMol &mol, const Compute2DCoordsData &data);
unsigned compute_2d_coords_mimic_distmat(
    RWMol &mol, rust::Slice<const double> dist_mat,
    const Compute2DCoordsMimicDistmatData &data);

// ref_patt is used only with has_ref_patt.
rust::Vec<IndexPair> generate_depiction_matching_2d_structure(
    RWMol &mol, const RWMol &reference, const RWMol &ref_patt,
    bool has_ref_patt, const ConstrainedDepictionParamsData &data);
void generate_depiction_matching_2d_structure_atom_map(
    RWMol &mol, const RWMol &reference, rust::Slice<const IndexPair> atom_map,
    const ConstrainedDepictionParamsData &data);
void generate_depiction_matching_3d_structure(
    RWMol &mol, const RWMol &reference, int conf_id, const RWMol &ref_patt,
    bool has_ref_patt, bool accept_failure, bool force_rdkit);

void straighten_depiction(RWMol &mol, int conf_id, bool minimize_rotation);
double normalize_depiction(RWMol &mol, int conf_id, int canonicalize,
                           double scale_factor);

}  // namespace rdkit_cxx
