/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for GraphMol/GaussianShape (rdkit.Chem.rdGaussianShape): Gaussian
// shape overlays of molecules' conformers, with pharmacophore features.

#pragma once

#include <cstdint>
#include <memory>
#include <mutex>

#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using GaussianShapeInput = RDKit::GaussianShape::ShapeInput;
// Scoped enums, whose underlying type is int, so cxx checks their values.
using StartMode = RDKit::GaussianShape::StartMode;
using OptimMode = RDKit::GaussianShape::OptimMode;

// ShapeInput.cpp fills getPh4Patterns' static list of feature patterns on
// first use, without a lock. Every construction of a shape, the only path
// to it, holds this mutex (CONVENTIONS.md U4), before control_c_mutex(),
// which the patterns' matches need; so a closure that holds that one must
// not wait for it.
std::mutex &gaussian_shape_mutex();

// The shapes' Gaussians multiply, add, and divide their atoms' coordinates:
// throws ValueErrorException unless every conformer's are finite and within
// 1e8 (gaussian_shape.cpp says why).
void check_shape_coordinates(const RDKit::ROMol &mol);

// Shared structs, defined by the bridge (src/ffi/gaussian_shape.rs).
struct GaussianShapeInputOptionsData;
struct ShapeOverlayOptionsData;

GaussianShapeInputOptionsData gaussian_shape_input_options_default();
ShapeOverlayOptionsData shape_overlay_options_default();

std::unique_ptr<GaussianShapeInput> gaussian_shape_input_new(
    const RWMol &mol, std::int32_t conf_id,
    const GaussianShapeInputOptionsData &shape_opt,
    const ShapeOverlayOptionsData &overlay_opts);
std::unique_ptr<GaussianShapeInput> gaussian_shape_input_clone(
    const GaussianShapeInput &shape);
std::unique_ptr<std::string> gaussian_shape_input_get_smiles(
    const GaussianShapeInput &shape);
void gaussian_shape_input_set_active_shape(GaussianShapeInput &shape,
                                           std::uint32_t shape_num);
std::uint32_t gaussian_shape_input_get_active_shape(
    const GaussianShapeInput &shape);
std::uint32_t gaussian_shape_input_num_atoms(const GaussianShapeInput &shape);
std::uint32_t gaussian_shape_input_num_features(
    const GaussianShapeInput &shape);
std::uint32_t gaussian_shape_input_num_shapes(const GaussianShapeInput &shape);
double gaussian_shape_input_shape_volume(const GaussianShapeInput &shape);
double gaussian_shape_input_color_volume(const GaussianShapeInput &shape);
void gaussian_shape_input_normalize_coords(GaussianShapeInput &shape);
std::unique_ptr<RWMol> gaussian_shape_input_shape_to_mol(
    const GaussianShapeInput &shape, bool include_colors, bool with_bonds);
// The three scores, this shape's and the fit's shape numbers (as doubles),
// then the 16 values of the transformation.
rust::Vec<double> gaussian_shape_input_best_similarity(
    GaussianShapeInput &shape, const GaussianShapeInput &fit_shape,
    double threshold, const ShapeOverlayOptionsData &overlay_opts);
double gaussian_shape_input_max_possible_similarity(
    const GaussianShapeInput &shape, const GaussianShapeInput &fit_shape,
    const ShapeOverlayOptionsData &overlay_opts);

// The three scores.
rust::Vec<double> gaussian_align_mol(
    const RWMol &ref, RWMol &fit, const GaussianShapeInputOptionsData &ref_opts,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t ref_conf_id,
    std::int32_t fit_conf_id);
rust::Vec<double> gaussian_align_mol_shape_input(
    const GaussianShapeInput &ref_shape, RWMol &fit,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t fit_conf_id);
// The three scores, then the 16 values of the transformation.
rust::Vec<double> gaussian_align_shapes(
    const GaussianShapeInput &ref_shape, GaussianShapeInput &fit_shape,
    const ShapeOverlayOptionsData &overlay_opts);
// The three scores and the two overlap volumes.
rust::Vec<double> gaussian_score_mol(
    const RWMol &ref, const RWMol &fit,
    const GaussianShapeInputOptionsData &ref_opts,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t ref_conf_id,
    std::int32_t fit_conf_id);
rust::Vec<double> gaussian_score_mol_shape_input(
    const GaussianShapeInput &ref_shape, const RWMol &fit,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t fit_conf_id);
rust::Vec<double> gaussian_score_shape(
    const GaussianShapeInput &ref_shape, const GaussianShapeInput &fit_shape,
    const ShapeOverlayOptionsData &overlay_opts);
// The numbers of reference and fit shapes, the scores by rows, the best
// reference and fit shapes, then the 16 values of the transformation.
rust::Vec<double> gaussian_score_molecule_all_conformers(
    const RWMol &ref, const RWMol &fit,
    const GaussianShapeInputOptionsData &ref_opts,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts);

}  // namespace rdkit_cxx
