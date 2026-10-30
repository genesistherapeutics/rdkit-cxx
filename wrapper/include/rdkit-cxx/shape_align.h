/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for External/pubchem_shape (rdkit.Chem.rdShapeAlign): PubChem's
// Gaussian shape alignment. The shapes are read-only from Rust, since the
// alignment trusts their vectors to agree with each other.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/PubChemShape.hpp>

#include "rdkit-cxx/conformer.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// PubChemShape.hpp declares its ShapeInput in the global namespace.
using PubChemShapeInput = ::ShapeInput;

// Shared structs, defined by the bridges (src/ffi/shape_align.rs and
// src/ffi/geometry.rs).
struct ShapeInputOptionsData;
struct Point3D;

ShapeInputOptionsData pubchem_shape_input_options_default();

std::unique_ptr<PubChemShapeInput> pubchem_prepare_conformer(
    const RWMol &mol, std::int32_t conf_id, const ShapeInputOptionsData &opts);
std::unique_ptr<PubChemShapeInput> pubchem_shape_input_clone(
    const PubChemShapeInput &shape);
rust::Vec<float> pubchem_shape_input_coord(const PubChemShapeInput &shape);
rust::Vec<double> pubchem_shape_input_alpha_vector(
    const PubChemShapeInput &shape);
rust::Vec<std::uint32_t> pubchem_shape_input_atom_type_vector(
    const PubChemShapeInput &shape);
rust::Vec<std::uint32_t> pubchem_shape_input_volume_atom_index_vector(
    const PubChemShapeInput &shape);
rust::Vec<double> pubchem_shape_input_shift(const PubChemShapeInput &shape);
rust::Vec<double> pubchem_shape_input_inertial_rot(
    const PubChemShapeInput &shape);
double pubchem_shape_input_sov(const PubChemShapeInput &shape);
double pubchem_shape_input_sof(const PubChemShapeInput &shape);

// The scores, shape then color.
rust::Vec<double> pubchem_align_mol(
    const RWMol &ref, RWMol &probe, const ShapeInputOptionsData &ref_opts,
    const ShapeInputOptionsData &probe_opts, std::int32_t ref_conf_id,
    std::int32_t probe_conf_id, double opt_param, std::uint32_t max_preiters,
    std::uint32_t max_postiters);
rust::Vec<double> pubchem_align_mol_shape_input(
    const PubChemShapeInput &ref_shape, RWMol &probe,
    std::int32_t probe_conf_id, bool use_colors, double opt_param,
    std::uint32_t max_preiters, std::uint32_t max_postiters,
    bool apply_ref_shift);
// The scores, then the 12 values of the matrix.
rust::Vec<double> pubchem_align_shapes(const PubChemShapeInput &ref_shape,
                                       PubChemShapeInput &probe_shape,
                                       double opt_param,
                                       std::uint32_t max_preiters,
                                       std::uint32_t max_postiters);
void pubchem_transform_conformer_to_shape(rust::Slice<const double> final_trans,
                                          rust::Slice<const double> final_rot,
                                          rust::Slice<const double> matrix,
                                          const PubChemShapeInput &probe_shape,
                                          Conformer &probe_conformer);
rust::Vec<double> pubchem_score_mol(const RWMol &mol1, const RWMol &mol2,
                                    const ShapeInputOptionsData &mol1_opts,
                                    const ShapeInputOptionsData &mol2_opts,
                                    std::int32_t mol1_conf_id,
                                    std::int32_t mol2_conf_id);
rust::Vec<double> pubchem_score_mol_shape_input(
    const PubChemShapeInput &shape, const RWMol &mol,
    const ShapeInputOptionsData &mol_opts, std::int32_t mol_conf_id);
rust::Vec<double> pubchem_score_shape(const PubChemShapeInput &shape1,
                                      const PubChemShapeInput &shape2,
                                      bool use_colors);

}  // namespace rdkit_cxx
