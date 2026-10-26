/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdDistGeom: conformer embedding, bounds matrices, and experimental
// torsions.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/DistGeomHelpers/Embedder.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using EmbedFailureCauses = RDKit::DGeomHelpers::EmbedFailureCauses;
using InitialEmbeddingMode = RDKit::DGeomHelpers::InitialEmbeddingMode;
using EmbedFF = RDKit::DGeomHelpers::EmbedFF;

// Shared types, defined by the bridges (src/ffi/embedder.rs, geometry.rs,
// chirality.rs).
struct Point3D;
struct IndexPair;
struct EmbedParametersData;
struct ExperimentalTorsionData;
enum class EmbedPreset : std::uint8_t;

// A default-constructed EmbedParameters, or one of RDKit's presets.
EmbedParametersData embed_parameters_preset(EmbedPreset preset);
// Embeds the molecule as DGeomHelpers::EmbedMultipleConfs, and writes back
// what RDKit changes in the parameters.
rust::Vec<std::int32_t> embed_multiple_confs(RWMol &mol,
                                             std::uint32_t num_confs,
                                             EmbedParametersData &params);
std::unique_ptr<std::string> embed_parameters_to_json(
    const EmbedParametersData &params);
rust::Vec<double> get_molecule_bounds_matrix(const RWMol &mol,
                                             const EmbedParametersData &params,
                                             bool do_triangle_smoothing,
                                             bool scale_vdw, bool set15bounds,
                                             bool set14bounds,
                                             bool set13bounds);
rust::Vec<ExperimentalTorsionData> get_experimental_torsions(
    const RWMol &mol, bool use_exp_torsion_angle_prefs,
    bool use_small_ring_torsions, bool use_macrocycle_torsions,
    bool use_basic_knowledge, std::uint32_t et_version, bool verbose);

// The checks embed_multiple_confs makes of the molecule (U1, U3), for RDKit
// code that embeds it itself with these parameters, which have no
// coordinate map, CPCI, or bounds matrix. Throws ValueErrorException.
void check_embed_mol(const RDKit::ROMol &mol,
                     const RDKit::DGeomHelpers::EmbedParameters &params);

// RDKit's computeInitialCoords (DistGeom/DistGeomUtils.cpp) casts
// sum(d_ij^2) / n, over the distances it picks between the bounds, to int,
// which is undefined above INT_MAX. Throws ValueErrorException unless
// (n - 1) * max_bound^2 / 2 fits (U1).
void check_initial_coords_cast(unsigned n, double max_bound);

}  // namespace rdkit_cxx
