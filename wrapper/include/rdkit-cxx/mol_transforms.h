/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolTransforms: centroids, principal axes, canonical orientations, and
// internal coordinates of conformers.

#pragma once

#include <cstdint>

#include "rdkit-cxx/conformer.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/geometry.rs).
struct Point3D;

Point3D compute_centroid(const Conformer &conf, bool ignore_hs,
                         bool has_weights, rust::Slice<const double> weights);
// The 4 x 4 transform, row-major.
rust::Vec<double> compute_canonical_transform(const Conformer &conf,
                                              bool has_center,
                                              const Point3D &center,
                                              bool normalize_covar,
                                              bool ignore_hs);
// False when the eigensolver fails; the axes are row-major.
bool compute_principal_axes_and_moments(const Conformer &conf,
                                        bool from_gyration_matrix,
                                        bool ignore_hs, bool has_weights,
                                        rust::Slice<const double> weights,
                                        rust::Vec<double> &axes,
                                        rust::Vec<double> &moments);
void transform_conformer(Conformer &conf, rust::Slice<const double> trans);
void canonicalize_conformer(Conformer &conf, bool has_center,
                            const Point3D &center, bool normalize_covar,
                            bool ignore_hs);
void canonicalize_mol(RWMol &mol, bool normalize_covar, bool ignore_hs);
double get_bond_length(const Conformer &conf, std::uint32_t i_atom_id,
                       std::uint32_t j_atom_id);
void set_bond_length(Conformer &conf, std::uint32_t i_atom_id,
                     std::uint32_t j_atom_id, double value);
double get_angle(const Conformer &conf, std::uint32_t i_atom_id,
                 std::uint32_t j_atom_id, std::uint32_t k_atom_id,
                 bool degrees);
void set_angle(Conformer &conf, std::uint32_t i_atom_id,
               std::uint32_t j_atom_id, std::uint32_t k_atom_id, double value,
               bool degrees);
double get_dihedral(const Conformer &conf, std::uint32_t i_atom_id,
                    std::uint32_t j_atom_id, std::uint32_t k_atom_id,
                    std::uint32_t l_atom_id, bool degrees);
void set_dihedral(Conformer &conf, std::uint32_t i_atom_id,
                  std::uint32_t j_atom_id, std::uint32_t k_atom_id,
                  std::uint32_t l_atom_id, double value, bool degrees);

}  // namespace rdkit_cxx
