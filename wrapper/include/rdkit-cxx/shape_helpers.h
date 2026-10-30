/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for GraphMol/ShapeHelpers (rdkit.Chem.rdShapeHelpers): shapes of
// conformers encoded on grids, and their boxes.

#pragma once

#include <cstdint>

#include "rdkit-cxx/conformer.h"
#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/grid.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A shared struct, defined by the bridges (src/ffi/geometry.rs).
struct Point3D;

void encode_shape(const RWMol &mol, UniformGrid3D &grid, std::int32_t conf_id,
                  bool has_trans, rust::Slice<const double> trans,
                  double vdw_scale, double step_size, std::int32_t max_layers,
                  bool ignore_hs);
// kind: 't' Tversky, 'd' Tanimoto distance, 'p' protrude distance.
double shape_compare(const RWMol &mol1, const RWMol &mol2, std::uint8_t kind,
                     double alpha, double beta, std::int32_t conf_id1,
                     std::int32_t conf_id2, double grid_spacing,
                     DiscreteValueType bits_per_point, double vdw_scale,
                     double step_size, std::int32_t max_layers, bool ignore_hs,
                     bool allow_reordering);
// The lower and upper corners, or with `dims` the dimensions and offset.
rust::Vec<Point3D> compute_conf_box(const Conformer &conf, bool dims,
                                    bool has_trans,
                                    rust::Slice<const double> trans,
                                    double padding);
rust::Vec<Point3D> compute_union_box(const Point3D &lower1,
                                     const Point3D &upper1,
                                     const Point3D &lower2,
                                     const Point3D &upper2);

}  // namespace rdkit_cxx
