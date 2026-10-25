/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// DistanceGeometry.DistGeom: triangle smoothing and embedding of bounds
// matrices.

#pragma once

#include <cstdint>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/distance_geometry.rs,
// geometry.rs, chirality.rs).
struct Point3D;
struct IndexPair;

// Smooths the n * n row-major bounds in place.
bool do_triangle_smoothing(rust::Vec<double> &bounds, std::uint32_t n,
                           double tol);
rust::Vec<Point3D> embed_bounds_matrix(
    rust::Slice<const double> bounds, std::uint32_t n, std::int32_t max_iters,
    bool randomize_on_failure, std::int32_t num_zero_fail,
    rust::Slice<const IndexPair> weight_pairs,
    rust::Slice<const double> weights, std::int32_t random_seed);

}  // namespace rdkit_cxx
