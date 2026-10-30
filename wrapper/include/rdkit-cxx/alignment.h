/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Numerics.rdAlignment: the optimal superposition of two sets of points.

#pragma once

#include <cstdint>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/geometry.rs).
struct Point3D;

// The sum of squared deviations; the 4 x 4 transform goes to transform,
// row-major.
double align_point_pairs(rust::Slice<const Point3D> ref_points,
                         rust::Slice<const Point3D> probe_points,
                         bool has_weights, rust::Slice<const double> weights,
                         bool reflect, std::uint32_t max_iterations,
                         rust::Vec<double> &transform);

}  // namespace rdkit_cxx
