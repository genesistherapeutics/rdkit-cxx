// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/distance_geometry.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/ChiralSet.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <ForceField/ForceField.h>
#include <Numerics/SymmMatrix.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/embedder.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/distance_geometry.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/force_field.h"
#include "rdkit-cxx/global_state.h"

namespace rdkit_cxx {

namespace {

// As the checks in Code/DistGeom/Wrap/DistGeom.cpp, whose nrows * nrows
// overflows int above 46340 rows in doTriangleSmoothing and wraps in unsigned
// int above 65535 in embedBoundsMatrix; matrices are capped where n * n fits
// in int (U3).
DistGeom::BoundsMatrix::DATA_SPTR copy_bounds(rust::Slice<const double> bounds,
                                              std::uint32_t n) {
  if (n == 0) {
    throw ValueErrorException("The array has to have a nonzero size");
  }
  if (n > 46340 || bounds.size() != static_cast<std::size_t>(n) * n) {
    throw ValueErrorException("The array has to be square");
  }
  DistGeom::BoundsMatrix::DATA_SPTR data(new double[bounds.size()]);
  for (std::size_t i = 0; i < bounds.size(); ++i) {
    data[i] = bounds[i];
  }
  return data;
}

}  // namespace

// As doTriangleSmoothing in DistGeom.cpp.
bool do_triangle_smoothing(rust::Vec<double> &bounds, std::uint32_t n,
                           double tol) {
  auto data = copy_bounds({bounds.data(), bounds.size()}, n);
  DistGeom::BoundsMatrix matrix(n, data);
  bool res = DistGeom::triangleSmoothBounds(&matrix, tol);
  for (std::size_t i = 0; i < bounds.size(); ++i) {
    bounds[i] = data[i];
  }
  return res;
}

// As embedBoundsMatrix in DistGeom.cpp, which draws from RDKit's global
// random generator, seeded or not (random_mutex); adds iter * 999 to the
// seed in int, which overflows, where this wraps; and minimizes in three
// dimensions. Its starting coordinates cast the summed squares of
// distances between the bounds to int (check_initial_coords_cast), and an
// infinite bound becomes an infinite distance.
rust::Vec<Point3D> embed_bounds_matrix(
    rust::Slice<const double> bounds, std::uint32_t n, std::int32_t max_iters,
    bool randomize_on_failure, std::int32_t num_zero_fail,
    rust::Slice<const IndexPair> weight_pairs,
    rust::Slice<const double> weights, std::int32_t random_seed) {
  auto data = copy_bounds(bounds, n);
  double max_bound = 0.0;
  for (std::size_t i = 0; i < bounds.size(); ++i) {
    if (!std::isfinite(bounds[i])) {
      throw ValueErrorException("the bounds matrix has non-finite values");
    }
    max_bound = std::max(max_bound, std::abs(bounds[i]));
  }
  check_initial_coords_cast(n, max_bound);
  check_minimize_size(n, 3);
  DistGeom::BoundsMatrix matrix(n, data);

  std::vector<RDGeom::Point3D> positions(n);
  RDGeom::PointPtrVect pointers;
  for (auto &position : positions) {
    pointers.push_back(&position);
  }
  RDNumeric::DoubleSymmMatrix dist_mat(n, 0.0);

  std::lock_guard<std::mutex> lock(random_mutex());
  bool got_coords = false;
  std::int32_t seed = random_seed;
  for (std::int32_t iter = 0; iter < max_iters && !got_coords; ++iter) {
    DistGeom::pickRandomDistMat(matrix, dist_mat, seed);
    got_coords = DistGeom::computeInitialCoords(
        dist_mat, pointers, randomize_on_failure,
        static_cast<unsigned>(num_zero_fail), seed);
    if (seed >= 0) {
      seed = static_cast<std::int32_t>(static_cast<std::int64_t>(seed) +
                                       static_cast<std::int64_t>(iter) * 999);
    }
  }
  if (!got_coords) {
    throw ValueErrorException("could not embed matrix");
  }
  std::map<std::pair<int, int>, double> weight_map;
  for (std::size_t i = 0; i < weight_pairs.size(); ++i) {
    weight_map[{static_cast<int>(weight_pairs[i].first),
                static_cast<int>(weight_pairs[i].second)}] = weights[i];
  }
  DistGeom::VECT_CHIRALSET csets;
  std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
      matrix, pointers, csets, 0.0, 0.0, &weight_map));
  field->initialize();
  if (field->calcEnergy() > 1e-5) {
    int need_more = 1;
    while (need_more) {
      need_more = field->minimize();
    }
  }
  rust::Vec<Point3D> result;
  for (const auto &position : positions) {
    result.push_back(Point3D{position.x, position.y, position.z});
  }
  return result;
}

}  // namespace rdkit_cxx
