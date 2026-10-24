/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Uniform 3D grids. Grid indices are checked (CONVENTIONS.md E3), and so are
// the parameters the in-place operators need to match. RDKit casts points,
// sizes, and radii in grid units to int without range checks, so the wrapper
// treats far-away points as outside the grid and rejects grids, spheres, and
// windows too large to index. A lock guards the global cache that the
// centroid functions write.

#pragma once

#include <cstdint>
#include <memory>

#include <Geometry/GridUtils.h>
#include <Geometry/UniformGrid3D.h>
#include <Geometry/UniformRealValueGrid3D.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/geometry.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using UniformGrid3D = RDGeom::UniformGrid3D;
using UniformRealValueGrid3D = RDGeom::UniformRealValueGrid3D;

// A shared struct, defined by the bridge (src/ffi/grid.rs).
struct GridIndices;

std::unique_ptr<UniformGrid3D> uniform_grid3d_new(double dim_x, double dim_y,
                                                  double dim_z, double spacing,
                                                  DiscreteValueType val_type,
                                                  bool has_offset,
                                                  const Point3D &offset);
std::unique_ptr<UniformGrid3D> uniform_grid3d_from_pkl(
    rust::Slice<const std::uint8_t> pkl);
std::unique_ptr<UniformGrid3D> uniform_grid3d_clone(const UniformGrid3D &grid);
rust::Vec<std::uint8_t> uniform_grid3d_to_binary(const UniformGrid3D &grid);
std::int32_t uniform_grid3d_get_grid_point_index(const UniformGrid3D &grid,
                                                 const Point3D &point);
std::int32_t uniform_grid3d_get_grid_index(const UniformGrid3D &grid,
                                           std::uint32_t xi, std::uint32_t yi,
                                           std::uint32_t zi);
GridIndices uniform_grid3d_get_grid_indices(const UniformGrid3D &grid,
                                            std::uint32_t idx);
std::int32_t uniform_grid3d_get_val_point(const UniformGrid3D &grid,
                                          const Point3D &pt);
std::uint32_t uniform_grid3d_get_val(const UniformGrid3D &grid,
                                     std::uint32_t id);
void uniform_grid3d_set_val(UniformGrid3D &grid, std::uint32_t id,
                            std::uint32_t val);
void uniform_grid3d_set_val_point(UniformGrid3D &grid, const Point3D &pt,
                                  std::uint32_t val);
Point3D uniform_grid3d_get_grid_point_loc(const UniformGrid3D &grid,
                                          std::uint32_t point_id);
std::uint32_t uniform_grid3d_get_size(const UniformGrid3D &grid);
std::uint32_t uniform_grid3d_get_num_x(const UniformGrid3D &grid);
std::uint32_t uniform_grid3d_get_num_y(const UniformGrid3D &grid);
std::uint32_t uniform_grid3d_get_num_z(const UniformGrid3D &grid);
Point3D uniform_grid3d_get_offset(const UniformGrid3D &grid);
double uniform_grid3d_get_spacing(const UniformGrid3D &grid);
std::unique_ptr<DiscreteValueVect> uniform_grid3d_get_occupancy_vect(
    const UniformGrid3D &grid);
bool uniform_grid3d_compare_params(const UniformGrid3D &grid,
                                   const UniformGrid3D &other);
void uniform_grid3d_set_sphere_occupancy(UniformGrid3D &grid,
                                         const Point3D &center, double radius,
                                         double step_size, int max_layers,
                                         bool ignore_out_of_bound);
// `op` is '+', '-', '&', or '|'.
void uniform_grid3d_assign_op(UniformGrid3D &grid, const UniformGrid3D &other,
                              std::uint8_t op);
void write_grid_to_file(const UniformGrid3D &grid,
                        rust::Slice<const std::uint8_t> filename);
double tversky_index(const UniformGrid3D &grid1, const UniformGrid3D &grid2,
                     double alpha, double beta);
double tanimoto_distance(const UniformGrid3D &grid1,
                         const UniformGrid3D &grid2);
double protrude_distance(const UniformGrid3D &grid1,
                         const UniformGrid3D &grid2);
// The centroid, with the weight sum in `weight_sum`.
Point3D compute_grid_centroid(const UniformGrid3D &grid, const Point3D &pt,
                              double window_radius, double &weight_sum);
rust::Vec<Point3D> find_grid_terminal_points(const UniformGrid3D &grid,
                                             double window_radius,
                                             double inclusion_fraction);

std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_default();
std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_new(
    double dim_x, double dim_y, double dim_z, double spacing, bool has_offset,
    const Point3D &offset);
std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_from_pkl(
    rust::Slice<const std::uint8_t> pkl);
std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_clone(
    const UniformRealValueGrid3D &grid);
rust::Vec<std::uint8_t> uniform_real_value_grid3d_to_binary(
    const UniformRealValueGrid3D &grid);
std::int32_t uniform_real_value_grid3d_get_grid_point_index(
    const UniformRealValueGrid3D &grid, const Point3D &point);
std::int32_t uniform_real_value_grid3d_get_grid_index(
    const UniformRealValueGrid3D &grid, std::uint32_t xi, std::uint32_t yi,
    std::uint32_t zi);
GridIndices uniform_real_value_grid3d_get_grid_indices(
    const UniformRealValueGrid3D &grid, std::uint32_t idx);
double uniform_real_value_grid3d_get_val_point(
    const UniformRealValueGrid3D &grid, const Point3D &pt);
double uniform_real_value_grid3d_get_val(const UniformRealValueGrid3D &grid,
                                         std::uint32_t id);
void uniform_real_value_grid3d_set_val(UniformRealValueGrid3D &grid,
                                       std::uint32_t id, double val);
void uniform_real_value_grid3d_set_val_point(UniformRealValueGrid3D &grid,
                                             const Point3D &pt, double val);
Point3D uniform_real_value_grid3d_get_grid_point_loc(
    const UniformRealValueGrid3D &grid, std::uint32_t point_id);
std::uint32_t uniform_real_value_grid3d_get_size(
    const UniformRealValueGrid3D &grid);
std::uint32_t uniform_real_value_grid3d_get_num_x(
    const UniformRealValueGrid3D &grid);
std::uint32_t uniform_real_value_grid3d_get_num_y(
    const UniformRealValueGrid3D &grid);
std::uint32_t uniform_real_value_grid3d_get_num_z(
    const UniformRealValueGrid3D &grid);
Point3D uniform_real_value_grid3d_get_offset(
    const UniformRealValueGrid3D &grid);
double uniform_real_value_grid3d_get_spacing(
    const UniformRealValueGrid3D &grid);
std::unique_ptr<RealValueVect> uniform_real_value_grid3d_get_occupancy_vect(
    const UniformRealValueGrid3D &grid);
// `what` is 'v' (vectors), 'p' (parameters), or 'g' (both).
bool uniform_real_value_grid3d_compare(const UniformRealValueGrid3D &grid,
                                       const UniformRealValueGrid3D &other,
                                       std::uint8_t what);
void uniform_real_value_grid3d_assign_op(UniformRealValueGrid3D &grid,
                                         const UniformRealValueGrid3D &other,
                                         std::uint8_t op);

}  // namespace rdkit_cxx
