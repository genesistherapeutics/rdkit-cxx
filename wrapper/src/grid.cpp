/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/grid.h"

#include <algorithm>
#include <cmath>
#include <mutex>

#include <Geometry/Grid3D.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/grid.rs.h"
#include "rdkit-cxx/pickle_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

RDGeom::Point3D rd(const Point3D &p) { return RDGeom::Point3D(p.x, p.y, p.z); }
Point3D ours(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

// RDKit casts grid coordinates to int and dimensions over the spacing to
// unsigned int without range checks, which is UB for NaN or a value out of
// range. The checks below keep every such cast in range, and keep grid point
// indices, which RDKit returns as int, below INT_MAX.
constexpr double kMaxInt = 2147483647.0;

// Whether `point` rounds to a grid point, computed as getGridPointIndex
// computes it but in double, so a point far outside the grid is not cast.
template <class Grid>
bool inside(const Grid &grid, const RDGeom::Point3D &point) {
  auto t = point;
  t -= grid.getOffset();
  t /= grid.getSpacing();
  auto in = [](double c, unsigned n) {
    return c + 0.5 >= 0.0 && c + 0.5 < static_cast<double>(n);
  };
  return in(t.x, grid.getNumX()) && in(t.y, grid.getNumY()) &&
         in(t.z, grid.getNumZ());
}

// Checks a new grid before initGrid casts its size. RDKit's own
// preconditions report nonpositive (and NaN) dimensions and spacings. A grid
// with no points along an axis can't be copied, since the copy's dimension
// is then zero, so it is rejected too.
void check_dimensions(double dim_x, double dim_y, double dim_z, double spacing,
                      bool has_offset, const Point3D &offset) {
  if (!(dim_x > 0.0 && dim_y > 0.0 && dim_z > 0.0 && spacing > 0.0)) {
    return;
  }
  if (!std::isfinite(spacing)) {
    throw ValueErrorException("the grid spacing must be finite");
  }
  if (has_offset && !(std::isfinite(offset.x) && std::isfinite(offset.y) &&
                      std::isfinite(offset.z))) {
    throw ValueErrorException("the grid offset must be finite");
  }
  double points = 1.0;
  for (double dim : {dim_x, dim_y, dim_z}) {
    double n = std::floor(dim / spacing + 0.5);
    if (!(n <= kMaxInt)) {
      throw ValueErrorException("the grid has too many points");
    }
    if (n < 1.0) {
      throw ValueErrorException("the grid has no points along an axis");
    }
    points *= n;
  }
  if (points > kMaxInt) {
    throw ValueErrorException("the grid has too many points");
  }
}

// computeGridIndices casts the window radius in grid points to int and sums
// the squares of three such ints, so it must stay below sqrt(INT_MAX / 3).
void check_window_radius(const UniformGrid3D &grid, double window_radius) {
  double radius = std::ceil(window_radius / grid.getSpacing());
  if (!(std::abs(radius) <= 26754.0)) {
    throw ValueErrorException("the window radius is too large for the grid");
  }
}

// computeGridIndices writes a global cache with no lock.
std::mutex grid_indices_mutex;

template <class Grid>
GridIndices indices(const Grid &grid, std::uint32_t idx) {
  check_index(idx, grid.getSize());
  unsigned xi, yi, zi;
  grid.getGridIndices(idx, xi, yi, zi);
  return GridIndices{xi, yi, zi};
}

template <class Grid>
void assign_op(Grid &grid, const Grid &other, std::uint8_t op) {
  if (!grid.compareParams(other)) {
    throw ValueErrorException("the grids' parameters differ");
  }
  switch (op) {
    case '+':
      grid += other;
      return;
    case '-':
      grid -= other;
      return;
    case '&':
      grid &= other;
      return;
    case '|':
      grid |= other;
      return;
  }
  throw ValueErrorException("unsupported grid operation");
}

}  // namespace

std::unique_ptr<UniformGrid3D> uniform_grid3d_new(double dim_x, double dim_y,
                                                  double dim_z, double spacing,
                                                  DiscreteValueType val_type,
                                                  bool has_offset,
                                                  const Point3D &offset) {
  check_dimensions(dim_x, dim_y, dim_z, spacing, has_offset, offset);
  auto off = rd(offset);
  return std::make_unique<UniformGrid3D>(dim_x, dim_y, dim_z, spacing, val_type,
                                         has_offset ? &off : nullptr);
}

std::unique_ptr<UniformGrid3D> uniform_grid3d_from_pkl(
    rust::Slice<const std::uint8_t> pkl) {
  auto text = to_std(pkl);
  check_uniform_grid3d_pickle(text);
  return std::make_unique<UniformGrid3D>(text);
}

std::unique_ptr<UniformGrid3D> uniform_grid3d_clone(const UniformGrid3D &grid) {
  return std::make_unique<UniformGrid3D>(grid);
}

rust::Vec<std::uint8_t> uniform_grid3d_to_binary(const UniformGrid3D &grid) {
  return to_bytes(grid.toString());
}

std::int32_t uniform_grid3d_get_grid_point_index(const UniformGrid3D &grid,
                                                 const Point3D &point) {
  return inside(grid, rd(point)) ? grid.getGridPointIndex(rd(point)) : -1;
}

std::int32_t uniform_grid3d_get_grid_index(const UniformGrid3D &grid,
                                           std::uint32_t xi, std::uint32_t yi,
                                           std::uint32_t zi) {
  return grid.getGridIndex(xi, yi, zi);
}

GridIndices uniform_grid3d_get_grid_indices(const UniformGrid3D &grid,
                                            std::uint32_t idx) {
  return indices(grid, idx);
}

std::int32_t uniform_grid3d_get_val_point(const UniformGrid3D &grid,
                                          const Point3D &pt) {
  return inside(grid, rd(pt)) ? grid.getVal(rd(pt)) : -1;
}

std::uint32_t uniform_grid3d_get_val(const UniformGrid3D &grid,
                                     std::uint32_t id) {
  check_index(id, grid.getSize());
  return grid.getVal(id);
}

void uniform_grid3d_set_val(UniformGrid3D &grid, std::uint32_t id,
                            std::uint32_t val) {
  check_index(id, grid.getSize());
  grid.setVal(id, val);
}

void uniform_grid3d_set_val_point(UniformGrid3D &grid, const Point3D &pt,
                                  std::uint32_t val) {
  if (inside(grid, rd(pt))) {
    grid.setVal(rd(pt), val);
  }
}

Point3D uniform_grid3d_get_grid_point_loc(const UniformGrid3D &grid,
                                          std::uint32_t point_id) {
  check_index(point_id, grid.getSize());
  return ours(grid.getGridPointLoc(point_id));
}

std::uint32_t uniform_grid3d_get_size(const UniformGrid3D &grid) {
  return grid.getSize();
}

std::uint32_t uniform_grid3d_get_num_x(const UniformGrid3D &grid) {
  return grid.getNumX();
}

std::uint32_t uniform_grid3d_get_num_y(const UniformGrid3D &grid) {
  return grid.getNumY();
}

std::uint32_t uniform_grid3d_get_num_z(const UniformGrid3D &grid) {
  return grid.getNumZ();
}

Point3D uniform_grid3d_get_offset(const UniformGrid3D &grid) {
  return ours(grid.getOffset());
}

double uniform_grid3d_get_spacing(const UniformGrid3D &grid) {
  return grid.getSpacing();
}

std::unique_ptr<DiscreteValueVect> uniform_grid3d_get_occupancy_vect(
    const UniformGrid3D &grid) {
  return std::make_unique<DiscreteValueVect>(*grid.getOccupancyVect());
}

bool uniform_grid3d_compare_params(const UniformGrid3D &grid,
                                   const UniformGrid3D &other) {
  return grid.compareParams(other);
}

void uniform_grid3d_set_sphere_occupancy(UniformGrid3D &grid,
                                         const Point3D &center, double radius,
                                         double step_size, int max_layers,
                                         bool ignore_out_of_bound) {
  if (!inside(grid, rd(center))) {
    if (ignore_out_of_bound) {
      return;
    }
    throw RDGeom::GridException("Center outside the grid boundary");
  }
  // RDKit casts the sphere's extent in grid points to int. A negative step
  // can make it cast a negative layer number to unsigned int.
  if (!(step_size >= 0.0)) {
    throw ValueErrorException("the step size must not be negative");
  }
  double max_val = (1u << grid.getOccupancyVect()->getNumBitsPerVal()) - 1;
  double extent = (std::abs(radius) + max_val * step_size) / grid.getSpacing();
  double max_num = std::max({grid.getNumX(), grid.getNumY(), grid.getNumZ()});
  if (!(extent + max_num < kMaxInt)) {
    throw ValueErrorException("the sphere is too large for the grid");
  }
  grid.setSphereOccupancy(rd(center), radius, step_size, max_layers,
                          ignore_out_of_bound);
}

void uniform_grid3d_assign_op(UniformGrid3D &grid, const UniformGrid3D &other,
                              std::uint8_t op) {
  assign_op(grid, other, op);
}

void write_grid_to_file(const UniformGrid3D &grid,
                        rust::Slice<const std::uint8_t> filename) {
  RDGeom::writeGridToFile(grid, to_std(filename));
}

double tversky_index(const UniformGrid3D &grid1, const UniformGrid3D &grid2,
                     double alpha, double beta) {
  return RDGeom::tverskyIndex(grid1, grid2, alpha, beta);
}

double tanimoto_distance(const UniformGrid3D &grid1,
                         const UniformGrid3D &grid2) {
  return RDGeom::tanimotoDistance(grid1, grid2);
}

double protrude_distance(const UniformGrid3D &grid1,
                         const UniformGrid3D &grid2) {
  return RDGeom::protrudeDistance(grid1, grid2);
}

Point3D compute_grid_centroid(const UniformGrid3D &grid, const Point3D &pt,
                              double window_radius, double &weight_sum) {
  check_window_radius(grid, window_radius);
  // RDKit uses `pt` only for its grid point index, so any point outside the
  // grid gives the same result.
  auto point = rd(pt);
  if (!inside(grid, point)) {
    point = grid.getOffset() - RDGeom::Point3D(grid.getSpacing(), 0.0, 0.0);
  }
  std::lock_guard<std::mutex> lock(grid_indices_mutex);
  return ours(
      RDGeom::computeGridCentroid(grid, point, window_radius, weight_sum));
}

rust::Vec<Point3D> find_grid_terminal_points(const UniformGrid3D &grid,
                                             double window_radius,
                                             double inclusion_fraction) {
  check_window_radius(grid, window_radius);
  std::vector<RDGeom::Point3D> points;
  {
    std::lock_guard<std::mutex> lock(grid_indices_mutex);
    points =
        RDGeom::findGridTerminalPoints(grid, window_radius, inclusion_fraction);
  }
  rust::Vec<Point3D> result;
  for (const auto &p : points) {
    result.push_back(ours(p));
  }
  return result;
}

std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_default() {
  return std::make_unique<UniformRealValueGrid3D>();
}

std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_new(
    double dim_x, double dim_y, double dim_z, double spacing, bool has_offset,
    const Point3D &offset) {
  check_dimensions(dim_x, dim_y, dim_z, spacing, has_offset, offset);
  auto off = rd(offset);
  return std::make_unique<UniformRealValueGrid3D>(dim_x, dim_y, dim_z, spacing,
                                                  has_offset ? &off : nullptr);
}

std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_from_pkl(
    rust::Slice<const std::uint8_t> pkl) {
  auto text = to_std(pkl);
  check_uniform_real_value_grid3d_pickle(text);
  return std::make_unique<UniformRealValueGrid3D>(text);
}

std::unique_ptr<UniformRealValueGrid3D> uniform_real_value_grid3d_clone(
    const UniformRealValueGrid3D &grid) {
  return std::make_unique<UniformRealValueGrid3D>(grid);
}

rust::Vec<std::uint8_t> uniform_real_value_grid3d_to_binary(
    const UniformRealValueGrid3D &grid) {
  return to_bytes(grid.toString());
}

std::int32_t uniform_real_value_grid3d_get_grid_point_index(
    const UniformRealValueGrid3D &grid, const Point3D &point) {
  return inside(grid, rd(point)) ? grid.getGridPointIndex(rd(point)) : -1;
}

std::int32_t uniform_real_value_grid3d_get_grid_index(
    const UniformRealValueGrid3D &grid, std::uint32_t xi, std::uint32_t yi,
    std::uint32_t zi) {
  return grid.getGridIndex(xi, yi, zi);
}

GridIndices uniform_real_value_grid3d_get_grid_indices(
    const UniformRealValueGrid3D &grid, std::uint32_t idx) {
  return indices(grid, idx);
}

double uniform_real_value_grid3d_get_val_point(
    const UniformRealValueGrid3D &grid, const Point3D &pt) {
  return inside(grid, rd(pt)) ? grid.getVal(rd(pt)) : -1.0;
}

double uniform_real_value_grid3d_get_val(const UniformRealValueGrid3D &grid,
                                         std::uint32_t id) {
  check_index(id, grid.getSize());
  return grid.getVal(id);
}

void uniform_real_value_grid3d_set_val(UniformRealValueGrid3D &grid,
                                       std::uint32_t id, double val) {
  check_index(id, grid.getSize());
  grid.setVal(id, val);
}

void uniform_real_value_grid3d_set_val_point(UniformRealValueGrid3D &grid,
                                             const Point3D &pt, double val) {
  if (inside(grid, rd(pt))) {
    grid.setVal(rd(pt), val);
  }
}

Point3D uniform_real_value_grid3d_get_grid_point_loc(
    const UniformRealValueGrid3D &grid, std::uint32_t point_id) {
  check_index(point_id, grid.getSize());
  return ours(grid.getGridPointLoc(point_id));
}

std::uint32_t uniform_real_value_grid3d_get_size(
    const UniformRealValueGrid3D &grid) {
  return grid.getSize();
}

std::uint32_t uniform_real_value_grid3d_get_num_x(
    const UniformRealValueGrid3D &grid) {
  return grid.getNumX();
}

std::uint32_t uniform_real_value_grid3d_get_num_y(
    const UniformRealValueGrid3D &grid) {
  return grid.getNumY();
}

std::uint32_t uniform_real_value_grid3d_get_num_z(
    const UniformRealValueGrid3D &grid) {
  return grid.getNumZ();
}

Point3D uniform_real_value_grid3d_get_offset(
    const UniformRealValueGrid3D &grid) {
  return ours(grid.getOffset());
}

double uniform_real_value_grid3d_get_spacing(
    const UniformRealValueGrid3D &grid) {
  return grid.getSpacing();
}

std::unique_ptr<RealValueVect> uniform_real_value_grid3d_get_occupancy_vect(
    const UniformRealValueGrid3D &grid) {
  return std::make_unique<RealValueVect>(*grid.getOccupancyVect());
}

bool uniform_real_value_grid3d_compare(const UniformRealValueGrid3D &grid,
                                       const UniformRealValueGrid3D &other,
                                       std::uint8_t what) {
  switch (what) {
    case 'v':
      return grid.compareVectors(other);
    case 'p':
      return grid.compareParams(other);
    case 'g':
      return grid.compareGrids(other);
  }
  throw ValueErrorException("unsupported grid comparison");
}

void uniform_real_value_grid3d_assign_op(UniformRealValueGrid3D &grid,
                                         const UniformRealValueGrid3D &other,
                                         std::uint8_t op) {
  assign_op(grid, other, op);
}

}  // namespace rdkit_cxx
