// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/shape_helpers.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <utility>
#include <vector>

#include <Geometry/GridUtils.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ShapeHelpers/ShapeEncoder.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/coordinate_checks.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

// As _copyTransform in Code/GraphMol/ShapeHelpers/Wrap/rdShapeHelpers.cpp.
std::unique_ptr<RDGeom::Transform3D> to_transform(
    bool has_trans, rust::Slice<const double> trans) {
  if (!has_trans) {
    return nullptr;
  }
  if (trans.size() != 16) {
    throw ValueErrorException(
        "The transform has to be square matrix, of size 4x4");
  }
  auto transform = std::make_unique<RDGeom::Transform3D>();
  std::memcpy(transform->getData(), trans.data(), 16 * sizeof(double));
  return transform;
}

struct Sphere {
  RDGeom::Point3D center;
  double radius;
};

// The spheres EncodeShape (ShapeEncoder.cpp) sets: one per atom, hydrogens
// left out with ignore_hs, at the transformed position, with the scaled van
// der Waals radius.
std::vector<Sphere> spheres(const Conformer &conf,
                            const RDGeom::Transform3D *trans, double vdw_scale,
                            bool ignore_hs) {
  std::vector<Sphere> result;
  const auto &mol = conf.getOwningMol();
  for (const auto atom : mol.atoms()) {
    auto anum = atom->getAtomicNum();
    if (anum == 1 && ignore_hs) {
      continue;
    }
    auto rad = RDKit::PeriodicTable::getTable()->getRvdw(anum);
    RDGeom::Point3D loc = conf.getAtomPos(atom->getIdx());
    if (trans) {
      trans->TransformPoint(loc);
    }
    result.push_back(Sphere{loc, vdw_scale * rad});
  }
  return result;
}

// What createShapeGrids (ShapeUtils.cpp) builds before encoding both
// conformers: grids over the union of their boxes, in the first's canonical
// frame. initGrid casts the grid's size to unsigned int and
// setSphereOccupancy each sphere's extent to int, so those are checked
// first; every center lies inside the padded box.
void check_shape_grids(const Conformer &conf1, const Conformer &conf2,
                       double spacing, DiscreteValueType bits_per_point,
                       double vdw_scale, double step_size, bool ignore_hs) {
  check_finite_coordinates(conf1, "the shape grids");
  check_finite_coordinates(conf2, "the shape grids");
  std::unique_ptr<RDGeom::Transform3D> trans(
      MolTransforms::computeCanonicalTransform(conf1));
  RDGeom::Point3D lower1, upper1, lower2, upper2, lower, upper;
  RDKit::MolShapes::computeConfBox(conf1, lower1, upper1, trans.get());
  RDKit::MolShapes::computeConfBox(conf2, lower2, upper2, trans.get());
  RDKit::MolShapes::computeUnionBox(lower1, upper1, lower2, upper2, lower,
                                    upper);
  upper -= lower;
  check_grid_dimensions(upper.x, upper.y, upper.z, spacing, true,
                        Point3D{lower.x, lower.y, lower.z});
  if (!(upper.x > 0.0 && upper.y > 0.0 && upper.z > 0.0 && spacing > 0.0)) {
    return;  // initGrid's preconditions throw.
  }
  double max_num = 0.0;
  for (double dim : {upper.x, upper.y, upper.z}) {
    max_num = std::max(max_num, std::floor(dim / spacing + 0.5));
  }
  unsigned bits = 1u << static_cast<unsigned>(bits_per_point);
  for (const auto *conf : {&conf1, &conf2}) {
    for (const auto &sphere :
         spheres(*conf, trans.get(), vdw_scale, ignore_hs)) {
      check_sphere_extent(spacing, static_cast<unsigned>(max_num), bits,
                          sphere.radius, step_size);
    }
  }
}

}  // namespace

// As MolShapes::EncodeShape, but here, so that a center RDKit would cast to
// int out of range is skipped as RDKit skips any center outside the grid.
// Every sphere is checked before the grid changes.
void encode_shape(const RWMol &mol, UniformGrid3D &grid, std::int32_t conf_id,
                  bool has_trans, rust::Slice<const double> trans,
                  double vdw_scale, double step_size, std::int32_t max_layers,
                  bool ignore_hs) {
  auto transform = to_transform(has_trans, trans);
  const auto &conf = mol.getConformer(conf_id);
  std::vector<Sphere> inside;
  for (const auto &sphere :
       spheres(conf, transform.get(), vdw_scale, ignore_hs)) {
    if (grid_contains(grid, sphere.center)) {
      check_sphere_extent(grid, sphere.radius, step_size);
      inside.push_back(sphere);
    }
  }
  for (const auto &sphere : inside) {
    grid.setSphereOccupancy(sphere.center, sphere.radius, step_size,
                            max_layers);
  }
}

double shape_compare(const RWMol &mol1, const RWMol &mol2, std::uint8_t kind,
                     double alpha, double beta, std::int32_t conf_id1,
                     std::int32_t conf_id2, double grid_spacing,
                     DiscreteValueType bits_per_point, double vdw_scale,
                     double step_size, std::int32_t max_layers, bool ignore_hs,
                     bool allow_reordering) {
  const auto &conf1 = mol1.getConformer(conf_id1);
  const auto &conf2 = mol2.getConformer(conf_id2);
  // protrudeDistance replaces the spacing with 0.5 (ShapeUtils.cpp).
  double spacing = kind == 'p' ? 0.5 : grid_spacing;
  check_shape_grids(conf1, conf2, spacing, bits_per_point, vdw_scale, step_size,
                    ignore_hs);
  switch (kind) {
    case 't':
      return RDKit::MolShapes::tverskyIndex(
          mol1, mol2, alpha, beta, conf_id1, conf_id2, grid_spacing,
          bits_per_point, vdw_scale, step_size, max_layers, ignore_hs);
    case 'd':
      return RDKit::MolShapes::tanimotoDistance(
          mol1, mol2, conf_id1, conf_id2, grid_spacing, bits_per_point,
          vdw_scale, step_size, max_layers, ignore_hs);
    case 'p':
      return RDKit::MolShapes::protrudeDistance(
          mol1, mol2, conf_id1, conf_id2, grid_spacing, bits_per_point,
          vdw_scale, step_size, max_layers, ignore_hs, allow_reordering);
  }
  throw ValueErrorException("unknown shape comparison");
}

rust::Vec<Point3D> compute_conf_box(const Conformer &conf, bool dims,
                                    bool has_trans,
                                    rust::Slice<const double> trans,
                                    double padding) {
  auto transform = to_transform(has_trans, trans);
  RDGeom::Point3D first, second;
  if (dims) {
    RDKit::MolShapes::computeConfDimsAndOffset(conf, first, second,
                                               transform.get(), padding);
  } else {
    RDKit::MolShapes::computeConfBox(conf, first, second, transform.get(),
                                     padding);
  }
  rust::Vec<Point3D> result;
  result.push_back(Point3D{first.x, first.y, first.z});
  result.push_back(Point3D{second.x, second.y, second.z});
  return result;
}

rust::Vec<Point3D> compute_union_box(const Point3D &lower1,
                                     const Point3D &upper1,
                                     const Point3D &lower2,
                                     const Point3D &upper2) {
  auto rd = [](const Point3D &p) { return RDGeom::Point3D(p.x, p.y, p.z); };
  RDGeom::Point3D lower, upper;
  RDKit::MolShapes::computeUnionBox(rd(lower1), rd(upper1), rd(lower2),
                                    rd(upper2), lower, upper);
  rust::Vec<Point3D> result;
  result.push_back(Point3D{lower.x, lower.y, lower.z});
  result.push_back(Point3D{upper.x, upper.y, upper.z});
  return result;
}

}  // namespace rdkit_cxx
