// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/geometry.h"

#include <Geometry/point.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

RDGeom::Point2D rd(const Point2D &p) { return RDGeom::Point2D(p.x, p.y); }
RDGeom::Point3D rd(const Point3D &p) { return RDGeom::Point3D(p.x, p.y, p.z); }
Point2D ours(const RDGeom::Point2D &p) { return Point2D{p.x, p.y}; }
Point3D ours(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

RDGeom::PointND rd(rust::Slice<const double> p) {
  RDGeom::PointND res(p.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    res[i] = p[i];
  }
  return res;
}

rust::Vec<double> ours(const RDGeom::PointND &p) {
  rust::Vec<double> res;
  for (unsigned i = 0; i < p.dimension(); ++i) {
    res.push_back(p[i]);
  }
  return res;
}

void check_same_dimension(rust::Slice<const double> p,
                          rust::Slice<const double> other) {
  if (p.size() != other.size()) {
    throw ValueErrorException("PointND dimension mismatch");
  }
}

}  // namespace

rust::Vec<IndexList> to_index_lists(
    const std::vector<std::vector<int>> &lists) {
  rust::Vec<IndexList> result;
  for (const auto &indices : lists) {
    IndexList list;
    for (auto idx : indices) {
      list.indices.push_back(static_cast<std::uint32_t>(idx));
    }
    result.push_back(std::move(list));
  }
  return result;
}

Point2D point2d_normalize(const Point2D &p) {
  auto res = rd(p);
  res.normalize();
  return ours(res);
}

double point2d_length(const Point2D &p) { return rd(p).length(); }

double point2d_length_sq(const Point2D &p) { return rd(p).lengthSq(); }

double point2d_dot_product(const Point2D &p, const Point2D &other) {
  return rd(p).dotProduct(rd(other));
}

double point2d_angle_to(const Point2D &p, const Point2D &other) {
  return rd(p).angleTo(rd(other));
}

double point2d_signed_angle_to(const Point2D &p, const Point2D &other) {
  return rd(p).signedAngleTo(rd(other));
}

Point2D point2d_direction_vector(const Point2D &p, const Point2D &other) {
  return ours(rd(p).directionVector(rd(other)));
}

Point3D point3d_normalize(const Point3D &p) {
  auto res = rd(p);
  res.normalize();
  return ours(res);
}

double point3d_length(const Point3D &p) { return rd(p).length(); }

double point3d_length_sq(const Point3D &p) { return rd(p).lengthSq(); }

double point3d_dot_product(const Point3D &p, const Point3D &other) {
  return rd(p).dotProduct(rd(other));
}

Point3D point3d_cross_product(const Point3D &p, const Point3D &other) {
  return ours(rd(p).crossProduct(rd(other)));
}

double point3d_angle_to(const Point3D &p, const Point3D &other) {
  return rd(p).angleTo(rd(other));
}

double point3d_signed_angle_to(const Point3D &p, const Point3D &other) {
  return rd(p).signedAngleTo(rd(other));
}

Point3D point3d_direction_vector(const Point3D &p, const Point3D &other) {
  return ours(rd(p).directionVector(rd(other)));
}

double compute_dihedral_angle(const Point3D &pt1, const Point3D &pt2,
                              const Point3D &pt3, const Point3D &pt4) {
  return RDGeom::computeDihedralAngle(rd(pt1), rd(pt2), rd(pt3), rd(pt4));
}

double compute_signed_dihedral_angle(const Point3D &pt1, const Point3D &pt2,
                                     const Point3D &pt3, const Point3D &pt4) {
  return RDGeom::computeSignedDihedralAngle(rd(pt1), rd(pt2), rd(pt3), rd(pt4));
}

rust::Vec<double> point_nd_normalize(rust::Slice<const double> p) {
  auto res = rd(p);
  res.normalize();
  return ours(res);
}

double point_nd_length(rust::Slice<const double> p) { return rd(p).length(); }

double point_nd_length_sq(rust::Slice<const double> p) {
  return rd(p).lengthSq();
}

double point_nd_dot_product(rust::Slice<const double> p,
                            rust::Slice<const double> other) {
  check_same_dimension(p, other);
  return rd(p).dotProduct(rd(other));
}

double point_nd_angle_to(rust::Slice<const double> p,
                         rust::Slice<const double> other) {
  check_same_dimension(p, other);
  return rd(p).angleTo(rd(other));
}

rust::Vec<double> point_nd_direction_vector(rust::Slice<const double> p,
                                            rust::Slice<const double> other) {
  check_same_dimension(p, other);
  return ours(rd(p).directionVector(rd(other)));
}

}  // namespace rdkit_cxx
