/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Value types shared across areas. The bridge (src/ffi/geometry.rs) defines
// them; wrappers that use them include "rdkit-cxx/ffi/geometry.rs.h".

#pragma once

#include <vector>

#include "rust/cxx.h"

namespace rdkit_cxx {

// rdGeometry.Point2D and Point3D as plain data. RDGeom's points have vtables,
// so they cannot cross the bridge themselves.
struct Point2D;
struct Point3D;

// A list of atom or bond indices. cxx has no Vec<Vec<T>>, so nested lists
// cross as Vec<IndexList>.
struct IndexList;

// Nested indices as IndexLists, for the wrappers that return them. RDKit's
// indices are never negative.
rust::Vec<IndexList> to_index_lists(const std::vector<std::vector<int>> &lists);

// As the methods of RDGeom::Point2D, Point3D, and PointND (Geometry/point.h).
// PointND crosses as its coordinates.
Point2D point2d_normalize(const Point2D &p);
double point2d_length(const Point2D &p);
double point2d_length_sq(const Point2D &p);
double point2d_dot_product(const Point2D &p, const Point2D &other);
double point2d_angle_to(const Point2D &p, const Point2D &other);
double point2d_signed_angle_to(const Point2D &p, const Point2D &other);
Point2D point2d_direction_vector(const Point2D &p, const Point2D &other);

Point3D point3d_normalize(const Point3D &p);
double point3d_length(const Point3D &p);
double point3d_length_sq(const Point3D &p);
double point3d_dot_product(const Point3D &p, const Point3D &other);
Point3D point3d_cross_product(const Point3D &p, const Point3D &other);
double point3d_angle_to(const Point3D &p, const Point3D &other);
double point3d_signed_angle_to(const Point3D &p, const Point3D &other);
Point3D point3d_direction_vector(const Point3D &p, const Point3D &other);
double compute_dihedral_angle(const Point3D &pt1, const Point3D &pt2,
                              const Point3D &pt3, const Point3D &pt4);
double compute_signed_dihedral_angle(const Point3D &pt1, const Point3D &pt2,
                                     const Point3D &pt3, const Point3D &pt4);

rust::Vec<double> point_nd_normalize(rust::Slice<const double> p);
double point_nd_length(rust::Slice<const double> p);
double point_nd_length_sq(rust::Slice<const double> p);
double point_nd_dot_product(rust::Slice<const double> p,
                            rust::Slice<const double> other);
double point_nd_angle_to(rust::Slice<const double> p,
                         rust::Slice<const double> other);
rust::Vec<double> point_nd_direction_vector(rust::Slice<const double> p,
                                            rust::Slice<const double> other);

}  // namespace rdkit_cxx
