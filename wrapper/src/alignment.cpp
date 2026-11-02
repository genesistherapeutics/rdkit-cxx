// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/alignment.h"

#include <memory>
#include <vector>

#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

// As AlignPointPairs in Code/Numerics/Alignment/Wrap/rdAlignment.cpp; RDKit
// checks that the weights are positive. Python treats an empty weights list
// as no weights and reports empty points as "Empty sequence passed in"; here
// empty weights must match the points like any others.
double align_point_pairs(rust::Slice<const Point3D> ref_points,
                         rust::Slice<const Point3D> probe_points,
                         bool has_weights, rust::Slice<const double> weights,
                         bool reflect, std::uint32_t max_iterations,
                         rust::Vec<double> &transform) {
  if (ref_points.empty() || probe_points.empty()) {
    throw ValueErrorException("expected a non-empty sequence of points");
  }
  if (ref_points.size() != probe_points.size()) {
    throw ValueErrorException("Mis-match in number of points");
  }
  std::vector<RDGeom::Point3D> ref_owned, probe_owned;
  for (const auto &p : ref_points) {
    ref_owned.emplace_back(p.x, p.y, p.z);
  }
  for (const auto &p : probe_points) {
    probe_owned.emplace_back(p.x, p.y, p.z);
  }
  RDGeom::Point3DConstPtrVect ref_ptrs, probe_ptrs;
  for (const auto &p : ref_owned) {
    ref_ptrs.push_back(&p);
  }
  for (const auto &p : probe_owned) {
    probe_ptrs.push_back(&p);
  }
  std::unique_ptr<RDNumeric::DoubleVector> weights_vec;
  if (has_weights) {
    if (weights.size() != ref_owned.size()) {
      throw ValueErrorException(
          "Number of weights supplied do not match the number of points");
    }
    weights_vec = std::make_unique<RDNumeric::DoubleVector>(weights.size());
    for (std::size_t i = 0; i < weights.size(); ++i) {
      weights_vec->setVal(i, weights[i]);
    }
  }
  RDGeom::Transform3D trans;
  double ssd = RDNumeric::Alignments::AlignPoints(
      ref_ptrs, probe_ptrs, trans, weights_vec.get(), reflect, max_iterations);
  transform.clear();
  const double *data = trans.getData();
  for (unsigned i = 0; i < 16; ++i) {
    transform.push_back(data[i]);
  }
  return ssd;
}

}  // namespace rdkit_cxx
