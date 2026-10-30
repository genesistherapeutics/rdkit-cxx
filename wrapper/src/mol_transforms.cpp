// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_transforms.h"

#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

#include <Geometry/Transform3D.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/coordinate_checks.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/mol_transforms.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

std::unique_ptr<RDGeom::Point3D> to_center(bool has_center,
                                           const Point3D &center) {
  if (!has_center) {
    return nullptr;
  }
  return std::make_unique<RDGeom::Point3D>(center.x, center.y, center.z);
}

rust::Vec<double> to_values(const RDGeom::Transform3D &trans) {
  rust::Vec<double> values;
  const double *data = trans.getData();
  for (unsigned i = 0; i < 16; ++i) {
    values.push_back(data[i]);
  }
  return values;
}

// The functions read and move atoms by index, through RDKit's range checks;
// the wrapper's come first (E3).
void check_atoms(const Conformer &conf,
                 std::initializer_list<std::uint32_t> indices) {
  for (auto idx : indices) {
    check_index(idx, conf.getNumAtoms());
  }
}

}  // namespace

// As ComputeCentroid in Code/GraphMol/MolTransforms/Wrap/rdMolTransforms.cpp,
// which binds computeCentroid directly; RDKit checks the number of weights.
Point3D compute_centroid(const Conformer &conf, bool ignore_hs,
                         bool has_weights, rust::Slice<const double> weights) {
  std::vector<double> weights_vec(weights.begin(), weights.end());
  auto p = MolTransforms::computeCentroid(conf, ignore_hs,
                                          has_weights ? &weights_vec : nullptr);
  return Point3D{p.x, p.y, p.z};
}

rust::Vec<double> compute_canonical_transform(const Conformer &conf,
                                              bool has_center,
                                              const Point3D &center,
                                              bool normalize_covar,
                                              bool ignore_hs) {
  check_finite_coordinates(conf, "the canonical transform");
  auto origin = to_center(has_center, center);
  std::unique_ptr<RDGeom::Transform3D> trans(
      MolTransforms::computeCanonicalTransform(conf, origin.get(),
                                               normalize_covar, ignore_hs));
  return to_values(*trans);
}

// As computePrincAxesMomentsHelper in rdMolTransforms.cpp, which recomputes
// rather than reading the molecule's cached axes (force). RDKit caches the
// axes in the molecule's properties, where the PMI and PBF descriptors read
// them (U8).
bool compute_principal_axes_and_moments(const Conformer &conf,
                                        bool from_gyration_matrix,
                                        bool ignore_hs, bool has_weights,
                                        rust::Slice<const double> weights,
                                        rust::Vec<double> &axes,
                                        rust::Vec<double> &moments) {
  std::vector<double> weights_vec(weights.begin(), weights.end());
  if (has_weights && weights_vec.size() != conf.getNumAtoms()) {
    throw ValueErrorException(
        "The Python container must have length equal to conf.GetNumAtoms()");
  }
  check_finite_coordinates(conf, "the principal axes");
  for (double w : weights_vec) {
    if (!(std::abs(w) <= 1e9)) {
      throw ValueErrorException(
          "the principal axes need finite weights within 1e9");
    }
  }
  DescriptorCacheGuard guard(conf.getOwningMol());
  Eigen::Matrix3d eigen_axes;
  Eigen::Vector3d eigen_moments;
  const auto *weights_ptr = has_weights ? &weights_vec : nullptr;
  bool success =
      from_gyration_matrix
          ? MolTransforms::computePrincipalAxesAndMomentsFromGyrationMatrix(
                conf, eigen_axes, eigen_moments, ignore_hs, true, weights_ptr)
          : MolTransforms::computePrincipalAxesAndMoments(
                conf, eigen_axes, eigen_moments, ignore_hs, true, weights_ptr);
  axes.clear();
  moments.clear();
  if (success) {
    for (unsigned y = 0; y < 3; ++y) {
      for (unsigned x = 0; x < 3; ++x) {
        axes.push_back(eigen_axes(y, x));
      }
      moments.push_back(eigen_moments(y));
    }
  }
  return success;
}

// transConformer in rdMolTransforms.cpp copies rows * rows values into the
// 4 x 4 transform whatever the array's shape, writing past its end for a
// larger one.
void transform_conformer(Conformer &conf, rust::Slice<const double> trans) {
  if (trans.size() != 16) {
    throw ValueErrorException("the transform must be a 4 x 4 matrix");
  }
  RDGeom::Transform3D transform;
  std::memcpy(transform.getData(), trans.data(), 16 * sizeof(double));
  MolTransforms::transformConformer(conf, transform);
}

void canonicalize_conformer(Conformer &conf, bool has_center,
                            const Point3D &center, bool normalize_covar,
                            bool ignore_hs) {
  check_finite_coordinates(conf, "canonicalization");
  auto origin = to_center(has_center, center);
  MolTransforms::canonicalizeConformer(conf, origin.get(), normalize_covar,
                                       ignore_hs);
}

void canonicalize_mol(RWMol &mol, bool normalize_covar, bool ignore_hs) {
  for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
    check_finite_coordinates(**it, "canonicalization");
  }
  MolTransforms::canonicalizeMol(mol, normalize_covar, ignore_hs);
}

double get_bond_length(const Conformer &conf, std::uint32_t i_atom_id,
                       std::uint32_t j_atom_id) {
  check_atoms(conf, {i_atom_id, j_atom_id});
  return MolTransforms::getBondLength(conf, i_atom_id, j_atom_id);
}

void set_bond_length(Conformer &conf, std::uint32_t i_atom_id,
                     std::uint32_t j_atom_id, double value) {
  check_atoms(conf, {i_atom_id, j_atom_id});
  MolTransforms::setBondLength(conf, i_atom_id, j_atom_id, value);
}

double get_angle(const Conformer &conf, std::uint32_t i_atom_id,
                 std::uint32_t j_atom_id, std::uint32_t k_atom_id,
                 bool degrees) {
  check_atoms(conf, {i_atom_id, j_atom_id, k_atom_id});
  return degrees
             ? MolTransforms::getAngleDeg(conf, i_atom_id, j_atom_id, k_atom_id)
             : MolTransforms::getAngleRad(conf, i_atom_id, j_atom_id,
                                          k_atom_id);
}

void set_angle(Conformer &conf, std::uint32_t i_atom_id,
               std::uint32_t j_atom_id, std::uint32_t k_atom_id, double value,
               bool degrees) {
  check_atoms(conf, {i_atom_id, j_atom_id, k_atom_id});
  if (degrees) {
    MolTransforms::setAngleDeg(conf, i_atom_id, j_atom_id, k_atom_id, value);
  } else {
    MolTransforms::setAngleRad(conf, i_atom_id, j_atom_id, k_atom_id, value);
  }
}

double get_dihedral(const Conformer &conf, std::uint32_t i_atom_id,
                    std::uint32_t j_atom_id, std::uint32_t k_atom_id,
                    std::uint32_t l_atom_id, bool degrees) {
  check_atoms(conf, {i_atom_id, j_atom_id, k_atom_id, l_atom_id});
  return degrees ? MolTransforms::getDihedralDeg(conf, i_atom_id, j_atom_id,
                                                 k_atom_id, l_atom_id)
                 : MolTransforms::getDihedralRad(conf, i_atom_id, j_atom_id,
                                                 k_atom_id, l_atom_id);
}

void set_dihedral(Conformer &conf, std::uint32_t i_atom_id,
                  std::uint32_t j_atom_id, std::uint32_t k_atom_id,
                  std::uint32_t l_atom_id, double value, bool degrees) {
  check_atoms(conf, {i_atom_id, j_atom_id, k_atom_id, l_atom_id});
  if (degrees) {
    MolTransforms::setDihedralDeg(conf, i_atom_id, j_atom_id, k_atom_id,
                                  l_atom_id, value);
  } else {
    MolTransforms::setDihedralRad(conf, i_atom_id, j_atom_id, k_atom_id,
                                  l_atom_id, value);
  }
}

}  // namespace rdkit_cxx
