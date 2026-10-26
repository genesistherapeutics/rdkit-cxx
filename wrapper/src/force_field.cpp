// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/force_field.h"

#include <string>

#include <ForceField/AngleConstraint.h>
#include <ForceField/DistanceConstraint.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/PositionConstraint.h>
#include <ForceField/TorsionConstraint.h>
#include <GraphMol/Trajectory/Snapshot.h>
#include <GraphMol/Trajectory/Trajectory.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/force_field.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// The contributions index the positions they were given, but the energy and
// gradient with explicit coordinates, and the gradient without them, size
// their arrays by numPoints, which only initialize() updates. A point added
// since would be read or written past those arrays' ends.
void check_initialized(const ForceFields::ForceField &field) {
  if (field.positions().size() != field.numPoints()) {
    throw ValueErrorException(
        "the force field has points added since it was initialized; call "
        "Initialize first");
  }
}

std::vector<double> check_pos(const ForceFields::ForceField &field,
                              rust::Slice<const double> pos) {
  check_initialized(field);
  // As PyForceField::calcEnergyWithPos in Code/ForceField/Wrap/ForceField.cpp.
  if (pos.size() !=
      static_cast<std::size_t>(field.dimension()) * field.numPoints()) {
    throw ValueErrorException(
        "The Python container must have length equal to Dimension() * "
        "NumPoints()");
  }
  return std::vector<double>(pos.begin(), pos.end());
}

std::size_t num_coords(const ForceFields::ForceField &field) {
  return static_cast<std::size_t>(field.dimension()) * field.numPoints();
}

TermParams not_found() { return TermParams{false, 0, {}}; }

TermParams found(unsigned type, std::initializer_list<double> values) {
  rust::Vec<double> result;
  for (auto value : values) {
    result.push_back(value);
  }
  return TermParams{true, type, std::move(result)};
}

}  // namespace

void check_minimize_size(std::size_t num_points, unsigned dimension) {
  if (num_points * dimension > 65535) {
    throw ValueErrorException(
        "RDKit's minimizer overflows above 65535 coordinates (points times "
        "dimension)");
  }
}

void check_mmff_indices(const MMFFMolProperties &props, const RWMol &mol,
                        std::initializer_list<unsigned> indices) {
  for (auto idx : indices) {
    check_index(idx, mol.getNumAtoms());
    check_index(idx, props.num_atoms());
  }
}

double force_field_calc_energy(const ForceField &ff) {
  return ff.field().calcEnergy();
}

double force_field_calc_energy_with_pos(ForceField &ff,
                                        rust::Slice<const double> pos) {
  auto coords = check_pos(ff.field(), pos);
  return ff.field().calcEnergy(coords.data());
}

rust::Vec<double> force_field_calc_grad(const ForceField &ff) {
  check_initialized(ff.field());
  std::vector<double> grad(num_coords(ff.field()), 0.0);
  ff.field().calcGrad(grad.data());
  return to_vec(grad);
}

rust::Vec<double> force_field_calc_grad_with_pos(
    ForceField &ff, rust::Slice<const double> pos) {
  auto coords = check_pos(ff.field(), pos);
  std::vector<double> grad(coords.size(), 0.0);
  ff.field().calcGrad(coords.data(), grad.data());
  return to_vec(grad);
}

// As PyForceField::positions; every point is a Point3D.
rust::Vec<double> force_field_positions(const ForceField &ff) {
  rust::Vec<double> coords;
  for (const auto *point : ff.field().positions()) {
    for (unsigned j = 0; j < 3; ++j) {
      coords.push_back((*point)[j]);
    }
  }
  return coords;
}

unsigned force_field_dimension(const ForceField &ff) {
  return ff.field().dimension();
}

unsigned force_field_num_points(const ForceField &ff) {
  return ff.field().numPoints();
}

int force_field_minimize(ForceField &ff, unsigned max_its, double force_tol,
                         double energy_tol) {
  check_minimize_size(ff.field().positions().size(), ff.field().dimension());
  return ff.field().minimize(max_its, force_tol, energy_tol);
}

// As PyForceField::minimizeTrajectory. A snapshot's coordinates are read
// through a Trajectory, which is what Snapshot::getPoint3D requires.
int force_field_minimize_trajectory(ForceField &ff, unsigned snapshot_freq,
                                    unsigned max_its, double force_tol,
                                    double energy_tol,
                                    rust::Vec<double> &energies,
                                    rust::Vec<double> &positions) {
  auto &field = ff.field();
  check_minimize_size(field.positions().size(), field.dimension());
  auto snapshots = std::make_unique<RDKit::SnapshotVect>();
  int res = field.minimize(snapshot_freq, snapshots.get(), max_its, force_tol,
                           energy_tol);
  unsigned num_points = field.numPoints();
  unsigned dimension = field.dimension();
  RDKit::Trajectory trajectory(dimension, num_points, snapshots.release());
  energies.clear();
  positions.clear();
  for (unsigned i = 0; i < trajectory.size(); ++i) {
    const auto &snapshot = trajectory.getSnapshot(i);
    energies.push_back(snapshot.getEnergy());
    for (unsigned p = 0; p < num_points; ++p) {
      auto point = snapshot.getPoint3D(p);
      for (unsigned j = 0; j < dimension && j < 3; ++j) {
        positions.push_back(point[j]);
      }
    }
  }
  return res;
}

// The constraints range-check their indices, and Python's IndexError comes
// from the wrapper's check first (E3).
void force_field_add_distance_constraint(ForceField &ff, unsigned idx1,
                                         unsigned idx2, double min_len,
                                         double max_len,
                                         double force_constant) {
  auto &field = ff.field();
  check_index(idx1, field.positions().size());
  check_index(idx2, field.positions().size());
  field.contribs().emplace_back(new ForceFields::DistanceConstraintContrib(
      &field, idx1, idx2, min_len, max_len, force_constant));
}

// RDKit keeps any index here and checks it when it computes a gradient.
void force_field_add_fixed_point(ForceField &ff, unsigned idx) {
  auto &field = ff.field();
  check_index(idx, field.positions().size());
  field.fixedPoints().push_back(static_cast<int>(idx));
}

// UFF's and MMFF's constraints are the same classes
// (ForceField/UFF/DistanceConstraint.h, MMFF/DistanceConstraint.h).
void force_field_add_relative_distance_constraint(ForceField &ff, unsigned idx1,
                                                  unsigned idx2, bool relative,
                                                  double min_len,
                                                  double max_len,
                                                  double force_constant) {
  auto &field = ff.field();
  check_index(idx1, field.positions().size());
  check_index(idx2, field.positions().size());
  field.contribs().emplace_back(new ForceFields::DistanceConstraintContrib(
      &field, idx1, idx2, relative, min_len, max_len, force_constant));
}

void force_field_add_angle_constraint(ForceField &ff, unsigned idx1,
                                      unsigned idx2, unsigned idx3,
                                      bool relative, double min_angle_deg,
                                      double max_angle_deg,
                                      double force_constant) {
  auto &field = ff.field();
  for (auto idx : {idx1, idx2, idx3}) {
    check_index(idx, field.positions().size());
  }
  field.contribs().emplace_back(new ForceFields::AngleConstraintContrib(
      &field, idx1, idx2, idx3, relative, min_angle_deg, max_angle_deg,
      force_constant));
}

void force_field_add_torsion_constraint(ForceField &ff, unsigned idx1,
                                        unsigned idx2, unsigned idx3,
                                        unsigned idx4, bool relative,
                                        double min_dihedral_deg,
                                        double max_dihedral_deg,
                                        double force_constant) {
  auto &field = ff.field();
  for (auto idx : {idx1, idx2, idx3, idx4}) {
    check_index(idx, field.positions().size());
  }
  field.contribs().emplace_back(new ForceFields::TorsionConstraintContrib(
      &field, idx1, idx2, idx3, idx4, relative, min_dihedral_deg,
      max_dihedral_deg, force_constant));
}

void force_field_add_position_constraint(ForceField &ff, unsigned idx,
                                         double max_displ,
                                         double force_constant) {
  auto &field = ff.field();
  check_index(idx, field.positions().size());
  field.contribs().emplace_back(new ForceFields::PositionConstraintContrib(
      &field, idx, max_displ, force_constant));
}

void force_field_initialize(ForceField &ff) { ff.field().initialize(); }

// As PyForceField::addExtraPoint, which returns the new number of points.
unsigned force_field_add_extra_point(ForceField &ff, double x, double y,
                                     double z, bool fixed) {
  auto &field = ff.field();
  ff.extra_points().push_back(std::make_unique<RDGeom::Point3D>(x, y, z));
  field.positions().push_back(ff.extra_points().back().get());
  auto count = static_cast<unsigned>(field.positions().size());
  if (fixed) {
    field.fixedPoints().push_back(static_cast<int>(count - 1));
  }
  return count;
}

Point3D force_field_get_extra_point_pos(const ForceField &ff, unsigned idx) {
  check_index(idx, ff.extra_points().size());
  const auto &point = *ff.extra_points()[idx];
  return Point3D{point.x, point.y, point.z};
}

std::unique_ptr<MMFFMolProperties> mmff_mol_properties_new(
    RWMol &mol, rust::Str variant, MMFFVerbosity verbosity) {
  return std::make_unique<MMFFMolProperties>(mol, std::string(variant),
                                             verbosity);
}

bool mmff_mol_properties_is_valid(const MMFFMolProperties &props) {
  return props.props().isValid();
}

// The getters range-check their indices with URANGE_CHECK; the wrapper's
// check comes first (E3).
unsigned mmff_mol_properties_get_atom_type(const MMFFMolProperties &props,
                                           unsigned idx) {
  check_index(idx, props.num_atoms());
  return props.props().getMMFFAtomType(idx);
}

double mmff_mol_properties_get_formal_charge(const MMFFMolProperties &props,
                                             unsigned idx) {
  check_index(idx, props.num_atoms());
  return props.props().getMMFFFormalCharge(idx);
}

double mmff_mol_properties_get_partial_charge(const MMFFMolProperties &props,
                                              unsigned idx) {
  check_index(idx, props.num_atoms());
  return props.props().getMMFFPartialCharge(idx);
}

// The parameter getters, as PyMMFFMolProperties' in
// Code/ForceField/Wrap/ForceField.cpp, take atom indices into the molecule
// passed and look their types up in the properties, so each index is checked
// against both.
TermParams mmff_mol_properties_get_bond_stretch_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2) {
  check_mmff_indices(props, mol, {idx1, idx2});
  unsigned bond_type;
  ForceFields::MMFF::MMFFBond params;
  if (!props.props().getMMFFBondStretchParams(mol, idx1, idx2, bond_type,
                                              params)) {
    return not_found();
  }
  return found(bond_type, {params.kb, params.r0});
}

TermParams mmff_mol_properties_get_angle_bend_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3) {
  check_mmff_indices(props, mol, {idx1, idx2, idx3});
  unsigned angle_type;
  ForceFields::MMFF::MMFFAngle params;
  if (!props.props().getMMFFAngleBendParams(mol, idx1, idx2, idx3, angle_type,
                                            params)) {
    return not_found();
  }
  return found(angle_type, {params.ka, params.theta0});
}

TermParams mmff_mol_properties_get_stretch_bend_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3) {
  check_mmff_indices(props, mol, {idx1, idx2, idx3});
  unsigned stretch_bend_type;
  ForceFields::MMFF::MMFFStbn params;
  ForceFields::MMFF::MMFFBond bond_params[2];
  ForceFields::MMFF::MMFFAngle angle_params;
  if (!props.props().getMMFFStretchBendParams(mol, idx1, idx2, idx3,
                                              stretch_bend_type, params,
                                              bond_params, angle_params)) {
    return not_found();
  }
  return found(stretch_bend_type, {params.kbaIJK, params.kbaKJI});
}

TermParams mmff_mol_properties_get_torsion_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3, unsigned idx4) {
  check_mmff_indices(props, mol, {idx1, idx2, idx3, idx4});
  unsigned torsion_type;
  ForceFields::MMFF::MMFFTor params;
  if (!props.props().getMMFFTorsionParams(mol, idx1, idx2, idx3, idx4,
                                          torsion_type, params)) {
    return not_found();
  }
  return found(torsion_type, {params.V1, params.V2, params.V3});
}

TermParams mmff_mol_properties_get_oop_bend_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3, unsigned idx4) {
  check_mmff_indices(props, mol, {idx1, idx2, idx3, idx4});
  ForceFields::MMFF::MMFFOop params;
  if (!props.props().getMMFFOopBendParams(mol, idx1, idx2, idx3, idx4,
                                          params)) {
    return not_found();
  }
  return found(0, {params.koop});
}

TermParams mmff_mol_properties_get_vdw_params(const MMFFMolProperties &props,
                                              unsigned idx1, unsigned idx2) {
  check_index(idx1, props.num_atoms());
  check_index(idx2, props.num_atoms());
  ForceFields::MMFF::MMFFVdWRijstarEps params;
  if (!props.props().getMMFFVdWParams(idx1, idx2, params)) {
    return not_found();
  }
  return found(0, {params.R_ij_starUnscaled, params.epsilonUnscaled,
                   params.R_ij_star, params.epsilon});
}

void mmff_mol_properties_set_dielectric_model(MMFFMolProperties &props,
                                              std::uint8_t model) {
  props.props().setMMFFDielectricModel(model);
}

std::uint8_t mmff_mol_properties_get_dielectric_model(
    const MMFFMolProperties &props) {
  return props.props().getMMFFDielectricModel();
}

void mmff_mol_properties_set_dielectric_constant(MMFFMolProperties &props,
                                                 double constant) {
  props.props().setMMFFDielectricConstant(constant);
}

double mmff_mol_properties_get_dielectric_constant(
    const MMFFMolProperties &props) {
  return props.props().getMMFFDielectricConstant();
}

void mmff_mol_properties_set_term(MMFFMolProperties &props, MMFFTerm term,
                                  bool state) {
  auto &p = props.props();
  switch (term) {
    case MMFFTerm::Bond:
      p.setMMFFBondTerm(state);
      break;
    case MMFFTerm::Angle:
      p.setMMFFAngleTerm(state);
      break;
    case MMFFTerm::StretchBend:
      p.setMMFFStretchBendTerm(state);
      break;
    case MMFFTerm::Oop:
      p.setMMFFOopTerm(state);
      break;
    case MMFFTerm::Torsion:
      p.setMMFFTorsionTerm(state);
      break;
    case MMFFTerm::VdW:
      p.setMMFFVdWTerm(state);
      break;
    case MMFFTerm::Ele:
      p.setMMFFEleTerm(state);
      break;
  }
}

bool mmff_mol_properties_get_term(const MMFFMolProperties &props,
                                  MMFFTerm term) {
  const auto &p = props.props();
  switch (term) {
    case MMFFTerm::Bond:
      return p.getMMFFBondTerm();
    case MMFFTerm::Angle:
      return p.getMMFFAngleTerm();
    case MMFFTerm::StretchBend:
      return p.getMMFFStretchBendTerm();
    case MMFFTerm::Oop:
      return p.getMMFFOopTerm();
    case MMFFTerm::Torsion:
      return p.getMMFFTorsionTerm();
    case MMFFTerm::VdW:
      return p.getMMFFVdWTerm();
    case MMFFTerm::Ele:
      return p.getMMFFEleTerm();
  }
  return false;
}

void mmff_mol_properties_set_variant(MMFFMolProperties &props,
                                     rust::Str variant) {
  props.props().setMMFFVariant(std::string(variant));
}

std::unique_ptr<std::string> mmff_mol_properties_get_variant(
    const MMFFMolProperties &props) {
  return std::make_unique<std::string>(props.props().getMMFFVariant());
}

void mmff_mol_properties_set_verbosity(MMFFMolProperties &props,
                                       MMFFVerbosity verbosity) {
  props.props().setMMFFVerbosity(verbosity);
}

}  // namespace rdkit_cxx
