// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mif.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <ForceField/MMFF/Params.h>
#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// RDKit's unit conversion for cube files (MIFDescriptors.cpp).
constexpr double kBohr = 0.529177249;

std::vector<RDGeom::Point3D> to_points(rust::Slice<const Point3D> points) {
  std::vector<RDGeom::Point3D> result;
  result.reserve(points.size());
  for (const auto &p : points) {
    result.emplace_back(p.x, p.y, p.z);
  }
  return result;
}

// As Python's extractChargesAndPositions; RDKit's own precondition is a
// RuntimeError.
void check_charges(rust::Slice<const double> charges,
                   rust::Slice<const Point3D> positions) {
  if (charges.size() != positions.size()) {
    throw ValueErrorException(
        "positions and charges must have the same length");
  }
}

// HBond's constructor stores each nitrogen's and oxygen's bond directions in
// arrays of 12 points (findSpecials, findAcceptors, findDonors, and their
// Unfixed forms in MIFDescriptors.cpp), and writes past the end for an atom
// with more neighbors.
void check_h_bond_degrees(const RWMol &mol) {
  for (const auto *atom : mol.atoms()) {
    auto num = atom->getAtomicNum();
    if ((num == 7 || num == 8) && atom->getDegree() > 12) {
      throw ValueErrorException(
          "HBond needs nitrogens and oxygens with at most 12 neighbors; atom " +
          std::to_string(atom->getIdx()) + " has more");
    }
  }
}

// readFromCubeStream (MIFDescriptors.cpp) reads a cube file with the
// stream extractors into variables it leaves uninitialized, so a value that
// does not parse, or a file that ends early, has it read indeterminate
// values; it makes a conformer of the atom count it reads, about 2^32 atoms
// for a negative count; and it makes the grid from the dimensions and
// spacing, which the grid casts to unsigned int (grid.h). This reads the
// file as RDKit does and refuses it unless every value RDKit reads parses,
// the atom count is not negative, and the grid fits in int, which also
// bounds the atoms and grid points by the file's size. Different spacings
// and dimensions that are not positive it leaves to RDKit, which throws.
void check_cube_text(const std::string &text) {
  std::istringstream in(text);
  std::string line;
  std::getline(in, line);
  std::getline(in, line);
  int n_atoms = 0, dim_x = 0, dim_y = 0, dim_z = 0;
  double x = 0.0, y = 0.0, z = 0.0, sx = 0.0, sy = 0.0, sz = 0.0, t = 0.0;
  in >> n_atoms >> x >> y >> z;
  in >> dim_x >> sx >> t >> t;
  in >> dim_y >> t >> sy >> t;
  in >> dim_z >> t >> t >> sz;
  if (!in) {
    throw ValueErrorException("the cube file's header does not parse");
  }
  if (std::fabs(sx - sy) > 0.0001 || std::fabs(sx - sz) > 0.0001) {
    return;
  }
  double spacing = sx * kBohr;
  double dims[3] = {spacing * dim_x, spacing * dim_y, spacing * dim_z};
  if (!(dims[0] > 0.0 && dims[1] > 0.0 && dims[2] > 0.0 && spacing > 0.0)) {
    return;
  }
  check_grid_dimensions(dims[0], dims[1], dims[2], spacing, true,
                        Point3D{x * kBohr, y * kBohr, z * kBohr});
  if (n_atoms < 0) {
    throw ValueErrorException("the cube file's atom count is negative");
  }
  for (int i = 0; i < n_atoms; ++i) {
    int num = 0;
    in >> num >> t >> x >> y >> z;
    if (!in) {
      throw ValueErrorException(
          "the cube file's atoms do not parse, or are fewer than it says");
    }
  }
  // dim_x, dim_y, and dim_z are positive here, and their product fits in
  // int.
  std::int64_t points = std::int64_t{dim_x} * dim_y * dim_z;
  for (std::int64_t i = 0; i < points; ++i) {
    in >> t;
    if (!in) {
      throw ValueErrorException(
          "the cube file's grid values do not parse, or are fewer than its "
          "grid has points");
    }
  }
}

}  // namespace

std::unique_ptr<UniformRealValueGrid3D> construct_grid(const RWMol &mol,
                                                       std::int32_t conf_id,
                                                       double margin,
                                                       double spacing) {
  // Without conformers, constructGrid's precondition throws.
  if (mol.getNumConformers()) {
    // The box constructGrid computes, whose size the grid casts to unsigned
    // int (grid.h).
    const auto &points = mol.getConformer(conf_id).getPositions();
    RDGeom::Point3D lo =
        points.empty() ? RDGeom::Point3D(0.0, 0.0, 0.0) : points.front();
    auto hi = lo;
    for (const auto &pt : points) {
      for (unsigned i = 0; i < 3; ++i) {
        lo[i] = std::min(lo[i], pt[i]);
        hi[i] = std::max(hi[i], pt[i]);
      }
    }
    RDGeom::Point3D widen(margin, margin, margin);
    lo -= widen;
    hi += widen;
    check_grid_dimensions(hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, spacing, true,
                          Point3D{lo.x, lo.y, lo.z});
  }
  return RDMIF::constructGrid(mol, conf_id, margin, spacing);
}

std::unique_ptr<Coulomb> coulomb_new(const RWMol &mol, std::int32_t conf_id,
                                     double probe_charge, bool abs_val,
                                     rust::Str charge_key,
                                     double softcore_param, double cutoff) {
  return std::make_unique<Coulomb>(mol, conf_id, probe_charge, abs_val,
                                   std::string(charge_key), softcore_param,
                                   cutoff);
}

std::unique_ptr<Coulomb> coulomb_from_charges(
    rust::Slice<const double> charges, rust::Slice<const Point3D> positions,
    double probe_charge, bool abs_val, double softcore_param, double cutoff) {
  check_charges(charges, positions);
  return std::make_unique<Coulomb>(
      std::vector<double>(charges.begin(), charges.end()), to_points(positions),
      probe_charge, abs_val, softcore_param, cutoff);
}

double coulomb_call(const Coulomb &functor, double x, double y, double z,
                    double threshold) {
  return functor(x, y, z, threshold);
}

std::unique_ptr<CoulombDielectric> coulomb_dielectric_new(
    const RWMol &mol, std::int32_t conf_id, double probe_charge, bool abs_val,
    rust::Str charge_key, double softcore_param, double cutoff, double epsilon,
    double xi) {
  // For a molecule without atoms, RDKit's loop over pairs of atoms runs to
  // the number of atoms minus 1, which wraps, and never ends. With no atoms
  // the charges constructor makes the same functor; the precondition and
  // the conformer's lookup still throw first.
  if (mol.getNumAtoms() == 0 && mol.getNumConformers() > 0) {
    mol.getConformer(conf_id);
    return std::make_unique<CoulombDielectric>(
        std::vector<double>(), std::vector<RDGeom::Point3D>(), probe_charge,
        abs_val, softcore_param, cutoff, epsilon, xi);
  }
  return std::make_unique<CoulombDielectric>(
      mol, conf_id, probe_charge, abs_val, std::string(charge_key),
      softcore_param, cutoff, epsilon, xi);
}

std::unique_ptr<CoulombDielectric> coulomb_dielectric_from_charges(
    rust::Slice<const double> charges, rust::Slice<const Point3D> positions,
    double probe_charge, bool abs_val, double softcore_param, double cutoff,
    double epsilon, double xi) {
  check_charges(charges, positions);
  return std::make_unique<CoulombDielectric>(
      std::vector<double>(charges.begin(), charges.end()), to_points(positions),
      probe_charge, abs_val, softcore_param, cutoff, epsilon, xi);
}

double coulomb_dielectric_call(const CoulombDielectric &functor, double x,
                               double y, double z, double threshold) {
  return functor(x, y, z, threshold);
}

std::unique_ptr<MMFFVdWaals> mmff_vd_waals_new(const RWMol &mol,
                                               std::int32_t conf_id,
                                               std::uint32_t probe_atom_type,
                                               bool scaling, double cutoff) {
  // RDKit looks the probe's parameters up without checking that there are
  // any, and dereferences the null it gets for an unknown type.
  if (!(*RDKit::MMFF::DefaultParameters::getMMFFVdW())(probe_atom_type)) {
    throw ValueErrorException(
        "MMFF has no van der Waals parameters for atom "
        "type " +
        std::to_string(probe_atom_type));
  }
  return std::make_unique<MMFFVdWaals>(mol, conf_id, probe_atom_type, scaling,
                                       cutoff);
}

double mmff_vd_waals_call(const MMFFVdWaals &functor, double x, double y,
                          double z, double threshold) {
  return functor(x, y, z, threshold);
}

std::unique_ptr<UFFVdWaals> uff_vd_waals_new(const RWMol &mol,
                                             std::int32_t conf_id,
                                             rust::Str probe_atom_type,
                                             double cutoff) {
  // As for MMFF: RDKit dereferences the null it gets for an unknown type.
  std::string type(probe_atom_type);
  if (!(*ForceFields::UFF::ParamCollection::getParams())(type)) {
    throw ValueErrorException("UFF has no parameters for atom type " + type);
  }
  return std::make_unique<UFFVdWaals>(mol, conf_id, type, cutoff);
}

double uff_vd_waals_call(const UFFVdWaals &functor, double x, double y,
                         double z, double threshold) {
  return functor(x, y, z, threshold);
}

std::unique_ptr<HBond> h_bond_new(const RWMol &mol, std::int32_t conf_id,
                                  rust::Str probe_atom_type, bool fixed,
                                  double cutoff) {
  check_h_bond_degrees(mol);
  // findSpecials matches a serine pattern (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return std::make_unique<HBond>(mol, conf_id, std::string(probe_atom_type),
                                 fixed, cutoff);
}

double h_bond_call(const HBond &functor, double x, double y, double z,
                   double threshold) {
  return functor(x, y, z, threshold);
}

std::unique_ptr<Hydrophilic> hydrophilic_new(const RWMol &mol,
                                             std::int32_t conf_id, bool fixed,
                                             double cutoff) {
  check_h_bond_degrees(mol);
  // Two HBond constructors, which match patterns (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return std::make_unique<Hydrophilic>(mol, conf_id, fixed, cutoff);
}

double hydrophilic_call(const Hydrophilic &functor, double x, double y,
                        double z, double threshold) {
  return functor(x, y, z, threshold);
}

#define RDKIT_CXX_CALCULATE_DESCRIPTORS(name, type)               \
  void name(UniformRealValueGrid3D &grid, const type &descriptor, \
            double threshold) {                                   \
    RDMIF::calculateDescriptors(grid, descriptor, threshold);     \
  }

RDKIT_CXX_CALCULATE_DESCRIPTORS(calculate_descriptors_coulomb, Coulomb)
RDKIT_CXX_CALCULATE_DESCRIPTORS(calculate_descriptors_coulomb_dielectric,
                                CoulombDielectric)
RDKIT_CXX_CALCULATE_DESCRIPTORS(calculate_descriptors_mmff_vd_waals,
                                MMFFVdWaals)
RDKIT_CXX_CALCULATE_DESCRIPTORS(calculate_descriptors_uff_vd_waals, UFFVdWaals)
RDKIT_CXX_CALCULATE_DESCRIPTORS(calculate_descriptors_h_bond, HBond)
RDKIT_CXX_CALCULATE_DESCRIPTORS(calculate_descriptors_hydrophilic, Hydrophilic)

#undef RDKIT_CXX_CALCULATE_DESCRIPTORS

// writeToCubeFile opens the file itself and ignores a stream that fails
// (U12), so the wrapper opens it and checks the stream after writing.
void write_to_cube_file(const UniformRealValueGrid3D &grid,
                        rust::Slice<const std::uint8_t> filename,
                        const RWMol *mol, std::int32_t conf_id) {
  auto name = to_path(filename);
  if (mol) {
    mol->getConformer(conf_id);  // throws before the file is made
  }
  std::ofstream out(name);
  if (!out) {
    throw RDKit::BadFileException("Bad output file " + name);
  }
  RDMIF::writeToCubeStream(grid, out, mol, conf_id);
  out.close();
  if (!out) {
    throw RDKit::BadFileException("could not write " + name);
  }
}

std::unique_ptr<RWMol> read_from_cube_file(
    rust::Slice<const std::uint8_t> filename, UniformRealValueGrid3D &grid) {
  auto text = read_file(to_path(filename));
  check_cube_text(text);
  std::istringstream in(text);
  return RDMIF::readFromCubeStream(grid, in);
}

}  // namespace rdkit_cxx
