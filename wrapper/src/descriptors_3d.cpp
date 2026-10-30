// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/descriptors_3d.h"

#include <algorithm>
#include <cmath>
#include <mutex>
#include <string>
#include <vector>

#include <GraphMol/Descriptors/MolDescriptors3D.h>
#include <GraphMol/Descriptors/USRDescriptor.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/coordinate_checks.h"
#include "rdkit-cxx/descriptors.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/descriptors_3d.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace Desc = RDKit::Descriptors;

// The checks below read the conformer, so they run only when there is one;
// otherwise RDKit's own precondition reports it.
bool has_conformers(const RWMol &mol) { return mol.getNumConformers() > 0; }

// AUTOCORR3D, WHIM, and GETAWAY take &v[0] of per-atom vectors, which is
// undefined for a molecule without atoms (GETAWAY crashes).
void check_has_atoms(const RWMol &mol) {
  if (mol.getNumAtoms() == 0) {
    throw ValueErrorException("these descriptors need at least one atom");
  }
}

// The descriptors that solve matrices built from the conformer's
// coordinates (coordinate_checks.h).
void check_coordinates(const RWMol &mol, std::int32_t conf_id,
                       const std::string &what) {
  if (has_conformers(mol)) {
    check_finite_coordinates(mol.getConformer(conf_id), what);
  }
}

// WHIM weights the coordinates' covariance matrix, which it solves, by a
// custom atom property's values.
void check_custom_weights(const RWMol &mol, const std::string &custom) {
  for (const auto *atom : mol.atoms()) {
    double v = 1.0;
    atom->getPropIfPresent(custom, v);
    if (!(std::abs(v) <= 1e9)) {
      throw ValueErrorException("WHIM needs the values of " + custom +
                                " to be finite and within 1e9");
    }
  }
}

// The checks the property-table descriptors need. All but WHIM also read
// the distance matrices, which they index in int (U3, U6).
template <class Call>
rust::Vec<double> table_descriptor(const RWMol &mol, const std::string &custom,
                                   bool needs_atoms, bool reads_matrices,
                                   Call call) {
  if (reads_matrices) {
    check_matrix_size(mol);
  }
  if (has_conformers(mol)) {
    if (needs_atoms) {
      check_has_atoms(mol);
    }
    if (custom.empty()) {
      check_property_table_atoms(mol);
    }
  }
  MatrixCacheGuard guard(mol);
  std::vector<double> result;
  call(result);
  return to_vec(result);
}

// EEM's parameter tables (A1 to B3 in EEM.cpp) have 36 entries, indexed by
// atomic number.
constexpr unsigned int eem_table_size = 36;

// EEM's atom type: the highest bond order at the atom, aromatic counting as
// 2 (getAtomtype in EEM.cpp). RDKit leaves an atom's entry of its b vector
// uninitialized for types above 3, and reads it.
unsigned int eem_atom_type(const RDKit::ROMol &mol, const RDKit::Atom *atom) {
  unsigned int type = 1;
  for (const auto *bond : mol.atomBonds(atom)) {
    double order = bond->getBondTypeAsDouble();
    if (order == 1.5) {
      order = 2.0;
    }
    type = std::max(type, static_cast<unsigned int>(order));
  }
  return type;
}

RDGeom::Point3D theirs(const Point3D &p) {
  return RDGeom::Point3D(p.x, p.y, p.z);
}

Point3D ours(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

rust::Vec<DoubleList> to_lists(
    const std::vector<std::vector<double>> &distributions) {
  rust::Vec<DoubleList> lists;
  for (const auto &values : distributions) {
    lists.push_back(DoubleList{to_vec(values)});
  }
  return lists;
}

// As GetUSR and GetUSRCAT in the Python wrappers.
void check_usr_molecule(const RWMol &mol) {
  if (mol.getNumConformers() == 0) {
    throw ValueErrorException("no conformers");
  }
  if (mol.getNumAtoms() < 3) {
    throw ValueErrorException("too few atoms (minimum three)");
  }
}

}  // namespace

// The moment descriptors cache their moments without the conformer (U8),
// and Python's force=True default is kept.
#define RDKIT_CXX_MOMENTS(name, function)                   \
  double name(const RWMol &mol, std::int32_t conf_id,       \
              bool use_atomic_masses) {                     \
    check_coordinates(mol, conf_id, "the moments");         \
    DescriptorCacheGuard guard(mol);                        \
    return function(mol, conf_id, use_atomic_masses, true); \
  }

RDKIT_CXX_MOMENTS(calc_pmi1, Desc::PMI1)
RDKIT_CXX_MOMENTS(calc_pmi2, Desc::PMI2)
RDKIT_CXX_MOMENTS(calc_pmi3, Desc::PMI3)
RDKIT_CXX_MOMENTS(calc_npr1, Desc::NPR1)
RDKIT_CXX_MOMENTS(calc_npr2, Desc::NPR2)
RDKIT_CXX_MOMENTS(calc_radius_of_gyration, Desc::radiusOfGyration)
RDKIT_CXX_MOMENTS(calc_inertial_shape_factor, Desc::inertialShapeFactor)
RDKIT_CXX_MOMENTS(calc_eccentricity, Desc::eccentricity)
RDKIT_CXX_MOMENTS(calc_asphericity, Desc::asphericity)

#undef RDKIT_CXX_MOMENTS

double calc_spherocity_index(const RWMol &mol, std::int32_t conf_id) {
  check_coordinates(mol, conf_id, "the spherocity index");
  DescriptorCacheGuard guard(mol);
  return Desc::spherocityIndex(mol, conf_id, true);
}

// PBF passes its weights pointer, null, as the force argument of
// computePrincipalAxesAndMomentsFromGyrationMatrix, so it reads the cached
// axes of whichever conformer came first (U8).
double calc_pbf(const RWMol &mol, std::int32_t conf_id) {
  check_coordinates(mol, conf_id, "PBF");
  DescriptorCacheGuard guard(mol);
  return Desc::PBF(mol, conf_id);
}

rust::Vec<double> calc_autocorr_3d(const RWMol &mol, std::int32_t conf_id,
                                   rust::Str custom_atom_property) {
  std::string custom(custom_atom_property);
  return table_descriptor(mol, custom, true, true,
                          [&](std::vector<double> &res) {
                            Desc::AUTOCORR3D(mol, res, conf_id, custom);
                          });
}

rust::Vec<double> calc_rdf(const RWMol &mol, std::int32_t conf_id,
                           rust::Str custom_atom_property) {
  std::string custom(custom_atom_property);
  return table_descriptor(
      mol, custom, false, true,
      [&](std::vector<double> &res) { Desc::RDF(mol, res, conf_id, custom); });
}

rust::Vec<double> calc_morse(const RWMol &mol, std::int32_t conf_id,
                             rust::Str custom_atom_property) {
  std::string custom(custom_atom_property);
  return table_descriptor(mol, custom, false, true,
                          [&](std::vector<double> &res) {
                            Desc::MORSE(mol, res, conf_id, custom);
                          });
}

rust::Vec<double> calc_whim(const RWMol &mol, std::int32_t conf_id,
                            double thresh, rust::Str custom_atom_property) {
  std::string custom(custom_atom_property);
  check_coordinates(mol, conf_id, "WHIM");
  if (!custom.empty()) {
    check_custom_weights(mol, custom);
  }
  return table_descriptor(mol, custom, true, false,
                          [&](std::vector<double> &res) {
                            Desc::WHIM(mol, res, conf_id, thresh, custom);
                          });
}

// GETAWAY converts the precision to unsigned int, which is undefined for a
// negative, huge, or NaN value, and prints each leverage with that many
// digits into a 32-byte buffer, which 25 or more overflow. It sorts the
// atoms' leverages with std::sort, which is undefined for NaN, and divides
// by the atoms' distances (coordinate_checks.h).
rust::Vec<double> calc_getaway(const RWMol &mol, std::int32_t conf_id,
                               double precision,
                               rust::Str custom_atom_property) {
  if (!(precision >= 0.0 && precision < 25.0)) {
    throw ValueErrorException("GETAWAY's precision must be in [0, 25)");
  }
  std::string custom(custom_atom_property);
  if (has_conformers(mol)) {
    const auto &conf = mol.getConformer(conf_id);
    check_finite_coordinates(conf, "GETAWAY");
    check_distinct_positions(conf, "GETAWAY");
  }
  return table_descriptor(
      mol, custom, true, true, [&](std::vector<double> &res) {
        Desc::GETAWAY(mol, res, conf_id, static_cast<unsigned int>(precision),
                      custom);
      });
}

rust::Vec<double> calc_coulomb_mat(const RWMol &mol, std::int32_t conf_id) {
  std::vector<std::vector<double>> matrix;
  Desc::CoulombMat(mol, matrix, conf_id);
  rust::Vec<double> values;
  values.reserve(matrix.size() * matrix.size());
  for (const auto &row : matrix) {
    for (double v : row) {
      values.push_back(v);
    }
  }
  return values;
}

// EEM takes a non-const molecule but only reads it, so it is given a
// kekulized copy, the molecule it computes on.
rust::Vec<double> calc_eem_charges(const RWMol &mol, std::int32_t conf_id) {
  RDKit::RWMol copy(mol);
  if (has_conformers(mol)) {
    // EEM sizes its (n + 1) * (n + 1) matrix in int.
    if (mol.getNumAtoms() > 46339) {
      throw ValueErrorException(
          "EEM sizes its matrix in int, which overflows above 46339 atoms");
    }
    // EEM solves a matrix of the inverse distances (coordinate_checks.h).
    const auto &conf = mol.getConformer(conf_id);
    check_finite_coordinates(conf, "EEM");
    check_distinct_positions(conf, "EEM");
    RDKit::MolOps::Kekulize(copy, true);
    for (const auto *atom : copy.atoms()) {
      if (static_cast<unsigned int>(atom->getAtomicNum()) >= eem_table_size) {
        throw ValueErrorException(
            "EEM's parameter tables cover atomic numbers below 36; atom " +
            std::to_string(atom->getIdx()) + " has atomic number " +
            std::to_string(atom->getAtomicNum()));
      }
      if (eem_atom_type(copy, atom) > 3) {
        throw ValueErrorException(
            "EEM has parameters for bond orders up to 3; atom " +
            std::to_string(atom->getIdx()) + " has a higher one");
      }
    }
  }
  std::vector<double> result;
  Desc::EEM(copy, result, conf_id);
  return to_vec(result);
}

rust::Vec<double> get_usr(const RWMol &mol, std::int32_t conf_id) {
  check_usr_molecule(mol);
  std::vector<double> descriptor(12);
  Desc::USR(mol, descriptor, conf_id);
  return to_vec(descriptor);
}

// USRCAT reads the conformer's coordinates at the selected atoms without
// checking them, so each atom number is checked first: 0, which Python's
// subtraction wraps, and numbers past the last atom are IndexErrors.
rust::Vec<double> get_usrcat(const RWMol &mol, bool has_atom_selections,
                             rust::Slice<const IndexList> atom_selections,
                             std::int32_t conf_id) {
  check_usr_molecule(mol);
  std::vector<std::vector<unsigned int>> atom_ids;
  if (has_atom_selections) {
    if (atom_selections.empty()) {
      throw ValueErrorException("empty atom selections");
    }
    for (const auto &selection : atom_selections) {
      std::vector<unsigned int> ids;
      for (auto number : selection.indices) {
        auto idx = number - 1;
        check_index(idx, mol.getNumAtoms());
        ids.push_back(idx);
      }
      atom_ids.push_back(std::move(ids));
    }
  }
  std::vector<double> descriptor(
      12 * ((has_atom_selections ? atom_ids.size() : 4) + 1));
  // Without selections, USRCAT finds its features by substructure matches
  // (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  Desc::USRCAT(mol, descriptor, atom_ids, conf_id);
  return to_vec(descriptor);
}

rust::Vec<DoubleList> get_usr_distributions(rust::Slice<const Point3D> coords,
                                            rust::Vec<Point3D> &points) {
  if (coords.empty()) {
    throw ValueErrorException("no coordinates");
  }
  std::vector<RDGeom::Point3D> stored;
  stored.reserve(coords.size());
  for (const auto &p : coords) {
    stored.push_back(theirs(p));
  }
  RDGeom::Point3DConstPtrVect ptrs;
  for (const auto &p : stored) {
    ptrs.push_back(&p);
  }
  std::vector<RDGeom::Point3D> refs(4);
  std::vector<std::vector<double>> distances(4);
  Desc::calcUSRDistributions(ptrs, distances, refs);
  for (const auto &p : refs) {
    points.push_back(ours(p));
  }
  return to_lists(distances);
}

rust::Vec<DoubleList> get_usr_distributions_from_points(
    rust::Slice<const Point3D> coords, rust::Slice<const Point3D> points) {
  if (coords.empty()) {
    throw ValueErrorException("no coordinates");
  }
  if (points.empty()) {
    throw ValueErrorException("no points");
  }
  std::vector<RDGeom::Point3D> stored;
  stored.reserve(coords.size());
  for (const auto &p : coords) {
    stored.push_back(theirs(p));
  }
  RDGeom::Point3DConstPtrVect ptrs;
  for (const auto &p : stored) {
    ptrs.push_back(&p);
  }
  std::vector<RDGeom::Point3D> refs;
  for (const auto &p : points) {
    refs.push_back(theirs(p));
  }
  std::vector<std::vector<double>> distances(refs.size());
  Desc::calcUSRDistributionsFromPoints(ptrs, refs, distances);
  return to_lists(distances);
}

rust::Vec<double> get_usr_from_distributions(
    rust::Slice<const DoubleList> distances) {
  if (distances.empty()) {
    throw ValueErrorException("no distances");
  }
  std::vector<std::vector<double>> dist;
  for (const auto &list : distances) {
    if (list.values.empty()) {
      throw ValueErrorException("distances missing");
    }
    dist.emplace_back(list.values.begin(), list.values.end());
  }
  std::vector<double> descriptor(12);
  Desc::calcUSRFromDistributions(dist, descriptor);
  return to_vec(descriptor);
}

// As GetUSRScore in the Python wrappers: no weights means 1 for each set of
// 12 values.
double get_usr_score(rust::Slice<const double> descriptor1,
                     rust::Slice<const double> descriptor2,
                     rust::Slice<const double> weights) {
  if (descriptor1.size() != descriptor2.size()) {
    throw ValueErrorException("descriptors must have the same length");
  }
  std::size_t num_weights = descriptor1.size() / 12;
  std::vector<double> w(num_weights, 1.0);
  if (!weights.empty()) {
    if (weights.size() != num_weights) {
      throw ValueErrorException("number of weights is not correct");
    }
    w.assign(weights.begin(), weights.end());
  }
  std::vector<double> d1(descriptor1.begin(), descriptor1.end());
  std::vector<double> d2(descriptor2.begin(), descriptor2.end());
  return Desc::calcUSRScore(d1, d2, w);
}

}  // namespace rdkit_cxx
