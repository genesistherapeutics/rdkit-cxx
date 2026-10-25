// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/descriptors_3d.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <GraphMol/Descriptors/MolDescriptors3D.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/descriptors.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/descriptors_3d.rs.h"
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

// GETAWAY sorts the atoms' leverages with std::sort, which is undefined for
// NaN, and non-finite or huge coordinates make NaN leverages. Within 1e9 they
// stay finite.
void check_coordinates(const RWMol &mol, std::int32_t conf_id) {
  for (const auto &pos : mol.getConformer(conf_id).getPositions()) {
    for (double v : {pos.x, pos.y, pos.z}) {
      if (!(std::abs(v) <= 1e9)) {
        throw ValueErrorException(
            "GETAWAY needs finite coordinates within 1e9");
      }
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

}  // namespace

// The moment descriptors cache their moments without the conformer (U8),
// and Python's force=True default is kept.
#define RDKIT_CXX_MOMENTS(name, function)                   \
  double name(const RWMol &mol, std::int32_t conf_id,       \
              bool use_atomic_masses) {                     \
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
  DescriptorCacheGuard guard(mol);
  return Desc::spherocityIndex(mol, conf_id, true);
}

// PBF passes its weights pointer, null, as the force argument of
// computePrincipalAxesAndMomentsFromGyrationMatrix, so it reads the cached
// axes of whichever conformer came first (U8).
double calc_pbf(const RWMol &mol, std::int32_t conf_id) {
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
  return table_descriptor(mol, custom, true, false,
                          [&](std::vector<double> &res) {
                            Desc::WHIM(mol, res, conf_id, thresh, custom);
                          });
}

// GETAWAY converts the precision to unsigned int, which is undefined for a
// negative, huge, or NaN value, and prints each leverage with that many
// digits into a 32-byte buffer, which 25 or more overflow.
rust::Vec<double> calc_getaway(const RWMol &mol, std::int32_t conf_id,
                               double precision,
                               rust::Str custom_atom_property) {
  if (!(precision >= 0.0 && precision < 25.0)) {
    throw ValueErrorException("GETAWAY's precision must be in [0, 25)");
  }
  std::string custom(custom_atom_property);
  if (has_conformers(mol)) {
    check_coordinates(mol, conf_id);
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

}  // namespace rdkit_cxx
