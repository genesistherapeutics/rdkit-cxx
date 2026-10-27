/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdForceField's ForceField and MMFFMolProperties.

#pragma once

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include <ForceField/ForceField.h>
#include <Geometry/point.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// RDKit's MMFF_VERBOSITY values, defined by the bridge
// (src/ffi/force_field.rs): RDKit 2025.03 declares them in an unnamed enum.
enum class MMFFVerbosity : std::uint8_t;

// Shared structs, defined by the bridges (src/ffi/force_field.rs,
// geometry.rs).
struct Point3D;
struct TermParams;
enum class MMFFTerm : std::uint8_t;

// A force field and the extra points added to it, as PyForceField in
// Code/ForceField/Wrap/PyForceField.h, but owned by unique_ptr where Python
// shares them. Its other positions point into a conformer, which the Rust
// ForceField borrows for as long as it lives.
class ForceField {
 public:
  explicit ForceField(std::unique_ptr<ForceFields::ForceField> field)
      : d_field(std::move(field)) {}
  ForceFields::ForceField &field() { return *d_field; }
  const ForceFields::ForceField &field() const { return *d_field; }
  // Each extra point lives in its own allocation, so the force field's
  // pointers to them survive later additions.
  std::vector<std::unique_ptr<RDGeom::Point3D>> &extra_points() {
    return d_extra_points;
  }
  const std::vector<std::unique_ptr<RDGeom::Point3D>> &extra_points() const {
    return d_extra_points;
  }

 private:
  std::unique_ptr<ForceFields::ForceField> d_field;
  std::vector<std::unique_ptr<RDGeom::Point3D>> d_extra_points;
};

// MMFF atom types and charges for a molecule, which remember how many atoms
// it had: RDKit's getters take atom indices for any molecule.
class MMFFMolProperties {
 public:
  MMFFMolProperties(RDKit::ROMol &mol, const std::string &variant,
                    std::uint8_t verbosity)
      : d_props(std::make_unique<RDKit::MMFF::MMFFMolProperties>(mol, variant,
                                                                 verbosity)),
        d_num_atoms(mol.getNumAtoms()) {}
  // RDKit 2025.03 declares MMFFMolProperties' getters non-const, though
  // they only read (RDKit 2026.03 declares them const), so the const
  // overload gives them the object, which is never const itself.
  RDKit::MMFF::MMFFMolProperties &props() const { return *d_props; }
  unsigned num_atoms() const { return d_num_atoms; }

 private:
  std::unique_ptr<RDKit::MMFF::MMFFMolProperties> d_props;
  unsigned d_num_atoms;
};

// BFGSOpt::minimize sizes its inverse Hessian as dim * dim in unsigned int,
// which wraps above 65535 coordinates, and then writes past its end. Throws
// ValueErrorException for more coordinates (U3).
void check_minimize_size(std::size_t num_points, unsigned dimension);

// Throws IndexError unless each of the molecule's atom indices is in range for
// both the molecule and the properties (E3).
void check_mmff_indices(const MMFFMolProperties &props, const RWMol &mol,
                        std::initializer_list<unsigned> indices);

double force_field_calc_energy(const ForceField &ff);
double force_field_calc_energy_with_pos(ForceField &ff,
                                        rust::Slice<const double> pos);
rust::Vec<double> force_field_calc_grad(const ForceField &ff);
rust::Vec<double> force_field_calc_grad_with_pos(ForceField &ff,
                                                 rust::Slice<const double> pos);
rust::Vec<double> force_field_positions(const ForceField &ff);
unsigned force_field_dimension(const ForceField &ff);
unsigned force_field_num_points(const ForceField &ff);
int force_field_minimize(ForceField &ff, unsigned max_its, double force_tol,
                         double energy_tol);
int force_field_minimize_trajectory(ForceField &ff, unsigned snapshot_freq,
                                    unsigned max_its, double force_tol,
                                    double energy_tol,
                                    rust::Vec<double> &energies,
                                    rust::Vec<double> &positions);
void force_field_add_distance_constraint(ForceField &ff, unsigned idx1,
                                         unsigned idx2, double min_len,
                                         double max_len, double force_constant);
void force_field_add_fixed_point(ForceField &ff, unsigned idx);
void force_field_add_relative_distance_constraint(ForceField &ff, unsigned idx1,
                                                  unsigned idx2, bool relative,
                                                  double min_len,
                                                  double max_len,
                                                  double force_constant);
void force_field_add_angle_constraint(ForceField &ff, unsigned idx1,
                                      unsigned idx2, unsigned idx3,
                                      bool relative, double min_angle_deg,
                                      double max_angle_deg,
                                      double force_constant);
void force_field_add_torsion_constraint(ForceField &ff, unsigned idx1,
                                        unsigned idx2, unsigned idx3,
                                        unsigned idx4, bool relative,
                                        double min_dihedral_deg,
                                        double max_dihedral_deg,
                                        double force_constant);
void force_field_add_position_constraint(ForceField &ff, unsigned idx,
                                         double max_displ,
                                         double force_constant);
void force_field_initialize(ForceField &ff);
unsigned force_field_add_extra_point(ForceField &ff, double x, double y,
                                     double z, bool fixed);
Point3D force_field_get_extra_point_pos(const ForceField &ff, unsigned idx);

std::unique_ptr<MMFFMolProperties> mmff_mol_properties_new(
    RWMol &mol, rust::Str variant, MMFFVerbosity verbosity);
bool mmff_mol_properties_is_valid(const MMFFMolProperties &props);
unsigned mmff_mol_properties_get_atom_type(const MMFFMolProperties &props,
                                           unsigned idx);
double mmff_mol_properties_get_formal_charge(const MMFFMolProperties &props,
                                             unsigned idx);
double mmff_mol_properties_get_partial_charge(const MMFFMolProperties &props,
                                              unsigned idx);
TermParams mmff_mol_properties_get_bond_stretch_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2);
TermParams mmff_mol_properties_get_angle_bend_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3);
TermParams mmff_mol_properties_get_stretch_bend_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3);
TermParams mmff_mol_properties_get_torsion_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3, unsigned idx4);
TermParams mmff_mol_properties_get_oop_bend_params(
    const MMFFMolProperties &props, const RWMol &mol, unsigned idx1,
    unsigned idx2, unsigned idx3, unsigned idx4);
TermParams mmff_mol_properties_get_vdw_params(const MMFFMolProperties &props,
                                              unsigned idx1, unsigned idx2);
void mmff_mol_properties_set_dielectric_model(MMFFMolProperties &props,
                                              std::uint8_t model);
void mmff_mol_properties_set_dielectric_constant(MMFFMolProperties &props,
                                                 double constant);
void mmff_mol_properties_set_term(MMFFMolProperties &props, MMFFTerm term,
                                  bool state);
void mmff_mol_properties_set_variant(MMFFMolProperties &props,
                                     rust::Str variant);
void mmff_mol_properties_set_verbosity(MMFFMolProperties &props,
                                       MMFFVerbosity verbosity);

}  // namespace rdkit_cxx
