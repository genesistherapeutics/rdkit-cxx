/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for GraphMol/MolInteractionFields (rdkit.Chem.rdMIF): molecular
// interaction fields, computed on grids by functors that RDKit builds from a
// molecule's conformer, and Gaussian cube files.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolInteractionFields/MIFDescriptors.h>

#include "rdkit-cxx/grid.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using Coulomb = RDMIF::Coulomb;
using CoulombDielectric = RDMIF::CoulombDielectric;
using MMFFVdWaals = RDMIF::MMFFVdWaals;
using UFFVdWaals = RDMIF::UFFVdWaals;
using HBond = RDMIF::HBond;
using Hydrophilic = RDMIF::Hydrophilic;

// A shared struct, defined by the bridges (src/ffi/geometry.rs).
struct Point3D;

std::unique_ptr<UniformRealValueGrid3D> construct_grid(const RWMol &mol,
                                                       std::int32_t conf_id,
                                                       double margin,
                                                       double spacing);

std::unique_ptr<Coulomb> coulomb_new(const RWMol &mol, std::int32_t conf_id,
                                     double probe_charge, bool abs_val,
                                     rust::Str charge_key,
                                     double softcore_param, double cutoff);
std::unique_ptr<Coulomb> coulomb_from_charges(
    rust::Slice<const double> charges, rust::Slice<const Point3D> positions,
    double probe_charge, bool abs_val, double softcore_param, double cutoff);
double coulomb_call(const Coulomb &functor, double x, double y, double z,
                    double threshold);

std::unique_ptr<CoulombDielectric> coulomb_dielectric_new(
    const RWMol &mol, std::int32_t conf_id, double probe_charge, bool abs_val,
    rust::Str charge_key, double softcore_param, double cutoff, double epsilon,
    double xi);
std::unique_ptr<CoulombDielectric> coulomb_dielectric_from_charges(
    rust::Slice<const double> charges, rust::Slice<const Point3D> positions,
    double probe_charge, bool abs_val, double softcore_param, double cutoff,
    double epsilon, double xi);
double coulomb_dielectric_call(const CoulombDielectric &functor, double x,
                               double y, double z, double threshold);

std::unique_ptr<MMFFVdWaals> mmff_vd_waals_new(const RWMol &mol,
                                               std::int32_t conf_id,
                                               std::uint32_t probe_atom_type,
                                               bool scaling, double cutoff);
double mmff_vd_waals_call(const MMFFVdWaals &functor, double x, double y,
                          double z, double threshold);

std::unique_ptr<UFFVdWaals> uff_vd_waals_new(const RWMol &mol,
                                             std::int32_t conf_id,
                                             rust::Str probe_atom_type,
                                             double cutoff);
double uff_vd_waals_call(const UFFVdWaals &functor, double x, double y,
                         double z, double threshold);

std::unique_ptr<HBond> h_bond_new(const RWMol &mol, std::int32_t conf_id,
                                  rust::Str probe_atom_type, bool fixed,
                                  double cutoff);
double h_bond_call(const HBond &functor, double x, double y, double z,
                   double threshold);

std::unique_ptr<Hydrophilic> hydrophilic_new(const RWMol &mol,
                                             std::int32_t conf_id, bool fixed,
                                             double cutoff);
double hydrophilic_call(const Hydrophilic &functor, double x, double y,
                        double z, double threshold);

void calculate_descriptors_coulomb(UniformRealValueGrid3D &grid,
                                   const Coulomb &descriptor, double threshold);
void calculate_descriptors_coulomb_dielectric(
    UniformRealValueGrid3D &grid, const CoulombDielectric &descriptor,
    double threshold);
void calculate_descriptors_mmff_vd_waals(UniformRealValueGrid3D &grid,
                                         const MMFFVdWaals &descriptor,
                                         double threshold);
void calculate_descriptors_uff_vd_waals(UniformRealValueGrid3D &grid,
                                        const UFFVdWaals &descriptor,
                                        double threshold);
void calculate_descriptors_h_bond(UniformRealValueGrid3D &grid,
                                  const HBond &descriptor, double threshold);
void calculate_descriptors_hydrophilic(UniformRealValueGrid3D &grid,
                                       const Hydrophilic &descriptor,
                                       double threshold);

// `mol` is null or a molecule.
void write_to_cube_file(const UniformRealValueGrid3D &grid,
                        rust::Slice<const std::uint8_t> filename,
                        const RWMol *mol, std::int32_t conf_id);
// Reads the grid into `grid`; returns the molecule, or null for none.
std::unique_ptr<RWMol> read_from_cube_file(
    rust::Slice<const std::uint8_t> filename, UniformRealValueGrid3D &grid);

}  // namespace rdkit_cxx
