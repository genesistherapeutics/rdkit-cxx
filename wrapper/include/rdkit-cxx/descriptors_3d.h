/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDescriptors's 3D descriptors. The moment descriptors read and write
// the descriptor caches, and the others the distance matrices, so their
// wrappers clear both around each call (CONVENTIONS.md U6, U8).

#pragma once

#include <cstdint>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/descriptors_3d.rs,
// geometry.rs).
struct DoubleList;
struct IndexList;
struct Point3D;

double calc_pmi1(const RWMol &mol, std::int32_t conf_id,
                 bool use_atomic_masses);
double calc_pmi2(const RWMol &mol, std::int32_t conf_id,
                 bool use_atomic_masses);
double calc_pmi3(const RWMol &mol, std::int32_t conf_id,
                 bool use_atomic_masses);
double calc_npr1(const RWMol &mol, std::int32_t conf_id,
                 bool use_atomic_masses);
double calc_npr2(const RWMol &mol, std::int32_t conf_id,
                 bool use_atomic_masses);
double calc_radius_of_gyration(const RWMol &mol, std::int32_t conf_id,
                               bool use_atomic_masses);
double calc_inertial_shape_factor(const RWMol &mol, std::int32_t conf_id,
                                  bool use_atomic_masses);
double calc_eccentricity(const RWMol &mol, std::int32_t conf_id,
                         bool use_atomic_masses);
double calc_asphericity(const RWMol &mol, std::int32_t conf_id,
                        bool use_atomic_masses);
double calc_spherocity_index(const RWMol &mol, std::int32_t conf_id);
double calc_pbf(const RWMol &mol, std::int32_t conf_id);

// An empty custom_atom_property means RDKit's property tables.
rust::Vec<double> calc_autocorr_3d(const RWMol &mol, std::int32_t conf_id,
                                   rust::Str custom_atom_property);
rust::Vec<double> calc_rdf(const RWMol &mol, std::int32_t conf_id,
                           rust::Str custom_atom_property);
rust::Vec<double> calc_morse(const RWMol &mol, std::int32_t conf_id,
                             rust::Str custom_atom_property);
rust::Vec<double> calc_whim(const RWMol &mol, std::int32_t conf_id,
                            double thresh, rust::Str custom_atom_property);
rust::Vec<double> calc_getaway(const RWMol &mol, std::int32_t conf_id,
                               double precision,
                               rust::Str custom_atom_property);
// Row-major, numAtoms by numAtoms.
rust::Vec<double> calc_coulomb_mat(const RWMol &mol, std::int32_t conf_id);
rust::Vec<double> calc_eem_charges(const RWMol &mol, std::int32_t conf_id);

// The USR descriptors. GetUSRCAT's atom selections number atoms from 1, as
// in Python, whose wrapper subtracts 1.
rust::Vec<double> get_usr(const RWMol &mol, std::int32_t conf_id);
rust::Vec<double> get_usrcat(const RWMol &mol, bool has_atom_selections,
                             rust::Slice<const IndexList> atom_selections,
                             std::int32_t conf_id);
// The four reference points go to `points`, as Python's points list.
rust::Vec<DoubleList> get_usr_distributions(rust::Slice<const Point3D> coords,
                                            rust::Vec<Point3D> &points);
rust::Vec<DoubleList> get_usr_distributions_from_points(
    rust::Slice<const Point3D> coords, rust::Slice<const Point3D> points);
rust::Vec<double> get_usr_from_distributions(
    rust::Slice<const DoubleList> distances);
double get_usr_score(rust::Slice<const double> descriptor1,
                     rust::Slice<const double> descriptor2,
                     rust::Slice<const double> weights);

}  // namespace rdkit_cxx
