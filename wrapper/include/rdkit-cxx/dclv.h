/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDescriptors's DoubleCubicLatticeVolume, which keeps a reference to
// its molecule: the Rust type borrows the molecule for its lifetime.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/Descriptors/DCLV.h>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using DoubleCubicLatticeVolume = RDKit::Descriptors::DoubleCubicLatticeVolume;

struct Point3D;

// Empty radii are RDKit's van der Waals radii, as Python's empty list is.
std::unique_ptr<DoubleCubicLatticeVolume> dclv_new(
    const RWMol &mol, rust::Slice<const double> radii, bool is_protein,
    bool include_ligand, double probe_radius, std::int32_t conf_id);
double dclv_get_surface_area(DoubleCubicLatticeVolume &dclv);
double dclv_get_atom_surface_area(DoubleCubicLatticeVolume &dclv,
                                  std::uint32_t atom_idx);
double dclv_get_polar_surface_area(DoubleCubicLatticeVolume &dclv,
                                   bool include_s_and_p, bool include_hs);
double dclv_get_partial_surface_area(
    DoubleCubicLatticeVolume &dclv,
    rust::Slice<const std::uint32_t> atom_indices);
// Each atom with points, its number of points, and the points in order.
void dclv_get_surface_points(DoubleCubicLatticeVolume &dclv,
                             rust::Vec<std::uint32_t> &atoms,
                             rust::Vec<std::uint32_t> &counts,
                             rust::Vec<Point3D> &points);
double dclv_get_volume(DoubleCubicLatticeVolume &dclv);
double dclv_get_vdw_volume(DoubleCubicLatticeVolume &dclv);
double dclv_get_atom_volume(DoubleCubicLatticeVolume &dclv,
                            std::uint32_t atom_idx, double solvent_radius);
double dclv_get_polar_volume(DoubleCubicLatticeVolume &dclv,
                             bool include_s_and_p, bool include_hs);
double dclv_get_partial_volume(DoubleCubicLatticeVolume &dclv,
                               rust::Slice<const std::uint32_t> atom_indices);
double dclv_get_compactness(DoubleCubicLatticeVolume &dclv);
double dclv_get_packing_density(DoubleCubicLatticeVolume &dclv);

}  // namespace rdkit_cxx
