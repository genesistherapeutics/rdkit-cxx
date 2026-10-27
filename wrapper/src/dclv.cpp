// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/dclv.h"

#include <cmath>
#include <string>

#include <GraphMol/Descriptors/DCLV.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/dclv.rs.h"

namespace rdkit_cxx {

namespace {

// findNeighbours converts each atom's position, less the grid's corner and
// the search radius, scaled by up to 160 voxels per unit, to int, and
// createVoxelGrid scaled positions to unsigned voxel indices: undefined for
// values that are not finite, or beyond int. An atom outside the grid's
// box is searched too (HETATM atoms that includeLigand leaves out of the
// grid), so every coordinate is bounded, and the probe radius with it:
// within 1e6, every scaled value stays under 1e9.
void check_finite_within(double value, double bound, const char *what) {
  if (!(std::abs(value) <= bound)) {
    throw ValueErrorException(std::string("DoubleCubicLatticeVolume needs ") +
                              what);
  }
}

// AtomRecord's constructor C-casts any monomer info to AtomPDBResidueInfo,
// which is undefined for a plain AtomMonomerInfo, and getAtomElem reads the
// atom name's second character, and with isProtein its third, without
// checking its length.
void check_monomer_infos(const RWMol &mol, bool is_protein) {
  const std::size_t min_name = is_protein ? 2 : 1;
  for (const auto *atom : mol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    if (!info) {
      continue;
    }
    if (!dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info)) {
      throw ValueErrorException(
          "DoubleCubicLatticeVolume needs PDB residue info for every atom "
          "with monomer info; atom " +
          std::to_string(atom->getIdx()) + " has other monomer info");
    }
    if (info->getName().size() < min_name) {
      throw ValueErrorException(
          "DoubleCubicLatticeVolume needs atom names of at least " +
          std::to_string(min_name) + " characters" +
          (is_protein ? " with isProtein" : "") + "; atom " +
          std::to_string(atom->getIdx()) + "'s is shorter");
    }
  }
}

}  // namespace

DclvResults dclv_compute(const RWMol &mol, bool is_protein, bool include_ligand,
                         double probe_radius, std::int32_t depth,
                         std::int32_t dot_density) {
  // With a dot density, generateElemPoints computes each element's dots
  // and discards them, and calculateVolume divides by the standard dots'
  // count, which only a dot density of 0 initializes; a negative density
  // converts a negative dot count to unsigned long.
  if (dot_density != 0) {
    throw ValueErrorException(
        "DoubleCubicLatticeVolume's dotDensity must be 0: RDKit 2025.03 "
        "reads an uninitialized count for any other");
  }
  // A negative probe radius takes RDKit's default; NaN is kept.
  if (!(probe_radius < 0.0)) {
    check_finite_within(probe_radius, 1e6, "a finite probe radius within 1e6");
  }
  // RDKit measures the default conformer, of the molecule without its
  // hydrogens, whose positions are the same.
  if (mol.getNumConformers()) {
    for (const auto &pos : mol.getConformer().getPositions()) {
      for (double v : {pos.x, pos.y, pos.z}) {
        check_finite_within(v, 1e6, "finite coordinates within 1e6");
      }
    }
  }
  check_monomer_infos(mol, is_protein);
  RDKit::Descriptors::DoubleCubicLatticeVolume dclv(
      mol, is_protein, include_ligand, probe_radius, depth, dot_density);
  DclvResults results;
  results.surface_area = dclv.getSurfaceArea();
  results.volume = dclv.getVolume();
  results.vdw_volume = dclv.getVDWVolume();
  results.compactness = dclv.getCompactness();
  results.packing_density = dclv.getPackingDensity();
  return results;
}

}  // namespace rdkit_cxx
