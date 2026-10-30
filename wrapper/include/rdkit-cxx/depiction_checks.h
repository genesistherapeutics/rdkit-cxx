/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of what RDKit's 2D depiction and CoordGen trust (CONVENTIONS.md
// U1). Each throws ValueErrorException. The callers hold lock_depiction(),
// so the depiction preference cannot change between a check and the call.

#pragma once

#include <GraphMol/Conformer.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/ROMol.h>

namespace rdkit_cxx {

// The largest coordinate magnitude a depiction function accepts. RDKit sorts
// angles and distances it derives from coordinates (std::sort in
// EmbeddedFrag.cpp, with a pair comparison that NaN makes inconsistent) and
// NormalizeDepiction casts bond lengths times 10 to int; coordinates within
// this bound keep every value it derives finite and in range.
constexpr double kMaxDepictionCoordinate = 1e5;

// CoordGen computes in float and squares sums of products of its
// coordinates (sketcherMinimizer::svd). The bound is on coordinates as
// CoordGen receives them, times coordgenScaling: kMaxDepictionCoordinate at
// the default scaling of 50.
constexpr double kMaxCoordGenCoordinate = 5e6;

// Throws unless v is finite and within kMaxDepictionCoordinate.
void check_depiction_coordinate(double v, const char *what);

// Throws unless every coordinate of the conformer passes
// check_depiction_coordinate.
void check_depiction_conformer(const RDKit::Conformer &conf, const char *what);

// A bond length the depictor scales every coordinate by: throws unless a
// positive one passes check_depiction_coordinate (zero and negative lengths
// keep RDKit's).
void check_bond_length(double bond_length);

// Python's Compute2DCoords, Compute2DCoordsMimicDistmat, and
// Compute2DCoordsForReaction set RDDepict::BOND_LEN for one call when
// bondLength > 0. This restores it even when RDKit throws. Its owner holds
// lock_depiction().
class BondLengthOverride {
 public:
  explicit BondLengthOverride(double bond_length)
      : d_old(RDDepict::BOND_LEN), d_set(bond_length > 0) {
    if (d_set) {
      RDDepict::BOND_LEN = bond_length;
    }
  }
  ~BondLengthOverride() {
    if (d_set) {
      RDDepict::BOND_LEN = d_old;
    }
  }
  BondLengthOverride(const BondLengthOverride &) = delete;
  BondLengthOverride &operator=(const BondLengthOverride &) = delete;

 private:
  double d_old;
  bool d_set;
};

// The conformer getConformer(id) returns, or null where it would throw.
const RDKit::Conformer *find_conformer(const RDKit::ROMol &mol, int id);

// What RDKit's own depictor (computeInitialCoords and EmbeddedFrag) trusts
// about a molecule.
void check_rdkit_depictor_mol(const RDKit::ROMol &mol);

// What CoordGen trusts about a molecule.
void check_coord_gen_mol(const RDKit::ROMol &mol,
                         bool treat_nonterminal_bonds_to_metal_as_zobs);

// Whether RDDepict::compute2DCoords would hand the molecule to CoordGen.
bool depicts_with_coord_gen(bool force_rdkit);

// check_coord_gen_mol or check_rdkit_depictor_mol, for whichever depictor
// RDDepict::compute2DCoords(mol, ..., force_rdkit) runs.
void check_compute_2d_coords_mol(const RDKit::ROMol &mol, bool force_rdkit);

}  // namespace rdkit_cxx
