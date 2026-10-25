/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of what RDKit's 2D depiction and CoordGen trust (CONVENTIONS.md
// U1). Each throws ValueErrorException. The callers hold lock_depiction(),
// so the depiction preference cannot change between a check and the call.

#pragma once

#include <GraphMol/Conformer.h>
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
