/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of whole molecules, for RDKit code that trusts them.

#pragma once

#include <GraphMol/ROMol.h>

namespace rdkit_cxx {

// RDKit caches a molecule's distance, adjacency, and path matrices as unsized
// arrays in computed properties (Matrices.cpp), and its graph edits do not
// all clear them (RWMol::addAtom, addBond, removeBond, insertMol), so a stale
// matrix is read past its end. No molecule Rust holds keeps one between calls
// (CONVENTIONS.md U6).

// Removes every cached matrix from the molecule's properties.
void clear_matrix_caches(const RDKit::ROMol &mol);

// Clears the molecule's cached matrices when it is made and again when it is
// destroyed, so an RDKit call made in its scope neither reads a stale matrix
// nor leaves one behind, even when it throws.
class MatrixCacheGuard {
 public:
  explicit MatrixCacheGuard(const RDKit::ROMol &mol) : d_mol(mol) {
    clear_matrix_caches(d_mol);
  }
  ~MatrixCacheGuard() { clear_matrix_caches(d_mol); }
  MatrixCacheGuard(const MatrixCacheGuard &) = delete;
  MatrixCacheGuard &operator=(const MatrixCacheGuard &) = delete;

 private:
  const RDKit::ROMol &d_mol;
};

// RDKit computes a molecule's n*n matrices in int, which overflows above
// 46340 atoms. Throws ValueErrorException for a larger molecule (U3).
void check_matrix_size(const RDKit::ROMol &mol);

}  // namespace rdkit_cxx
