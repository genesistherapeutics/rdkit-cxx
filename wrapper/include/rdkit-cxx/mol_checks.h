/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of whole molecules, for RDKit code that trusts them.

#pragma once

#include <cstdint>

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

// RDKit stores each atom's formal charge in 8 bits but sums charges over
// atoms in int, which only 2^24 or more atoms can overflow (CONVENTIONS.md
// U9). Throws ValueErrorException unless the charges' magnitudes, plus
// headroom, sum to at most INT_MAX: then every sum of the charges over a
// subset of the atoms, its negation, and its difference with a count of at
// most headroom fit in int.
void check_formal_charge_sum(const RDKit::ROMol &mol,
                             std::int64_t headroom = 0);

// RDKit's descriptors cache their results, and per-atom contributions, in a
// molecule's computed properties, and read them back unless forced, after
// edits and for other conformers; the Chi descriptors index a cached per-atom
// vector by atom without checking its length. No molecule Rust holds keeps
// one between calls (CONVENTIONS.md U8).

// Removes every descriptor cache from the molecule's properties.
void clear_descriptor_caches(const RDKit::ROMol &mol);

// As MatrixCacheGuard, for the descriptor caches.
class DescriptorCacheGuard {
 public:
  explicit DescriptorCacheGuard(const RDKit::ROMol &mol) : d_mol(mol) {
    clear_descriptor_caches(d_mol);
  }
  ~DescriptorCacheGuard() { clear_descriptor_caches(d_mol); }
  DescriptorCacheGuard(const DescriptorCacheGuard &) = delete;
  DescriptorCacheGuard &operator=(const DescriptorCacheGuard &) = delete;

 private:
  const RDKit::ROMol &d_mol;
};

}  // namespace rdkit_cxx
