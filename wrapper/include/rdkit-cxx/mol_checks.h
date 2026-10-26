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

// RDKit takes a monomer info's type for its class and casts an info of type
// PDBRESIDUE to AtomPDBResidueInfo, so only one of those may have it
// (CONVENTIONS.md U10). A pickle can give an atom a plain info of that type,
// or a type AtomMonomerType lacks; throws ValueErrorException for either,
// in the molecule or in the molecules of its recursive queries, which RDKit
// pickles as molecules.
void check_monomer_types(const RDKit::ROMol &mol);

// RWMol::removeAtom and batchRemoveAtoms (so commitBatchEdit too) rewrite
// every bond's _MolFileBondEndPts, a property any caller can set: they read
// the value's front() without checking that it is not empty, and for one in
// parentheses they read the first space-separated token without checking
// that there is one and convert each token's std::stod to unsigned int,
// which is undefined for a value that does not fit. Throws
// ValueErrorException for a value they would read so, or one with a token
// std::stod rejects, which they would throw for halfway through the edit.
void check_bond_end_points(const RDKit::ROMol &mol);

// When RWMol removes an atom or bond it drops the substance groups that hold
// it, and those whose PARENT chain leads to a dropped one: it follows each
// group's PARENT to the group with that "index", both properties any caller
// can set, recursively and without noting where it has been, so a cycle never
// ends. Throws ValueErrorException for a cycle, or a chain longer than
// MAX_PARENT_CHAIN, and for a PARENT or index that is not an unsigned int as
// RDKit throws midway through the edit.
constexpr unsigned MAX_PARENT_CHAIN = 4000;
void check_substance_group_parents(const RDKit::ROMol &mol);

// What RWMol trusts in a molecule when a wrapper's RDKit call removes atoms or
// bonds from it (or a copy): check_bond_end_points and
// check_substance_group_parents.
void check_removals(const RDKit::ROMol &mol);

}  // namespace rdkit_cxx
