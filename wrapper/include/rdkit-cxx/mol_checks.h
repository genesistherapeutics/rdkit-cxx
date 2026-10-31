/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of whole molecules, for RDKit code that trusts them.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

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

// What a molecule pickle, which MolPickler reads without these checks, can
// give a molecule that RDKit or Rust trusts:
// - RDKit takes a monomer info's type for its class and casts an info of
//   type PDBRESIDUE to AtomPDBResidueInfo, so only one of those may have it
//   (CONVENTIONS.md U10). A pickle can give an atom a plain info of that
//   type, or a type AtomMonomerType lacks.
// - MolPickler casts the bytes it reads for atoms' chiral tags and
//   hybridizations, bonds' types, directions, and stereo, and stereo
//   groups' types to the enums unchecked, and the Rust enums have no
//   variant for a value RDKit's lacks (U2).
// Throws ValueErrorException for any of these, in the molecule or in the
// molecules of its recursive queries, which RDKit pickles as molecules.
// Every wrapper that reads a molecule pickle, alone or inside another
// format, calls it on the molecules it read.
void check_pickled_mol(const RDKit::ROMol &mol);

// MolPickler writes some counts in one byte and truncates larger ones without
// a word, so it writes a pickle it cannot read back: the reader takes the
// rest of the list for the pickle's next fields, which can be bytes the
// caller chose (a bracket's coordinates), and parses them unchecked (U2):
// - every query node's number of children (pickleQuery), for atom and bond
//   queries in a molecule of any size;
// - in a molecule with at most 255 atoms and at most 255 bonds, which it
//   pickles with one-byte indices: each substance group's numbers of atoms,
//   parent atoms, bonds, brackets, CStates, and attachment points
//   (_pickleSubstanceGroup), and the number of stereo groups and each one's
//   numbers of atoms and bonds (_pickleStereo).
// It also counts the properties of the molecule, each conformer, atom, bond,
// and substance group in two bytes (streamWriteProps), wrapping past 65535.
// Throws ValueErrorException for a count past 255 (or 65535) in any of these,
// in the molecule or in the molecules of its recursive queries, which
// MolPickler pickles as molecules with their own one-byte rule. Every wrapper
// that hands a molecule to MolPickler, alone or inside another format, calls
// it on the molecules it will pickle (CONVENTIONS.md U2).
void check_picklable(const RDKit::ROMol &mol);

// check_picklable on a reaction's reactant, product, and agent templates,
// which ReactionPickler pickles as molecules.
template <class Reaction>
void check_picklable_templates(const Reaction &rxn) {
  for (const auto *templates :
       {&rxn.getReactants(), &rxn.getProducts(), &rxn.getAgents()}) {
    for (const auto &mol : *templates) {
      check_picklable(*mol);
    }
  }
}

// As check_picklable, for the parts FragCatGenerator::addFragsFromMol makes
// of the molecule and pickles in its catalog's entries: subgraphs of at most
// max_part_bonds bonds, made with copyMolSubset, which keeps a substance
// group only when the part holds all its atoms and bonds, and, from RDKit
// 2025.09 on (copySelectedStereoGroups in Subset.cpp), every stereo group
// that shares an atom (and a bond, if it has bonds) with the part, with the
// members the part holds. So a part has at most 255 atoms and bonds whenever
// max_part_bonds is below 255, and a small part can come from a molecule of
// any size: the one-byte limits apply to the molecule's own groups whatever
// its size, except to substance groups whose distinct atoms or bonds no part
// has room for; every stereo group counts, with at most a part's atoms and
// bonds. Recursive queries' molecules are pickled whole, and checked as
// check_picklable checks them.
void check_picklable_parts(const RDKit::ROMol &mol, unsigned max_part_bonds);

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

// As check_removals, for an RDKit call that combines the molecules
// (insertMol, combineMols) and then removes atoms or bonds: RDKit does not
// renumber the substance groups' indices, so the combination's PARENT
// chains can run through the groups of several molecules. Which group an
// index finds depends on which molecules RDKit combines and in what order
// (molzip leaves out fragments whose SMILES repeat), so this follows every
// group with the PARENT's index, and refuses a cycle or a chain longer than
// MAX_PARENT_CHAIN through any of them.
void check_removals(const std::vector<const RDKit::ROMol *> &mols);

// RWMol::replaceBond deletes the bond it replaces but, unlike replaceAtom,
// leaves the molecule's stereo groups (an atropisomer's group holds its bond)
// pointing at it, a pointer that copying, pickling, stereo perception, and
// the CXSMILES and mol-block writers then read (CONVENTIONS.md U13).
//
// RDKit's own calls of it that the wrappers reach, where a group can hold
// the bond and be read after, are guarded: adjustQueryProperties
// (adjust_query.cpp, generalized_substruct.cpp, reactions' adjust_templates,
// and the depictor's allowRGroups) and TautomerQuery::fromMol (replaced by
// tautomer_query.cpp's own). MolInterchange's reader is refused such a
// document (mol_interchange.cpp). The others leave no pointer that is read:
// MolFragmenter replaces bonds it has just added; ReactionRunner's products'
// groups hold no bonds; the mol-file parser replaces coordinate bonds before
// stereo perception puts bonds in groups; the adjustments of R-group
// decomposition's prepareCore and of the abbreviations replace no bonds; and
// the copies that symmetrizeTerminalAtoms (in the embedder's pruning) and
// R-group decomposition (its matching core, and its tautomer query's
// template) replace bonds in are only matched, which reads groups' atoms.

// As RWMol::replaceBond, and then points the stereo groups that held the old
// bond at the new one, as replaceAtom does for atoms.
void replace_bond(RDKit::RWMol &mol, unsigned idx, RDKit::Bond *bond,
                  bool preserve_props, bool keep_s_groups);

// Whether any of the molecule's stereo groups holds a bond.
bool has_stereo_group_bonds(const RDKit::ROMol &mol);

// Takes the bonds out of the molecule's stereo groups, for a copy that an
// RDKit call replaces bonds in and that nothing returns afterwards.
void clear_stereo_group_bonds(RDKit::ROMol &mol);

// Takes the bonds out of the molecule's stereo groups when it is made, for
// an RDKit call that replaces bonds without renumbering bonds or stereo
// groups (adjustQueryProperties), and when it is destroyed, even when the
// call throws, puts back into each group the bonds now at the indices it
// held, as replace_bond does.
class StereoGroupBondGuard {
 public:
  explicit StereoGroupBondGuard(RDKit::ROMol &mol);
  ~StereoGroupBondGuard();
  StereoGroupBondGuard(const StereoGroupBondGuard &) = delete;
  StereoGroupBondGuard &operator=(const StereoGroupBondGuard &) = delete;

 private:
  RDKit::ROMol &d_mol;
  // Each group's bond indices; empty when no group held a bond.
  std::vector<std::vector<unsigned>> d_bonds;
};

}  // namespace rdkit_cxx
