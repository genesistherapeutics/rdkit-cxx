// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/depiction_checks.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/Chirality.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDConfig.h>
#include <coordgen/sketcherMinimizerAtom.h>

#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"

namespace rdkit_cxx {

namespace {

constexpr std::int64_t kIntMax = std::numeric_limits<int>::max();

// CoordGen::addCoords gives ZERO bonds order 0, DOUBLE 2, TRIPLE 3, and every
// other bond 1.
int coord_gen_bond_order(const RDKit::Bond *bond) {
  switch (bond->getBondType()) {
    case RDKit::Bond::ZERO:
      return 0;
    case RDKit::Bond::DOUBLE:
      return 2;
    case RDKit::Bond::TRIPLE:
      return 3;
    default:
      return 1;
  }
}

// sketcherMinimizer::canonicalOrdering orders atoms by Morgan scores it sums
// in int (sketcherMinimizer::morganScores), which overflow for long chains:
// a chain of 40 carbons. This repeats the sums in 64 bits on the graph
// CoordGen sees, which sketcherMinimizer::initialize gets by dropping
// zero-order bonds.
void check_coord_gen_morgan_scores(const RDKit::ROMol &mol, bool zobs) {
  const auto n = mol.getNumAtoms();
  if (n < 2) {
    return;
  }
  std::vector<std::pair<unsigned, unsigned>> bonds;
  for (const auto bond : mol.bonds()) {
    auto order = coord_gen_bond_order(bond);
    const auto begin = bond->getBeginAtom();
    const auto end = bond->getEndAtom();
    // With treatNonterminalBondsToMetalAsZOBs, initialize makes a single or
    // double bond to a metal zero-order unless either atom has one bond.
    if (zobs && (order == 1 || order == 2) && begin->getDegree() != 1 &&
        end->getDegree() != 1 &&
        (sketcherMinimizerAtom::isMetal(begin->getAtomicNum()) ||
         sketcherMinimizerAtom::isMetal(end->getAtomicNum()))) {
      order = 0;
    }
    if (order != 0) {
      bonds.emplace_back(begin->getIdx(), end->getIdx());
    }
  }
  const auto overflow = [] {
    return ValueErrorException(
        "CoordGen sums Morgan scores in int, which overflow for this molecule "
        "(a chain of about 40 atoms or more); use RDKit's depictor");
  };
  std::vector<std::int64_t> old_scores(n, 1);
  std::vector<std::int64_t> new_scores(n, 0);
  std::size_t old_ties = n;
  while (true) {
    for (const auto &[a, b] : bonds) {
      new_scores[a] += old_scores[b];
      new_scores[b] += old_scores[a];
      if (new_scores[a] > kIntMax || new_scores[b] > kIntMax) {
        throw overflow();
      }
    }
    auto ordered = new_scores;
    std::sort(ordered.begin(), ordered.end());
    std::size_t ties = 0;
    for (std::size_t j = 1; j < ordered.size(); ++j) {
      if (ordered[j] == ordered[j - 1]) {
        ++ties;
      }
    }
    if (ties >= old_ties) {
      break;
    }
    old_ties = ties;
    old_scores = new_scores;
  }
  // canonicalOrdering then ranks each atom by score * 100 + atomic number.
  for (const auto atom : mol.atoms()) {
    if (old_scores[atom->getIdx()] * 100 + atom->getAtomicNum() > kIntMax) {
      throw overflow();
    }
  }
}

}  // namespace

void check_depiction_coordinate(double v, const char *what) {
  if (!(std::abs(v) <= kMaxDepictionCoordinate)) {
    throw ValueErrorException(std::string(what) +
                              " must be finite and within 1e5 for 2D "
                              "depiction");
  }
}

void check_depiction_conformer(const RDKit::Conformer &conf, const char *what) {
  for (const auto &p : conf.getPositions()) {
    for (double v : {p.x, p.y, p.z}) {
      check_depiction_coordinate(v, what);
    }
  }
}

void check_bond_length(double bond_length) {
  if (bond_length > 0) {
    check_depiction_coordinate(bond_length, "bondLength");
  }
}

const RDKit::Conformer *find_conformer(const RDKit::ROMol &mol, int id) {
  if (mol.getNumConformers() == 0) {
    return nullptr;
  }
  if (id < 0) {
    return &mol.getConformer();
  }
  for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
    if ((*it)->getId() == static_cast<unsigned>(id)) {
      return it->get();
    }
  }
  return nullptr;
}

void check_rdkit_depictor_mol(const RDKit::ROMol &mol) {
  // The depictor removes collisions with RDKit's distance matrix
  // (MolOps::getDistanceMat), which it computes in int.
  check_matrix_size(mol);
  // computeInitialCoords symmetrizes the SSSR with dative bonds.
  check_relevant_cycles(mol, true);
  const std::int64_t n = mol.getNumAtoms();
  for (const auto atom : mol.atoms()) {
    // computeInitialCoords ranks each atom it starts a fragment from as
    // getAtomDepictRank * numAtoms + index, in int, and dereferences the
    // atom of the lowest rank below INT_MAX, past the end if there is none.
    if (RDDepict::getAtomDepictRank(atom) * n + (n - 1) >= kIntMax) {
      throw ValueErrorException(
          "RDKit's depictor ranks atoms in int, which overflows for this "
          "many atoms with explicit hydrogens");
    }
    switch (atom->getChiralTag()) {
      case RDKit::Atom::CHI_SQUAREPLANAR:
        // embedSquarePlanar reads the first neighbor unchecked.
        if (atom->getDegree() == 0) {
          throw ValueErrorException(
              "RDKit's depictor needs a neighbor for each square-planar "
              "atom (atom " +
              std::to_string(atom->getIdx()) + ")");
        }
        break;
      case RDKit::Atom::CHI_TRIGONALBIPYRAMIDAL: {
        // embedTBP places the atoms that are not axial at three equatorial
        // positions, read from a table it indexes unchecked.
        const auto axial1 =
            RDKit::Chirality::getTrigonalBipyramidalAxialAtom(atom);
        const auto axial2 =
            RDKit::Chirality::getTrigonalBipyramidalAxialAtom(atom, -1);
        unsigned equatorial = 0;
        for (const auto nbr : mol.atomNeighbors(atom)) {
          if (nbr != axial1 && nbr != axial2) {
            ++equatorial;
          }
        }
        if (equatorial > 3) {
          throw ValueErrorException(
              "RDKit's depictor places at most three equatorial neighbors "
              "around a trigonal-bipyramidal atom (atom " +
              std::to_string(atom->getIdx()) + ")");
        }
        break;
      }
      default:
        break;
    }
  }
}

void check_coord_gen_mol(const RDKit::ROMol &mol,
                         bool treat_nonterminal_bonds_to_metal_as_zobs) {
  // sketcherMinimizerMolecule::totalCharge sums each fragment's charges in
  // int (U9).
  check_formal_charge_sum(mol);
  check_coord_gen_morgan_scores(mol, treat_nonterminal_bonds_to_metal_as_zobs);
}

bool depicts_with_coord_gen(bool force_rdkit) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return !force_rdkit && RDDepict::preferCoordGen;
#else
  (void)force_rdkit;
  return false;
#endif
}

void check_compute_2d_coords_mol(const RDKit::ROMol &mol, bool force_rdkit) {
  if (depicts_with_coord_gen(force_rdkit)) {
    // compute2DCoords gives CoordGen its default parameters.
    check_coord_gen_mol(mol, false);
  } else {
    check_rdkit_depictor_mol(mol);
  }
}

}  // namespace rdkit_cxx
