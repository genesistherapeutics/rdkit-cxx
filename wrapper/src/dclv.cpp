// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/dclv.h"

#include <cmath>
#include <memory>
#include <string>

#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/dclv.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

// The largest coordinate and probe radius the checks below allow. RDKit
// 2024.03's voxel grid (DCLV.cpp, State::findNeighbours) converts
// voxX * (x - range - minX) to int, with voxX at most 16 / 0.1 and range the
// probe's diameter plus two atomic radii, and createVoxelGrid converts the
// scaled coordinates to the unsigned indices of its 16x16x16 grid. Within
// these bounds every value fits, and the grid's spread keeps its +0.1.
constexpr double kMaxCoordinate = 1e6;
constexpr double kMaxProbeRadius = 1e5;

void check_finite_within(double value, double bound, const char *what) {
  if (!(std::abs(value) <= bound)) {
    throw ValueErrorException(std::string("DoubleCubicLatticeVolume needs ") +
                              what);
  }
}

// AtomRecord::isSolvent: atoms of these residues are left out.
bool is_solvent(const std::string &res_name) {
  for (const char *name :
       {"DOD", "D20", "HOH", "H20", "SOL", "SO4", "SUL", "WAT", "TIP", "P04"}) {
    if (res_name == name) {
      return true;
    }
  }
  return false;
}

// The 20 residues AtomRecord does not flag as HETATM.
bool is_standard_residue(const std::string &res_name) {
  for (const char *name :
       {"ALA", "GLY", "LEU", "SER", "VAL", "THR", "LYS", "ASP", "ILE", "ASN",
        "GLU", "PRO", "ARG", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP"}) {
    if (res_name == name) {
      return true;
    }
  }
  return false;
}

// What DCLV.cpp trusts of the molecule it computes on, the input without
// hydrogens:
// - AtomRecord C-casts any atom's monomer info to AtomPDBResidueInfo, which
//   is undefined for a plain AtomMonomerInfo;
// - getAtomElem reads characters 1 (and, for a protein, 2) of each kept
//   atom's PDB name through c_str(), past the end of a shorter name;
// - every coordinate is converted to an integer grid index (above);
// - with isProtein and not includeLigand, the grid holds only the atoms of
//   the 20 standard residues, and when there are none createVoxelGrid
//   returns before it sets the grid's scale, which calculateVolume then
//   reads for every atom.
void check_molecule(const RDKit::ROMol &mol, bool is_protein,
                    bool include_ligand) {
  if (mol.getNumConformers() == 0) {
    return;  // getConformer throws
  }
  const auto &conf = mol.getConformer();
  bool kept = false;
  bool gridded = false;
  for (const auto *atom : mol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    bool het = true;
    if (info) {
      const auto *pdb = dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info);
      if (!pdb) {
        throw ValueErrorException(
            "DoubleCubicLatticeVolume needs PDB residue info for every atom "
            "with monomer info; atom " +
            std::to_string(atom->getIdx()) + " has other monomer info");
      }
      if (is_solvent(pdb->getResidueName())) {
        continue;
      }
      std::size_t needed = is_protein ? 2 : 1;
      if (info->getName().size() < needed) {
        throw ValueErrorException(
            "DoubleCubicLatticeVolume needs PDB atom names of at least " +
            std::to_string(needed) + " characters; atom " +
            std::to_string(atom->getIdx()) + "'s is shorter");
      }
      het = !is_standard_residue(pdb->getResidueName());
    }
    kept = true;
    gridded = gridded || !is_protein || include_ligand || !het;
    const auto &pos = conf.getAtomPos(atom->getIdx());
    for (double v : {pos.x, pos.y, pos.z}) {
      check_finite_within(v, kMaxCoordinate, "finite coordinates within 1e6");
    }
  }
  if (kept && !gridded) {
    throw ValueErrorException(
        "with isProtein and not includeLigand, DoubleCubicLatticeVolume needs "
        "an atom of a standard residue");
  }
}

}  // namespace

DoubleCubicLatticeVolumeData dclv_compute(const RWMol &mol, bool is_protein,
                                          bool include_ligand,
                                          double probe_radius,
                                          std::int32_t depth,
                                          std::int32_t dot_density) {
  // With a dot density, generateSurfacePoints never sizes the standard dots,
  // and calculateVolume divides by their count, which is then indeterminate.
  if (dot_density != 0) {
    throw ValueErrorException(
        "DoubleCubicLatticeVolume in RDKit 2024.03 reads an uninitialized "
        "dot count for a nonzero dotDensity; use 0");
  }
  check_finite_within(probe_radius, kMaxProbeRadius,
                      "a finite probe radius within 1e5");
  // As the constructor does first, so the checks see the atoms it keeps.
  check_removals(mol);
  std::unique_ptr<RDKit::ROMol> heavy(RDKit::MolOps::removeAllHs(mol, false));
  check_molecule(*heavy, is_protein, !is_protein || include_ligand);
  RDKit::Descriptors::DoubleCubicLatticeVolume dclv(
      *heavy, is_protein, include_ligand, probe_radius, depth, dot_density);
  return DoubleCubicLatticeVolumeData{
      dclv.getSurfaceArea(), dclv.getVolume(), dclv.getVDWVolume(),
      dclv.getCompactness(), dclv.getPackingDensity()};
}

}  // namespace rdkit_cxx
