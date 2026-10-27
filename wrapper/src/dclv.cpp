// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/dclv.h"

#include <cmath>
#include <memory>
#include <string>

#include <GraphMol/Descriptors/DCLV.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/dclv.rs.h"

namespace rdkit_cxx {

namespace {

// The residues DCLV.cpp treats as protein (residueCheck), and those it
// leaves out as solvent (AtomRecord::isSolvent).
bool is_standard_residue(const std::string &name) {
  static const char *const names[] = {
      "ALA", "GLY", "LEU", "SER", "VAL", "THR", "LYS", "ASP", "ILE", "ASN",
      "GLU", "PRO", "ARG", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP"};
  for (const char *n : names) {
    if (name == n) {
      return true;
    }
  }
  return false;
}

bool is_solvent(const std::string &name) {
  static const char *const names[] = {"DOD", "D20", "HOH", "H20", "SOL",
                                      "SO4", "SUL", "WAT", "TIP", "P04"};
  for (const char *n : names) {
    if (name == n) {
      return true;
    }
  }
  return false;
}

void check_finite_within(double value, double bound, const char *what) {
  if (!(std::abs(value) <= bound)) {
    throw ValueErrorException(std::string("DoubleCubicLatticeVolume needs ") +
                              what);
  }
}

// What RDKit 2024.09's constructor (DCLV.cpp) trusts, checked on the
// molecule it computes on, the input without hydrogens:
// - Any monomer info is C-cast to AtomPDBResidueInfo, protein or not, which
//   is undefined for a plain AtomMonomerInfo.
// - getAtomElem reads characters 1 and, for a protein, 2 of each atom's
//   name, past the end of a shorter PDB atom name.
// - The voxel grid converts scaled coordinates and the probe radius to int
//   and unsigned: undefined for values that are not finite or are huge.
//   Within 1e6 the products stay below 1e9.
// - For a protein without its ligand, the grid holds only standard
//   residues; with none, RDKit reads the grid's uninitialized scales for
//   every other atom.
void check_molecule(const RDKit::ROMol &nmol, bool is_protein,
                    bool include_ligand) {
  bool any_member = false;
  bool any_gridded = false;
  for (const auto *atom : nmol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    std::string res_name = "UNK";
    if (info) {
      const auto *pdb = dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info);
      if (!pdb) {
        throw ValueErrorException(
            "DoubleCubicLatticeVolume needs PDB residue info for every atom "
            "with monomer info; atom " +
            std::to_string(atom->getIdx()) + " has other monomer info");
      }
      if (info->getName().size() < (is_protein ? 2u : 1u)) {
        throw ValueErrorException(
            "DoubleCubicLatticeVolume needs PDB atom names of at least " +
            std::string(is_protein ? "two characters" : "one character") +
            "; atom " + std::to_string(atom->getIdx()) + "'s is shorter");
      }
      res_name = pdb->getResidueName();
    }
    if (is_solvent(res_name)) {
      continue;
    }
    any_member = true;
    if (info && is_standard_residue(res_name)) {
      any_gridded = true;
    }
  }
  if (is_protein && !include_ligand && any_member && !any_gridded) {
    throw ValueErrorException(
        "DoubleCubicLatticeVolume without the ligand needs an atom of a "
        "standard residue");
  }
  if (nmol.getNumConformers()) {
    for (const auto &pos : nmol.getConformer().getPositions()) {
      for (double v : {pos.x, pos.y, pos.z}) {
        check_finite_within(v, 1e6, "finite coordinates within 1e6");
      }
    }
  }
}

}  // namespace

DclvResults dclv_compute(const RWMol &mol, bool is_protein, bool include_ligand,
                         double probe_radius, std::int32_t depth,
                         std::int32_t dot_density) {
  // With a dotDensity, RDKit 2024.09 never stores the dots it generates and
  // divides the volumes by a count it never set.
  if (dot_density != 0) {
    throw ValueErrorException(
        "DoubleCubicLatticeVolume's dotDensity must be 0: RDKit 2024.09 "
        "computes from uninitialized values with any other");
  }
  check_finite_within(probe_radius, 1e6, "a finite probe radius within 1e6");
  // RDKit removes the hydrogens the same way first.
  std::unique_ptr<RDKit::ROMol> nmol(RDKit::MolOps::removeAllHs(mol, false));
  // As RDKit's protein default with include_ligand forced on otherwise.
  check_molecule(*nmol, is_protein, is_protein ? include_ligand : true);
  RDKit::Descriptors::DoubleCubicLatticeVolume dclv(
      mol, is_protein, include_ligand, probe_radius, depth, dot_density);
  return DclvResults{dclv.getSurfaceArea(), dclv.getVolume(),
                     dclv.getVDWVolume(), dclv.getCompactness(),
                     dclv.getPackingDensity()};
}

}  // namespace rdkit_cxx
