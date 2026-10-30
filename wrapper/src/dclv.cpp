// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/dclv.h"

#include <cmath>
#include <mutex>
#include <string>
#include <vector>

#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/Exceptions.h>
#include <boost/dynamic_bitset.hpp>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/dclv.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

// DCLV.cpp keeps the atom its last surface test hit in a global int,
// recordCache, which every object reads and writes without a lock (U4), and
// which one object's call can leave pointing past another's atoms. Each call
// that tests points holds this lock; each resets recordCache first.
std::mutex &record_cache_mutex() {
  // A leaf: held for one RDKit call that takes no other lock, so a closure
  // running under control_c_mutex (U4) can take it without deadlock. A
  // wrapper that takes another lock while holding it must first call
  // check_not_in_control_c_callback(), as rgroup_mutex() does.
  static std::mutex mutex;
  return mutex;
}

// The voxel grid casts scaled coordinates, radii, and the probe radius to
// int and unsigned and indexes a 16x16x16 grid with them: undefined, or past
// the grid, for values that are not finite or are huge. Coordinates within
// 1e8 and radii within 1e6 keep every cast in range.
void check_finite_within(double value, double bound, const char *what) {
  if (!(std::abs(value) <= bound)) {
    throw ValueErrorException(std::string("DoubleCubicLatticeVolume needs ") +
                              what);
  }
}

// checkExcludedAtoms static_casts any monomer info to AtomPDBResidueInfo,
// which is undefined for a plain AtomMonomerInfo.
void check_monomer_infos(const RWMol &mol) {
  for (const auto *atom : mol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    if (info && !dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info)) {
      throw ValueErrorException(
          "with isProtein, DoubleCubicLatticeVolume needs PDB residue info "
          "for every atom with monomer info; atom " +
          std::to_string(atom->getIdx()) + " has other monomer info");
    }
  }
}

// RDKit lists each atom's neighbors only for atoms with a nonzero radius,
// but finds an atom's list at its atom index, so after an atom with radius
// 0 every later atom reads another's list, and the last read past the end.
// Atoms with radius 0 (or excluded as solvent) must therefore come last.
void check_zero_radii_last(const DoubleCubicLatticeVolume &dclv) {
  bool seen_zero = false;
  for (std::size_t i = 0; i < dclv.radii_.size(); ++i) {
    if (dclv.radii_[i] == 0.0) {
      seen_zero = true;
    } else if (seen_zero) {
      throw ValueErrorException(
          "DoubleCubicLatticeVolume needs atoms with radius 0, and excluded "
          "atoms, after all the others; atom " +
          std::to_string(i) + " follows one");
    }
  }
}

// RDKit 2026.03's includeAsPolar (DCLV.cpp), with includeHs, asks whether
// each of a hydrogen's neighbors is polar by calling itself with includeHs,
// so two bonded hydrogens recurse until the stack overflows (a segfault in
// Python for H2; rdkit/rdkit#9687, fixed after 2026.09.1). A hydrogen
// bonded to a hydrogen is refused for the polar descriptors with hydrogens.
void check_polar_hydrogens(const DoubleCubicLatticeVolume &dclv,
                           bool include_hs) {
  if (!include_hs) {
    return;
  }
  for (const auto *bond : dclv.mol.bonds()) {
    if (bond->getBeginAtom()->getAtomicNum() == 1 &&
        bond->getEndAtom()->getAtomicNum() == 1) {
      throw ValueErrorException(
          "with includeHs, DoubleCubicLatticeVolume's polar descriptors need "
          "no hydrogen bonded to a hydrogen, which RDKit 2026.03 follows "
          "back and forth until the stack overflows; bond " +
          std::to_string(bond->getIdx()) + " joins two");
    }
  }
}

boost::dynamic_bitset<> atom_set(const DoubleCubicLatticeVolume &dclv,
                                 rust::Slice<const std::uint32_t> indices) {
  boost::dynamic_bitset<> atoms(dclv.mol.getNumAtoms());
  for (auto idx : indices) {
    check_index(idx, dclv.mol.getNumAtoms());
    atoms.set(idx);
  }
  return atoms;
}

}  // namespace

std::unique_ptr<DoubleCubicLatticeVolume> dclv_new(
    const RWMol &mol, rust::Slice<const double> radii, bool is_protein,
    bool include_ligand, double probe_radius, std::int32_t conf_id) {
  // The constructor takes the largest radius with *std::max_element, which
  // is undefined for no atoms.
  if (mol.getNumAtoms() == 0) {
    throw ValueErrorException("DoubleCubicLatticeVolume needs atoms");
  }
  // RDKit indexes the radii by atom without checking their number.
  if (!radii.empty() && radii.size() != mol.getNumAtoms()) {
    throw ValueErrorException(
        "DoubleCubicLatticeVolume needs a radius for each atom");
  }
  for (double radius : radii) {
    check_finite_within(radius, 1e6, "finite radii within 1e6");
  }
  check_finite_within(probe_radius, 1e6, "a finite probe radius within 1e6");
  for (const auto &pos : mol.getConformer(conf_id).getPositions()) {
    for (double v : {pos.x, pos.y, pos.z}) {
      check_finite_within(v, 1e8, "finite coordinates within 1e8");
    }
  }
  if (is_protein) {
    check_monomer_infos(mol);
  }
  auto dclv = std::make_unique<DoubleCubicLatticeVolume>(
      mol, std::vector<double>(radii.begin(), radii.end()), is_protein,
      include_ligand, probe_radius, conf_id);
  check_zero_radii_last(*dclv);
  return dclv;
}

#define RDKIT_CXX_DCLV(name, method)                        \
  double name(DoubleCubicLatticeVolume &dclv) {             \
    std::lock_guard<std::mutex> lock(record_cache_mutex()); \
    return dclv.method();                                   \
  }

RDKIT_CXX_DCLV(dclv_get_surface_area, getSurfaceArea)
RDKIT_CXX_DCLV(dclv_get_volume, getVolume)
RDKIT_CXX_DCLV(dclv_get_vdw_volume, getVDWVolume)
RDKIT_CXX_DCLV(dclv_get_compactness, getCompactness)
RDKIT_CXX_DCLV(dclv_get_packing_density, getPackingDensity)

#undef RDKIT_CXX_DCLV

double dclv_get_atom_surface_area(DoubleCubicLatticeVolume &dclv,
                                  std::uint32_t atom_idx) {
  check_index(atom_idx, dclv.mol.getNumAtoms());
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  return dclv.getAtomSurfaceArea(atom_idx);
}

double dclv_get_polar_surface_area(DoubleCubicLatticeVolume &dclv,
                                   bool include_s_and_p, bool include_hs) {
  check_polar_hydrogens(dclv, include_hs);
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  return dclv.getPolarSurfaceArea(include_s_and_p, include_hs);
}

double dclv_get_partial_surface_area(
    DoubleCubicLatticeVolume &dclv,
    rust::Slice<const std::uint32_t> atom_indices) {
  auto atoms = atom_set(dclv, atom_indices);
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  return dclv.getPartialSurfaceArea(atoms);
}

void dclv_get_surface_points(DoubleCubicLatticeVolume &dclv,
                             rust::Vec<std::uint32_t> &atoms,
                             rust::Vec<std::uint32_t> &counts,
                             rust::Vec<Point3D> &points) {
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  for (const auto &[atom, atom_points] : dclv.getSurfacePoints()) {
    atoms.push_back(atom);
    counts.push_back(static_cast<std::uint32_t>(atom_points.size()));
    for (const auto &p : atom_points) {
      points.push_back(Point3D{p.x, p.y, p.z});
    }
  }
}

double dclv_get_atom_volume(DoubleCubicLatticeVolume &dclv,
                            std::uint32_t atom_idx, double solvent_radius) {
  check_index(atom_idx, dclv.mol.getNumAtoms());
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  return dclv.getAtomVolume(atom_idx, solvent_radius);
}

double dclv_get_polar_volume(DoubleCubicLatticeVolume &dclv,
                             bool include_s_and_p, bool include_hs) {
  check_polar_hydrogens(dclv, include_hs);
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  return dclv.getPolarVolume(include_s_and_p, include_hs);
}

double dclv_get_partial_volume(DoubleCubicLatticeVolume &dclv,
                               rust::Slice<const std::uint32_t> atom_indices) {
  auto atoms = atom_set(dclv, atom_indices);
  std::lock_guard<std::mutex> lock(record_cache_mutex());
  return dclv.getPartialVolume(atoms);
}

}  // namespace rdkit_cxx
