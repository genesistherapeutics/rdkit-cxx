// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/conformer.h"

#include <iterator>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

Point3D to_point(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

RDGeom::Point3D from_point(const Point3D &p) {
  return RDGeom::Point3D(p.x, p.y, p.z);
}

}  // namespace

std::unique_ptr<Conformer> conformer_new() {
  return std::make_unique<Conformer>();
}

std::unique_ptr<Conformer> conformer_with_num_atoms(unsigned num_atoms) {
  return std::make_unique<Conformer>(num_atoms);
}

// Conformer's copy constructor keeps the source's owning molecule, so a
// copy taken from a molecule would point at it after it is gone, and every
// RDKit function that calls getOwningMol on the copy would read freed memory
// (Python's Conformer(other) does). The copy here has no owner: its
// positions, id, dimensionality, and properties.
std::unique_ptr<Conformer> conformer_clone(const Conformer &conf) {
  auto copy = std::make_unique<Conformer>(conf.getNumAtoms());
  static_cast<RDKit::RDProps &>(*copy) =
      static_cast<const RDKit::RDProps &>(conf);
  copy->getPositions() = conf.getPositions();
  copy->setId(conf.getId());
  copy->set3D(conf.is3D());
  return copy;
}

unsigned conformer_get_id(const Conformer &conf) { return conf.getId(); }

unsigned conformer_get_num_atoms(const Conformer &conf) {
  return conf.getNumAtoms();
}

bool conformer_is_3d(const Conformer &conf) { return conf.is3D(); }

bool conformer_has_owning_mol(const Conformer &conf) {
  return conf.hasOwningMol();
}

Point3D conformer_get_atom_position(const Conformer &conf, unsigned aid) {
  check_index(aid, conf.getNumAtoms());
  return to_point(conf.getAtomPos(aid));
}

rust::Vec<Point3D> conformer_get_positions(const Conformer &conf) {
  rust::Vec<Point3D> positions;
  for (const auto &p : conf.getPositions()) {
    positions.push_back(to_point(p));
  }
  return positions;
}

void conformer_set_atom_position(Conformer &conf, unsigned aid,
                                 const Point3D &loc) {
  check_index(aid, conf.getNumAtoms());
  conf.setAtomPos(aid, from_point(loc));
}

// As SetPos in Code/GraphMol/nbWrap/Conformer.cpp.
void conformer_set_positions(Conformer &conf,
                             rust::Slice<const Point3D> positions) {
  if (positions.size() != conf.getNumAtoms()) {
    throw ValueErrorException(
        "Position array shape doesn't equal the number of atoms in the "
        "conformer");
  }
  auto &pos = conf.getPositions();
  for (std::size_t i = 0; i < positions.size(); ++i) {
    pos[i] = from_point(positions[i]);
  }
}

void conformer_set_3d(Conformer &conf, bool v) { conf.set3D(v); }

void conformer_set_id(Conformer &conf, unsigned id) { conf.setId(id); }

unsigned mol_get_num_conformers(const RWMol &mol) {
  return mol.getNumConformers();
}

const Conformer &mol_get_conformer(const RWMol &mol, int id) {
  return mol.getConformer(id);
}

Conformer &mol_get_conformer_mut(RWMol &mol, int id) {
  return mol.getConformer(id);
}

const Conformer &mol_get_conformer_at(const RWMol &mol, unsigned position) {
  check_index(position, mol.getNumConformers());
  auto it = mol.beginConformers();
  std::advance(it, position);
  return **it;
}

// As Mol.AddConformer in Code/GraphMol/nbWrap/Mol.cpp, which copies.
unsigned mol_add_conformer(RWMol &mol, const Conformer &conf, bool assign_id) {
  return mol.addConformer(new Conformer(conf), assign_id);
}

void mol_remove_conformer(RWMol &mol, unsigned id) { mol.removeConformer(id); }

void mol_remove_all_conformers(RWMol &mol) { mol.clearConformers(); }

}  // namespace rdkit_cxx
