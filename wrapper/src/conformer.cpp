// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/conformer.h"

#include <iterator>
#include <limits>
#include <string>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

Point3D to_point(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

RDGeom::Point3D from_point(const Point3D &p) {
  return RDGeom::Point3D(p.x, p.y, p.z);
}

constexpr unsigned kMaxId = std::numeric_limits<int>::max();

// getConformer takes an int, so a larger id could not be looked up, and
// would read as a negative id, which means the default conformer.
void check_id(unsigned id) {
  if (id > kMaxId) {
    throw ValueErrorException("conformer ids are at most " +
                              std::to_string(kMaxId));
  }
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

// As SetPos in Code/GraphMol/Wrap/Conformer.cpp.
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

void conformer_set_id(Conformer &conf, unsigned id) {
  check_id(id);
  conf.setId(id);
}

unsigned mol_get_num_conformers(const RWMol &mol) {
  return mol.getNumConformers();
}

const Conformer &mol_get_conformer(const RWMol &mol, int id) {
  return mol.getConformer(id);
}

// The first conformer with exactly `id`: unlike getConformer, no id stands
// for the default conformer.
const Conformer &mol_find_conformer(const RWMol &mol, unsigned id) {
  for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
    if ((*it)->getId() == id) {
      return **it;
    }
  }
  throw ValueErrorException("no conformer has id " + std::to_string(id));
}

Conformer &mol_find_conformer_mut(RWMol &mol, unsigned id) {
  return const_cast<Conformer &>(mol_find_conformer(mol, id));
}

// Lookups by id find the first conformer with it, so an id another
// conformer has would make that one unreachable.
void mol_set_conformer_id(RWMol &mol, unsigned current, unsigned id) {
  check_id(id);
  if (id == current) {
    return;
  }
  for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
    if ((*it)->getId() == id) {
      throw ValueErrorException("another conformer has id " +
                                std::to_string(id));
    }
  }
  mol_find_conformer_mut(mol, current).setId(id);
}

const Conformer &mol_get_conformer_at(const RWMol &mol, unsigned position) {
  check_index(position, mol.getNumConformers());
  auto it = mol.beginConformers();
  std::advance(it, position);
  return **it;
}

// As AddMolConformer in Code/GraphMol/Wrap/Mol.cpp, which copies.
// addConformer assigns the largest id, as an int, plus one, which overflows
// when that is INT_MAX.
unsigned mol_add_conformer(RWMol &mol, const Conformer &conf, bool assign_id) {
  if (assign_id) {
    for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
      if ((*it)->getId() == kMaxId) {
        throw ValueErrorException("no conformer id is left to assign");
      }
    }
  } else {
    check_id(conf.getId());
  }
  return mol.addConformer(new Conformer(conf), assign_id);
}

void mol_remove_conformer(RWMol &mol, unsigned id) { mol.removeConformer(id); }

void mol_remove_all_conformers(RWMol &mol) { mol.clearConformers(); }

}  // namespace rdkit_cxx
