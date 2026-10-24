// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/atom.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>

namespace rdkit_cxx {

std::unique_ptr<Atom> atom_new(rust::Str what) {
  return std::make_unique<Atom>(std::string(what));
}

std::unique_ptr<Atom> atom_from_num(unsigned num) {
  return std::make_unique<Atom>(num);
}

std::unique_ptr<Atom> atom_clone(const Atom &atom) {
  return std::make_unique<Atom>(atom);
}

unsigned atom_get_idx(const Atom &atom) { return atom.getIdx(); }

int atom_get_atomic_num(const Atom &atom) { return atom.getAtomicNum(); }

std::unique_ptr<std::string> atom_get_symbol(const Atom &atom) {
  return std::make_unique<std::string>(atom.getSymbol());
}

int atom_get_formal_charge(const Atom &atom) { return atom.getFormalCharge(); }

bool atom_get_is_aromatic(const Atom &atom) { return atom.getIsAromatic(); }

unsigned atom_get_isotope(const Atom &atom) { return atom.getIsotope(); }

unsigned atom_get_num_explicit_hs(const Atom &atom) {
  return atom.getNumExplicitHs();
}

HybridizationType atom_get_hybridization(const Atom &atom) {
  return atom.getHybridization();
}

ChiralType atom_get_chiral_tag(const Atom &atom) {
  return atom.getChiralTag();
}

void atom_set_atomic_num(Atom &atom, int new_num) {
  atom.setAtomicNum(new_num);
}

void atom_set_formal_charge(Atom &atom, int what) {
  atom.setFormalCharge(what);
}

void atom_set_is_aromatic(Atom &atom, bool what) { atom.setIsAromatic(what); }

void atom_set_isotope(Atom &atom, unsigned what) { atom.setIsotope(what); }

void atom_set_num_explicit_hs(Atom &atom, unsigned what) {
  atom.setNumExplicitHs(what);
}

void atom_set_no_implicit(Atom &atom, bool what) { atom.setNoImplicit(what); }

void atom_set_hybridization(Atom &atom, HybridizationType what) {
  atom.setHybridization(what);
}

void atom_set_chiral_tag(Atom &atom, ChiralType what) {
  atom.setChiralTag(what);
}

unsigned atom_get_degree(const Atom &atom) { return atom.getDegree(); }

unsigned atom_get_total_degree(const Atom &atom) {
  return atom.getTotalDegree();
}

unsigned atom_get_total_num_hs(const Atom &atom, bool include_neighbors) {
  return atom.getTotalNumHs(include_neighbors);
}

unsigned atom_get_num_implicit_hs(const Atom &atom) {
  return atom.getNumImplicitHs();
}

unsigned atom_get_total_valence(const Atom &atom) {
  return atom.getTotalValence();
}

// As AtomIsInRing in Code/GraphMol/nbWrap/Atom.cpp, which finds the SSSR on
// demand.
bool atom_is_in_ring(const Atom &atom) {
  auto &mol = atom.getOwningMol();
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  return mol.getRingInfo()->numAtomRings(atom.getIdx()) != 0;
}

rust::Vec<std::uint32_t> atom_get_neighbors(const Atom &atom) {
  rust::Vec<std::uint32_t> indices;
  for (const auto *neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
    indices.push_back(neighbor->getIdx());
  }
  return indices;
}

rust::Vec<std::uint32_t> atom_get_bonds(const Atom &atom) {
  rust::Vec<std::uint32_t> indices;
  for (const auto *bond : atom.getOwningMol().atomBonds(&atom)) {
    indices.push_back(bond->getIdx());
  }
  return indices;
}

void atom_update_property_cache(Atom &atom, bool strict) {
  atom.updatePropertyCache(strict);
}

}  // namespace rdkit_cxx
