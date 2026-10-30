// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/atom.h"

#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace rdkit_cxx {

std::unique_ptr<Atom> atom_new(rust::Str what) {
  return std::make_unique<Atom>(std::string(what));
}

std::unique_ptr<Atom> atom_from_num(unsigned num) {
  return std::make_unique<Atom>(num);
}

// copy() keeps the dynamic type, so a query atom's copy keeps its query.
std::unique_ptr<Atom> atom_clone(const Atom &atom) {
  return std::unique_ptr<Atom>(atom.copy());
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

ChiralType atom_get_chiral_tag(const Atom &atom) { return atom.getChiralTag(); }

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

unsigned atom_noatom() { return RDKit::Chirality::StereoInfo::NOATOM; }

int atom_get_atom_map_num(const Atom &atom) { return atom.getAtomMapNum(); }

double atom_get_mass(const Atom &atom) { return atom.getMass(); }

bool atom_get_no_implicit(const Atom &atom) { return atom.getNoImplicit(); }

unsigned atom_get_num_radical_electrons(const Atom &atom) {
  return atom.getNumRadicalElectrons();
}

bool atom_has_owning_mol(const Atom &atom) { return atom.hasOwningMol(); }

bool atom_has_query(const Atom &atom) { return atom.hasQuery(); }

bool atom_needs_update_property_cache(const Atom &atom) {
  return atom.needsUpdatePropertyCache();
}

unsigned atom_get_valence(const Atom &atom, ValenceType which) {
  return atom.getValence(which);
}

std::unique_ptr<std::string> atom_get_smarts(const Atom &atom, bool do_kekule,
                                             bool all_hs_explicit,
                                             bool isomeric_smiles) {
  if (atom.hasQuery()) {
    return std::make_unique<std::string>(RDKit::SmartsWrite::GetAtomSmarts(
        static_cast<const RDKit::QueryAtom *>(&atom)));
  }
  return std::make_unique<std::string>(RDKit::SmilesWrite::GetAtomSmiles(
      &atom, do_kekule, nullptr, all_hs_explicit, isomeric_smiles));
}

void atom_set_atom_map_num(Atom &atom, int mapno, bool strict) {
  atom.setAtomMapNum(mapno, strict);
}

void atom_set_num_radical_electrons(Atom &atom, unsigned num) {
  atom.setNumRadicalElectrons(num);
}

void atom_clear_property_cache(Atom &atom) { atom.clearPropertyCache(); }

bool atom_invert_chirality(Atom &atom) { return atom.invertChirality(); }

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

// As AtomIsInRing in Code/GraphMol/Wrap/Atom.cpp, which finds the SSSR on
// demand.
bool atom_is_in_ring(const Atom &atom) {
  auto &mol = atom.getOwningMol();
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  return mol.getRingInfo()->numAtomRings(atom.getIdx()) != 0;
}

// As AtomIsInRingSize in Code/GraphMol/Wrap/Atom.cpp.
bool atom_is_in_ring_size(const Atom &atom, int size) {
  auto &mol = atom.getOwningMol();
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  return mol.getRingInfo()->isAtomInRingOfSize(atom.getIdx(), size);
}

bool atom_has_valence_violation(const Atom &atom) {
  return atom.hasValenceViolation();
}

std::unique_ptr<std::string> atom_get_alias(const Atom &atom) {
  return std::make_unique<std::string>(RDKit::getAtomAlias(&atom));
}

void atom_set_alias(Atom &atom, rust::Str rlabel) {
  RDKit::setAtomAlias(&atom, std::string(rlabel));
}

int atom_get_r_label(const Atom &atom) { return RDKit::getAtomRLabel(&atom); }

void atom_set_r_label(Atom &atom, int rlabel) {
  RDKit::setAtomRLabel(&atom, rlabel);
}

std::unique_ptr<std::string> atom_get_value(const Atom &atom) {
  return std::make_unique<std::string>(RDKit::getAtomValue(&atom));
}

void atom_set_value(Atom &atom, rust::Str rlabel) {
  RDKit::setAtomValue(&atom, std::string(rlabel));
}

std::unique_ptr<std::string> atom_get_supplemental_smiles_label(
    const Atom &atom) {
  return std::make_unique<std::string>(
      RDKit::getSupplementalSmilesLabel(&atom));
}

void atom_set_supplemental_smiles_label(Atom &atom, rust::Str label) {
  RDKit::setSupplementalSmilesLabel(&atom, std::string(label));
}

unsigned atom_get_num_pi_electrons(const Atom &atom) {
  return RDKit::numPiElectrons(atom);
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
