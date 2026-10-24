// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/bond.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

unsigned bond_get_idx(const Bond &bond) { return bond.getIdx(); }

BondType bond_get_bond_type(const Bond &bond) { return bond.getBondType(); }

double bond_get_bond_type_as_double(const Bond &bond) {
  return bond.getBondTypeAsDouble();
}

unsigned bond_get_begin_atom_idx(const Bond &bond) {
  return bond.getBeginAtomIdx();
}

unsigned bond_get_end_atom_idx(const Bond &bond) {
  return bond.getEndAtomIdx();
}

unsigned bond_get_other_atom_idx(const Bond &bond, unsigned this_idx) {
  return bond.getOtherAtomIdx(this_idx);
}

bool bond_get_is_aromatic(const Bond &bond) { return bond.getIsAromatic(); }

bool bond_get_is_conjugated(const Bond &bond) { return bond.getIsConjugated(); }

BondDir bond_get_bond_dir(const Bond &bond) { return bond.getBondDir(); }

BondStereo bond_get_stereo(const Bond &bond) { return bond.getStereo(); }

rust::Vec<std::int32_t> bond_get_stereo_atoms(const Bond &bond) {
  rust::Vec<std::int32_t> atoms;
  for (auto idx : bond.getStereoAtoms()) {
    atoms.push_back(idx);
  }
  return atoms;
}

bool bond_has_owning_mol(const Bond &bond) { return bond.hasOwningMol(); }

bool bond_has_query(const Bond &bond) { return bond.hasQuery(); }

// As BondIsInRing and BondIsInRingSize in Code/GraphMol/nbWrap/Bond.cpp,
// which find the SSSR on demand.
bool bond_is_in_ring(const Bond &bond) {
  auto &mol = bond.getOwningMol();
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  return mol.getRingInfo()->numBondRings(bond.getIdx()) != 0;
}

bool bond_is_in_ring_size(const Bond &bond, int size) {
  auto &mol = bond.getOwningMol();
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  return mol.getRingInfo()->isBondInRingOfSize(bond.getIdx(), size);
}

double bond_get_valence_contrib(const Bond &bond, const RDKit::Atom &at) {
  return bond.getValenceContrib(&at);
}

std::unique_ptr<std::string> bond_get_smarts(const Bond &bond,
                                             bool all_bonds_explicit) {
  if (bond.hasQuery()) {
    return std::make_unique<std::string>(RDKit::SmartsWrite::GetBondSmarts(
        static_cast<const RDKit::QueryBond *>(&bond)));
  }
  return std::make_unique<std::string>(
      RDKit::SmilesWrite::GetBondSmiles(&bond, -1, false, all_bonds_explicit));
}

void bond_set_bond_type(Bond &bond, BondType b_t) { bond.setBondType(b_t); }

void bond_set_is_aromatic(Bond &bond, bool what) { bond.setIsAromatic(what); }

void bond_set_is_conjugated(Bond &bond, bool what) {
  bond.setIsConjugated(what);
}

void bond_set_bond_dir(Bond &bond, BondDir what) { bond.setBondDir(what); }

void bond_set_stereo(Bond &bond, BondStereo what) { bond.setStereo(what); }

void bond_set_stereo_atoms(Bond &bond, unsigned bgn_idx, unsigned end_idx) {
  check_index(bgn_idx, bond.getOwningMol().getNumAtoms());
  check_index(end_idx, bond.getOwningMol().getNumAtoms());
  bond.setStereoAtoms(bgn_idx, end_idx);
}

bool bond_invert_chirality(Bond &bond) { return bond.invertChirality(); }

// copy() keeps the dynamic type, so a query bond's copy keeps its query.
std::unique_ptr<Bond> bond_clone(const Bond &bond) {
  return std::unique_ptr<Bond>(bond.copy());
}

}  // namespace rdkit_cxx
