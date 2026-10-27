// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/cip_labeler.h"

#include <mutex>

#include <GraphMol/Atropisomers.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <boost/dynamic_bitset.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

namespace {

// As pythonObjectToDynBitset in Code/RDBoost/Wrap.cpp, which
// assignCIPLabelsWrapHelper in Code/GraphMol/CIPLabeler/Wrap/rdCIPLabeler.cpp
// uses: each index is checked, since dynamic_bitset::set writes past the end
// for one beyond the atoms or bonds (its bound is only an assert). Python
// treats empty atom and bond lists as None and labels everything; here they
// are empty sets.
boost::dynamic_bitset<> to_bitset(rust::Slice<const std::uint32_t> indices,
                                  unsigned n) {
  boost::dynamic_bitset<> result(n);
  for (auto idx : indices) {
    check_index(idx, n);
    result.set(idx);
  }
  return result;
}

// RDKit 2024.09's AtropisomerBond leaves its carriers empty for a bond with
// an end of no other bonds (getAtropisomerAtomsAndBonds fails), and its label
// reads carriers[0] or carriers[1] when an end has more than one edge past the
// bond (one per other bond at most, rdkit/rdkit#9577). The wrapper leaves
// such a bond out unless neither end has two other bonds, keeping its labels
// as RDKit keeps those of a bond it cannot label. 2024.09's Sp2Bond reads no
// neighbor of its own, so double bonds are left to RDKit.
bool labelable(const Bond &bond) {
  if (bond.getStereo() != Bond::STEREOATROPCW &&
      bond.getStereo() != Bond::STEREOATROPCCW) {
    return true;
  }
  RDKit::Atropisomers::AtropAtomAndBondVec atoms_and_bonds[2];
  return RDKit::Atropisomers::getAtropisomerAtomsAndBonds(
             &bond, atoms_and_bonds, bond.getOwningMol()) ||
         (bond.getBeginAtom()->getDegree() <= 2 &&
          bond.getEndAtom()->getDegree() <= 2);
}

}  // namespace

void assign_cip_labels(RWMol &mol, bool has_atoms,
                       rust::Slice<const std::uint32_t> atoms, bool has_bonds,
                       rust::Slice<const std::uint32_t> bonds,
                       std::uint32_t max_recursive_iterations) {
  auto atom_set = to_bitset(atoms, mol.getNumAtoms());
  auto bond_set = to_bitset(bonds, mol.getNumBonds());
  if (!has_atoms && !has_bonds) {
    atom_set.set();
    bond_set.set();
  }
  for (auto idx = bond_set.find_first(); idx != bond_set.npos;
       idx = bond_set.find_next(idx)) {
    if (!labelable(*mol.getBondWithIdx(idx))) {
      bond_set.reset(idx);
    }
  }
  // RDKit 2026.03's assignCIPLabels resets ControlCHandler; 2024.09's does
  // not, and the lock costs it only concurrency (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  RDKit::CIPLabeler::assignCIPLabels(mol, atom_set, bond_set,
                                     max_recursive_iterations);
}

}  // namespace rdkit_cxx
