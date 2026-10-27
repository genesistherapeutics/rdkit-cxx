// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/cip_labeler.h"

#include <cstddef>
#include <mutex>
#include <string>

#include <GraphMol/Atropisomers.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
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

// RDKit 2025.03's AtropisomerBond constructor returns without carriers for
// a bond getAtropisomerAtomsAndBonds rejects (an end of no other bonds, say),
// and its label then reads both carriers from the empty vector. Such bonds
// are left unlabeled, as RDKit 2026.09 does (rdkit/rdkit#9577), clearing a
// stale label as RDKit clears each bond's before labeling it. (Its Sp2Bond
// takes findStereoAtoms's atoms as carriers, which it checks.)
bool labelable(const Bond &bond) {
  switch (bond.getStereo()) {
    case Bond::STEREOATROPCW:
    case Bond::STEREOATROPCCW: {
      RDKit::Atropisomers::AtropAtomAndBondVec atoms_and_bonds[2];
      return RDKit::Atropisomers::getAtropisomerAtomsAndBonds(
          &bond, atoms_and_bonds, bond.getOwningMol());
    }
    default:
      return true;
  }
}

// As findConfigs in CIPLabeler.cpp: the atoms and bonds RDKit labels.
bool is_config(const Atom &atom) {
  return atom.getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
         atom.getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW;
}

bool is_config(const Bond &bond) {
  switch (bond.getStereo()) {
    case Bond::STEREOE:
    case Bond::STEREOZ:
    case Bond::STEREOCIS:
    case Bond::STEREOTRANS:
    case Bond::STEREOATROPCW:
    case Bond::STEREOATROPCCW:
      return true;
    default:
      return false;
  }
}

// The most stereo atoms and bonds labeled in one call (assign_cip_labels).
constexpr std::size_t kMaxCipConfigs = 64;

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
    auto *bond = mol.getBondWithIdx(idx);
    if (!labelable(*bond)) {
      bond_set.reset(idx);
      bond->clearProp(RDKit::common_properties::_CIPCode);
    }
  }
  // RDKit 2025.03's PairList::addAndPair (rules/Pairlist.h), which Rules 4b
  // and 5 fill with the descriptors of a branch's stereo atoms and bonds,
  // sets bit 63 - size of a 64-bit pairing for each, which is a negative
  // shift once a list holds more than 64 (2026.03 removed it). A list holds
  // its reference descriptor and one for each node of a branch of the
  // digraph that has one: in a tree or a single ring each other stereo atom
  // or bond is at most one node of a branch, so at most kMaxCipConfigs of
  // them keep every list within 64. (In fused ring systems one atom can be
  // several nodes of a branch; a sound bound there would refuse common
  // polycycles, so it is not checked: README, Known limitations.)
  std::size_t configs = 0;
  for (auto idx = atom_set.find_first(); idx != atom_set.npos;
       idx = atom_set.find_next(idx)) {
    configs += is_config(*mol.getAtomWithIdx(idx));
  }
  for (auto idx = bond_set.find_first(); idx != bond_set.npos;
       idx = bond_set.find_next(idx)) {
    configs += is_config(*mol.getBondWithIdx(idx));
  }
  if (configs > kMaxCipConfigs) {
    throw ValueErrorException(
        "RDKit 2025.03's CIP labeler shifts by the number of stereo "
        "descriptors it pairs, which is undefined past 64, so it labels at "
        "most " +
        std::to_string(kMaxCipConfigs) + " stereo atoms and bonds at a time; " +
        std::to_string(configs) + " were asked for");
  }
  // assignCIPLabels calls ControlCHandler::reset, which installs its SIGINT
  // handler (CIPLabeler.cpp).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  RDKit::CIPLabeler::assignCIPLabels(mol, atom_set, bond_set,
                                     max_recursive_iterations);
}

}  // namespace rdkit_cxx
