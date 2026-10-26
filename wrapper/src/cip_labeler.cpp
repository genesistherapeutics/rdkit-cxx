// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/cip_labeler.h"

#include <mutex>

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <boost/dynamic_bitset.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

namespace {

// As pyObjToBitset in Code/GraphMol/CIPLabeler/nbWrap/rdCIPLabeler.cpp,
// which sets each index unchecked: dynamic_bitset::set writes past the end
// for an index beyond the atoms or bonds (its bound is only an assert).
boost::dynamic_bitset<> to_bitset(rust::Slice<const std::uint32_t> indices,
                                  unsigned n) {
  boost::dynamic_bitset<> result(n);
  for (auto idx : indices) {
    check_index(idx, n);
    result.set(idx);
  }
  return result;
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
  // assignCIPLabels constructs a ControlCHandler (CIPLabeler.cpp).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  RDKit::CIPLabeler::assignCIPLabels(mol, atom_set, bond_set,
                                     max_recursive_iterations);
}

}  // namespace rdkit_cxx
