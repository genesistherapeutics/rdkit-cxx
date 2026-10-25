// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_checks.h"

#include <string>
#include <vector>

#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDValue.h>
#include <boost/shared_array.hpp>

namespace rdkit_cxx {

// The matrices are the only properties RDKit stores as shared_arrays, under
// names a caller can prefix, so they are found by type rather than by name.
void clear_matrix_caches(const RDKit::ROMol &mol) {
  std::vector<std::string> keys;
  for (const auto &pair : mol.getDict()) {
    if (RDKit::rdvalue_is<boost::shared_array<double>>(pair.val) ||
        RDKit::rdvalue_is<boost::shared_array<int>>(pair.val)) {
      keys.push_back(pair.key);
    }
  }
  for (const auto &key : keys) {
    mol.clearProp(key);
  }
}

// The names RDKit caches descriptors under: Crippen.cpp, MolSurf.cpp (the
// TPSA names end in its includeSandP flag), ConnectivityDescriptors.cpp,
// PMI.cpp, and the principal axes and moments in MolTransforms.cpp, which
// PMI.cpp and PBF.cpp read.
void clear_descriptor_caches(const RDKit::ROMol &mol) {
  static const char *const names[] = {"_connectivityHKDeltas",
                                      "_connectivityNVals",
                                      "_crippenLogP",
                                      "_crippenMR",
                                      "_crippenLogPContribs",
                                      "_crippenMRContribs",
                                      "_labuteASA",
                                      "_labuteAtomContribs",
                                      "_labuteAtomHContrib",
                                      "_tpsa-0",
                                      "_tpsa-1",
                                      "_tpsaAtomContribs-0",
                                      "_tpsaAtomContribs-1",
                                      "_PMI1",
                                      "_PMI2",
                                      "_PMI3",
                                      "_PMI1_mass",
                                      "_PMI2_mass",
                                      "_PMI3_mass",
                                      "_PMI1_cov",
                                      "_PMI2_cov",
                                      "_PMI3_cov",
                                      "_PMI1_mass_cov",
                                      "_PMI2_mass_cov",
                                      "_PMI3_mass_cov",
                                      "_principalAxes",
                                      "_principalMoments",
                                      "_principalAxes_noH",
                                      "_principalMoments_noH",
                                      "_principalAxes_cov",
                                      "_principalMoments_cov",
                                      "_principalAxes_noH_cov",
                                      "_principalMoments_noH_cov"};
  for (const char *name : names) {
    if (mol.hasProp(name)) {
      mol.clearProp(name);
    }
  }
}

void check_matrix_size(const RDKit::ROMol &mol) {
  if (mol.getNumAtoms() > 46340) {
    throw ValueErrorException(
        "RDKit computes atom matrices in int, which overflows above 46340 "
        "atoms");
  }
}

void check_monomer_types(const RDKit::ROMol &mol) {
  using RDKit::AtomMonomerInfo;
  for (const auto *atom : mol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    if (!info) {
      continue;
    }
    int type = info->getMonomerType();
    if (type < AtomMonomerInfo::UNKNOWN || type > AtomMonomerInfo::OTHER) {
      throw ValueErrorException("atom " + std::to_string(atom->getIdx()) +
                                " has unknown AtomMonomerType " +
                                std::to_string(type));
    }
    if (type == AtomMonomerInfo::PDBRESIDUE &&
        !dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info)) {
      throw ValueErrorException(
          "atom " + std::to_string(atom->getIdx()) +
          " has type PDBRESIDUE, which only AtomPDBResidueInfo can have");
    }
  }
}

}  // namespace rdkit_cxx
