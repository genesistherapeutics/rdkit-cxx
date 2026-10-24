// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_checks.h"

#include <string>
#include <vector>

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

void check_matrix_size(const RDKit::ROMol &mol) {
  if (mol.getNumAtoms() > 46340) {
    throw ValueErrorException(
        "RDKit computes atom matrices in int, which overflows above 46340 "
        "atoms");
  }
}

}  // namespace rdkit_cxx
