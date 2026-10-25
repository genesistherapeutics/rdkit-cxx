// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/partial_charges.h"

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

void compute_gasteiger_charges(RWMol &mol, std::int32_t n_iter,
                               bool throw_on_param_failure) {
  // RDKit counts the iterations as unsigned, so a negative count would run
  // for about four billion.
  if (n_iter < 0) {
    throw ValueErrorException("nIter must not be negative");
  }
  RDKit::computeGasteigerCharges(mol, n_iter, throw_on_param_failure);
}

}  // namespace rdkit_cxx
