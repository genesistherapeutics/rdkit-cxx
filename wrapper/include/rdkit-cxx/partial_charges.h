/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdPartialCharges: Gasteiger charges.

#pragma once

#include <cstdint>

#include "rdkit-cxx/mol.h"

namespace rdkit_cxx {

void compute_gasteiger_charges(RWMol &mol, std::int32_t n_iter,
                               bool throw_on_param_failure);

}  // namespace rdkit_cxx
