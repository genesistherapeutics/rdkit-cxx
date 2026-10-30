/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdCIPLabeler: CIP labels from the new labeler.

#pragma once

#include <cstdint>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Labels the given atoms and bonds; with neither list, the whole molecule.
void assign_cip_labels(RWMol &mol, bool has_atoms,
                       rust::Slice<const std::uint32_t> atoms, bool has_bonds,
                       rust::Slice<const std::uint32_t> bonds,
                       std::uint32_t max_recursive_iterations);

}  // namespace rdkit_cxx
