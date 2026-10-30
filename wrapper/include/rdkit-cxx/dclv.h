/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDescriptors's DoubleCubicLatticeVolume. RDKit 2024.09's computes
// everything in its constructor, so the wrapper constructs one and returns
// its results.

#pragma once

#include <cstdint>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct DclvResults;

DclvResults dclv_compute(const RWMol &mol, bool is_protein, bool include_ligand,
                         double probe_radius, std::int32_t depth,
                         std::int32_t dot_density);

}  // namespace rdkit_cxx
