/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDescriptors's DoubleCubicLatticeVolume, which RDKit 2024.03 computes
// in full when it is constructed: the wrapper returns its results.

#pragma once

#include <cstdint>

#include <GraphMol/Descriptors/DCLV.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct DoubleCubicLatticeVolumeData;

DoubleCubicLatticeVolumeData dclv_compute(const RWMol &mol, bool is_protein,
                                          bool include_ligand,
                                          double probe_radius,
                                          std::int32_t depth,
                                          std::int32_t dot_density);

}  // namespace rdkit_cxx
