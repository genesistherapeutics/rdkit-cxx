/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdConformerParser: conformers from Amber trajectories.

#pragma once

#include <cstdint>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

rust::Vec<std::uint32_t> add_conformers_from_amber_trajectory(
    RWMol &mol, rust::Slice<const std::uint8_t> path, std::int32_t num_confs,
    bool clear_confs);

}  // namespace rdkit_cxx
