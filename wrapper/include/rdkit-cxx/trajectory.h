/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdtrajectory: RDKit's AMBER and GROMOS trajectory readers.

#pragma once

#include <cstdint>

#include "rust/cxx.h"

namespace rdkit_cxx {

std::uint32_t read_amber_trajectory(rust::Slice<const std::uint8_t> path,
                                    std::uint32_t num_points,
                                    rust::Vec<double> &coords);

std::uint32_t read_gromos_trajectory(rust::Slice<const std::uint8_t> path,
                                     std::uint32_t num_points,
                                     rust::Vec<double> &coords);

}  // namespace rdkit_cxx
