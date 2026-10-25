/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit's 2D depiction keeps its settings in globals without a lock:
// RDDepict's bond length (BOND_LEN, which Python's Compute2DCoords changes
// for one call), preferCoordGen, and the ring-system templates, all read by
// every compute2DCoords; and CoordGen's template list
// (sketcherMinimizer::m_templates), filled on first use, and its template
// directory, which a CoordGen call can write. Every wrapper whose RDKit call
// may compute 2D coordinates or change those settings holds this lock
// (CONVENTIONS.md U4), which serializes 2D depiction across threads.

#pragma once

#include <mutex>

namespace rdkit_cxx {

std::unique_lock<std::mutex> lock_depiction();

}  // namespace rdkit_cxx
