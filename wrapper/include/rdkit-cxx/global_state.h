/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Locks for process-wide state RDKit writes without one (CONVENTIONS.md U4).
// A wrapper that needs both takes control_c_mutex first.

#pragma once

#include <mutex>

namespace rdkit_cxx {

// RDKit::ControlCHandler (RDGeneral/ControlCHandler.h) swaps the SIGINT
// handler and keeps the previous one, and its flag, in plain static
// variables, so two calls that construct one at once race. The embedder,
// DetermineBondOrders, and the CIP labeler construct one, and so do the
// fingerprint generators that run the CIP labeler for chirality.
std::mutex &control_c_mutex();

// With the new stereo perception on, the fingerprint code runs the CIP
// labeler for chirality (MorganGenerator.cpp, and getAtomCode in
// FingerprintUtil.cpp, which the atom-pair code uses), so a call that may
// run it holds control_c_mutex then. The legacy perception, RDKit's
// default, never labels there; the lock returned then holds nothing.
std::unique_lock<std::mutex> lock_control_c_for_cip();

// RDKit's random generator (getRandomGenerator and getDoubleRandomSource in
// RDGeneral/utils.cpp) is a global that callers seed and draw from with no
// lock: the embedder without a seed, DistGeom's EmbedBoundsMatrix, and
// MolAlign's randomTransform use it.
std::mutex &random_mutex();

}  // namespace rdkit_cxx
