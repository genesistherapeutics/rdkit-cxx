/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit::ControlCHandler (RDGeneral/ControlCHandler.h) swaps the process's
// SIGINT handler and keeps the previous one, and its flag, in plain static
// variables, so two calls that construct one at once race. The embedder,
// DetermineBondOrders, and the CIP labeler construct one, and so does the
// fingerprint code that runs the CIP labeler for chirality. Every wrapper
// call that can construct one holds this mutex for the whole RDKit call
// (CONVENTIONS.md U4). A wrapper that also needs random_mutex()
// (global_state.h) takes this one first.

#pragma once

#include <mutex>

namespace rdkit_cxx {

// Defined inline so every wrapper source shares one mutex.
inline std::mutex &control_c_mutex() {
  static std::mutex mutex;
  return mutex;
}

// With the new stereo perception on, the fingerprint code runs the CIP
// labeler for chirality (MorganGenerator.cpp, and getAtomCode in
// FingerprintUtil.cpp, which the atom-pair code uses), so a call that may
// run it holds control_c_mutex then. The legacy perception, RDKit's
// default, never labels there; the lock returned then holds nothing.
std::unique_lock<std::mutex> lock_control_c_for_cip();

}  // namespace rdkit_cxx
