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
#include <stdexcept>

namespace rdkit_cxx {

// Set while this thread runs a Rust closure that RDKit calls from inside a
// substructure match (FMCS's final match checks run inside vf2), so with a
// ControlCHandler alive: the thread holds control_c_mutex, which it cannot
// release for the closure. A wrapper the closure calls that would take the
// mutex again (undefined for a std::mutex the thread owns), or a lock that
// comes before it in U4's order (which could deadlock with a thread waiting
// for this one), fails with this check instead.
inline thread_local bool in_control_c_callback = false;

inline void check_not_in_control_c_callback() {
  if (in_control_c_callback) {
    throw std::runtime_error(
        "a closure RDKit calls inside a substructure match cannot call "
        "wrappers that match substructures, depict, or take another lock "
        "that comes before the match's (CONVENTIONS.md U4)");
  }
}

// Defined inline so every wrapper source shares one mutex. Every lock site
// gets the mutex here, so each checks check_not_in_control_c_callback.
inline std::mutex &control_c_mutex() {
  check_not_in_control_c_callback();
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
