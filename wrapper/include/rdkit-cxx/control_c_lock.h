/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit::ControlCHandler (RDGeneral/ControlCHandler.h) swaps the process's
// SIGINT handler and keeps the previous one, and its flag, in plain static
// variables, so two calls that construct one at once race, and can leave
// RDKit's handler installed after both. Every substructure match constructs one
// (boost::vf2 and vf2_all in Substruct/vf2.hpp), so it is reached from the
// match methods and the transforms that match, from descriptors and
// fingerprints that type atoms by patterns (Crippen, the SMARTS counts, MACCS,
// pattern fingerprints, feature invariants), from reactions, alignment,
// force-field builders, depiction with templates, and standardization; the
// embedder, DetermineBondOrders, and the CIP labeler construct one directly,
// and so does the fingerprint code that runs the CIP labeler for chirality.
// Every wrapper call that can construct one holds this mutex for the whole
// RDKit call, and runs RDKit's own threads, which the mutex cannot order, on
// the calling thread instead (CONVENTIONS.md U4). A wrapper that also needs
// random_mutex() (global_state.h) takes this one first.

#pragma once

#include <mutex>

namespace rdkit_cxx {

// Defined inline so every wrapper source shares one mutex.
inline std::mutex &control_c_mutex() {
  static std::mutex mutex;
  return mutex;
}

// Holds control_c_mutex() for an RDKit call that runs code that must not
// hold it: Rust closures, which may call wrappers that take it, and steps
// that take a lock ordered before it. A ReleaseControlC around such code
// releases the lock the thread's innermost HoldControlC holds, if it holds
// it, and takes it back after. RDKit runs such code between matches, when
// no ControlCHandler is alive.
class HoldControlC {
 public:
  HoldControlC();
  ~HoldControlC();
  HoldControlC(const HoldControlC &) = delete;
  HoldControlC &operator=(const HoldControlC &) = delete;

 private:
  std::unique_lock<std::mutex> d_lock;
  std::unique_lock<std::mutex> *d_outer;
};

class ReleaseControlC {
 public:
  ReleaseControlC();
  ~ReleaseControlC();
  ReleaseControlC(const ReleaseControlC &) = delete;
  ReleaseControlC &operator=(const ReleaseControlC &) = delete;

 private:
  std::unique_lock<std::mutex> *d_held;
};

// With the new stereo perception on, the fingerprint code runs the CIP
// labeler for chirality (MorganGenerator.cpp, and getAtomCode in
// FingerprintUtil.cpp, which the atom-pair code uses), so a call that may
// run it holds control_c_mutex then. The legacy perception, RDKit's
// default, never labels there; the lock returned then holds nothing.
std::unique_lock<std::mutex> lock_control_c_for_cip();

}  // namespace rdkit_cxx
