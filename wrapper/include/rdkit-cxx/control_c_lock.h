/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit::ControlCHandler (RDGeneral/ControlCHandler.h) swaps the process's
// SIGINT handler and keeps the previous one, and its flag, in plain static
// variables, so two threads that use it at once race. In RDKit 2025.09 the
// embedder, DetermineBondOrders, the CIP labeler (which the fingerprint code
// runs for chirality), and the synthon search call the static
// ControlCHandler::reset(), which installs RDKit's handler and keeps nothing
// to restore; RDKit's handler swallows the next SIGINT and installs whatever
// handler a ControlCHandler last saved. Substructure matches do not use it in
// 2025.09 (RDKit 2026.03's boost::vf2 and vf2_all in Substruct/vf2.hpp call
// reset() too), but the wrappers hold this mutex around every call that
// matches substructures as well, as the rdkit-2026.03 branch does: the
// match methods and the transforms that match, descriptors and fingerprints
// that type atoms by patterns (Crippen, the SMARTS counts, MACCS, pattern
// fingerprints, feature invariants), reactions, alignment, force-field
// builders, depiction with templates, standardization, the catalogs, and
// MCS. Every wrapper call that can reach a handler holds this mutex for the
// whole RDKit call, and runs RDKit's own threads, which the mutex cannot
// order, on the calling thread instead (CONVENTIONS.md U4). A wrapper that
// also needs random_mutex() (global_state.h) takes this one first.
//
// The mutex saves SIGINT's disposition when it is taken and restores it when
// it is released, so the process's own handler is back after every call
// that holds it, as a scoped ControlCHandler (RDKit 2026.09,
// rdkit/rdkit#9393) leaves it. A handler another thread installs while such
// a call runs is replaced by the one before.

#pragma once

#include <mutex>
#include <stdexcept>

#include <signal.h>

namespace rdkit_cxx {

// Set while this thread runs a Rust closure that RDKit calls from inside a
// substructure match (FMCS's final match checks run inside vf2): the thread
// holds control_c_mutex, which it cannot release for the closure. A wrapper the
// closure calls that would take the mutex again (undefined for a std::mutex the
// thread owns), or a lock that comes before it in U4's order (which could
// deadlock with a thread waiting for this one), fails with this check instead.
inline thread_local bool in_control_c_callback = false;

inline void check_not_in_control_c_callback() {
  if (in_control_c_callback) {
    throw std::runtime_error(
        "a closure RDKit calls inside a substructure match cannot call "
        "wrappers that match substructures, depict, or take another lock "
        "that comes before the match's (CONVENTIONS.md U4)");
  }
}

// Whether this thread holds control_c_mutex().
inline thread_local bool control_c_held = false;

// control_c_mutex()'s type: a std::mutex that records whether this thread
// holds it, so ReleaseControlC can tell when the thread holds it through a
// plain lock, which only that lock can release.
class ControlCMutex {
 public:
  void lock() {
    d_mutex.lock();
    taken();
  }
  bool try_lock() {
    if (!d_mutex.try_lock()) {
      return false;
    }
    taken();
    return true;
  }
  void unlock() {
    sigaction(SIGINT, &d_sigint, nullptr);
    control_c_held = false;
    d_mutex.unlock();
  }

 private:
  void taken() {
    control_c_held = true;
    sigaction(SIGINT, nullptr, &d_sigint);
  }

  std::mutex d_mutex;
  // SIGINT's disposition when the mutex was taken.
  struct sigaction d_sigint{};
};

// Defined inline so every wrapper source shares one mutex. Every lock site
// gets the mutex here, so each checks check_not_in_control_c_callback.
inline ControlCMutex &control_c_mutex() {
  check_not_in_control_c_callback();
  static ControlCMutex mutex;
  return mutex;
}

// Holds control_c_mutex() for an RDKit call that runs code that must not
// hold it: Rust closures, which may call wrappers that take it, and steps
// that take a lock ordered before it. A ReleaseControlC around such code
// releases the lock the thread's innermost HoldControlC holds, if it holds
// it, and takes it back after. RDKit runs such code between matches, when
// no ControlCHandler is alive. When the thread holds the mutex through
// another lock, which a ReleaseControlC cannot release, the code inside runs
// as in_control_c_callback, so a wrapper it calls that would take the mutex
// again fails instead.
class HoldControlC {
 public:
  HoldControlC();
  ~HoldControlC();
  HoldControlC(const HoldControlC &) = delete;
  HoldControlC &operator=(const HoldControlC &) = delete;

 private:
  std::unique_lock<ControlCMutex> d_lock;
  std::unique_lock<ControlCMutex> *d_outer;
};

class ReleaseControlC {
 public:
  ReleaseControlC();
  ~ReleaseControlC();
  ReleaseControlC(const ReleaseControlC &) = delete;
  ReleaseControlC &operator=(const ReleaseControlC &) = delete;

 private:
  std::unique_lock<ControlCMutex> *d_held;
  // Whether this set in_control_c_callback.
  bool d_marked = false;
};

// With the new stereo perception on, the fingerprint code runs the CIP
// labeler for chirality (MorganGenerator.cpp, and getAtomCode in
// FingerprintUtil.cpp, which the atom-pair code uses), so a call that may
// run it holds control_c_mutex then. The legacy perception, RDKit's
// default, never labels there; the lock returned then holds nothing.
std::unique_lock<ControlCMutex> lock_control_c_for_cip();

}  // namespace rdkit_cxx
