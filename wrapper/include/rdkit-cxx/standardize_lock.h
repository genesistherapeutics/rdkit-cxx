/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// How the rdMolStandardize wrappers hold control_c_mutex.
//
// Every substructure match constructs a ControlCHandler (control_c_lock.h),
// and most rdMolStandardize calls match: the normalizer's transformations
// (runReactant), the reionizer's and uncharger's patterns, the fragment
// remover's fragments, the metal disconnector's queries, FragmentValidation,
// tautomer enumeration and the substructure scores, and so cleanup, the
// parents, the pipeline, and standardizeSmiles. Their wrappers hold the
// mutex for the whole RDKit call through a StandardizeMatching
// (CONVENTIONS.md U4).
//
// RDKit calls the Rust callbacks (validation methods, the tautomer
// enumerator's callback, and score functions) only between its matches: a
// ControlCHandler lives only inside a SubstructMatch (vf2 and vf2_all in
// Substruct/vf2.hpp), which calls none of them. So a StandardizeUnlocked
// releases the thread's lock while a closure runs and takes it back after,
// even when the closure throws, and the closure may call wrappers that take
// the lock themselves.

#pragma once

#include <mutex>
#include <stdexcept>

#include "rdkit-cxx/control_c_lock.h"

namespace rdkit_cxx {

// The lock the thread's innermost StandardizeMatching took, if any. It is
// not held while a StandardizeUnlocked has released it.
inline std::unique_lock<std::mutex> *&standardize_control_c() {
  thread_local std::unique_lock<std::mutex> *held = nullptr;
  return held;
}

// Holds control_c_mutex until destroyed when the call matches. A thread
// already holding it through another StandardizeMatching, outside a
// callback that released it, would deadlock on the mutex, which is
// undefined for std::mutex, so that is refused. Callbacks always release it,
// so no call path does that.
class StandardizeMatching {
 public:
  explicit StandardizeMatching(bool matches = true)
      : d_outer(standardize_control_c()) {
    if (!matches) {
      return;
    }
    if (d_outer && d_outer->owns_lock()) {
      throw std::logic_error(
          "an rdMolStandardize call that matches substructures was made "
          "while this thread holds control_c_mutex");
    }
    d_lock = std::unique_lock<std::mutex>(control_c_mutex());
    standardize_control_c() = &d_lock;
    d_pushed = true;
  }
  ~StandardizeMatching() {
    if (d_pushed) {
      standardize_control_c() = d_outer;
    }
  }
  StandardizeMatching(const StandardizeMatching &) = delete;
  StandardizeMatching &operator=(const StandardizeMatching &) = delete;

 private:
  std::unique_lock<std::mutex> *d_outer;
  std::unique_lock<std::mutex> d_lock;
  bool d_pushed = false;
};

// Releases the thread's StandardizeMatching lock, if it holds it, until
// destroyed. A callback can run without it: one called by a wrapper that
// does not match (a score function given to pickCanonical, or a callback
// called directly from Rust).
class StandardizeUnlocked {
 public:
  StandardizeUnlocked() {
    auto *held = standardize_control_c();
    if (held && held->owns_lock()) {
      held->unlock();
      d_held = held;
    }
  }
  ~StandardizeUnlocked() {
    if (d_held) {
      d_held->lock();
    }
  }
  StandardizeUnlocked(const StandardizeUnlocked &) = delete;
  StandardizeUnlocked &operator=(const StandardizeUnlocked &) = delete;

 private:
  std::unique_lock<std::mutex> *d_held = nullptr;
};

}  // namespace rdkit_cxx
