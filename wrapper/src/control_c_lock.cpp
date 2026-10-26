// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/control_c_lock.h"

#include <GraphMol/Chirality.h>

namespace rdkit_cxx {

namespace {

// The lock the thread's innermost HoldControlC holds, if any.
thread_local std::unique_lock<std::mutex> *held_control_c = nullptr;

}  // namespace

HoldControlC::HoldControlC()
    : d_lock(control_c_mutex()), d_outer(held_control_c) {
  held_control_c = &d_lock;
}

HoldControlC::~HoldControlC() { held_control_c = d_outer; }

ReleaseControlC::ReleaseControlC() : d_held(held_control_c) {
  if (d_held && d_held->owns_lock()) {
    d_held->unlock();
  } else {
    d_held = nullptr;
  }
}

ReleaseControlC::~ReleaseControlC() {
  if (d_held) {
    d_held->lock();
  }
}

std::unique_lock<std::mutex> lock_control_c_for_cip() {
  if (RDKit::Chirality::getUseLegacyStereoPerception()) {
    return {};
  }
  return std::unique_lock<std::mutex>(control_c_mutex());
}

}  // namespace rdkit_cxx
