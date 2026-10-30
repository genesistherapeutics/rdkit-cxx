/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Locks for process-wide state RDKit shares between modules and touches with
// no lock of its own (CONVENTIONS.md U4).

#pragma once

#include <mutex>

namespace rdkit_cxx {

// RDKit's random number generator (getRandomGenerator, getRandomVal, and
// getDoubleRandomSource in RDGeneral/utils.cpp) is a plain global that RDKit
// seeds and draws from with no lock, so two threads that reach it at once
// race. Every wrapper call that can reach it holds this one mutex for the
// whole RDKit call. Defined inline so every wrapper source shares one mutex.
inline std::mutex &random_mutex() {
  static std::mutex mutex;
  return mutex;
}

}  // namespace rdkit_cxx
