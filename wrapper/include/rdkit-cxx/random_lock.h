/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit draws random numbers from one process-wide generator,
// RDKit::getRandomGenerator() in RDGeneral/utils.cpp, which it seeds and
// advances without a lock: random SMILES (doRandom) and the depictor's random
// bond flips both use it. Every wrapper whose RDKit call may draw from it
// holds this lock (CONVENTIONS.md U4). A wrapper that also holds
// lock_depiction() takes that first.

#pragma once

#include <mutex>

namespace rdkit_cxx {

std::unique_lock<std::mutex> lock_random_generator();

}  // namespace rdkit_cxx
