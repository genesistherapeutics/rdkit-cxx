/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Threads for the wrappers that parallelize RDKit work themselves.

#pragma once

#include <functional>

namespace rdkit_cxx {

// The number of threads RDKit's getNumThreadsToUse (RDGeneral/RDThreads.h)
// picks: target when positive, else the hardware's concurrency plus target,
// and at least 1. RDKit negates target in int, which is undefined for
// INT_MIN; this computes in 64 bits.
unsigned num_threads_to_use(int target);

// The most writer threads a caller may ask RDKit's MultithreadedMolSupplier
// to start. It starts every one at its first next(), so a count the system
// cannot start fails there, and one this large would only cost memory.
inline constexpr unsigned kMaxSupplierThreads = 1024;

// Calls work(t) for each t in [0, n), each on its own thread (the calling
// thread runs t = 0), and once every call has returned rethrows the first
// exception any of them threw. An exception that leaves a std::thread calls
// std::terminate, which is why RDKit's own threads are not used where its
// work can throw.
void run_on_threads(unsigned n, const std::function<void(unsigned)> &work);

}  // namespace rdkit_cxx
