/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Murtagh's hierarchical clustering (hc_ and distdriver_, f2c translations
// in Code/ML/Cluster/Murtagh/hc.c and hcdriver.c), which
// rdkit.ML.Cluster.Clustering and HierarchicalClusterPicker call. What the
// C trusts is checked here, for both (CONVENTIONS.md U1).

#pragma once

#include <cstdint>
#include <mutex>
#include <vector>

namespace rdkit_cxx {

// hc_ keeps its loop indices and merge candidates in static variables, and
// distdriver_ its work arrays, so two calls at once race. Every call holds
// this mutex for the whole RDKit call (U4). Defined inline so every wrapper
// source shares one mutex.
inline std::mutex &murtagh_mutex() {
  // A leaf: held for one RDKit call that takes no other lock, so a closure
  // running under control_c_mutex (U4) can take it without deadlock. A
  // wrapper that takes another lock while holding it must first call
  // check_not_in_control_c_callback(), as rgroup_mutex() does.
  static std::mutex mutex;
  return mutex;
}

// The clustering criteria hc_ implements (iopt), 1 (Ward's) through 7
// (centroid). hc_ applies no update for another value, so the merges mean
// nothing; ValueErrorException.
void check_murtagh_method(std::int64_t method);

// Readies a condensed distance matrix of n items (n * (n - 1) / 2 values,
// the lower triangle by rows) for hc_, and returns the power of two by
// which the merge distances it then returns must be scaled (ldexp) to be
// those of the given matrix. Throws ValueErrorException for fewer than 2
// items, a matrix of another length, or a value that is not finite.
//
// hc_ reads and writes out of bounds for fewer than 2 items (heap
// corruption in Python). It chooses each merge among the clusters whose
// distance is below its sentinel 1e20, and when none is it merges the
// clusters its static variables kept from an earlier merge or call (0 at
// first), out of bounds (a segfault or heap corruption in Python for
// distances of 1e21, or infinity). Every criterion keeps the distances
// within 8 n^2 times the largest magnitude given: the single, complete,
// average, and McQuitty updates pick or average distances, and Ward's, the
// centroid, and the median updates give distances between points in the
// convex hull of a pseudo-Euclidean embedding of the matrix (each point's
// squared norm is at most 2 n times the largest magnitude, so each distance
// at most 8 n times it), times at most n / 2 for Ward's. Infinite values
// are refused, and so is NaN, after which HierarchicalClusterPicker::pick
// reads a representative it never set. A matrix whose bound reaches 2^56 is
// scaled down by a power of two first. The scaling is exact (save for
// values so much smaller than the largest that they become subnormal) and
// every update is linear, so the merges are the same.
int prepare_murtagh_distances(std::vector<double> &dists, std::uint64_t n);

}  // namespace rdkit_cxx
