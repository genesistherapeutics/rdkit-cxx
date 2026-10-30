// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/murtagh_checks.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

// f2c's integer, which hc_ and distdriver_ index with, is long: n * (n - 1)
// / 2 in 32 bits would overflow for 65537 items.
static_assert(sizeof(long) >= 8, "hc_ indexes the distances in long");

void check_murtagh_method(std::int64_t method) {
  if (method < 1 || method > 7) {
    throw ValueErrorException(
        "the clustering method must be 1 (Ward's) through 7 (centroid)");
  }
}

int prepare_murtagh_distances(std::vector<double> &dists, std::uint64_t n) {
  if (n < 2) {
    throw ValueErrorException("Murtagh's clustering needs at least 2 items");
  }
  if (dists.size() != n * (n - 1) / 2) {
    throw ValueErrorException(
        "the distance matrix must have n * (n - 1) / 2 values");
  }
  double largest = 0.0;
  for (auto d : dists) {
    if (!std::isfinite(d)) {
      throw ValueErrorException("Murtagh's clustering needs finite distances");
    }
    largest = std::max(largest, std::abs(d));
  }
  // The bound, 8 n^2 times the largest magnitude, is below
  // 2^(largest_exp + count_exp), which is computed without overflowing.
  constexpr int kLimitExponent = 56;
  int largest_exp = 0;
  int count_exp = 0;
  std::frexp(largest, &largest_exp);
  std::frexp(8.0 * static_cast<double>(n) * static_cast<double>(n), &count_exp);
  int shift = largest_exp + count_exp - kLimitExponent;
  if (shift <= 0) {
    return 0;
  }
  for (auto &d : dists) {
    d = std::ldexp(d, -shift);
  }
  return shift;
}

}  // namespace rdkit_cxx
