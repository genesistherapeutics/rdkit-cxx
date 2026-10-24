/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks for the binary formats RDKit parses without validating them. RDKit
// reads these pickles' sizes and counts as given, so malformed bytes can make
// it divide by zero, read past its buffers, or build an object whose later
// use does. Each check throws ValueErrorException for bytes RDKit would
// misread and leaves the rest, which RDKit rejects itself, to RDKit.

#pragma once

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string_view>
#include <utility>

#include <DataStructs/SparseIntVect.h>
#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

// ExplicitBitVect and SparseBitVect pickles.
void check_bit_vect_pickle(std::string_view pkl);
// DiscreteValueVect pickles.
void check_discrete_value_vect_pickle(std::string_view pkl);
// RealValueVect pickles.
void check_real_value_vect_pickle(std::string_view pkl);
// UniformGrid3D pickles, including the DiscreteValueVect inside.
void check_uniform_grid3d_pickle(std::string_view pkl);
// UniformRealValueGrid3D pickles, including the RealValueVect inside.
void check_uniform_real_value_grid3d_pickle(std::string_view pkl);

// RDKit sums a SparseIntVect's values, and their absolute values, in int, so
// the absolute values must sum to at most INT_MAX (which also rules out
// abs(INT_MIN)).
template <class IndexType>
void check_sparse_int_vect_sums(const RDKit::SparseIntVect<IndexType> &siv) {
  long long sum = 0;
  for (const auto &[idx, val] : siv.getNonzeroElements()) {
    sum += std::llabs(static_cast<long long>(val));
    if (sum > std::numeric_limits<int>::max()) {
      throw ValueErrorException("SparseIntVect values sum past INT_MAX");
    }
  }
}

// A SparseIntVect read from a pickle: RDKit never checks its length or its
// indices, which the wrapper uses to size and index its own buffers.
template <class IndexType>
void check_sparse_int_vect(const RDKit::SparseIntVect<IndexType> &siv) {
  if (std::cmp_less(siv.getLength(), 0)) {
    throw ValueErrorException("negative length in SparseIntVect pickle");
  }
  for (const auto &[idx, val] : siv.getNonzeroElements()) {
    if (std::cmp_less(idx, 0) || !std::cmp_less(idx, siv.getLength())) {
      throw ValueErrorException("index out of range in SparseIntVect pickle");
    }
  }
}

}  // namespace rdkit_cxx
