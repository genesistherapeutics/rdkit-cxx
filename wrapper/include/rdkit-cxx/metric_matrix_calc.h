/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMetricMatrixCalc: the lower triangles of distance and similarity
// matrices (DataManip/MetricMatrixCalc), and the borrowed lists of bit
// vectors that it and the pickers take.

#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "rdkit-cxx/data_structs.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Bit vectors Rust lends for one call: pointers to vectors Rust owns, so a
// list never outlives the call it is made for (with_bit_vect_list in
// src/data_manip/metric/rd_metric_matrix_calc.rs).
class BitVectList {
 public:
  std::vector<const BitVect *> items;
};

std::unique_ptr<BitVectList> bit_vect_list_new();
void bit_vect_list_push(BitVectList &list, const BitVect &bv);

// RDKit compares two bit vectors of different sizes by folding the larger
// (SimilarityWrapper in DataStructs/BitOps.h), dividing their sizes, which
// is a division by zero (SIGFPE in Python) when one has no bits, and
// FoldFingerprint computes the sizes in int, which a SparseBitVect of more
// than INT_MAX bits overflows. Throws ValueErrorException for a list whose
// sizes differ unless every size is between 1 and INT_MAX.
void check_foldable(const BitVectList &list);

// The values are the rows, nrows of ncols each.
rust::Vec<double> get_euclidean_dist_mat(rust::Slice<const double> values,
                                         std::uint32_t nrows,
                                         std::uint32_t ncols);
rust::Vec<double> get_tanimoto_dist_mat(const BitVectList &list);
rust::Vec<double> get_tanimoto_sim_mat(const BitVectList &list);

}  // namespace rdkit_cxx
