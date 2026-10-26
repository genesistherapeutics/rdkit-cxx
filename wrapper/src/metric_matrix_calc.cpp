// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/metric_matrix_calc.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <variant>
#include <vector>

#include <DataManip/MetricMatrixCalc/MetricFuncs.h>
#include <DataManip/MetricMatrixCalc/MetricMatrixCalc.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

namespace rdkit_cxx {

namespace {

using RDDataManip::MetricMatrixCalc;

// calcMetricMatrix (MetricMatrixCalc.h) indexes the triangle with
// i * (i - 1) / 2 in unsigned int, which wraps for a row past 65536, so
// RDKit writes the distances of more rows over those of earlier ones.
constexpr std::uint32_t kMaxItems = 65537;

void check_num_items(std::uint32_t n) {
  if (n > kMaxItems) {
    throw ValueErrorException(
        "at most 65537 items: RDKit indexes the matrix in unsigned int");
  }
}

// The triangle's storage. calcMetricMatrix refuses a null pointer, which an
// empty vector may give, so it always has room for one value.
std::vector<double> triangle(std::uint32_t n) {
  std::size_t len = static_cast<std::size_t>(n) * (n - 1) / 2;
  std::vector<double> values(std::max<std::size_t>(len, 1));
  values.resize(len);
  return values;
}

rust::Vec<double> to_rust(const std::vector<double> &values) {
  rust::Vec<double> result;
  result.reserve(values.size());
  for (auto v : values) {
    result.push_back(v);
  }
  return result;
}

// The list's vectors as RDKit's, for calcMetricMatrix, which reads them
// with operator[] as Python's NbSequenceHolder gives them.
template <class T>
class BitVectRefs {
 public:
  explicit BitVectRefs(const BitVectList &list) : d_list(list) {}
  const T &operator[](unsigned int i) const {
    return std::get<T>(d_list.items[i]->v);
  }

 private:
  const BitVectList &d_list;
};

// As getTanimotoDistMat and getTanimotoSimMat in
// Code/DataManip/MetricMatrixCalc/Wrap/rdMetricMatrixCalc.cpp.
template <bool distance>
rust::Vec<double> tanimoto_mat(const BitVectList &list) {
  auto n = static_cast<std::uint32_t>(list.items.size());
  CHECK_INVARIANT(n > 1, "");
  check_num_items(n);
  check_foldable(list);
  auto values = triangle(n);
  std::visit(
      [&](const auto &first) {
        using T = std::decay_t<decltype(first)>;
        BitVectRefs<T> refs(list);
        MetricMatrixCalc<BitVectRefs<T>, T> calc;
        if constexpr (distance) {
          calc.setMetricFunc(&RDDataManip::TanimotoDistanceMetric<T, T>);
        } else {
          calc.setMetricFunc(&RDDataManip::TanimotoSimilarityMetric<T, T>);
        }
        calc.calcMetricMatrix(refs, n, 0, values.data());
      },
      list.items[0]->v);
  return to_rust(values);
}

}  // namespace

std::unique_ptr<BitVectList> bit_vect_list_new() {
  return std::make_unique<BitVectList>();
}

void bit_vect_list_push(BitVectList &list, const BitVect &bv) {
  list.items.push_back(&bv);
}

void check_foldable(const BitVectList &list) {
  if (list.items.empty()) {
    return;
  }
  auto size = [](const BitVect *bv) {
    return std::visit([](const auto &x) { return x.getNumBits(); }, bv->v);
  };
  auto [lo, hi] = std::minmax_element(
      list.items.begin(), list.items.end(),
      [&](const BitVect *a, const BitVect *b) { return size(a) < size(b); });
  auto smallest = size(*lo);
  auto largest = size(*hi);
  if (smallest == largest) {
    return;
  }
  if (smallest == 0 ||
      largest > static_cast<unsigned>(std::numeric_limits<int>::max())) {
    throw ValueErrorException(
        "bit vectors of different sizes must have between 1 and 2^31 - 1 "
        "bits: RDKit folds the larger by the ratio of their sizes");
  }
}

// As getEuclideanDistMat's sequence path in
// Code/DataManip/MetricMatrixCalc/Wrap/rdMetricMatrixCalc.cpp; its NumPy
// path converts each value to double the same way.
rust::Vec<double> get_euclidean_dist_mat(rust::Slice<const double> values,
                                         std::uint32_t nrows,
                                         std::uint32_t ncols) {
  CHECK_INVARIANT(nrows > 0, "Empty list passed in");
  check_num_items(nrows);
  if (values.size() != static_cast<std::size_t>(nrows) * ncols) {
    throw ValueErrorException("All subsequences must be the same length");
  }
  std::vector<const double *> rows(nrows);
  for (std::uint32_t i = 0; i < nrows; ++i) {
    rows[i] = values.data() + static_cast<std::size_t>(i) * ncols;
  }
  auto result = triangle(nrows);
  MetricMatrixCalc<std::vector<const double *>, const double *> calc;
  calc.setMetricFunc(
      &RDDataManip::EuclideanDistanceMetric<const double *, const double *>);
  calc.calcMetricMatrix(rows, nrows, ncols, result.data());
  return to_rust(result);
}

rust::Vec<double> get_tanimoto_dist_mat(const BitVectList &list) {
  return tanimoto_mat<true>(list);
}

rust::Vec<double> get_tanimoto_sim_mat(const BitVectList &list) {
  return tanimoto_mat<false>(list);
}

}  // namespace rdkit_cxx
