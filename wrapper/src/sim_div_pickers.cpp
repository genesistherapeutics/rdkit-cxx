// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/sim_div_pickers.h"

#include <cstdint>
#include <limits>
#include <mutex>
#include <vector>

// LeaderPicker.h chooses between two LeaderPickerStates, whose picks differ
// for a negative threshold, by the platform and RDK_BUILD_THREADSAFE_SSS,
// which RDConfig.h defines; included first, it chooses as for Python.
// clang-format off
#include <RDGeneral/RDConfig.h>
#include <SimDivPickers/LeaderPicker.h>
// clang-format on
#include <DataStructs/BitOps.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>
#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/sim_div_pickers.rs.h"
#include "rdkit-cxx/murtagh_checks.h"

// Python's module on Linux, a threaded build (S2), uses the threaded state.
#if defined(__linux__) && !defined(USE_THREADED_LEADERPICKER)
#error "LeaderPicker.h chose another LeaderPickerState than Python's"
#endif

namespace rdkit_cxx {

namespace {

// getDistFromLTM (DistPicker.cpp), through which the pickers read a
// distance matrix, indexes it with i * (i - 1) / 2 in unsigned int, which
// wraps for an item past 65536, so RDKit reads another pair's distance.
constexpr std::uint32_t kMaxMatrixPool = 65537;

void check_matrix(rust::Slice<const double> dist_mat, std::uint32_t pool_size) {
  if (pool_size > kMaxMatrixPool) {
    throw ValueErrorException(
        "at most 65537 items: RDKit indexes the matrix in unsigned int");
  }
  // Python's wrappers read poolSize * (poolSize - 1) / 2 values whatever
  // the array's length (a segfault in Python for a shorter one).
  std::uint64_t pool = pool_size;
  if (pool > 1 && dist_mat.size() < pool * (pool - 1) / 2) {
    throw ValueErrorException(
        "the distance matrix has fewer than poolSize * (poolSize - 1) / 2 "
        "values");
  }
}

// RDKit takes the first picks as ints and casts each back to unsigned, so
// an index past INT_MAX comes back as itself and is refused as any index
// past the pool is.
RDKit::INT_VECT first_pick_vect(rust::Slice<const std::uint32_t> picks) {
  RDKit::INT_VECT result;
  for (auto pick : picks) {
    result.push_back(static_cast<int>(pick));
  }
  return result;
}

rust::Vec<std::uint32_t> to_rust(const RDKit::INT_VECT &picks) {
  rust::Vec<std::uint32_t> result;
  for (auto pick : picks) {
    result.push_back(static_cast<std::uint32_t>(pick));
  }
  return result;
}

// The Rust function. A panic in it is the rust::Error call_distance_function
// throws, which the lazy pickers, which run on this thread, let through.
class RustDistance {
 public:
  explicit RustDistance(DistanceFunction &func) : d_func(func) {}
  double operator()(unsigned int i, unsigned int j) {
    return call_distance_function(d_func, i, j);
  }

 private:
  DistanceFunction &d_func;
};

// One minus the Tanimoto similarity, as pyBVFunctor in
// Code/SimDivPickers/Wrap/PickerHelpers.h computes it.
class BitVectDistance {
 public:
  explicit BitVectDistance(std::vector<const ExplicitBitVect *> bvs)
      : d_bvs(std::move(bvs)) {}
  double operator()(unsigned int i, unsigned int j) {
    return 1. - TanimotoSimilarity(*d_bvs[i], *d_bvs[j]);
  }

 private:
  std::vector<const ExplicitBitVect *> d_bvs;
};

// The pool's vectors. Python's wrappers read objects[i] for every i below
// poolSize, which raises IndexError past the end.
BitVectDistance pool_vectors(const BitVectList &objects,
                             std::uint32_t pool_size) {
  auto len = static_cast<unsigned>(objects.items.size());
  if (pool_size > objects.items.size()) {
    throw IndexError(len, len);
  }
  std::vector<const ExplicitBitVect *> bvs;
  bvs.reserve(pool_size);
  for (std::uint32_t i = 0; i < pool_size; ++i) {
    bvs.push_back(&as_explicit(*objects.items[i]));
  }
  return BitVectDistance(std::move(bvs));
}

// MaxMinPicker::lazyPick picks, in each round, the item farthest from the
// picks so far, starting from a distance of -1; when no item's distance is
// above -1 (or every one is NaN) it picks none and dereferences null (a
// segfault in Python for a distance function that returns -2, or NaN).
template <class F>
class MaxMinDistance {
 public:
  explicit MaxMinDistance(F &func) : d_func(func) {}
  double operator()(unsigned int i, unsigned int j) {
    double d = d_func(i, j);
    if (!(d > -1.0)) {
      throw ValueErrorException(
          "MaxMinPicker needs distances above -1 (and not NaN)");
    }
    return d;
  }

 private:
  F &d_func;
};

// MaxMinPicker::lazyPick sizes its per-item array as poolSize *
// sizeof(MaxMinPickInfo) in unsigned int, which wraps for a pool of 2^28
// items, and then writes every item's entry past the small array.
void check_max_min_pool(std::uint32_t pool_size) {
  if (std::uint64_t{pool_size} * sizeof(RDPickers::MaxMinPickInfo) >
      std::numeric_limits<unsigned>::max()) {
    throw ValueErrorException(
        "MaxMinPicker picks from at most 2^28 - 1 items: RDKit sizes its "
        "work array in unsigned int");
  }
}

template <class F>
PicksWithThreshold max_min_lazy(F &func, std::uint32_t pool_size,
                                std::uint32_t pick_size,
                                rust::Slice<const std::uint32_t> first_picks,
                                std::int32_t seed, double threshold) {
  check_max_min_pool(pool_size);
  MaxMinDistance<F> checked(func);
  RDPickers::MaxMinPicker picker;
  auto picks = picker.lazyPick(checked, pool_size, pick_size,
                               first_pick_vect(first_picks), seed, threshold);
  PicksWithThreshold result;
  result.picks = to_rust(picks);
  result.threshold = threshold;
  return result;
}

// LeaderPicker::lazyPick numbers the items in int and counts its blocks of
// items as (poolSize + 32767) / 32768 in unsigned int, which wraps for a
// pool past INT_MAX.
//
// With more than one thread RDKit compacts the pool on its own pthreads,
// created without checking for failure, where an exception from the
// distance function (a Rust panic, or bit vectors of different sizes) calls
// std::terminate (an abort in Python, and Python's own functions fail
// there without the GIL). The picks do not depend on the thread count, so
// one thread, this one, makes them.
template <class F>
rust::Vec<std::uint32_t> leader_lazy(
    F &func, std::uint32_t pool_size, double threshold, std::uint32_t pick_size,
    rust::Slice<const std::uint32_t> first_picks) {
  if (pool_size > static_cast<unsigned>(std::numeric_limits<int>::max())) {
    throw ValueErrorException(
        "LeaderPicker picks from at most 2^31 - 1 items: RDKit numbers them "
        "in int");
  }
  RDPickers::LeaderPicker picker;
  return to_rust(picker.lazyPick(func, pool_size, pick_size,
                                 first_pick_vect(first_picks), threshold, 1));
}

// A copy of the matrix's first poolSize * (poolSize - 1) / 2 values, which
// hc_ changes, checked and scaled for it (murtagh_checks.h). No distance
// comes back, so the scale is not needed.
std::vector<double> murtagh_matrix(rust::Slice<const double> dist_mat,
                                   std::uint32_t pool_size) {
  if (pool_size < 2) {
    throw ValueErrorException("Murtagh's clustering needs at least 2 items");
  }
  std::uint64_t pool = pool_size;
  auto len = static_cast<std::size_t>(pool * (pool - 1) / 2);
  if (dist_mat.size() < len) {
    throw ValueErrorException(
        "the distance matrix has fewer than poolSize * (poolSize - 1) / 2 "
        "values");
  }
  std::vector<double> dists(dist_mat.begin(), dist_mat.begin() + len);
  prepare_murtagh_distances(dists, pool);
  return dists;
}

// HierarchicalClusterPicker::cluster merges poolSize - pickSize times, and
// hc_ gives poolSize - 1 merges, so with no picks RDKit reads a merge of
// cluster 0, which it indexes as -1 (a segfault in Python).
void check_cluster_picks(std::uint32_t pick_size) {
  if (pick_size == 0) {
    throw ValueErrorException("pickSize must be at least 1");
  }
}

}  // namespace

// As MaxMinPicks in Code/SimDivPickers/Wrap/MaxMinPicker.cpp, with
// MaxMinPicker::pick's lazyPick of the matrix.
rust::Vec<std::uint32_t> max_min_picker_pick(
    rust::Slice<const double> dist_mat, std::uint32_t pool_size,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks,
    std::int32_t seed) {
  if (pick_size >= pool_size) {
    throw ValueErrorException("pickSize must be less than poolSize");
  }
  check_matrix(dist_mat, pool_size);
  RDPickers::distmatFunctor matrix(dist_mat.data());
  MaxMinDistance<RDPickers::distmatFunctor> checked(matrix);
  RDPickers::MaxMinPicker picker;
  return to_rust(picker.lazyPick(checked, pool_size, pick_size,
                                 first_pick_vect(first_picks), seed));
}

PicksWithThreshold max_min_picker_lazy_pick(
    DistanceFunction &func, std::uint32_t pool_size, std::uint32_t pick_size,
    rust::Slice<const std::uint32_t> first_picks, std::int32_t seed,
    double threshold) {
  RustDistance distance(func);
  return max_min_lazy(distance, pool_size, pick_size, first_picks, seed,
                      threshold);
}

PicksWithThreshold max_min_picker_lazy_bit_vector_pick(
    const BitVectList &objects, std::uint32_t pool_size,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks,
    std::int32_t seed, double threshold) {
  auto distance = pool_vectors(objects, pool_size);
  return max_min_lazy(distance, pool_size, pick_size, first_picks, seed,
                      threshold);
}

rust::Vec<std::uint32_t> leader_picker_lazy_pick(
    DistanceFunction &func, std::uint32_t pool_size, double threshold,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks) {
  RustDistance distance(func);
  return leader_lazy(distance, pool_size, threshold, pick_size, first_picks);
}

rust::Vec<std::uint32_t> leader_picker_lazy_bit_vector_pick(
    const BitVectList &objects, std::uint32_t pool_size, double threshold,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks) {
  auto distance = pool_vectors(objects, pool_size);
  return leader_lazy(distance, pool_size, threshold, pick_size, first_picks);
}

// As HierarchicalPicks in
// Code/SimDivPickers/Wrap/HierarchicalClusterPicker.cpp. RDKit's pick
// reads the matrix hc_ changed through getDistFromLTM.
rust::Vec<std::uint32_t> hierarchical_cluster_picker_pick(
    ClusterMethod method, rust::Slice<const double> dist_mat,
    std::uint32_t pool_size, std::uint32_t pick_size) {
  if (pick_size >= pool_size) {
    throw ValueErrorException("pickSize must be less than poolSize");
  }
  check_cluster_picks(pick_size);
  check_matrix(dist_mat, pool_size);
  auto dists = murtagh_matrix(dist_mat, pool_size);
  RDPickers::HierarchicalClusterPicker picker(method);
  std::lock_guard<std::mutex> lock(murtagh_mutex());
  return to_rust(picker.pick(dists.data(), pool_size, pick_size));
}

// As HierarchicalClusters, which leaves the sizes to cluster's precondition.
rust::Vec<IndexList> hierarchical_cluster_picker_cluster(
    ClusterMethod method, rust::Slice<const double> dist_mat,
    std::uint32_t pool_size, std::uint32_t pick_size) {
  PRECONDITION(pool_size >= pick_size,
               "pickSize cannot be larger than the poolSize");
  check_cluster_picks(pick_size);
  auto dists = murtagh_matrix(dist_mat, pool_size);
  RDPickers::HierarchicalClusterPicker picker(method);
  std::lock_guard<std::mutex> lock(murtagh_mutex());
  return to_index_lists(picker.cluster(dists.data(), pool_size, pick_size));
}

}  // namespace rdkit_cxx
