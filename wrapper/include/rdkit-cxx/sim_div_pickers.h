/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdSimDivPickers: the MaxMin, Leader, and hierarchical-clustering
// diversity pickers (Code/SimDivPickers). The pickers keep no state RDKit
// uses between calls, so each call makes its own.

#pragma once

#include <cstdint>

#include <SimDivPickers/HierarchicalClusterPicker.h>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/metric_matrix_calc.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using ClusterMethod = RDPickers::HierarchicalClusterPicker::ClusterMethod;

// Defined by the bridge (src/ffi/sim_div_pickers.rs): a Rust distance
// function lent for one call, and a pick list with MaxMin's last distance.
struct DistanceFunction;
struct PicksWithThreshold;

// MaxMinPicker. RDKit picks the first item with its own generator, seeded
// from std::random_device for a negative seed, not RDKit's shared one.
rust::Vec<std::uint32_t> max_min_picker_pick(
    rust::Slice<const double> dist_mat, std::uint32_t pool_size,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks,
    std::int32_t seed);
PicksWithThreshold max_min_picker_lazy_pick(
    DistanceFunction &func, std::uint32_t pool_size, std::uint32_t pick_size,
    rust::Slice<const std::uint32_t> first_picks, std::int32_t seed,
    double threshold);
PicksWithThreshold max_min_picker_lazy_bit_vector_pick(
    const BitVectList &objects, std::uint32_t pool_size,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks,
    std::int32_t seed, double threshold);

// LeaderPicker.
rust::Vec<std::uint32_t> leader_picker_lazy_pick(
    DistanceFunction &func, std::uint32_t pool_size, double threshold,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks);
rust::Vec<std::uint32_t> leader_picker_lazy_bit_vector_pick(
    const BitVectList &objects, std::uint32_t pool_size, double threshold,
    std::uint32_t pick_size, rust::Slice<const std::uint32_t> first_picks);

// HierarchicalClusterPicker.
rust::Vec<std::uint32_t> hierarchical_cluster_picker_pick(
    ClusterMethod method, rust::Slice<const double> dist_mat,
    std::uint32_t pool_size, std::uint32_t pick_size);
rust::Vec<IndexList> hierarchical_cluster_picker_cluster(
    ClusterMethod method, rust::Slice<const double> dist_mat,
    std::uint32_t pool_size, std::uint32_t pick_size);

}  // namespace rdkit_cxx
