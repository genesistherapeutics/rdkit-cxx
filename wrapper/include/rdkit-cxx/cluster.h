/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdkit.ML.Cluster.Clustering: Murtagh's hierarchical clustering
// (Code/ML/Cluster/Murtagh), as its wrapper nbClustering.cpp calls it.

#pragma once

#include <cstdint>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared struct, defined by the bridge (src/ffi/cluster.rs).
struct MurtaghClusters;

// The rows are nrows of ncols values each.
MurtaghClusters murtagh_cluster(rust::Slice<const double> data,
                                std::uint64_t nrows, std::uint64_t ncols,
                                std::int32_t n_pts, std::int32_t sz,
                                std::int32_t option);
MurtaghClusters murtagh_dist_cluster(rust::Slice<const double> data,
                                     std::int32_t n_pts, std::int32_t option);

}  // namespace rdkit_cxx
