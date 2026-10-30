// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/cluster.h"

#include <cmath>
#include <cstdint>
#include <mutex>
#include <vector>

#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

#include "rdkit-cxx/ffi/cluster.rs.h"
#include "rdkit-cxx/murtagh_checks.h"

// hcdriver.c, in RDKit's hc library, which installs no header for it.
extern "C" int distdriver_(long *n, long *len, double *dists, long *toggle,
                           long *ia, long *ib, double *crit);

namespace rdkit_cxx {

namespace {

// As clusterit and distclusterit in Code/ML/Cluster/Murtagh/nbClustering.cpp,
// with the distances checked and scaled first (murtagh_checks.h). The
// distances are a copy: hc_ changes them, which Python's does in the
// caller's array.
MurtaghClusters cluster_distances(std::vector<double> dists, long n,
                                  std::int32_t option) {
  check_murtagh_method(option);
  int shift = prepare_murtagh_distances(dists, static_cast<std::uint64_t>(n));
  std::vector<long> ia(n, 0);
  std::vector<long> ib(n, 0);
  std::vector<double> crit(n, 0.0);
  long len = n * (n - 1) / 2;
  long iopt = option;
  {
    std::lock_guard<std::mutex> lock(murtagh_mutex());
    distdriver_(&n, &len, dists.data(), &iopt, ia.data(), ib.data(),
                crit.data());
  }
  MurtaghClusters result;
  for (long i = 0; i < n; ++i) {
    result.ia.push_back(ia[i]);
    result.ib.push_back(ib[i]);
    result.crit.push_back(std::ldexp(crit[i], shift));
  }
  return result;
}

}  // namespace

// As Clustering_MurtaghCluster: the distances are the rows' squared
// Euclidean distances, summed in clusterit's order.
MurtaghClusters murtagh_cluster(rust::Slice<const double> data,
                                std::uint64_t nrows, std::uint64_t ncols,
                                std::int32_t n_pts, std::int32_t sz,
                                std::int32_t option) {
  CHECK_INVARIANT(nrows == static_cast<std::uint64_t>(n_pts),
                  "data row count does not match nPts");
  CHECK_INVARIANT(ncols == static_cast<std::uint64_t>(sz),
                  "data column count does not match sz");
  if (n_pts < 2) {
    throw ValueErrorException("Murtagh's clustering needs at least 2 items");
  }
  long n = n_pts;
  long m = sz;
  std::vector<double> dists(static_cast<std::size_t>(n * (n - 1) / 2), 0.0);
  std::size_t pos = 0;
  for (long i = 1; i < n; ++i) {
    for (long j = 0; j < i; ++j) {
      for (long k = 0; k < m; ++k) {
        double tmp = data[i * m + k] - data[j * m + k];
        dists[pos] += tmp * tmp;
      }
      ++pos;
    }
  }
  return cluster_distances(std::move(dists), n, option);
}

// As Clustering_MurtaghDistCluster, which reads n_pts * (n_pts - 1) / 2
// values whatever the array's length (a segfault in Python for a shorter
// one).
MurtaghClusters murtagh_dist_cluster(rust::Slice<const double> data,
                                     std::int32_t n_pts, std::int32_t option) {
  CHECK_INVARIANT(data.size() > 0, "distance matrix must not be empty");
  if (n_pts < 2) {
    throw ValueErrorException("Murtagh's clustering needs at least 2 items");
  }
  long n = n_pts;
  auto len = static_cast<std::size_t>(n * (n - 1) / 2);
  if (data.size() < len) {
    throw ValueErrorException(
        "the distance matrix has fewer than nPts * (nPts - 1) / 2 values");
  }
  return cluster_distances(
      std::vector<double>(data.begin(), data.begin() + len), n, option);
}

}  // namespace rdkit_cxx
