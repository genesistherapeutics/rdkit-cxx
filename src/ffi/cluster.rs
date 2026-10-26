// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// The merges `distdriver_` returns: clusters `ia[i]` and `ib[i]`
    /// (1-based) merge at `crit[i]`, for `i` below `n - 1`.
    struct MurtaghClusters {
        ia: Vec<i64>,
        ib: Vec<i64>,
        crit: Vec<f64>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/cluster.h");

        fn murtagh_cluster(
            data: &[f64],
            nrows: u64,
            ncols: u64,
            n_pts: i32,
            sz: i32,
            option: i32,
        ) -> Result<MurtaghClusters>;
        fn murtagh_dist_cluster(data: &[f64], n_pts: i32, option: i32) -> Result<MurtaghClusters>;
    }
}
