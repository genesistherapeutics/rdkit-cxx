// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.ML.Cluster.Clustering`: Murtagh's hierarchical
//! clustering.
//!
//! `option` picks the criterion, as `rdkit.ML.Cluster.Murtagh` names them:
//! 1 Ward's minimum variance, 2 single link, 3 complete link, 4 average
//! link (UPGMA), 5 McQuitty's, 6 median (Gower's), and 7 centroid. Each
//! function returns the merges as `(ia, ib, crit)`, each of `n_pts`
//! entries: for `i` below `n_pts - 1`, clusters `ia[i]` and `ib[i]`
//! (numbered from 1, the merged cluster keeping `ia[i]`'s number) merge at
//! `crit[i]`, and the last entries are 0.
//!
//! The C code (`hc.c`, which Python calls unguarded) reads and writes out
//! of bounds for fewer than 2 items, for a distance matrix shorter than
//! `n_pts`' (Python's reads past the array), and when no distance is below
//! its 1e20 sentinel; these functions refuse the first two and distances
//! that are not finite with `Error::Value`, and scale distances that could
//! reach the sentinel down by an exact power of two, scaling `crit` back.
//! They also refuse an `option` outside 1 to 7, for which the C code
//! applies no criterion. Calls hold a lock, since the C code keeps its
//! state in static variables.
//!
//! ```no_run
//! use rdkit_cxx::ml::cluster::clustering::murtagh_dist_cluster;
//!
//! let (ia, ib, crit) = murtagh_dist_cluster(&[1.0, 4.0, 3.5], 3, 2).unwrap();
//! assert_eq!((ia, ib, crit), (vec![1, 1, 0], vec![2, 3, 0], vec![1.0, 3.5, 0.0]));
//! ```

use crate::ffi::cluster::ffi;
use crate::{Error, Result};

fn merges(clusters: ffi::MurtaghClusters) -> (Vec<i64>, Vec<i64>, Vec<f64>) {
    (clusters.ia, clusters.ib, clusters.crit)
}

/// Clusters points by their squared Euclidean distances: the first `n_pts *
/// sz` values of `data`, row by row, as `n_pts` points of `sz` coordinates,
/// as Python's binding reads them. When the rows have `sz` values, those are
/// the first `n_pts` rows. Mirrors `Clustering.MurtaghCluster`.
///
/// `Error::Value` for fewer values, which Python's binding reads past, for
/// rows of different lengths, which NumPy cannot make an array of, and as
/// the module describes.
pub fn murtagh_cluster(
    data: &[Vec<f64>],
    n_pts: i32,
    sz: i32,
    option: i32,
) -> Result<(Vec<i64>, Vec<i64>, Vec<f64>)> {
    let ncols = data.first().map_or(0, |row| row.len());
    if data.iter().any(|row| row.len() != ncols) {
        return Err(Error::Value {
            message: "the rows must all have the same length".to_owned(),
        });
    }
    let (Ok(rows), Ok(cols)) = (u64::try_from(n_pts), u64::try_from(sz)) else {
        return Err(Error::Value {
            message: "nPts and sz must not be negative".to_owned(),
        });
    };
    let values = data.concat();
    let needed = rows * cols;
    if (values.len() as u64) < needed {
        return Err(Error::Value {
            message: format!(
                "data holds {} values, fewer than the nPts * sz = {needed} it reads",
                values.len()
            ),
        });
    }
    Ok(merges(ffi::murtagh_cluster(
        &values[..needed as usize],
        rows,
        cols,
        n_pts,
        sz,
        option,
    )?))
}

/// Clusters items by a condensed distance matrix: the lower triangle by
/// rows, the distance between items `i > j` at `i * (i - 1) / 2 + j`, of
/// which the first `n_pts * (n_pts - 1) / 2` values are read. Mirrors
/// `Clustering.MurtaghDistCluster`, which also changes the caller's array.
///
/// `Error::Runtime` for an empty `data`, as Python raises; `Error::Value`
/// as the module describes.
pub fn murtagh_dist_cluster(
    data: &[f64],
    n_pts: i32,
    option: i32,
) -> Result<(Vec<i64>, Vec<i64>, Vec<f64>)> {
    Ok(merges(ffi::murtagh_dist_cluster(data, n_pts, option)?))
}
