// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.Murtagh`, a Python module: Murtagh's
//! hierarchical clustering ([`clustering`](super::clustering)) returning
//! cluster trees.
//!
//! ```no_run
//! use rdkit_cxx::ml::cluster::murtagh::{cluster_data_dists, upgma};
//!
//! let tree = cluster_data_dists(&[1.0, 4.0, 3.5], 3, upgma()).unwrap();
//! ```

use crate::ml::cluster::clustering::{murtagh_cluster, murtagh_dist_cluster};
use crate::ml::cluster::clusters::{Cluster, ClusterParams};
use crate::{Error, Result};

/// Ward's minimum variance. Mirrors `Murtagh.WARDS`.
pub fn wards() -> i32 {
    1
}

/// Single linkage. Mirrors `Murtagh.SLINK`.
pub fn slink() -> i32 {
    2
}

/// Complete linkage. Mirrors `Murtagh.CLINK`.
pub fn clink() -> i32 {
    3
}

/// Group average linkage. Mirrors `Murtagh.UPGMA`.
pub fn upgma() -> i32 {
    4
}

/// McQuitty's method. Mirrors `Murtagh.MCQUITTY`.
pub fn mcquitty() -> i32 {
    5
}

/// Gower's median method. Mirrors `Murtagh.GOWER`.
pub fn gower() -> i32 {
    6
}

/// The centroid method. Mirrors `Murtagh.CENTROID`.
pub fn centroid() -> i32 {
    7
}

/// The methods offered to users, as `(name, method, description)`. Mirrors
/// `Murtagh.methods`, which leaves out McQuitty's and Gower's.
pub fn methods() -> Vec<(&'static str, i32, &'static str)> {
    vec![
        (
            "Ward's Minimum Variance",
            wards(),
            "Ward's Minimum Variance",
        ),
        ("Average Linkage", upgma(), "Group Average Linkage (UPGMA)"),
        ("Single Linkage", slink(), "Single Linkage (SLINK)"),
        ("Complete Linkage", clink(), "Complete Linkage (CLINK)"),
        ("Centroid", centroid(), "Centroid method"),
    ]
}

/// The tree of Murtagh's merges, as Python's `_ToClusters` builds it:
/// leaf `i` has index `i + 1` and data `i`, and merge `i` makes a cluster
/// with index `n_pts + i + 1`.
fn to_clusters(n_pts: usize, ia: &[i64], ib: &[i64], crit: &[f64]) -> Result<Vec<Cluster>> {
    let mut cs: Vec<Option<Cluster>> = (0..n_pts)
        .map(|i| {
            Some(Cluster::from_params(ClusterParams {
                metric: 0.0,
                data: Some(i as f64),
                index: i as i64 + 1,
                ..Default::default()
            }))
        })
        .collect();
    fn take(cs: &mut [Option<Cluster>], idx: i64) -> Result<Cluster> {
        usize::try_from(idx - 1)
            .ok()
            .and_then(|i| cs.get_mut(i))
            .and_then(Option::take)
            .ok_or_else(|| Error::Runtime {
                message: format!("the clustering merged cluster {idx} twice"),
            })
    }
    let mut last = None;
    for i in 0..ia.len().saturating_sub(1) {
        let c1 = take(&mut cs, ia[i])?;
        let c2 = take(&mut cs, ib[i])?;
        // Python puts the merged cluster where the first one was, and
        // returns the last one made.
        let at = (ia[i] - 1) as usize;
        cs[at] = Some(Cluster::from_params(ClusterParams {
            metric: crit[i],
            children: vec![c1, c2],
            index: (n_pts + i + 1) as i64,
            ..Default::default()
        }));
        last = Some(at);
    }
    last.and_then(|at| cs[at].take())
        .map(|c| vec![c])
        .ok_or_else(|| Error::Value {
            message: "there were no merges".to_owned(),
        })
}

/// Clusters the first `n_pts` rows of `data` (points of equal length) with
/// `method` (one of [`wards`] to [`centroid`]), returning a single-item
/// list holding the cluster tree; the rows after them are ignored, as in
/// Python. Mirrors `Murtagh.ClusterData(data, nPts, method)`.
///
/// Errors as [`murtagh_cluster`].
pub fn cluster_data(data: &[Vec<f64>], n_pts: i32, method: i32) -> Result<Vec<Cluster>> {
    let sz = data.first().map_or(0, |row| row.len());
    let sz = i32::try_from(sz).map_err(|_| Error::Value {
        message: "too many coordinates".to_owned(),
    })?;
    let (ia, ib, crit) = murtagh_cluster(data, n_pts, sz, method)?;
    to_clusters(n_pts.max(0) as usize, &ia, &ib, &crit)
}

/// Clusters `n_pts` items by a condensed distance matrix, the distance
/// between items `i < j` at `j * (j - 1) / 2 + i`. Mirrors
/// `Murtagh.ClusterData(data, nPts, method, isDistData=1)`.
///
/// Errors as [`murtagh_dist_cluster`].
pub fn cluster_data_dists(data: &[f64], n_pts: i32, method: i32) -> Result<Vec<Cluster>> {
    let (ia, ib, crit) = murtagh_dist_cluster(data, n_pts, method)?;
    to_clusters(n_pts.max(0) as usize, &ia, &ib, &crit)
}
