// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.Butina`, a Python module: Butina's
//! sphere-exclusion clustering (J. Chem. Inf. Comput. Sci. 39, 747-750,
//! 1999).
//!
//! Each function returns the clusters as lists of point indices, the
//! centroid first, in Python's order. Distances are compared with the
//! threshold as Python compares them, so the clusters are Python's.
//!
//! ```no_run
//! use rdkit_cxx::ml::cluster::butina::cluster_data_flat;
//!
//! // The lower triangle of the distance matrix, row by row.
//! let dists = [0.1, 0.9, 0.8, 0.2, 0.9, 0.95];
//! let clusters = cluster_data_flat(&dists, 4, 0.3).unwrap();
//! ```

use crate::py::np_sum;
use crate::{Error, Result};

fn value_error(message: String) -> Error {
    Error::Value { message }
}

/// The Euclidean distance between two points. Mirrors `Butina.EuclideanDist`.
///
/// As NumPy broadcasts, a point of one coordinate is compared with every
/// coordinate of the other; `Error::Value` for other points of different
/// lengths, as NumPy raises.
pub fn euclidean_dist(pi: &[f64], pj: &[f64]) -> Result<f64> {
    let n = pi.len().max(pj.len());
    if pi.len() != pj.len() && pi.len() != 1 && pj.len() != 1 {
        return Err(value_error(format!(
            "points of {} and {} coordinates cannot be compared",
            pi.len(),
            pj.len()
        )));
    }
    let at = |p: &[f64], k: usize| if p.len() == 1 { p[0] } else { p[k] };
    let squares: Vec<f64> = (0..n)
        .map(|k| {
            let d = at(pi, k) - at(pj, k);
            d.powi(2)
        })
        .collect();
    Ok(np_sum(&squares).sqrt())
}

/// The `n_pts` by `n_pts` matrix of distances between the first `n_pts`
/// items, with zeros on the diagonal. Mirrors
/// `Butina.compute_distance_matrix`, which calls `dist_func(data[i],
/// data[j])` for `j < i` only.
///
/// `Error::Index` when `data` has fewer than `n_pts` items, as Python
/// raises; errors from `dist_func` pass through.
pub fn compute_distance_matrix<T>(
    data: &[T],
    n_pts: usize,
    mut dist_func: impl FnMut(&T, &T) -> Result<f64>,
) -> Result<Vec<Vec<f64>>> {
    if n_pts > data.len() {
        return Err(Error::index(data.len(), data.len()));
    }
    let mut matrix = vec![vec![0.0; n_pts]; n_pts];
    for i in 0..n_pts {
        for j in 0..i {
            let d = dist_func(&data[i], &data[j])?;
            matrix[i][j] = d;
            matrix[j][i] = d;
        }
    }
    Ok(matrix)
}

/// Parameters for [`cluster_data_with_params`], with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClusterDataParams {
    /// `data` is the `n_pts` by `n_pts` distance matrix, not points.
    pub is_dist_data: bool,
    /// After each cluster, recounts the unassigned neighbors of the points
    /// it touched, so the next centroid is always the point with the most
    /// unassigned neighbors.
    pub reordering: bool,
}

/// Clusters the first `n_pts` points by their Euclidean distances; points
/// within `dist_thresh` of each other are neighbors. Mirrors
/// `Butina.ClusterData(data, nPts, distThresh)`.
///
/// Errors as [`compute_distance_matrix`] and [`euclidean_dist`].
pub fn cluster_data(data: &[Vec<f64>], n_pts: usize, dist_thresh: f64) -> Result<Vec<Vec<usize>>> {
    cluster_data_with_params(data, n_pts, dist_thresh, &ClusterDataParams::default())
}

/// Mirrors `Butina.ClusterData(data, nPts, distThresh, isDistData,
/// reordering=reordering)` for points or an `n_pts` by `n_pts` distance
/// matrix, which, as in Python, is used as it is, row `i` giving point
/// `i`'s neighbors.
///
/// `Error::Value` for a matrix of another shape, as Python raises; for
/// points, errors as [`cluster_data`].
pub fn cluster_data_with_params(
    data: &[Vec<f64>],
    n_pts: usize,
    dist_thresh: f64,
    params: &ClusterDataParams,
) -> Result<Vec<Vec<usize>>> {
    if params.is_dist_data {
        if data.len() != n_pts || data.iter().any(|row| row.len() != n_pts) {
            return Err(value_error(format!(
                "the distance matrix is not a matrix of the required shape ({n_pts}, {n_pts})"
            )));
        }
        Ok(butina(data.to_vec(), n_pts, dist_thresh, params.reordering))
    } else {
        let matrix = compute_distance_matrix(data, n_pts, |a, b| euclidean_dist(a, b))?;
        Ok(butina(matrix, n_pts, dist_thresh, params.reordering))
    }
}

/// Clusters `n_pts` points by a distance matrix stored as its lower
/// triangle, row by row: the distance between points `i > j` at
/// `i * (i - 1) / 2 + j`. Mirrors `Butina.ClusterData(dists, nPts,
/// distThresh, isDistData=True)`.
///
/// `Error::Value` when `data` does not hold `n_pts * (n_pts - 1) / 2`
/// distances, as Python raises.
pub fn cluster_data_flat(data: &[f64], n_pts: usize, dist_thresh: f64) -> Result<Vec<Vec<usize>>> {
    cluster_data_flat_with_reordering(data, n_pts, dist_thresh, false)
}

/// As [`cluster_data_flat`], with `reordering` as in
/// [`ClusterDataParams`]. Mirrors `Butina.ClusterData(dists, nPts,
/// distThresh, isDistData=True, reordering=reordering)`.
#[allow(clippy::needless_range_loop)]
pub fn cluster_data_flat_with_reordering(
    data: &[f64],
    n_pts: usize,
    dist_thresh: f64,
    reordering: bool,
) -> Result<Vec<Vec<usize>>> {
    if Some(data.len()) != n_pts.checked_mul(n_pts.saturating_sub(1)).map(|n| n / 2) {
        return Err(value_error(
            "Mismatched input data dimension and nPts".to_owned(),
        ));
    }
    let mut matrix = vec![vec![0.0; n_pts]; n_pts];
    let mut k = 0;
    for i in 0..n_pts {
        for j in 0..i {
            // NumPy adds the transpose, so each distance is 0 + d.
            matrix[i][j] = data[k] + 0.0;
            matrix[j][i] = 0.0 + data[k];
            k += 1;
        }
    }
    Ok(butina(matrix, n_pts, dist_thresh, reordering))
}

/// Clusters the first `n_pts` items of `data` by the distances
/// `dist_func` gives. Mirrors `Butina.ClusterData(data, nPts, distThresh,
/// distFunc=distFunc, reordering=reordering)`.
///
/// Errors as [`compute_distance_matrix`].
pub fn cluster_data_with_dist_func<T>(
    data: &[T],
    n_pts: usize,
    dist_thresh: f64,
    dist_func: impl FnMut(&T, &T) -> Result<f64>,
    reordering: bool,
) -> Result<Vec<Vec<usize>>> {
    let matrix = compute_distance_matrix(data, n_pts, dist_func)?;
    Ok(butina(matrix, n_pts, dist_thresh, reordering))
}

/// The clustering itself, from the full distance matrix.
fn butina(
    matrix: Vec<Vec<f64>>,
    n_pts: usize,
    dist_thresh: f64,
    reordering: bool,
) -> Vec<Vec<usize>> {
    let mut neighbor_lists: Vec<Vec<usize>> = matrix
        .iter()
        .map(|row| (0..n_pts).filter(|&j| row[j] <= dist_thresh).collect())
        .collect();
    // (number of neighbors, index), largest first, as Python sorts the
    // tuples in reverse.
    let sort = |indices: &mut Vec<(usize, usize)>| indices.sort_by(|a, b| b.cmp(a));
    let mut sorted_indices: Vec<(usize, usize)> = neighbor_lists
        .iter()
        .enumerate()
        .map(|(idx, neighbors)| (neighbors.len(), idx))
        .collect();
    sort(&mut sorted_indices);
    // Python pops from the front; `front` marks the next item instead.
    let mut clusters = Vec::new();
    let mut seen = vec![false; n_pts];
    let mut front = 0;
    while front < sorted_indices.len() && sorted_indices[front].0 > 1 {
        let idx = sorted_indices[front].1;
        front += 1;
        if seen[idx] {
            continue;
        }
        let mut cluster = vec![idx];
        seen[idx] = true;
        for &neighbor in &neighbor_lists[idx] {
            if !seen[neighbor] {
                cluster.push(neighbor);
                seen[neighbor] = true;
            }
        }
        if reordering {
            sorted_indices.drain(..front);
            front = 0;
            let mut affected = vec![false; n_pts];
            for &point in &cluster {
                for &neighbor in &neighbor_lists[point] {
                    if !seen[neighbor] {
                        affected[neighbor] = true;
                    }
                }
            }
            for element in &mut sorted_indices {
                let point = element.1;
                if affected[point] {
                    let new_neighbors: Vec<usize> = neighbor_lists[point]
                        .iter()
                        .copied()
                        .filter(|&nbr| !seen[nbr])
                        .collect();
                    *element = (new_neighbors.len(), point);
                    neighbor_lists[point] = new_neighbors;
                }
            }
            sort(&mut sorted_indices);
        }
        clusters.push(cluster);
    }
    for &(_, idx) in &sorted_indices[front..] {
        if !seen[idx] {
            clusters.push(vec![idx]);
        }
    }
    clusters
}
