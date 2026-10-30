// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.ClusterUtils`, a Python module: functions
//! that walk and split cluster trees.
//!
//! Like Python's, these sort the children of the clusters they visit in
//! place (see [`clusters`](super::clusters)). Python also marks the
//! clusters they visit with private attributes (`_aboveCentroid`,
//! `_clustCenter`, `_distToCenter`), which no public code reads; the ports
//! leave them out.

use std::cmp::Ordering;

use crate::ml::cluster::clusters::Cluster;
use crate::{Error, Result};

fn by_size_descending(a: &Cluster, b: &Cluster) -> Ordering {
    b.len().cmp(&a.len())
}

fn by_metric_descending(a: &Cluster, b: &Cluster) -> Ordering {
    b.get_metric()
        .partial_cmp(&a.get_metric())
        .unwrap_or(Ordering::Equal)
}

/// Sorts every cluster's children as `GetNodeList` does: by metric, then
/// by size, largest first.
fn sort_for_node_list(cluster: &mut Cluster) {
    if cluster.len() == 1 {
        return;
    }
    cluster.get_children();
    let children = cluster.children_mut();
    children.sort_by(by_size_descending);
    for child in children.iter_mut() {
        sort_for_node_list(child);
    }
}

fn collect_nodes<'a>(cluster: &'a Cluster, out: &mut Vec<&'a Cluster>) {
    if cluster.len() == 1 {
        out.push(cluster);
        return;
    }
    for child in cluster.children() {
        collect_nodes(child, out);
    }
    out.push(cluster);
}

/// Every node under the cluster, children before their parent, the
/// children of each node taken largest first. Mirrors
/// `ClusterUtils.GetNodeList`.
pub fn get_node_list(cluster: &mut Cluster) -> Vec<&Cluster> {
    sort_for_node_list(cluster);
    let mut nodes = Vec::new();
    collect_nodes(cluster, &mut nodes);
    nodes
}

/// Every node under the cluster, as [`get_node_list`] gives them. Mirrors
/// `ClusterUtils.GetNodesDownToCentroids(cluster)`, whose `above` argument
/// only sets the nodes' private `_aboveCentroid` attributes, and which
/// stops at no centroid, as no RDKit code marks one.
pub fn get_nodes_down_to_centroids(cluster: &mut Cluster) -> Vec<&Cluster> {
    get_node_list(cluster)
}

/// The data of the point in the cluster with the smallest summed distance
/// to the others, by the condensed distance matrix `dists` (the distance
/// between points `i < j` at `j * (j - 1) / 2 + i`), or -1 when every sum
/// reaches 1e24. Each point's data is its index into `dists`. Mirrors
/// `ClusterUtils.FindClusterCentroidFromDists`, which in RDKit 2025.09.6
/// raises `TypeError` for a cluster of two or more points (it indexes
/// `dists` with `/`, a float; RDKit 2026.03 uses `//`); this returns what
/// the fixed function does.
///
/// `Error::Value` for a point whose data is not a whole number from zero up
/// (Python indexes from the end for negative numbers, and fails for others);
/// `Error::Index` for a point past the end of `dists`.
pub fn find_cluster_centroid_from_dists(cluster: &mut Cluster, dists: &[f64]) -> Result<i64> {
    let pts = cluster
        .get_points()
        .into_iter()
        .map(|c| match c.get_data() {
            Some(data) if data >= 0.0 && data.fract() == 0.0 && data < 2f64.powi(53) => {
                Ok(data as i64)
            }
            _ => Err(Error::Value {
                message: "each point needs data, its index into the distance matrix".to_owned(),
            }),
        })
        .collect::<Result<Vec<i64>>>()?;
    let dist = |a: i64, b: i64| -> Result<f64> {
        let (row, col) = if b > a { (a, b) } else { (b, a) };
        let idx = col
            .checked_mul(col - 1)
            .and_then(|x| (x / 2).checked_add(row))
            .and_then(|x| usize::try_from(x).ok())
            .unwrap_or(usize::MAX);
        dists
            .get(idx)
            .copied()
            .ok_or(Error::index(idx, dists.len()))
    };
    let mut best = 1e24;
    let mut best_idx = -1;
    for &pt in &pts {
        let mut accum = 0.0;
        for &other in &pts {
            if other != pt {
                accum += dist(pt, other)?;
                if accum >= best {
                    break;
                }
            }
        }
        if accum < best {
            best = accum;
            best_idx = pt;
        }
    }
    // Python then reads each point's distance to the centroid, which fails
    // for a point past the end of `dists`.
    for &pt in &pts {
        if pt != best_idx && best_idx >= 0 {
            dist(pt, best_idx)?;
        }
    }
    Ok(best_idx)
}

/// A path from the root to a node: child indices, in the children's order
/// at the time.
type Path = Vec<usize>;

fn node_at<'a>(root: &'a Cluster, path: &[usize]) -> &'a Cluster {
    path.iter().fold(root, |node, &i| &node.children()[i])
}

fn node_at_mut<'a>(root: &'a mut Cluster, path: &[usize]) -> &'a mut Cluster {
    path.iter()
        .fold(root, |node, &i| &mut node.children_mut()[i])
}

fn check_split(cluster: &Cluster, n: usize) -> Result<()> {
    if cluster.len() < n {
        return Err(Error::Value {
            message: format!(
                "Cannot split cluster of length {} into {n} pieces",
                cluster.len()
            ),
        });
    }
    Ok(())
}

fn assertion_error() -> Error {
    Error::Value {
        message: "no cluster left to split".to_owned(),
    }
}

/// Python's `_BreadthFirstSplit`, as paths.
fn breadth_first_split(cluster: &mut Cluster, n: usize) -> Result<Vec<Path>> {
    let mut clusters: Vec<Path> = vec![Vec::new()];
    let mut nxt_idx = 0;
    for _ in 0..n.saturating_sub(1) {
        while nxt_idx < clusters.len() && node_at(cluster, &clusters[nxt_idx]).len() == 1 {
            nxt_idx += 1;
        }
        if nxt_idx >= clusters.len() {
            return Err(assertion_error());
        }
        let path = clusters.remove(nxt_idx);
        let node = node_at_mut(cluster, &path);
        node.get_children();
        node.children_mut().sort_by(by_metric_descending);
        for i in 0..node.children().len() {
            let mut child = path.clone();
            child.push(i);
            clusters.push(child);
        }
    }
    Ok(clusters)
}

/// Python's `_HeightFirstSplit`, as paths.
fn height_first_split(cluster: &mut Cluster, n: usize) -> Result<Vec<Path>> {
    let mut clusters: Vec<Path> = vec![Vec::new()];
    for _ in 0..n.saturating_sub(1) {
        let mut nxt_idx = 0;
        while nxt_idx < clusters.len() && node_at(cluster, &clusters[nxt_idx]).len() == 1 {
            nxt_idx += 1;
        }
        if nxt_idx >= clusters.len() {
            return Err(assertion_error());
        }
        let path = clusters[nxt_idx].clone();
        let node = node_at_mut(cluster, &path);
        node.get_children();
        for i in 0..node.children().len() {
            let mut child = path.clone();
            child.push(i);
            clusters.push(child);
        }
        clusters.remove(nxt_idx);
        let root: &Cluster = cluster;
        clusters.sort_by(|a, b| by_metric_descending(node_at(root, a), node_at(root, b)));
    }
    Ok(clusters)
}

/// Splits a cluster tree into `n` branches, splitting breadth first: the
/// branches come from splitting, in turn, each non-leaf branch in the
/// order they were found, children with larger metrics first. Mirrors
/// `ClusterUtils.SplitIntoNClusters(cluster, n)`.
///
/// `Error::Value` when the tree has fewer than `n` nodes, as Python raises,
/// or when it runs out of branches to split (Python's assertion).
pub fn split_into_n_clusters(cluster: &mut Cluster, n: usize) -> Result<Vec<&Cluster>> {
    split_into_n_clusters_with_breadth_first(cluster, n, true)
}

/// As [`split_into_n_clusters`]; with `breadth_first` false, each split
/// takes the branch with the largest metric that is not a leaf. Mirrors
/// `ClusterUtils.SplitIntoNClusters(cluster, n, breadthFirst)`.
pub fn split_into_n_clusters_with_breadth_first(
    cluster: &mut Cluster,
    n: usize,
    breadth_first: bool,
) -> Result<Vec<&Cluster>> {
    check_split(cluster, n)?;
    if cluster.len() == n {
        return Ok(cluster.get_points());
    }
    let paths = if breadth_first {
        breadth_first_split(cluster, n)?
    } else {
        height_first_split(cluster, n)?
    };
    let root: &Cluster = cluster;
    Ok(paths.iter().map(|path| node_at(root, path)).collect())
}
