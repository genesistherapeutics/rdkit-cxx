// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.Clusters`, a Python module: the `Cluster`
//! class for hierarchical cluster trees.
//!
//! A [`Cluster`] owns its children, so a tree changes only through its
//! root. Python sorts a cluster's children in place whenever they are read
//! (by metric) and whenever its points are gathered (by size), and later
//! reads see that order; the methods here sort the same way, which is why
//! they take `&mut self`. Python caches a cluster's points, which then go
//! stale when the tree below changes (a cluster given a single child by
//! `AddChild` keeps itself as its only point); [`Cluster::get_points`]
//! gathers them afresh.
//!
//! ```no_run
//! use rdkit_cxx::ml::cluster::clusters::{Cluster, ClusterParams};
//!
//! let mut root = Cluster::new_with_params(&ClusterParams {
//!     index: 1,
//!     metric: 1000.0,
//!     ..Default::default()
//! });
//! root.add_child(Cluster::new_with_params(&ClusterParams {
//!     index: 10,
//!     ..Default::default()
//! }));
//! ```

use std::cmp::Ordering;

use crate::py::format_f;
use crate::{Error, Result};

/// Python 2's `cmp`: -1, 0, or 1 as `t1` is less than, equal to, or greater
/// than `t2` (0 when they are unordered). Mirrors `Clusters.cmp`.
pub fn cmp<T: PartialOrd + ?Sized>(t1: &T, t2: &T) -> i32 {
    if t1 < t2 {
        -1
    } else if t1 > t2 {
        1
    } else {
        0
    }
}

/// The tolerance within which [`Cluster::compare_with_ignore_extras`] counts
/// metrics equal. Mirrors `Clusters.CMPTOL`.
pub fn cmptol() -> f64 {
    1e-6
}

/// Parameters for [`Cluster::new_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ClusterParams {
    /// How far apart the cluster's children are.
    pub metric: f64,
    /// The cluster's children.
    pub children: Vec<Cluster>,
    /// The cluster's location.
    pub position: Vec<f64>,
    /// The order in which the cluster was made.
    pub index: i64,
    /// The cluster's name; `None` names it by its index.
    pub name: Option<String>,
    /// The point's data. Python's can be any object; RDKit's code stores
    /// point indices, and numbers are what `Print` can show.
    pub data: Option<f64>,
}

impl Default for ClusterParams {
    fn default() -> ClusterParams {
        ClusterParams {
            metric: 0.0,
            children: Vec::new(),
            position: Vec::new(),
            index: -1,
            name: None,
            data: None,
        }
    }
}

/// Parameters for [`Cluster::print_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct PrintParams {
    /// The depth of the cluster, which indents it by two spaces a level.
    pub level: usize,
    /// Shows each cluster's data, when it has some.
    pub show_data: bool,
    /// The text between a cluster's name and its metric.
    pub offset: String,
}

impl Default for PrintParams {
    fn default() -> PrintParams {
        PrintParams {
            level: 0,
            show_data: false,
            offset: "\t".to_owned(),
        }
    }
}

/// A node of a hierarchical cluster tree, whose leaves are the clustered
/// points. Mirrors `Clusters.Cluster`.
#[derive(Clone, Debug, PartialEq)]
pub struct Cluster {
    metric: f64,
    children: Vec<Cluster>,
    position: Vec<f64>,
    index: i64,
    name: Option<String>,
    data: Option<f64>,
    len: usize,
}

impl Default for Cluster {
    /// A cluster with Python's defaults. Mirrors `Clusters.Cluster()`.
    fn default() -> Cluster {
        Cluster::new_with_params(&ClusterParams::default())
    }
}

fn by_metric(a: &Cluster, b: &Cluster) -> Ordering {
    a.metric.partial_cmp(&b.metric).unwrap_or(Ordering::Equal)
}

impl Cluster {
    /// Mirrors `Clusters.Cluster(metric, children, position, index, name,
    /// data)`.
    pub fn new_with_params(params: &ClusterParams) -> Cluster {
        Cluster::from_params(params.clone())
    }

    /// As [`new_with_params`](Cluster::new_with_params), taking the
    /// children rather than copying them.
    pub(crate) fn from_params(params: ClusterParams) -> Cluster {
        let mut cluster = Cluster {
            metric: params.metric,
            children: params.children,
            position: params.position,
            index: params.index,
            name: params.name,
            data: params.data,
            len: 1,
        };
        cluster.update_length();
        cluster
    }

    fn update_length(&mut self) {
        self.len = self.children.iter().map(|c| c.len).sum::<usize>() + 1;
    }

    /// Sets the metric. Mirrors `Cluster.SetMetric`.
    pub fn set_metric(&mut self, metric: f64) {
        self.metric = metric;
    }

    /// How far apart the cluster's children are. Mirrors `Cluster.GetMetric`.
    pub fn get_metric(&self) -> f64 {
        self.metric
    }

    /// Sets the index. Mirrors `Cluster.SetIndex`.
    pub fn set_index(&mut self, index: i64) {
        self.index = index;
    }

    /// The order in which the cluster was made. Mirrors `Cluster.GetIndex`.
    pub fn get_index(&self) -> i64 {
        self.index
    }

    /// Sets the position. Mirrors `Cluster.SetPosition`.
    pub fn set_position(&mut self, pos: Vec<f64>) {
        self.position = pos;
    }

    /// The cluster's location. Mirrors `Cluster.GetPosition`.
    pub fn get_position(&self) -> &[f64] {
        &self.position
    }

    /// Sets the data. Mirrors `Cluster.SetData`.
    pub fn set_data(&mut self, data: Option<f64>) {
        self.data = data;
    }

    /// The point's data. Mirrors `Cluster.GetData`.
    pub fn get_data(&self) -> Option<f64> {
        self.data
    }

    /// Sets the name. Mirrors `Cluster.SetName`.
    pub fn set_name(&mut self, name: Option<String>) {
        self.name = name;
    }

    /// The cluster's name, or `Cluster(<index>)` when it has none. Mirrors
    /// `Cluster.GetName`.
    pub fn get_name(&self) -> String {
        match &self.name {
            Some(name) => name.clone(),
            None => format!("Cluster({})", self.index),
        }
    }

    /// The number of nodes in the tree under the cluster, the cluster
    /// included. Mirrors `len(cluster)`.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Always `false`: a cluster counts itself.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Whether the cluster is a leaf, a single point. Mirrors
    /// `Cluster.IsTerminal`.
    pub fn is_terminal(&self) -> bool {
        self.len <= 1
    }

    /// The children, sorted in place by metric (a stable sort). Mirrors
    /// `Cluster.GetChildren`.
    pub fn get_children(&mut self) -> &[Cluster] {
        self.children.sort_by(by_metric);
        &self.children
    }

    /// The children, in their current order.
    pub(crate) fn children(&self) -> &[Cluster] {
        &self.children
    }

    /// The children, for reordering only: their sizes must not change.
    pub(crate) fn children_mut(&mut self) -> &mut [Cluster] {
        &mut self.children
    }

    /// Adds a child. Mirrors `Cluster.AddChild`.
    pub fn add_child(&mut self, child: Cluster) {
        self.children.push(child);
        self.gen_points();
        self.update_length();
    }

    /// Adds children. Mirrors `Cluster.AddChildren`.
    pub fn add_children(&mut self, children: Vec<Cluster>) {
        self.children.extend(children);
        self.gen_points();
        self.update_length();
    }

    /// Removes and returns the child at `which`, in the children's current
    /// order (the order [`get_children`](Self::get_children) last gave, or
    /// the order they were added in). Mirrors `Cluster.RemoveChild`, which
    /// takes the child itself.
    ///
    /// `Error::Index` for an index past the children.
    pub fn remove_child(&mut self, which: usize) -> Result<Cluster> {
        if which >= self.children.len() {
            return Err(Error::index(which, self.children.len()));
        }
        let child = self.children.remove(which);
        self.update_length();
        Ok(child)
    }

    /// Sorts the children as Python's `_GenPoints` does before gathering
    /// their points: by metric, then by size, largest first (stable sorts),
    /// throughout the tree. The (Python-)terminal clusters are left alone.
    pub(crate) fn gen_points(&mut self) {
        if self.len == 1 {
            return;
        }
        self.children.sort_by(by_metric);
        self.children.sort_by_key(|c| std::cmp::Reverse(c.len));
        for child in &mut self.children {
            child.gen_points();
        }
    }

    fn collect_points<'a>(&'a self, out: &mut Vec<&'a Cluster>) {
        if self.len == 1 {
            out.push(self);
        } else {
            for child in &self.children {
                child.collect_points(out);
            }
        }
    }

    /// The leaves under the cluster (itself, for a leaf), the children of
    /// each cluster taken largest first. Mirrors `Cluster.GetPoints`.
    pub fn get_points(&mut self) -> Vec<&Cluster> {
        self.gen_points();
        let mut points = Vec::new();
        self.collect_points(&mut points);
        points
    }

    /// The positions of [`get_points`](Self::get_points)' leaves. Mirrors
    /// `Cluster.GetPointsPositions`.
    pub fn get_points_positions(&mut self) -> Vec<&[f64]> {
        self.get_points()
            .into_iter()
            .map(|c| c.position.as_slice())
            .collect()
    }

    /// The first cluster in the tree, searched depth first, with the index
    /// `index`. Mirrors `Cluster.FindSubtree`.
    pub fn find_subtree(&self, index: i64) -> Option<&Cluster> {
        if index == self.index {
            return Some(self);
        }
        self.children.iter().find_map(|c| c.find_subtree(index))
    }

    /// The tree, one line per cluster (its name and metric), as
    /// `Cluster.Print()` prints it. Mirrors `Cluster.Print`, which prints
    /// the text this returns.
    pub fn print(&mut self) -> String {
        self.print_with_params(&PrintParams::default())
    }

    /// Mirrors `Cluster.Print(level, showData, offset)`.
    pub fn print_with_params(&mut self, params: &PrintParams) -> String {
        let mut out = String::new();
        self.print_to(params.level, params.show_data, &params.offset, &mut out);
        out
    }

    fn print_to(&mut self, level: usize, show_data: bool, offset: &str, out: &mut String) {
        let indent = "  ".repeat(level);
        match self.data.filter(|_| show_data) {
            Some(data) => out.push_str(&format!(
                "{indent}{}{offset} Data: {}\t Metric: {}\n",
                self.get_name(),
                format_f(data),
                format_f(self.metric)
            )),
            None => out.push_str(&format!(
                "{indent}{}{offset} Metric: {}\n",
                self.get_name(),
                format_f(self.metric)
            )),
        }
        self.children.sort_by(by_metric);
        for child in &mut self.children {
            child.print_to(level + 1, show_data, offset, out);
        }
    }

    /// Compares two trees by shape: -1, 0, or 1, as Python 2's `cmp`.
    /// Mirrors `Cluster.Compare(other)`.
    pub fn compare(&mut self, other: &mut Cluster) -> i32 {
        self.compare_with_ignore_extras(other, true)
    }

    /// Compares two trees: by size, then (unless `ignore_extras`) by metric
    /// (beyond [`cmptol`]), name, and position, then child by child, the
    /// children sorted by metric. Mirrors `Cluster.Compare(other,
    /// ignoreExtras)`.
    pub fn compare_with_ignore_extras(&mut self, other: &mut Cluster, ignore_extras: bool) -> i32 {
        let tv = cmp(&self.len, &other.len);
        if tv != 0 {
            return tv;
        }
        if !ignore_extras {
            let (m1, m2) = (self.metric, other.metric);
            if (m1 - m2).abs() > cmptol() {
                return cmp(&m1, &m2);
            }
            let names = cmp(&self.get_name(), &other.get_name());
            if names != 0 {
                return names;
            }
            let r = cmp(&self.position.len(), &other.position.len());
            if r != 0 {
                return r;
            }
            let r = list_cmp(&self.position, &other.position);
            if r != 0 {
                return r;
            }
        }
        self.children.sort_by(by_metric);
        other.children.sort_by(by_metric);
        let r = cmp(&self.children.len(), &other.children.len());
        if r != 0 {
            return r;
        }
        for (c1, c2) in self.children.iter_mut().zip(other.children.iter_mut()) {
            let t = c1.compare_with_ignore_extras(c2, ignore_extras);
            if t != 0 {
                return t;
            }
        }
        0
    }
}

/// Python's comparison of two lists of floats: by the first items that
/// differ (by `==`), else by length.
fn list_cmp(a: &[f64], b: &[f64]) -> i32 {
    for (x, y) in a.iter().zip(b) {
        if x != y {
            return cmp(x, y);
        }
    }
    cmp(&a.len(), &b.len())
}
