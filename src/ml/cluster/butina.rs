// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.Butina`, a Python module: Butina's
//! sphere-exclusion clustering (J. Chem. Inf. Comput. Sci. 39, 747-750,
//! 1999).
//!
//! Each function returns the clusters as lists of point indices, the
//! centroid first, in Python's order. Distances are compared with the
//! threshold as Python compares them, so the clusters are Python's. With
//! reordering, Python turns the neighbor lists it recounts into sets, and
//! the members of later clusters come in the order CPython's sets iterate
//! integers in; the port keeps its sets the way CPython lays them out.
//!
//! ```no_run
//! use rdkit_cxx::ml::cluster::butina::cluster_data_flat;
//!
//! // The lower triangle of the distance matrix, row by row.
//! let dists = [0.1, 0.9, 0.8, 0.2, 0.9, 0.95];
//! let clusters = cluster_data_flat(&dists, 4, 0.3).unwrap();
//! ```

use crate::py::np_sum;
use crate::rd_logger::Logger;
use crate::{Error, Result};

/// The Euclidean distance between two points. Mirrors `Butina.EuclideanDist`.
///
/// RDKit 2024.03.6's Python returns the array of the coordinates' absolute
/// differences (`sqrt(dv * dv)`), which `ClusterData` cannot compare with
/// its threshold for points of two or more coordinates (NumPy raises
/// `ValueError`); this computes what it means to, the square root of the
/// sum of the squared differences.
///
/// As NumPy broadcasts, a point of one coordinate is compared with every
/// coordinate of the other; `Error::Value` for other points of different
/// lengths, as NumPy raises.
pub fn euclidean_dist(pi: &[f64], pj: &[f64]) -> Result<f64> {
    let n = pi.len().max(pj.len());
    if pi.len() != pj.len() && pi.len() != 1 && pj.len() != 1 {
        return Err(Error::Value {
            message: format!(
                "points of {} and {} coordinates cannot be compared",
                pi.len(),
                pj.len()
            ),
        });
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

/// Clusters the first `n_pts` points by their Euclidean distances
/// ([`euclidean_dist`]); points within `dist_thresh` of each other are
/// neighbors. Mirrors `Butina.ClusterData(data, nPts, distThresh)`.
///
/// RDKit 2024.03.6's Python raises `ValueError` for two or more points of
/// two or more coordinates (see [`euclidean_dist`]); this clusters them as
/// it means to.
///
/// `Error::Index` when a pair of points to compare is past the end of
/// `data`, as Python raises; errors from [`euclidean_dist`] pass through.
pub fn cluster_data(data: &[Vec<f64>], n_pts: usize, dist_thresh: f64) -> Result<Vec<Vec<usize>>> {
    cluster_data_with_reordering(data, n_pts, dist_thresh, false)
}

/// As [`cluster_data`], recounting after each cluster the unassigned
/// neighbors of the points it touched when `reordering`, so the next
/// centroid is always the point with the most unassigned neighbors.
/// Mirrors `Butina.ClusterData(data, nPts, distThresh,
/// reordering=reordering)`.
pub fn cluster_data_with_reordering(
    data: &[Vec<f64>],
    n_pts: usize,
    dist_thresh: f64,
    reordering: bool,
) -> Result<Vec<Vec<usize>>> {
    cluster_data_with_dist_func(
        data,
        n_pts,
        dist_thresh,
        |a, b| euclidean_dist(a, b),
        reordering,
    )
}

/// Clusters `n_pts` points by a distance matrix stored as its lower
/// triangle, row by row: the distance between points `i > j` at
/// `i * (i - 1) / 2 + j`. Mirrors `Butina.ClusterData(dists, nPts,
/// distThresh, isDistData=True)`.
///
/// As Python, logs `Distance matrix is too long` to `rdApp.warning` when
/// `data` holds more than `n_pts * (n_pts - 1) / 2` distances, and uses the
/// first ones; `Error::Index` when it holds fewer, as Python raises.
pub fn cluster_data_flat(data: &[f64], n_pts: usize, dist_thresh: f64) -> Result<Vec<Vec<usize>>> {
    cluster_data_flat_with_reordering(data, n_pts, dist_thresh, false)
}

/// As [`cluster_data_flat`], with `reordering` as in
/// [`cluster_data_with_reordering`]. Mirrors `Butina.ClusterData(dists,
/// nPts, distThresh, isDistData=True, reordering=reordering)`.
pub fn cluster_data_flat_with_reordering(
    data: &[f64],
    n_pts: usize,
    dist_thresh: f64,
    reordering: bool,
) -> Result<Vec<Vec<usize>>> {
    // Python compares with its float division `nPts * (nPts - 1) / 2`.
    if (data.len() as f64) > (n_pts as f64) * (n_pts as f64 - 1.0) / 2.0 {
        Logger.warning("Distance matrix is too long");
    }
    let mut dm_idx = 0;
    butina(n_pts, dist_thresh, reordering, |_, _| {
        let dij = *data.get(dm_idx).ok_or(Error::index(dm_idx, data.len()))?;
        dm_idx += 1;
        Ok(dij)
    })
}

/// Clusters the first `n_pts` items of `data` by the distances
/// `dist_func` gives, called for each pair `(data[i], data[j])` with
/// `j < i`, in order. Mirrors `Butina.ClusterData(data, nPts, distThresh,
/// distFunc=distFunc, reordering=reordering)`.
///
/// `Error::Index` when a pair to compare is past the end of `data`, as
/// Python raises; errors from `dist_func` pass through.
pub fn cluster_data_with_dist_func<T>(
    data: &[T],
    n_pts: usize,
    dist_thresh: f64,
    mut dist_func: impl FnMut(&T, &T) -> Result<f64>,
    reordering: bool,
) -> Result<Vec<Vec<usize>>> {
    butina(n_pts, dist_thresh, reordering, |i, j| {
        let item = |k: usize| data.get(k).ok_or(Error::index(k, data.len()));
        dist_func(item(i)?, item(j)?)
    })
}

/// The clustering itself, from the distance `dist(i, j)` of each pair
/// `j < i`, asked for in Python's order.
fn butina(
    n_pts: usize,
    dist_thresh: f64,
    reordering: bool,
    mut dist: impl FnMut(usize, usize) -> Result<f64>,
) -> Result<Vec<Vec<usize>>> {
    let mut nbr_lists: Vec<Neighbors> = vec![Neighbors::List(Vec::new()); n_pts];
    for i in 0..n_pts {
        for j in 0..i {
            if dist(i, j)? <= dist_thresh {
                nbr_lists[i].push(j);
                nbr_lists[j].push(i);
            }
        }
    }
    // (number of neighbors, index), largest first, as Python sorts the
    // tuples in reverse.
    let sort = |indices: &mut [(usize, usize)]| indices.sort_by(|a, b| b.cmp(a));
    let mut t_lists: Vec<(usize, usize)> = nbr_lists
        .iter()
        .enumerate()
        .map(|(idx, nbrs)| (nbrs.len(), idx))
        .collect();
    sort(&mut t_lists);
    // Python pops from the front; `front` marks the next item instead.
    let mut res = Vec::new();
    let mut seen = vec![false; n_pts];
    let mut front = 0;
    while front < t_lists.len() {
        let idx = t_lists[front].1;
        front += 1;
        if seen[idx] {
            continue;
        }
        // As in Python, the centroid is not marked seen; nothing that could
        // still see it is left.
        let mut t_res = vec![idx];
        for nbr in nbr_lists[idx].iter() {
            if !seen[nbr] {
                t_res.push(nbr);
                seen[nbr] = true;
            }
        }
        if reordering {
            let mut nbr_nbr = vec![false; n_pts];
            for &t in &t_res {
                for nbr in nbr_lists[t].iter() {
                    nbr_nbr[nbr] = true;
                }
            }
            for entry in &mut t_lists[front..] {
                let y1 = entry.1;
                if seen[y1] || !nbr_nbr[y1] {
                    continue;
                }
                // `set(nbrLists[y1]).difference(tRes)`.
                let mut nbrs = match &nbr_lists[y1] {
                    Neighbors::List(list) => IntSet::from_list(list),
                    Neighbors::Set(set) => set.copy(),
                }
                .copy();
                nbrs.difference_update(&t_res);
                *entry = (nbrs.used, y1);
                nbr_lists[y1] = Neighbors::Set(nbrs);
            }
            sort(&mut t_lists[front..]);
        }
        res.push(t_res);
    }
    Ok(res)
}

/// A point's neighbors: Python's list, until reordering makes it a set.
#[derive(Clone)]
enum Neighbors {
    List(Vec<usize>),
    Set(IntSet),
}

impl Neighbors {
    fn push(&mut self, nbr: usize) {
        match self {
            Neighbors::List(list) => list.push(nbr),
            Neighbors::Set(_) => unreachable!("sets are made after the lists are"),
        }
    }

    fn len(&self) -> usize {
        match self {
            Neighbors::List(list) => list.len(),
            Neighbors::Set(set) => set.used,
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = usize> + '_> {
        match self {
            Neighbors::List(list) => Box::new(list.iter().copied()),
            Neighbors::Set(set) => Box::new(set.iter()),
        }
    }
}

/// A slot of an [`IntSet`]'s hash table.
#[derive(Clone, Copy, PartialEq)]
enum Slot {
    Unused,
    Dummy,
    Key(usize),
}

/// A set of non-negative integers laid out as CPython's `set` lays them
/// out (`Objects/setobject.c`), so it iterates in CPython's order: an
/// integer's hash is itself, and the table starts with 8 slots, probes 9
/// slots linearly before perturbing, and resizes as CPython's does.
#[derive(Clone)]
struct IntSet {
    table: Vec<Slot>,
    mask: usize,
    /// Active and dummy slots.
    fill: usize,
    /// Active slots.
    used: usize,
}

const MIN_SIZE: usize = 8;
const LINEAR_PROBES: usize = 9;
const PERTURB_SHIFT: u32 = 5;

impl IntSet {
    fn new() -> IntSet {
        IntSet {
            table: vec![Slot::Unused; MIN_SIZE],
            mask: MIN_SIZE - 1,
            fill: 0,
            used: 0,
        }
    }

    /// `set(list)`: adds the keys in order.
    fn from_list(keys: &[usize]) -> IntSet {
        let mut set = IntSet::new();
        for &key in keys {
            set.add(key);
        }
        set
    }

    fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.table.iter().filter_map(|slot| match slot {
            Slot::Key(key) => Some(*key),
            _ => None,
        })
    }

    /// `set_insert_clean`: puts a key known to be absent in a table with
    /// no dummies.
    fn insert_clean(table: &mut [Slot], mask: usize, key: usize) {
        let mut perturb = key;
        let mut i = key & mask;
        loop {
            if table[i] == Slot::Unused {
                table[i] = Slot::Key(key);
                return;
            }
            if i + LINEAR_PROBES <= mask {
                for j in 1..=LINEAR_PROBES {
                    if table[i + j] == Slot::Unused {
                        table[i + j] = Slot::Key(key);
                        return;
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb) & mask;
        }
    }

    /// `set_table_resize`: the smallest table of more than `min_used`
    /// slots, with the active keys reinserted in slot order.
    fn resize(&mut self, min_used: usize) {
        let mut size = MIN_SIZE;
        while size <= min_used {
            size <<= 1;
        }
        let old = std::mem::replace(&mut self.table, vec![Slot::Unused; size]);
        self.mask = size - 1;
        for slot in old {
            if let Slot::Key(key) = slot {
                IntSet::insert_clean(&mut self.table, self.mask, key);
            }
        }
        self.fill = self.used;
    }

    /// `set_add_entry`.
    fn add(&mut self, key: usize) {
        let mask = self.mask;
        let mut i = key & mask;
        let mut perturb = key;
        let mut free_slot = None;
        loop {
            let probes = if i + LINEAR_PROBES <= mask {
                LINEAR_PROBES
            } else {
                0
            };
            for j in i..=i + probes {
                match self.table[j] {
                    Slot::Unused => {
                        if let Some(free) = free_slot {
                            self.table[free] = Slot::Key(key);
                            self.used += 1;
                            return;
                        }
                        self.table[j] = Slot::Key(key);
                        self.fill += 1;
                        self.used += 1;
                        if self.fill * 5 >= mask * 3 {
                            self.resize(self.used * if self.used > 50000 { 2 } else { 4 });
                        }
                        return;
                    }
                    Slot::Key(k) if k == key => return,
                    Slot::Dummy if free_slot.is_none() => free_slot = Some(j),
                    _ => {}
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb) & mask;
        }
    }

    /// `set.copy()` (`set_merge` into a new set).
    fn copy(&self) -> IntSet {
        let mut set = IntSet::new();
        if self.used == 0 {
            return set;
        }
        if (set.fill + self.used) * 5 >= set.mask * 3 {
            set.resize((set.used + self.used) * 2);
        }
        if set.mask == self.mask && self.fill == self.used {
            set.table = self.table.clone();
        } else {
            for key in self.iter() {
                IntSet::insert_clean(&mut set.table, set.mask, key);
            }
        }
        set.fill = self.used;
        set.used = self.used;
        set
    }

    /// `set_difference_update_internal` for a list: discards each key, then
    /// resizes the dummies away when more than a quarter of the table holds
    /// them.
    fn difference_update(&mut self, keys: &[usize]) {
        for &key in keys {
            if let Some(slot) = self.table.iter_mut().find(|s| **s == Slot::Key(key)) {
                *slot = Slot::Dummy;
                self.used -= 1;
            }
        }
        if self.fill - self.used > self.mask / 4 {
            self.resize(self.used * if self.used > 50000 { 2 } else { 4 });
        }
    }
}
