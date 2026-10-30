// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.SimDivFilters.rdSimDivPickers`: pickers of diverse
//! subsets.
//!
//! Each picks from a pool of `pool_size` items, numbered from 0, by their
//! distances: a condensed distance matrix (the lower triangle by rows, the
//! distance between items `i > j` at `i * (i - 1) / 2 + j`, of which the
//! first `pool_size * (pool_size - 1) / 2` values are read, as
//! `rd_metric_matrix_calc` computes them), a distance function of two
//! items' indices, or the Tanimoto distance (1 - similarity) between the
//! first `pool_size` of some `ExplicitBitVect`s. The distance functions run
//! on the calling thread, and a panic in one ends the pick with
//! `Error::Runtime`. Bit vectors of different sizes are `Error::Value` when
//! compared, as in Python, and a pool larger than the bit vectors given is
//! `Error::Index`.
//!
//! ```no_run
//! use rdkit_cxx::sim_div_filters::MaxMinPicker;
//!
//! let points: [f64; 5] = [0.0, 1.0, 5.0, 6.0, 10.0];
//! let picks = MaxMinPicker::default()
//!     .lazy_pick(|i, j| (points[i as usize] - points[j as usize]).abs(), 5, 3)
//!     .unwrap();
//! assert_eq!(picks.len(), 3);
//! ```

use crate::Result;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::data_manip::metric::rd_metric_matrix_calc::with_bit_vect_list;
use crate::data_structs::ExplicitBitVect;
use crate::ffi::sim_div_pickers::{DistanceFunction, ffi};

/// The optional arguments of `MaxMinPicker`'s methods, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MaxMinPickParams {
    /// Items to pick first, in order. With none, RDKit picks the first item
    /// at random.
    pub first_picks: Vec<u32>,
    /// Seeds the random first pick. A negative seed seeds it from the
    /// system's random device, so the picks differ from call to call.
    pub seed: i32,
}

impl Default for MaxMinPickParams {
    fn default() -> MaxMinPickParams {
        MaxMinPickParams {
            first_picks: Vec::new(),
            seed: -1,
        }
    }
}

/// Picks a diverse subset by the MaxMin algorithm (Ashton et al., Quant.
/// Struct.-Act. Relat. 21 (2002), 598-604): each pick is the item farthest
/// from those already picked. Mirrors `rdSimDivPickers.MaxMinPicker`.
///
/// `Error::Value`, as RDKit raises, for an empty pool, a `pick_size` above
/// `pool_size`, and a first pick past the pool. RDKit picks by the largest
/// distance above -1 and dereferences null when no item's is (a segfault in
/// Python), so every distance a pick reads must be above -1 and not NaN,
/// else `Error::Value`; so is a pool of 2^28 items or more, for which RDKit
/// sizes its work array in `unsigned int`, which wraps. `pick_size` 0 picks
/// one item, as RDKit does.
#[derive(Debug, Default)]
pub struct MaxMinPicker {
    _private: (),
}

impl MaxMinPicker {
    /// Picks by a distance matrix, with no first picks and a random first
    /// item. Mirrors `MaxMinPicker.Pick`.
    ///
    /// `Error::Value` for a `pick_size` that is not below `pool_size`, as
    /// Python raises; for a matrix shorter than the pool's (read past its
    /// end in Python); and for more than 65537 items, past which RDKit's
    /// `unsigned int` index into the matrix wraps. Otherwise as the type
    /// describes.
    pub fn pick(&self, dist_mat: &[f64], pool_size: u32, pick_size: u32) -> Result<Vec<u32>> {
        self.pick_with_params(dist_mat, pool_size, pick_size, &MaxMinPickParams::default())
    }

    /// As [`pick`](Self::pick), with first picks and a seed.
    pub fn pick_with_params(
        &self,
        dist_mat: &[f64],
        pool_size: u32,
        pick_size: u32,
        params: &MaxMinPickParams,
    ) -> Result<Vec<u32>> {
        Ok(ffi::max_min_picker_pick(
            dist_mat,
            pool_size,
            pick_size,
            &params.first_picks,
            params.seed,
        )?)
    }

    /// Picks by a distance function of two items' indices, which RDKit
    /// calls as it needs distances and does not cache. Mirrors
    /// `MaxMinPicker.LazyPick`.
    pub fn lazy_pick<F>(&self, dist_func: F, pool_size: u32, pick_size: u32) -> Result<Vec<u32>>
    where
        F: FnMut(u32, u32) -> f64,
    {
        self.lazy_pick_with_params(
            dist_func,
            pool_size,
            pick_size,
            &MaxMinPickParams::default(),
        )
    }

    /// As [`lazy_pick`](Self::lazy_pick), with first picks and a seed.
    pub fn lazy_pick_with_params<F>(
        &self,
        dist_func: F,
        pool_size: u32,
        pick_size: u32,
        params: &MaxMinPickParams,
    ) -> Result<Vec<u32>>
    where
        F: FnMut(u32, u32) -> f64,
    {
        Ok(lazy(dist_func, pool_size, pick_size, params, -1.0)?.0)
    }

    /// Picks by the Tanimoto distances between the first `pool_size`
    /// vectors. Mirrors `MaxMinPicker.LazyBitVectorPick`.
    pub fn lazy_bit_vector_pick(
        &self,
        objects: &[ExplicitBitVect],
        pool_size: u32,
        pick_size: u32,
    ) -> Result<Vec<u32>> {
        self.lazy_bit_vector_pick_with_params(
            objects,
            pool_size,
            pick_size,
            &MaxMinPickParams::default(),
        )
    }

    /// As [`lazy_bit_vector_pick`](Self::lazy_bit_vector_pick), with first
    /// picks and a seed.
    pub fn lazy_bit_vector_pick_with_params(
        &self,
        objects: &[ExplicitBitVect],
        pool_size: u32,
        pick_size: u32,
        params: &MaxMinPickParams,
    ) -> Result<Vec<u32>> {
        Ok(lazy_bit_vectors(objects, pool_size, pick_size, params, -1.0)?.0)
    }

    /// As [`lazy_pick`](Self::lazy_pick), stopping before a pick whose
    /// distance to the picks is at most `threshold`, when `threshold` is
    /// not negative. Returns the picks and the last pick's distance (-1
    /// when the first picks fill `pick_size`). Mirrors
    /// `MaxMinPicker.LazyPickWithThreshold`.
    pub fn lazy_pick_with_threshold<F>(
        &self,
        dist_func: F,
        pool_size: u32,
        pick_size: u32,
        threshold: f64,
    ) -> Result<(Vec<u32>, f64)>
    where
        F: FnMut(u32, u32) -> f64,
    {
        self.lazy_pick_with_threshold_with_params(
            dist_func,
            pool_size,
            pick_size,
            threshold,
            &MaxMinPickParams::default(),
        )
    }

    /// As [`lazy_pick_with_threshold`](Self::lazy_pick_with_threshold),
    /// with first picks and a seed.
    pub fn lazy_pick_with_threshold_with_params<F>(
        &self,
        dist_func: F,
        pool_size: u32,
        pick_size: u32,
        threshold: f64,
        params: &MaxMinPickParams,
    ) -> Result<(Vec<u32>, f64)>
    where
        F: FnMut(u32, u32) -> f64,
    {
        lazy(dist_func, pool_size, pick_size, params, threshold)
    }

    /// As [`lazy_bit_vector_pick`](Self::lazy_bit_vector_pick), with a
    /// threshold as [`lazy_pick_with_threshold`](Self::lazy_pick_with_threshold).
    /// Mirrors `MaxMinPicker.LazyBitVectorPickWithThreshold`.
    pub fn lazy_bit_vector_pick_with_threshold(
        &self,
        objects: &[ExplicitBitVect],
        pool_size: u32,
        pick_size: u32,
        threshold: f64,
    ) -> Result<(Vec<u32>, f64)> {
        self.lazy_bit_vector_pick_with_threshold_with_params(
            objects,
            pool_size,
            pick_size,
            threshold,
            &MaxMinPickParams::default(),
        )
    }

    /// As
    /// [`lazy_bit_vector_pick_with_threshold`](Self::lazy_bit_vector_pick_with_threshold),
    /// with first picks and a seed.
    pub fn lazy_bit_vector_pick_with_threshold_with_params(
        &self,
        objects: &[ExplicitBitVect],
        pool_size: u32,
        pick_size: u32,
        threshold: f64,
        params: &MaxMinPickParams,
    ) -> Result<(Vec<u32>, f64)> {
        lazy_bit_vectors(objects, pool_size, pick_size, params, threshold)
    }
}

fn lazy<F>(
    mut dist_func: F,
    pool_size: u32,
    pick_size: u32,
    params: &MaxMinPickParams,
    threshold: f64,
) -> Result<(Vec<u32>, f64)>
where
    F: FnMut(u32, u32) -> f64,
{
    let mut func = DistanceFunction(&mut dist_func);
    // SAFETY: the wrapper calls the function only during this call.
    let result = unsafe {
        ffi::max_min_picker_lazy_pick(
            &mut func,
            pool_size,
            pick_size,
            &params.first_picks,
            params.seed,
            threshold,
        )
    }?;
    Ok((result.picks, result.threshold))
}

fn lazy_bit_vectors(
    objects: &[ExplicitBitVect],
    pool_size: u32,
    pick_size: u32,
    params: &MaxMinPickParams,
    threshold: f64,
) -> Result<(Vec<u32>, f64)> {
    let result = with_bit_vect_list(objects, |list| {
        ffi::max_min_picker_lazy_bit_vector_pick(
            list,
            pool_size,
            pick_size,
            &params.first_picks,
            params.seed,
            threshold,
        )
    })?;
    Ok((result.picks, result.threshold))
}

/// The optional arguments of `LeaderPicker`'s methods, with Python's
/// defaults (A2). Python's `numThreads` is not bound: the picks do not
/// depend on it, and RDKit's threads call `std::terminate` when the distance
/// function fails on them, so the picks are made on the calling thread.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct LeaderPickParams {
    /// The most items to pick; 0 picks until no item is left.
    pub pick_size: u32,
    /// Items to pick first, in order.
    pub first_picks: Vec<u32>,
}

/// Picks a diverse subset by Roger Sayle's Leader algorithm (a sphere
/// exclusion): each pick is the first item left, and removes from the pool
/// every item within `threshold` of it, so no two picks are within
/// `threshold` of each other. Mirrors `rdSimDivPickers.LeaderPicker`.
///
/// `Error::Value`, as RDKit raises, for an empty pool, a `pick_size` above
/// `pool_size`, and a first pick past the pool; and for a pool past
/// `i32::MAX` items, which RDKit numbers in `int`.
#[derive(Debug, Default)]
pub struct LeaderPicker {
    _private: (),
}

impl LeaderPicker {
    /// Picks by the Tanimoto distances between the first `pool_size`
    /// vectors. Mirrors `LeaderPicker.LazyBitVectorPick`.
    pub fn lazy_bit_vector_pick(
        &self,
        objects: &[ExplicitBitVect],
        pool_size: u32,
        threshold: f64,
    ) -> Result<Vec<u32>> {
        self.lazy_bit_vector_pick_with_params(
            objects,
            pool_size,
            threshold,
            &LeaderPickParams::default(),
        )
    }

    /// As [`lazy_bit_vector_pick`](Self::lazy_bit_vector_pick), with a pick
    /// size and first picks.
    pub fn lazy_bit_vector_pick_with_params(
        &self,
        objects: &[ExplicitBitVect],
        pool_size: u32,
        threshold: f64,
        params: &LeaderPickParams,
    ) -> Result<Vec<u32>> {
        Ok(with_bit_vect_list(objects, |list| {
            ffi::leader_picker_lazy_bit_vector_pick(
                list,
                pool_size,
                threshold,
                params.pick_size,
                &params.first_picks,
            )
        })?)
    }

    /// Picks by a distance function of two items' indices. Mirrors
    /// `LeaderPicker.LazyPick`.
    pub fn lazy_pick<F>(&self, dist_func: F, pool_size: u32, threshold: f64) -> Result<Vec<u32>>
    where
        F: FnMut(u32, u32) -> f64,
    {
        self.lazy_pick_with_params(
            dist_func,
            pool_size,
            threshold,
            &LeaderPickParams::default(),
        )
    }

    /// As [`lazy_pick`](Self::lazy_pick), with a pick size and first picks.
    pub fn lazy_pick_with_params<F>(
        &self,
        mut dist_func: F,
        pool_size: u32,
        threshold: f64,
        params: &LeaderPickParams,
    ) -> Result<Vec<u32>>
    where
        F: FnMut(u32, u32) -> f64,
    {
        let mut func = DistanceFunction(&mut dist_func);
        // SAFETY: the wrapper calls the function only during this call.
        Ok(unsafe {
            ffi::leader_picker_lazy_pick(
                &mut func,
                pool_size,
                threshold,
                params.pick_size,
                &params.first_picks,
            )
        }?)
    }
}

rdkit_enum! {
    /// A hierarchical clustering criterion. Mirrors
    /// `rdSimDivPickers.ClusterMethod`.
    ClusterMethod: ffi::ClusterMethod {
        Ward = WARD,
        Slink = SLINK,
        Clink = CLINK,
        Upgma = UPGMA,
        Mcquitty = MCQUITTY,
        Gower = GOWER,
        Centroid = CENTROID,
    }
}

/// Picks a diverse subset by clustering the pool hierarchically, with
/// Murtagh's code (as `rdkit_cxx::ml::cluster::clustering`, whose lock
/// it shares), into `pick_size` clusters, and picking from each the item
/// with the least sum of squared distances to the rest. Mirrors
/// `rdSimDivPickers.HierarchicalClusterPicker`.
///
/// Both methods take a distance matrix and return `Error::Value` for fewer
/// than 2 items, a `pick_size` of 0 (for which RDKit reads cluster -1, a
/// segfault in Python), a matrix shorter than the pool's (read past its end
/// in Python), and distances that are not finite. Distances large enough to
/// reach Murtagh's 1e20 sentinel, past which its merges index out of
/// bounds, are scaled down by an exact power of two first, which picks the
/// same items.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HierarchicalClusterPicker {
    cluster_method: ClusterMethod,
}

impl HierarchicalClusterPicker {
    /// A picker clustering by a criterion. Mirrors
    /// `HierarchicalClusterPicker(clusterMethod)`.
    pub fn new(cluster_method: ClusterMethod) -> HierarchicalClusterPicker {
        HierarchicalClusterPicker { cluster_method }
    }

    /// The picks, one per cluster. Mirrors `HierarchicalClusterPicker.Pick`.
    ///
    /// `Error::Value` for a `pick_size` that is not below `pool_size`, as
    /// Python raises, and for more than 65537 items, past which RDKit's
    /// `unsigned int` index into the matrix wraps; otherwise as the type
    /// describes.
    pub fn pick(&self, dist_mat: &[f64], pool_size: u32, pick_size: u32) -> Result<Vec<u32>> {
        Ok(ffi::hierarchical_cluster_picker_pick(
            self.cluster_method.to_ffi(),
            dist_mat,
            pool_size,
            pick_size,
        )?)
    }

    /// The `pick_size` clusters, each a list of items. Mirrors
    /// `HierarchicalClusterPicker.Cluster`.
    ///
    /// `Error::Runtime` for a `pick_size` above `pool_size`, as Python
    /// raises; otherwise as the type describes.
    pub fn cluster(
        &self,
        dist_mat: &[f64],
        pool_size: u32,
        pick_size: u32,
    ) -> Result<Vec<Vec<u32>>> {
        Ok(ffi::hierarchical_cluster_picker_cluster(
            self.cluster_method.to_ffi(),
            dist_mat,
            pool_size,
            pick_size,
        )?
        .into_iter()
        .map(|list| list.indices)
        .collect())
    }
}
