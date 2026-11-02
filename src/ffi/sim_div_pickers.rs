// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::panic::{AssertUnwindSafe, catch_unwind};

use crate::ffi::standardize_validate::panic_message;

/// A distance function, lent to C++ for one call.
pub(crate) struct DistanceFunction<'a>(pub(crate) &'a mut dyn FnMut(u32, u32) -> f64);

/// Called by the pickers, on the calling thread, with two items' indices.
/// A panic is caught here, at the boundary, and becomes the exception C++
/// sees, which ends the pick.
///
/// # Safety
///
/// cxx makes a function whose signature names a lifetime `unsafe` for C++:
/// the C++ caller must not keep `func` past the call that lent it, which
/// the wrappers do not.
unsafe fn call_distance_function(
    func: &mut DistanceFunction<'_>,
    i: u32,
    j: u32,
) -> Result<f64, String> {
    catch_unwind(AssertUnwindSafe(|| (func.0)(i, j)))
        .map_err(|payload| panic_message(payload, "a distance function"))
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDPickers::HierarchicalClusterPicker::ClusterMethod`; cxx checks the
    /// values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum ClusterMethod {
        WARD = 1,
        SLINK,
        CLINK,
        UPGMA,
        MCQUITTY,
        GOWER,
        CENTROID,
    }

    /// MaxMin's picks and the distance of its last pick.
    struct PicksWithThreshold {
        picks: Vec<u32>,
        threshold: f64,
    }

    extern "Rust" {
        type DistanceFunction<'a>;

        unsafe fn call_distance_function<'a>(
            func: &mut DistanceFunction<'a>,
            i: u32,
            j: u32,
        ) -> Result<f64>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/sim_div_pickers.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type ClusterMethod;
        type BitVectList = crate::ffi::metric_matrix_calc::ffi::BitVectList;
        type IndexList = crate::ffi::geometry::ffi::IndexList;

        fn max_min_picker_pick(
            dist_mat: &[f64],
            pool_size: u32,
            pick_size: u32,
            first_picks: &[u32],
            seed: i32,
        ) -> Result<Vec<u32>>;
        unsafe fn max_min_picker_lazy_pick<'a>(
            func: &mut DistanceFunction<'a>,
            pool_size: u32,
            pick_size: u32,
            first_picks: &[u32],
            seed: i32,
            threshold: f64,
        ) -> Result<PicksWithThreshold>;
        fn max_min_picker_lazy_bit_vector_pick(
            objects: &BitVectList,
            pool_size: u32,
            pick_size: u32,
            first_picks: &[u32],
            seed: i32,
            threshold: f64,
        ) -> Result<PicksWithThreshold>;

        unsafe fn leader_picker_lazy_pick<'a>(
            func: &mut DistanceFunction<'a>,
            pool_size: u32,
            threshold: f64,
            pick_size: u32,
            first_picks: &[u32],
        ) -> Result<Vec<u32>>;
        fn leader_picker_lazy_bit_vector_pick(
            objects: &BitVectList,
            pool_size: u32,
            threshold: f64,
            pick_size: u32,
            first_picks: &[u32],
        ) -> Result<Vec<u32>>;

        fn hierarchical_cluster_picker_pick(
            method: ClusterMethod,
            dist_mat: &[f64],
            pool_size: u32,
            pick_size: u32,
        ) -> Result<Vec<u32>>;
        fn hierarchical_cluster_picker_cluster(
            method: ClusterMethod,
            dist_mat: &[f64],
            pool_size: u32,
            pick_size: u32,
        ) -> Result<Vec<IndexList>>;
    }
}
