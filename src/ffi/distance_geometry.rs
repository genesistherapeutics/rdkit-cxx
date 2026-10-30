// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/distance_geometry.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;

        fn do_triangle_smoothing(bounds: &mut Vec<f64>, n: u32, tol: f64) -> Result<bool>;
        #[allow(clippy::too_many_arguments)]
        fn embed_bounds_matrix(
            bounds: &[f64],
            n: u32,
            max_iters: i32,
            randomize_on_failure: bool,
            num_zero_fail: i32,
            weight_pairs: &[IndexPair],
            weights: &[f64],
            random_seed: i32,
        ) -> Result<Vec<Point3D>>;
    }
}
