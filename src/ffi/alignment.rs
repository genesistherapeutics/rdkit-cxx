// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/alignment.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type Point3D = crate::ffi::geometry::ffi::Point3D;

        #[allow(clippy::too_many_arguments)]
        fn align_point_pairs(
            ref_points: &[Point3D],
            probe_points: &[Point3D],
            has_weights: bool,
            weights: &[f64],
            reflect: bool,
            max_iterations: u32,
            transform: &mut Vec<f64>,
        ) -> Result<f64>;
    }
}
