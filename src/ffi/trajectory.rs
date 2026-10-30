// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/trajectory.h");

        fn read_amber_trajectory(
            path: &[u8],
            num_points: u32,
            coords: &mut Vec<f64>,
        ) -> Result<u32>;

        fn read_gromos_trajectory(
            path: &[u8],
            num_points: u32,
            coords: &mut Vec<f64>,
        ) -> Result<u32>;
    }
}
