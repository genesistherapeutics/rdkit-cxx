// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/conformer_parser.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn add_conformers_from_amber_trajectory(
            mol: Pin<&mut RWMol>,
            path: &[u8],
            num_confs: i32,
            clear_confs: bool,
        ) -> Result<Vec<u32>>;
    }
}
