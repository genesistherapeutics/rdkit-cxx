// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/reduced_graphs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn generate_mol_extended_reduced_graph(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn generate_erg_fingerprint_for_reduced_graph(
            mol: &RWMol,
            fuzz_increment: f64,
            min_path: i32,
            max_path: i32,
        ) -> Result<Vec<f64>>;
        fn get_erg_fingerprint(
            mol: &RWMol,
            fuzz_increment: f64,
            min_path: i32,
            max_path: i32,
        ) -> Result<Vec<f64>>;
    }
}
