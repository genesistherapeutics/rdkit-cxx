// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// What `RDKit::Descriptors::DoubleCubicLatticeVolume` computes.
    #[derive(Clone, Copy, Debug, PartialEq)]
    struct DclvResults {
        surface_area: f64,
        volume: f64,
        vdw_volume: f64,
        compactness: f64,
        packing_density: f64,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/dclv.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn dclv_compute(
            mol: &RWMol,
            is_protein: bool,
            include_ligand: bool,
            probe_radius: f64,
            depth: i32,
            dot_density: i32,
        ) -> Result<DclvResults>;
    }
}
