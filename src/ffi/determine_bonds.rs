// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/determine_bonds.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn determine_connectivity(
            mol: Pin<&mut RWMol>,
            use_hueckel: bool,
            charge: i32,
            cov_factor: f64,
            use_vdw: bool,
        ) -> Result<()>;
        fn determine_bond_orders(
            mol: Pin<&mut RWMol>,
            charge: i32,
            allow_charged_fragments: bool,
            embed_chiral: bool,
            use_atom_map: bool,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn determine_bonds(
            mol: Pin<&mut RWMol>,
            use_hueckel: bool,
            charge: i32,
            cov_factor: f64,
            allow_charged_fragments: bool,
            embed_chiral: bool,
            use_atom_map: bool,
            use_vdw: bool,
        ) -> Result<()>;
        fn hueckel_enabled() -> Result<bool>;
    }
}
