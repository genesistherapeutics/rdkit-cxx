// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `FreeSASA::SASAOpts` (A7). Its enums have no fixed underlying type,
    /// which cxx cannot check, so they cross as their values.
    struct SASAOptsData {
        algorithm: i32,
        classifier: i32,
        probe_radius: f64,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/free_sasa.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;

        fn sasa_opts_default() -> Result<SASAOptsData>;
        fn classify_atoms(mol: Pin<&mut RWMol>, opts: &SASAOptsData) -> Result<Vec<f64>>;
        /// # Safety
        ///
        /// The query is null or valid.
        unsafe fn calc_sasa(
            mol: &RWMol,
            radii: &[f64],
            conf_idx: i32,
            query: *const Atom,
            opts: &SASAOptsData,
        ) -> Result<f64>;
        fn make_free_sasa_a_polar_atom_query() -> Result<UniquePtr<Atom>>;
        fn make_free_sasa_polar_atom_query() -> Result<UniquePtr<Atom>>;
    }
}
