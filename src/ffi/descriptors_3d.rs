// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/descriptors_3d.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn calc_pmi1(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_pmi2(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_pmi3(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_npr1(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_npr2(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_radius_of_gyration(
            mol: &RWMol,
            conf_id: i32,
            use_atomic_masses: bool,
        ) -> Result<f64>;
        fn calc_inertial_shape_factor(
            mol: &RWMol,
            conf_id: i32,
            use_atomic_masses: bool,
        ) -> Result<f64>;
        fn calc_eccentricity(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_asphericity(mol: &RWMol, conf_id: i32, use_atomic_masses: bool) -> Result<f64>;
        fn calc_spherocity_index(mol: &RWMol, conf_id: i32) -> Result<f64>;
        fn calc_pbf(mol: &RWMol, conf_id: i32) -> Result<f64>;

        fn calc_autocorr_3d(
            mol: &RWMol,
            conf_id: i32,
            custom_atom_property: &str,
        ) -> Result<Vec<f64>>;
        fn calc_rdf(mol: &RWMol, conf_id: i32, custom_atom_property: &str) -> Result<Vec<f64>>;
        fn calc_morse(mol: &RWMol, conf_id: i32, custom_atom_property: &str) -> Result<Vec<f64>>;
        fn calc_whim(
            mol: &RWMol,
            conf_id: i32,
            thresh: f64,
            custom_atom_property: &str,
        ) -> Result<Vec<f64>>;
        fn calc_getaway(
            mol: &RWMol,
            conf_id: i32,
            precision: f64,
            custom_atom_property: &str,
        ) -> Result<Vec<f64>>;
        fn calc_coulomb_mat(mol: &RWMol, conf_id: i32) -> Result<Vec<f64>>;
        fn calc_eem_charges(mol: &RWMol, conf_id: i32) -> Result<Vec<f64>>;
    }
}
