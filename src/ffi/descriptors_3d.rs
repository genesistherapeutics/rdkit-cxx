// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A list of `f64`s; nested lists cross as `Vec<DoubleList>`.
    struct DoubleList {
        values: Vec<f64>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/descriptors_3d.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type Point3D = crate::ffi::geometry::ffi::Point3D;

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
        fn calc_coulomb_mat(mol: &RWMol, conf_id: i32) -> Result<Vec<f64>>;
        fn calc_eem_charges(mol: &RWMol, conf_id: i32) -> Result<Vec<f64>>;

        fn get_usr(mol: &RWMol, conf_id: i32) -> Result<Vec<f64>>;
        fn get_usrcat(
            mol: &RWMol,
            has_atom_selections: bool,
            atom_selections: &[IndexList],
            conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn get_usr_distributions(
            coords: &[Point3D],
            points: &mut Vec<Point3D>,
        ) -> Result<Vec<DoubleList>>;
        fn get_usr_distributions_from_points(
            coords: &[Point3D],
            points: &[Point3D],
        ) -> Result<Vec<DoubleList>>;
        fn get_usr_from_distributions(distances: &[DoubleList]) -> Result<Vec<f64>>;
        fn get_usr_score(descriptor1: &[f64], descriptor2: &[f64], weights: &[f64]) -> Result<f64>;
    }
}
