// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_transforms.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Conformer = crate::ffi::conformer::ffi::Conformer;
        type Point3D = crate::ffi::geometry::ffi::Point3D;

        fn compute_centroid(
            conf: &Conformer,
            ignore_hs: bool,
            has_weights: bool,
            weights: &[f64],
        ) -> Result<Point3D>;
        fn compute_canonical_transform(
            conf: &Conformer,
            has_center: bool,
            center: &Point3D,
            normalize_covar: bool,
            ignore_hs: bool,
        ) -> Result<Vec<f64>>;
        #[allow(clippy::too_many_arguments)]
        fn compute_principal_axes_and_moments(
            conf: &Conformer,
            from_gyration_matrix: bool,
            ignore_hs: bool,
            has_weights: bool,
            weights: &[f64],
            axes: &mut Vec<f64>,
            moments: &mut Vec<f64>,
        ) -> Result<bool>;
        fn transform_conformer(conf: Pin<&mut Conformer>, trans: &[f64]) -> Result<()>;
        fn canonicalize_conformer(
            conf: Pin<&mut Conformer>,
            has_center: bool,
            center: &Point3D,
            normalize_covar: bool,
            ignore_hs: bool,
        ) -> Result<()>;
        fn canonicalize_mol(
            mol: Pin<&mut RWMol>,
            normalize_covar: bool,
            ignore_hs: bool,
        ) -> Result<()>;
        fn get_bond_length(conf: &Conformer, i_atom_id: u32, j_atom_id: u32) -> Result<f64>;
        fn set_bond_length(
            conf: Pin<&mut Conformer>,
            i_atom_id: u32,
            j_atom_id: u32,
            value: f64,
        ) -> Result<()>;
        fn get_angle(
            conf: &Conformer,
            i_atom_id: u32,
            j_atom_id: u32,
            k_atom_id: u32,
            degrees: bool,
        ) -> Result<f64>;
        fn set_angle(
            conf: Pin<&mut Conformer>,
            i_atom_id: u32,
            j_atom_id: u32,
            k_atom_id: u32,
            value: f64,
            degrees: bool,
        ) -> Result<()>;
        fn get_dihedral(
            conf: &Conformer,
            i_atom_id: u32,
            j_atom_id: u32,
            k_atom_id: u32,
            l_atom_id: u32,
            degrees: bool,
        ) -> Result<f64>;
        #[allow(clippy::too_many_arguments)]
        fn set_dihedral(
            conf: Pin<&mut Conformer>,
            i_atom_id: u32,
            j_atom_id: u32,
            k_atom_id: u32,
            l_atom_id: u32,
            value: f64,
            degrees: bool,
        ) -> Result<()>;
    }
}
