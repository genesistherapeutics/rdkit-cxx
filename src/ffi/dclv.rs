// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/dclv.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type DoubleCubicLatticeVolume;

        fn dclv_new(
            mol: &RWMol,
            radii: &[f64],
            is_protein: bool,
            include_ligand: bool,
            probe_radius: f64,
            conf_id: i32,
        ) -> Result<UniquePtr<DoubleCubicLatticeVolume>>;
        fn dclv_get_surface_area(dclv: Pin<&mut DoubleCubicLatticeVolume>) -> Result<f64>;
        fn dclv_get_atom_surface_area(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            atom_idx: u32,
        ) -> Result<f64>;
        fn dclv_get_polar_surface_area(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            include_s_and_p: bool,
            include_hs: bool,
        ) -> Result<f64>;
        fn dclv_get_partial_surface_area(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            atom_indices: &[u32],
        ) -> Result<f64>;
        fn dclv_get_surface_points(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            all_points: bool,
            atoms: &mut Vec<u32>,
            counts: &mut Vec<u32>,
            points: &mut Vec<Point3D>,
        ) -> Result<()>;
        fn dclv_get_volume(dclv: Pin<&mut DoubleCubicLatticeVolume>) -> Result<f64>;
        fn dclv_get_vdw_volume(dclv: Pin<&mut DoubleCubicLatticeVolume>) -> Result<f64>;
        fn dclv_get_atom_volume(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            atom_idx: u32,
            solvent_radius: f64,
        ) -> Result<f64>;
        fn dclv_get_polar_volume(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            include_s_and_p: bool,
            include_hs: bool,
        ) -> Result<f64>;
        fn dclv_get_partial_volume(
            dclv: Pin<&mut DoubleCubicLatticeVolume>,
            atom_indices: &[u32],
        ) -> Result<f64>;
        fn dclv_get_compactness(dclv: Pin<&mut DoubleCubicLatticeVolume>) -> Result<f64>;
        fn dclv_get_packing_density(dclv: Pin<&mut DoubleCubicLatticeVolume>) -> Result<f64>;
    }
}
