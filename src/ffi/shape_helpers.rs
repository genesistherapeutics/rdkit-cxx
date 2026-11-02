// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/shape_helpers.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Conformer = crate::ffi::conformer::ffi::Conformer;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type UniformGrid3D = crate::ffi::grid::ffi::UniformGrid3D;
        type DiscreteValueType = crate::ffi::data_structs::ffi::DiscreteValueType;

        #[allow(clippy::too_many_arguments)]
        fn encode_shape(
            mol: &RWMol,
            grid: Pin<&mut UniformGrid3D>,
            conf_id: i32,
            has_trans: bool,
            trans: &[f64],
            vdw_scale: f64,
            step_size: f64,
            max_layers: i32,
            ignore_hs: bool,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn shape_compare(
            mol1: &RWMol,
            mol2: &RWMol,
            kind: u8,
            alpha: f64,
            beta: f64,
            conf_id1: i32,
            conf_id2: i32,
            grid_spacing: f64,
            bits_per_point: DiscreteValueType,
            vdw_scale: f64,
            step_size: f64,
            max_layers: i32,
            ignore_hs: bool,
            allow_reordering: bool,
        ) -> Result<f64>;
        fn compute_conf_box(
            conf: &Conformer,
            dims: bool,
            has_trans: bool,
            trans: &[f64],
            padding: f64,
        ) -> Result<Vec<Point3D>>;
        fn compute_union_box(
            lower1: &Point3D,
            upper1: &Point3D,
            lower2: &Point3D,
            upper2: &Point3D,
        ) -> Result<Vec<Point3D>>;
    }
}
