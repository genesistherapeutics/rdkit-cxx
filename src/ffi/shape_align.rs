// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// pubchem_shape's `ShapeInputOptions`, with its pairs and triples as
    /// parallel vectors (A7).
    struct ShapeInputOptionsData {
        use_colors: bool,
        include_dummies: bool,
        dummy_radius: f64,
        atom_subset: Vec<u32>,
        not_color_atoms: Vec<u32>,
        atom_radii_atoms: Vec<u32>,
        atom_radii_values: Vec<f64>,
        custom_feature_types: Vec<u32>,
        custom_feature_positions: Vec<Point3D>,
        custom_feature_radii: Vec<f64>,
        normalize: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/shape_align.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Conformer = crate::ffi::conformer::ffi::Conformer;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type PubChemShapeInput;

        fn pubchem_shape_input_options_default() -> Result<ShapeInputOptionsData>;
        fn pubchem_prepare_conformer(
            mol: &RWMol,
            conf_id: i32,
            opts: &ShapeInputOptionsData,
        ) -> Result<UniquePtr<PubChemShapeInput>>;
        fn pubchem_shape_input_clone(
            shape: &PubChemShapeInput,
        ) -> Result<UniquePtr<PubChemShapeInput>>;
        fn pubchem_shape_input_coord(shape: &PubChemShapeInput) -> Result<Vec<f32>>;
        fn pubchem_shape_input_alpha_vector(shape: &PubChemShapeInput) -> Result<Vec<f64>>;
        fn pubchem_shape_input_atom_type_vector(shape: &PubChemShapeInput) -> Result<Vec<u32>>;
        fn pubchem_shape_input_volume_atom_index_vector(
            shape: &PubChemShapeInput,
        ) -> Result<Vec<u32>>;
        fn pubchem_shape_input_shift(shape: &PubChemShapeInput) -> Result<Vec<f64>>;
        fn pubchem_shape_input_inertial_rot(shape: &PubChemShapeInput) -> Result<Vec<f64>>;
        fn pubchem_shape_input_sov(shape: &PubChemShapeInput) -> Result<f64>;
        fn pubchem_shape_input_sof(shape: &PubChemShapeInput) -> Result<f64>;

        #[allow(clippy::too_many_arguments)]
        fn pubchem_align_mol(
            reference: &RWMol,
            probe: Pin<&mut RWMol>,
            ref_opts: &ShapeInputOptionsData,
            probe_opts: &ShapeInputOptionsData,
            ref_conf_id: i32,
            probe_conf_id: i32,
            opt_param: f64,
            max_preiters: u32,
            max_postiters: u32,
        ) -> Result<Vec<f64>>;
        #[allow(clippy::too_many_arguments)]
        fn pubchem_align_mol_shape_input(
            ref_shape: &PubChemShapeInput,
            probe: Pin<&mut RWMol>,
            probe_conf_id: i32,
            use_colors: bool,
            opt_param: f64,
            max_preiters: u32,
            max_postiters: u32,
            apply_ref_shift: bool,
        ) -> Result<Vec<f64>>;
        fn pubchem_align_shapes(
            ref_shape: &PubChemShapeInput,
            probe_shape: Pin<&mut PubChemShapeInput>,
            opt_param: f64,
            max_preiters: u32,
            max_postiters: u32,
        ) -> Result<Vec<f64>>;
        fn pubchem_transform_conformer_to_shape(
            final_trans: &[f64],
            final_rot: &[f64],
            matrix: &[f64],
            probe_shape: &PubChemShapeInput,
            probe_conformer: Pin<&mut Conformer>,
        ) -> Result<()>;
        fn pubchem_score_mol(
            mol1: &RWMol,
            mol2: &RWMol,
            mol1_opts: &ShapeInputOptionsData,
            mol2_opts: &ShapeInputOptionsData,
            mol1_conf_id: i32,
            mol2_conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn pubchem_score_mol_shape_input(
            shape: &PubChemShapeInput,
            mol: &RWMol,
            mol_opts: &ShapeInputOptionsData,
            mol_conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn pubchem_score_shape(
            shape1: &PubChemShapeInput,
            shape2: &PubChemShapeInput,
            use_colors: bool,
        ) -> Result<Vec<f64>>;
    }
}

// SAFETY: a ShapeInput is plain data (vectors and numbers), which the
// wrappers only read through const references or change through an
// exclusive one.
unsafe impl Send for ffi::PubChemShapeInput {}
unsafe impl Sync for ffi::PubChemShapeInput {}
