// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// pubchem_shape's `ShapeInputOptions`, with its pairs as parallel
    /// vectors (A7).
    struct ShapeInputOptionsData {
        use_colors: bool,
        include_dummies: bool,
        dummy_radius: f64,
        atom_subset: Vec<u32>,
        not_color_atoms: Vec<u32>,
        atom_radii_atoms: Vec<u32>,
        atom_radii_values: Vec<f64>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/shape_align.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Conformer = crate::ffi::conformer::ffi::Conformer;
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
        fn pubchem_shape_input_sov(shape: &PubChemShapeInput) -> Result<f64>;
        fn pubchem_shape_input_sof(shape: &PubChemShapeInput) -> Result<f64>;

        #[allow(clippy::too_many_arguments)]
        fn pubchem_align_mol(
            reference: &RWMol,
            probe: Pin<&mut RWMol>,
            ref_conf_id: i32,
            probe_conf_id: i32,
            use_colors: bool,
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
            matrix: &[f64],
            probe_shape: &PubChemShapeInput,
            probe_conformer: Pin<&mut Conformer>,
        ) -> Result<()>;
    }
}

// SAFETY: a ShapeInput is plain data (vectors and numbers), which the
// wrappers only read through const references or change through an
// exclusive one.
unsafe impl Send for ffi::PubChemShapeInput {}
unsafe impl Sync for ffi::PubChemShapeInput {}
