// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/shape_align.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type PubChemShapeInput;

        fn pubchem_prepare_conformer(
            mol: &RWMol,
            conf_id: i32,
            use_colors: bool,
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
            ref_mol: &RWMol,
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
        ) -> Result<Vec<f64>>;
    }
}

// SAFETY: a ShapeInput is plain data (vectors, a map and numbers), which the
// wrappers only read through const references or change through an
// exclusive one.
unsafe impl Send for ffi::PubChemShapeInput {}
unsafe impl Sync for ffi::PubChemShapeInput {}
