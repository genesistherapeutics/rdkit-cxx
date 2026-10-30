// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::GaussianShape::StartMode`; cxx checks the values.
    #[repr(i32)]
    #[allow(non_camel_case_types, clippy::upper_case_acronyms)]
    enum StartMode {
        ROTATE_0,
        ROTATE_180,
        ROTATE_180_WIGGLE,
        ROTATE_45,
        ROTATE_0_FRAGMENT,
        ROTATE_180_FRAGMENT,
        ROTATE_45_FRAGMENT,
        A_LA_PUBCHEM,
    }

    /// `RDKit::GaussianShape::OptimMode`; cxx checks the values.
    #[repr(i32)]
    #[allow(non_camel_case_types, clippy::upper_case_acronyms)]
    enum OptimMode {
        SHAPE_ONLY,
        SHAPE_PLUS_COLOR_SCORE,
        SHAPE_PLUS_COLOR,
    }

    /// `RDKit::GaussianShape::ShapeInputOptions` (A7). The custom features,
    /// one list per conformer, are flattened: `custom_feature_counts` gives
    /// each conformer's number of features, and `custom_feature_atom_counts`
    /// each feature's number of atoms.
    struct GaussianShapeInputOptionsData {
        use_colors: bool,
        custom_feature_counts: Vec<u32>,
        custom_feature_types: Vec<u32>,
        custom_feature_positions: Vec<Point3D>,
        custom_feature_radii: Vec<f64>,
        custom_feature_atom_counts: Vec<u32>,
        custom_feature_atoms: Vec<u32>,
        atom_subset: Vec<u32>,
        atom_radii_atoms: Vec<u32>,
        atom_radii_values: Vec<f64>,
        all_carbon_radii: bool,
        shape_prune_threshold: f64,
        sort_shapes: bool,
        include_dummies: bool,
    }

    /// `RDKit::GaussianShape::ShapeOverlayOptions` (A7).
    struct ShapeOverlayOptionsData {
        start_mode: StartMode,
        optim_mode: OptimMode,
        sim_alpha: f64,
        sim_beta: f64,
        opt_param: f64,
        n_steps: i32,
        normalize: bool,
        use_dist_cutoff: bool,
        dist_cutoff: f64,
        shape_convergence_criterion: f64,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/gaussian_shape.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type StartMode;
        type OptimMode;
        type GaussianShapeInput;

        fn gaussian_shape_input_options_default() -> Result<GaussianShapeInputOptionsData>;
        fn shape_overlay_options_default() -> Result<ShapeOverlayOptionsData>;

        fn gaussian_shape_input_new(
            mol: &RWMol,
            conf_id: i32,
            shape_opt: &GaussianShapeInputOptionsData,
            overlay_opts: &ShapeOverlayOptionsData,
        ) -> Result<UniquePtr<GaussianShapeInput>>;
        fn gaussian_shape_input_clone(
            shape: &GaussianShapeInput,
        ) -> Result<UniquePtr<GaussianShapeInput>>;
        fn gaussian_shape_input_get_smiles(
            shape: &GaussianShapeInput,
        ) -> Result<UniquePtr<CxxString>>;
        fn gaussian_shape_input_set_active_shape(
            shape: Pin<&mut GaussianShapeInput>,
            shape_num: u32,
        ) -> Result<()>;
        fn gaussian_shape_input_get_active_shape(shape: &GaussianShapeInput) -> Result<u32>;
        fn gaussian_shape_input_num_atoms(shape: &GaussianShapeInput) -> Result<u32>;
        fn gaussian_shape_input_num_features(shape: &GaussianShapeInput) -> Result<u32>;
        fn gaussian_shape_input_num_shapes(shape: &GaussianShapeInput) -> Result<u32>;
        fn gaussian_shape_input_shape_volume(shape: &GaussianShapeInput) -> Result<f64>;
        fn gaussian_shape_input_color_volume(shape: &GaussianShapeInput) -> Result<f64>;
        fn gaussian_shape_input_normalize_coords(shape: Pin<&mut GaussianShapeInput>)
        -> Result<()>;
        fn gaussian_shape_input_shape_to_mol(
            shape: &GaussianShapeInput,
            include_colors: bool,
            with_bonds: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn gaussian_shape_input_best_similarity(
            shape: Pin<&mut GaussianShapeInput>,
            fit_shape: &GaussianShapeInput,
            threshold: f64,
            overlay_opts: &ShapeOverlayOptionsData,
        ) -> Result<Vec<f64>>;
        fn gaussian_shape_input_max_possible_similarity(
            shape: &GaussianShapeInput,
            fit_shape: &GaussianShapeInput,
            overlay_opts: &ShapeOverlayOptionsData,
        ) -> Result<f64>;

        #[allow(clippy::too_many_arguments)]
        fn gaussian_align_mol(
            reference: &RWMol,
            fit: Pin<&mut RWMol>,
            ref_opts: &GaussianShapeInputOptionsData,
            fit_opts: &GaussianShapeInputOptionsData,
            overlay_opts: &ShapeOverlayOptionsData,
            ref_conf_id: i32,
            fit_conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn gaussian_align_mol_shape_input(
            ref_shape: &GaussianShapeInput,
            fit: Pin<&mut RWMol>,
            fit_opts: &GaussianShapeInputOptionsData,
            overlay_opts: &ShapeOverlayOptionsData,
            fit_conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn gaussian_align_shapes(
            ref_shape: &GaussianShapeInput,
            fit_shape: Pin<&mut GaussianShapeInput>,
            overlay_opts: &ShapeOverlayOptionsData,
        ) -> Result<Vec<f64>>;
        #[allow(clippy::too_many_arguments)]
        fn gaussian_score_mol(
            reference: &RWMol,
            fit: &RWMol,
            ref_opts: &GaussianShapeInputOptionsData,
            fit_opts: &GaussianShapeInputOptionsData,
            overlay_opts: &ShapeOverlayOptionsData,
            ref_conf_id: i32,
            fit_conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn gaussian_score_mol_shape_input(
            ref_shape: &GaussianShapeInput,
            fit: &RWMol,
            fit_opts: &GaussianShapeInputOptionsData,
            overlay_opts: &ShapeOverlayOptionsData,
            fit_conf_id: i32,
        ) -> Result<Vec<f64>>;
        fn gaussian_score_shape(
            ref_shape: &GaussianShapeInput,
            fit_shape: &GaussianShapeInput,
            overlay_opts: &ShapeOverlayOptionsData,
        ) -> Result<Vec<f64>>;
        fn gaussian_score_molecule_all_conformers(
            reference: &RWMol,
            fit: &RWMol,
            ref_opts: &GaussianShapeInputOptionsData,
            fit_opts: &GaussianShapeInputOptionsData,
            overlay_opts: &ShapeOverlayOptionsData,
        ) -> Result<Vec<f64>>;
    }
}

// SAFETY: a ShapeInput owns its data (vectors, a bitset, and strings); its
// const methods only read it, and the methods that cache or move take it
// exclusively.
unsafe impl Send for ffi::GaussianShapeInput {}
unsafe impl Sync for ffi::GaussianShapeInput {}
