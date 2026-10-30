// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::DGeomHelpers::EmbedFailureCauses`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum EmbedFailureCauses {
        INITIAL_COORDS = 0,
        FIRST_MINIMIZATION = 1,
        CHECK_TETRAHEDRAL_CENTERS = 2,
        CHECK_CHIRAL_CENTERS = 3,
        MINIMIZE_FOURTH_DIMENSION = 4,
        ETK_MINIMIZATION = 5,
        FINAL_CHIRAL_BOUNDS = 6,
        FINAL_CENTER_IN_VOLUME = 7,
        LINEAR_DOUBLE_BOND = 8,
        BAD_DOUBLE_BOND_STEREO = 9,
        CHECK_CHIRAL_CENTERS2 = 10,
        EXCEEDED_TIMEOUT = 11,
        END_OF_ENUM = 12,
    }

    /// A default `EmbedParameters`, or one of RDKit's presets.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum EmbedPreset {
        Default,
        ETKDG,
        ETKDGv2,
        ETKDGv3,
        SrETKDGv3,
        ETDG,
        KDG,
    }

    /// `RDKit::DGeomHelpers::EmbedParameters`, with its coordinate map and
    /// CPCI as parallel vectors and its bounds matrix flattened.
    struct EmbedParametersData {
        max_iterations: u32,
        num_threads: i32,
        timeout: u32,
        random_seed: i32,
        clear_confs: bool,
        use_random_coords: bool,
        box_size_mult: f64,
        rand_neg_eig: bool,
        num_zero_fail: u32,
        optimizer_force_tol: f64,
        basin_thresh: f64,
        ignore_smoothing_failures: bool,
        enforce_chirality: bool,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        et_version: u32,
        verbose: bool,
        prune_rms_thresh: f64,
        only_heavy_atoms_for_rms: bool,
        embed_fragments_separately: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle14config: bool,
        bounds_mat_force_scaling: f64,
        use_symmetry_for_pruning: bool,
        force_trans_amides: bool,
        track_failures: bool,
        enable_sequential_random_seeds: bool,
        symmetrize_conjugated_terminal_groups_for_pruning: bool,
        coord_map_atoms: Vec<u32>,
        coord_map_points: Vec<Point3D>,
        cpci_pairs: Vec<IndexPair>,
        cpci_values: Vec<f64>,
        has_bounds_mat: bool,
        bounds_mat_size: u32,
        bounds_mat: Vec<f64>,
        failures: Vec<u32>,
    }

    /// One of `GetExperimentalTorsions`' dictionaries.
    struct ExperimentalTorsionData {
        bond_index: u32,
        torsion_index: u64,
        smarts: UniquePtr<CxxString>,
        v: Vec<f64>,
        signs: Vec<i32>,
        atom_indices: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/embedder.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type EmbedFailureCauses;

        fn embed_parameters_preset(preset: EmbedPreset) -> Result<EmbedParametersData>;
        fn embed_multiple_confs(
            mol: Pin<&mut RWMol>,
            num_confs: u32,
            params: &mut EmbedParametersData,
        ) -> Result<Vec<i32>>;
        fn get_molecule_bounds_matrix(
            mol: &RWMol,
            set15bounds: bool,
            scale_vdw: bool,
            do_triangle_smoothing: bool,
            use_macrocycle14config: bool,
        ) -> Result<Vec<f64>>;
        #[allow(clippy::too_many_arguments)]
        fn get_experimental_torsions(
            mol: &RWMol,
            use_exp_torsion_angle_prefs: bool,
            use_small_ring_torsions: bool,
            use_macrocycle_torsions: bool,
            use_basic_knowledge: bool,
            et_version: u32,
            verbose: bool,
        ) -> Result<Vec<ExperimentalTorsionData>>;
    }
}
