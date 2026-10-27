// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `rdDepictor.Compute2DCoords`'s arguments; the coordinate map as
    /// parallel keys and values (A7).
    struct Compute2DCoordsData {
        canon_orient: bool,
        clear_confs: bool,
        coord_map_keys: Vec<u32>,
        coord_map_values: Vec<Point2D>,
        n_flips_per_sample: u32,
        n_sample: u32,
        sample_seed: i32,
        permute_deg4_nodes: i32,
        bond_length: f64,
        force_rdkit: bool,
        use_ring_templates: bool,
    }

    /// `rdDepictor.Compute2DCoordsMimicDistmat`'s optional arguments.
    struct Compute2DCoordsMimicDistmatData {
        canon_orient: bool,
        clear_confs: bool,
        weight_dist_mat: f64,
        n_flips_per_sample: u32,
        n_sample: u32,
        sample_seed: i32,
        permute_deg4_nodes: bool,
        bond_length: f64,
        force_rdkit: bool,
    }

    /// `RDDepict::ConstrainedDepictionParams`, plus the reference's
    /// conformer.
    struct ConstrainedDepictionParamsData {
        accept_failure: bool,
        force_rdkit: bool,
        allow_r_groups: bool,
        align_only: bool,
        adjust_mol_block_wedging: bool,
        existing_conf_id: i32,
        conf_id: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/depictor.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type Point2D = crate::ffi::geometry::ffi::Point2D;

        fn is_coord_gen_support_available() -> Result<bool>;
        fn set_prefer_coord_gen(val: bool) -> Result<()>;
        fn get_prefer_coord_gen() -> Result<bool>;
        fn replace_prefer_coord_gen(val: bool) -> Result<bool>;

        fn set_ring_system_templates(path: &[u8]) -> Result<()>;
        fn add_ring_system_templates(path: &[u8]) -> Result<()>;
        fn load_default_ring_system_templates() -> Result<()>;

        fn constrained_depiction_params_default() -> Result<ConstrainedDepictionParamsData>;

        fn compute_2d_coords(mol: Pin<&mut RWMol>, data: &Compute2DCoordsData) -> Result<u32>;
        fn compute_2d_coords_mimic_distmat(
            mol: Pin<&mut RWMol>,
            dist_mat: &[f64],
            data: &Compute2DCoordsMimicDistmatData,
        ) -> Result<u32>;

        fn generate_depiction_matching_2d_structure(
            mol: Pin<&mut RWMol>,
            reference: &RWMol,
            ref_patt: &RWMol,
            has_ref_patt: bool,
            data: &ConstrainedDepictionParamsData,
        ) -> Result<Vec<IndexPair>>;
        fn generate_depiction_matching_2d_structure_atom_map(
            mol: Pin<&mut RWMol>,
            reference: &RWMol,
            atom_map: &[IndexPair],
            data: &ConstrainedDepictionParamsData,
        ) -> Result<()>;
        fn generate_depiction_matching_3d_structure(
            mol: Pin<&mut RWMol>,
            reference: &RWMol,
            conf_id: i32,
            ref_patt: &RWMol,
            has_ref_patt: bool,
            accept_failure: bool,
            force_rdkit: bool,
        ) -> Result<()>;

        fn straighten_depiction(
            mol: Pin<&mut RWMol>,
            conf_id: i32,
            minimize_rotation: bool,
        ) -> Result<()>;
        fn normalize_depiction(
            mol: Pin<&mut RWMol>,
            conf_id: i32,
            canonicalize: i32,
            scale_factor: f64,
        ) -> Result<f64>;
    }
}
