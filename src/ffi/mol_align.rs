// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// The best-alignment functions' keyword arguments, with their atom maps
    /// flattened (each map's length in `map_counts`).
    struct BestAlignmentParamsData {
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        num_threads: i32,
        map_pairs: Vec<IndexPair>,
        map_counts: Vec<u32>,
        weights: Vec<f64>,
        prb_cid: i32,
        ref_cid: i32,
        reflect: bool,
        max_iters: u32,
    }

    /// The options every O3A function takes.
    struct O3AOptionsData {
        prb_cid: i32,
        ref_cid: i32,
        reflect: bool,
        max_iters: u32,
        options: u32,
        num_threads: i32,
        constraint_map: Vec<IndexPair>,
        has_constraint_weights: bool,
        constraint_weights: Vec<f64>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_align.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type MMFFMolProperties = crate::ffi::force_field::ffi::MMFFMolProperties;
        type O3AHolder;
        type O3AList;

        #[allow(clippy::too_many_arguments)]
        fn get_alignment_transform(
            prb_mol: &RWMol,
            ref_mol: &RWMol,
            prb_cid: i32,
            ref_cid: i32,
            atom_map: &[IndexPair],
            weights: &[f64],
            reflect: bool,
            max_iters: u32,
            trans: &mut Vec<f64>,
        ) -> Result<f64>;
        #[allow(clippy::too_many_arguments)]
        fn align_mol(
            prb_mol: Pin<&mut RWMol>,
            ref_mol: &RWMol,
            prb_cid: i32,
            ref_cid: i32,
            atom_map: &[IndexPair],
            weights: &[f64],
            reflect: bool,
            max_iters: u32,
        ) -> Result<f64>;
        fn get_best_alignment_transform(
            prb_mol: &RWMol,
            ref_mol: &RWMol,
            params: &BestAlignmentParamsData,
            trans: &mut Vec<f64>,
            best_match: &mut Vec<IndexPair>,
        ) -> Result<f64>;
        fn get_best_rms(
            prb_mol: Pin<&mut RWMol>,
            ref_mol: &RWMol,
            params: &BestAlignmentParamsData,
        ) -> Result<f64>;
        fn get_all_conformer_best_rms(
            mol: &RWMol,
            params: &BestAlignmentParamsData,
        ) -> Result<Vec<f64>>;
        fn calc_rms(
            prb_mol: &RWMol,
            ref_mol: &RWMol,
            params: &BestAlignmentParamsData,
        ) -> Result<f64>;
        #[allow(clippy::too_many_arguments)]
        fn align_mol_conformers(
            mol: Pin<&mut RWMol>,
            atom_ids: &[u32],
            conf_ids: &[u32],
            weights: &[f64],
            reflect: bool,
            max_iters: u32,
        ) -> Result<Vec<f64>>;
        fn random_transform(mol: Pin<&mut RWMol>, cid: i32, seed: i32) -> Result<()>;

        /// # Safety
        ///
        /// The properties are null or valid.
        unsafe fn get_mmff_o3a(
            prb_mol: Pin<&mut RWMol>,
            ref_mol: &RWMol,
            prb_props: *const MMFFMolProperties,
            ref_props: *const MMFFMolProperties,
            for_probe_confs: bool,
            options: &O3AOptionsData,
        ) -> Result<UniquePtr<O3AList>>;
        #[allow(clippy::too_many_arguments)]
        fn get_crippen_o3a(
            prb_mol: Pin<&mut RWMol>,
            ref_mol: &RWMol,
            has_prb_contribs: bool,
            prb_contribs: &[f64],
            has_ref_contribs: bool,
            ref_contribs: &[f64],
            for_probe_confs: bool,
            options: &O3AOptionsData,
        ) -> Result<UniquePtr<O3AList>>;
        fn o3a_list_size(list: &O3AList) -> Result<usize>;
        fn o3a_list_take(list: Pin<&mut O3AList>, i: usize) -> Result<UniquePtr<O3AHolder>>;
        fn o3a_align(o3a: Pin<&mut O3AHolder>) -> Result<f64>;
        fn o3a_trans(o3a: &O3AHolder, trans: &mut Vec<f64>) -> Result<f64>;
        fn o3a_score(o3a: &O3AHolder) -> Result<f64>;
        fn o3a_matches(o3a: &O3AHolder) -> Result<Vec<IndexPair>>;
        fn o3a_weights(o3a: &O3AHolder) -> Result<Vec<f64>>;
    }
}
