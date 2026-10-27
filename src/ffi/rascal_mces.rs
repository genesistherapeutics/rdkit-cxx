// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RascalMCES::RascalOptions`, without `doEquivBondPruning`, which
    /// Python does not expose.
    struct RascalOptionsData {
        similarity_threshold: f64,
        complete_aromatic_rings: bool,
        ring_matches_ring_only: bool,
        exact_connections_match: bool,
        single_largest_frag: bool,
        min_frag_size: i32,
        max_frag_separation: i32,
        all_best_mcess: bool,
        timeout: i32,
        return_empty_mces: bool,
        max_bond_match_pairs: u32,
        equivalent_atoms: String,
        ignore_bond_orders: bool,
        ignore_atom_aromaticity: bool,
        min_clique_size: u32,
    }

    /// `RascalMCES::RascalClusterOptions`.
    struct RascalClusterOptionsData {
        similarity_cutoff: f64,
        a: f64,
        b: f64,
        min_frag_size: u32,
        min_intra_cluster_sim: f64,
        cluster_merge_sim: f64,
        max_num_frags: u32,
        num_threads: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/rascal_mces.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type RascalResult;
        type RascalResultList;

        fn rascal_options_default() -> Result<RascalOptionsData>;
        fn rascal_cluster_options_default() -> Result<RascalClusterOptionsData>;

        fn find_mces(
            mol1: &RWMol,
            mol2: &RWMol,
            opts: &RascalOptionsData,
        ) -> Result<UniquePtr<RascalResultList>>;
        fn rascal_result_list_size(list: &RascalResultList) -> Result<usize>;
        fn rascal_result_list_take(
            list: Pin<&mut RascalResultList>,
            i: usize,
        ) -> Result<UniquePtr<RascalResult>>;

        fn rascal_result_clone(result: &RascalResult) -> Result<UniquePtr<RascalResult>>;
        fn rascal_result_smarts_string(result: &RascalResult) -> Result<UniquePtr<CxxString>>;
        fn rascal_result_bond_matches(result: &RascalResult) -> Result<Vec<IndexPair>>;
        fn rascal_result_atom_matches(result: &RascalResult) -> Result<Vec<IndexPair>>;
        fn rascal_result_largest_fragment_only(result: Pin<&mut RascalResult>) -> Result<()>;
        fn rascal_result_similarity(result: &RascalResult) -> Result<f64>;
        fn rascal_result_num_fragments(result: &RascalResult) -> Result<i32>;
        fn rascal_result_largest_fragment_size(result: &RascalResult) -> Result<u32>;
        fn rascal_result_tier1_sim(result: &RascalResult) -> Result<f64>;
        fn rascal_result_tier2_sim(result: &RascalResult) -> Result<f64>;
        fn rascal_result_timed_out(result: &RascalResult) -> Result<bool>;

        fn rascal_cluster(
            mols: Pin<&mut BorrowedMols>,
            opts: &RascalClusterOptionsData,
        ) -> Result<Vec<IndexList>>;
        fn rascal_butina_cluster(
            mols: Pin<&mut BorrowedMols>,
            opts: &RascalClusterOptionsData,
        ) -> Result<Vec<IndexList>>;
    }
}

// SAFETY: a RascalResult owns copies of its molecules (its constructors and
// copy constructor copy them). It is not Sync: its getters cache the SMARTS,
// fragment counts, and scores, and its molecules' distance matrices, in
// mutable members.
unsafe impl Send for ffi::RascalResult {}
