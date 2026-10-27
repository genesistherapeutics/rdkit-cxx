// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `SynthonSpaceSearch::SynthonSpaceSearchParams`, with the fields
    /// Python binds.
    struct SynthonSpaceSearchParamsData {
        max_hits: i64,
        max_num_frags: u64,
        hit_start: i64,
        random_sample: bool,
        random_seed: i32,
        build_hits: bool,
        num_random_sweeps: i32,
        similarity_cutoff: f64,
        frag_similarity_adjuster: f64,
        approx_similarity_adjuster: f64,
        min_hit_heavy_atoms: u32,
        max_hit_heavy_atoms: i32,
        min_hit_mol_wt: f64,
        max_hit_mol_wt: f64,
        min_hit_chiral_atoms: u32,
        max_hit_chiral_atoms: i32,
        time_out: u64,
        num_threads: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/synthon_space_search.h");
        include!("rdkit-cxx/ffi/substruct.rs.h");
        include!("rdkit-cxx/ffi/rascal_mces.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type FpGenerator = crate::ffi::fingerprint_generator::ffi::FpGenerator;
        type ExtendedQueryMol = crate::ffi::generalized_substruct::ffi::ExtendedQueryMol;
        type SubstructMatchParametersData =
            crate::ffi::substruct::ffi::SubstructMatchParametersData;
        type RascalOptionsData = crate::ffi::rascal_mces::ffi::RascalOptionsData;
        type SynthonSpace;
        type SearchResults;

        fn synthon_space_search_params_default() -> Result<SynthonSpaceSearchParamsData>;

        fn synthon_space_new() -> Result<UniquePtr<SynthonSpace>>;
        fn synthon_space_read_text_file(space: Pin<&mut SynthonSpace>, path: &[u8]) -> Result<()>;
        fn synthon_space_read_db_file(space: Pin<&mut SynthonSpace>, path: &[u8]) -> Result<()>;
        fn synthon_space_write_db_file(space: &SynthonSpace, path: &[u8]) -> Result<()>;
        fn synthon_space_write_enumerated_file(space: &SynthonSpace, path: &[u8]) -> Result<()>;
        fn synthon_space_get_num_reactions(space: &SynthonSpace) -> Result<u64>;
        fn synthon_space_get_num_products(space: &SynthonSpace) -> Result<u64>;
        fn synthon_space_summarise(space: &SynthonSpace) -> Result<()>;
        fn synthon_space_get_synthon_fingerprint_type(
            space: &SynthonSpace,
        ) -> Result<UniquePtr<CxxString>>;
        fn synthon_space_build_synthon_fingerprints(
            space: Pin<&mut SynthonSpace>,
            generator: &FpGenerator,
        ) -> Result<()>;

        fn synthon_space_substructure_search(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            match_params: &SubstructMatchParametersData,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<UniquePtr<SearchResults>>;
        fn synthon_space_substructure_search_extended(
            space: Pin<&mut SynthonSpace>,
            query: &ExtendedQueryMol,
            match_params: &SubstructMatchParametersData,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<UniquePtr<SearchResults>>;
        fn synthon_space_fingerprint_search(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            generator: &FpGenerator,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<UniquePtr<SearchResults>>;
        fn synthon_space_rascal_search(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            rascal_options: &RascalOptionsData,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<UniquePtr<SearchResults>>;

        fn substructure_result_get_hit_molecules(
            results: &SearchResults,
        ) -> Result<UniquePtr<MolList>>;
        fn substructure_result_get_max_num_results(results: &SearchResults) -> Result<u64>;
        fn substructure_result_get_timed_out(results: &SearchResults) -> Result<bool>;
        fn substructure_result_get_cancelled(results: &SearchResults) -> Result<bool>;

        fn formatted_integer_string(value: i64) -> Result<UniquePtr<CxxString>>;
    }
}

// SAFETY: a space owns its synthons and reactions, and nothing else holds
// them. Its molecules are RDKit's, which its const methods may write
// through (computed properties, CONVENTIONS.md S1), so it is not Sync; the
// Rust handle calls its searches through &mut.
unsafe impl Send for ffi::SynthonSpace {}

// SAFETY: a result owns its hit molecules, and nothing else holds them; as
// molecules, they are not Sync (S1).
unsafe impl Send for ffi::SearchResults {}
