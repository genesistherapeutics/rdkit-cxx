// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::panic::{AssertUnwindSafe, catch_unwind};

use cxx::{CxxString, UniquePtr};

use crate::chem::Mol;
use crate::chem::rdmolops::mol_list_into_vec;
use crate::error::Error;
use crate::ffi::standardize_validate::panic_message;

/// An incremental search's closure, lent to C++ for one call, with the
/// first error it returned (or the panic it raised), which the search then
/// returns in place of the exception C++ saw.
pub(crate) struct SearchCallback<'a> {
    pub(crate) callback: &'a mut dyn FnMut(Vec<Mol>) -> crate::Result<bool>,
    pub(crate) error: Option<Error>,
}

/// Called by the search, on the calling thread, with each batch of hits. An
/// error or a panic is caught here, at the boundary, and becomes the
/// exception C++ sees, which ends the search.
///
/// # Safety
///
/// cxx makes a function whose signature names a lifetime `unsafe` for C++:
/// the C++ caller must not keep `callback` past the call that lent it, which
/// the wrappers do not.
unsafe fn call_search_callback(
    callback: &mut SearchCallback<'_>,
    hits: UniquePtr<ffi::MolList>,
) -> Result<bool, String> {
    let hits = mol_list_into_vec(hits);
    match catch_unwind(AssertUnwindSafe(|| (callback.callback)(hits))) {
        Ok(Ok(stop)) => Ok(stop),
        Ok(Err(error)) => {
            let message = error.to_string();
            callback.error.get_or_insert(error);
            Err(message)
        }
        Err(payload) => {
            let message = panic_message(payload, "a search callback");
            callback.error.get_or_insert(Error::Runtime {
                message: message.clone(),
            });
            Err(message)
        }
    }
}

/// A user conformer generator, as `ShapeBuildParams` holds it.
pub(crate) type ConformerFn = dyn Fn(&str, u32) -> crate::Result<Option<Mol>> + Send + Sync;

/// A shape build's conformer generator, lent to C++ for one call, with the
/// first error it returned (or the panic it raised), which the build then
/// returns in place of the exception C++ saw.
pub(crate) struct ConformerCallback<'a> {
    pub(crate) generator: Option<&'a ConformerFn>,
    pub(crate) error: Option<Error>,
}

/// Called by the shape build, on the calling thread, with an isomer's
/// SMILES; a null molecule for `None`. An error or a panic is caught here,
/// at the boundary, and becomes the exception C++ sees.
///
/// # Safety
///
/// As `call_search_callback`'s.
unsafe fn call_conformer_generator(
    callback: &mut ConformerCallback<'_>,
    smiles: &CxxString,
    num_conformers: u32,
) -> Result<UniquePtr<ffi::RWMol>, String> {
    let Some(generator) = callback.generator else {
        return Err("RDKit called a conformer generator Rust did not set".to_string());
    };
    let call = || -> crate::Result<Option<Mol>> {
        let smiles = smiles.to_str().map_err(|err| Error::Value {
            message: format!("RDKit's SMILES is not UTF-8: {err}"),
        })?;
        generator(smiles, num_conformers)
    };
    match catch_unwind(AssertUnwindSafe(call)) {
        Ok(Ok(Some(mol))) => Ok(mol.into_ffi()),
        Ok(Ok(None)) => Ok(UniquePtr::null()),
        Ok(Err(error)) => {
            let message = error.to_string();
            callback.error.get_or_insert(error);
            Err(message)
        }
        Err(payload) => {
            let message = panic_message(payload, "a conformer generator");
            callback.error.get_or_insert(Error::Runtime {
                message: message.clone(),
            });
            Err(message)
        }
    }
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `SynthonSpaceSearch::SynthonSpaceSearchParams`, without the fields
    /// the bindings leave out (the manifest says why). The possible hits
    /// file is a path's bytes (A8).
    struct SynthonSpaceSearchParamsData {
        max_hits: i64,
        max_num_frags: u64,
        hit_start: i64,
        random_sample: bool,
        random_seed: i32,
        build_hits: bool,
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
        to_try_chunk_size: u64,
        num_threads: i32,
        possible_hits_file: Vec<u8>,
        max_possible_hits_to_write: u64,
        write_possible_hits_and_stop: bool,
    }

    /// `SynthonSpaceSearch::ShapeBuildParams`, with its
    /// `StereoEnumerationOptions` as plain fields and a flag for the user
    /// conformer generator, which crosses as a `ConformerCallback`.
    struct ShapeBuildParamsData {
        num_confs: u32,
        rms_threshold: f64,
        shape_sim_threshold: f64,
        num_threads: i32,
        random_seed: i32,
        stereo_try_embedding: bool,
        stereo_only_unassigned: bool,
        stereo_only_stereo_groups: bool,
        stereo_unique: bool,
        stereo_max_isomers: u64,
        stereo_random_seed: i32,
        max_synthon_atoms: u32,
        max_embed_attempts: u32,
        time_out: u32,
        interim_file: Vec<u8>,
        interim_writes: u64,
        has_user_conformer_generator: bool,
    }

    extern "Rust" {
        type SearchCallback<'a>;
        type ConformerCallback<'a>;

        unsafe fn call_conformer_generator<'a>(
            callback: &mut ConformerCallback<'a>,
            smiles: &CxxString,
            num_conformers: u32,
        ) -> Result<UniquePtr<RWMol>>;

        unsafe fn call_search_callback<'a>(
            callback: &mut SearchCallback<'a>,
            hits: UniquePtr<MolList>,
        ) -> Result<bool>;
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
        fn shape_build_params_default() -> Result<ShapeBuildParamsData>;

        fn synthon_space_new() -> Result<UniquePtr<SynthonSpace>>;
        fn synthon_space_read_text_file(space: Pin<&mut SynthonSpace>, path: &[u8]) -> Result<()>;
        fn synthon_space_read_db_file(space: Pin<&mut SynthonSpace>, path: &[u8]) -> Result<()>;
        fn synthon_space_write_db_file(space: &SynthonSpace, path: &[u8]) -> Result<()>;
        fn synthon_space_write_enumerated_file(space: &SynthonSpace, path: &[u8]) -> Result<()>;
        fn synthon_space_get_num_reactions(space: &SynthonSpace) -> Result<u64>;
        fn synthon_space_get_num_products(space: &SynthonSpace) -> Result<u64>;
        fn synthon_space_get_num_synthons(space: &SynthonSpace) -> Result<u64>;
        fn synthon_space_get_num_synthons_with_shapes(space: &SynthonSpace) -> Result<u64>;
        fn synthon_space_summarise(space: &SynthonSpace) -> Result<()>;
        fn synthon_space_report_synthon_usage(space: &SynthonSpace) -> Result<()>;
        fn synthon_space_get_synthon_fingerprint_type(
            space: &SynthonSpace,
        ) -> Result<UniquePtr<CxxString>>;
        fn synthon_space_build_synthon_fingerprints(
            space: Pin<&mut SynthonSpace>,
            generator: &FpGenerator,
        ) -> Result<()>;

        unsafe fn synthon_space_build_synthon_shapes<'a>(
            space: Pin<&mut SynthonSpace>,
            params: &ShapeBuildParamsData,
            generator: &mut ConformerCallback<'a>,
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
        fn synthon_space_substructure_search_lines(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            match_params: &SubstructMatchParametersData,
            params: &SynthonSpaceSearchParamsData,
            start_line: u64,
            finish_line: u64,
        ) -> Result<UniquePtr<SearchResults>>;
        fn synthon_space_substructure_search_extended_lines(
            space: Pin<&mut SynthonSpace>,
            query: &ExtendedQueryMol,
            match_params: &SubstructMatchParametersData,
            params: &SynthonSpaceSearchParamsData,
            start_line: u64,
            finish_line: u64,
        ) -> Result<UniquePtr<SearchResults>>;
        unsafe fn synthon_space_substructure_search_incremental<'a>(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            callback: &mut SearchCallback<'a>,
            match_params: &SubstructMatchParametersData,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<()>;

        fn synthon_space_fingerprint_search(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            generator: &FpGenerator,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<UniquePtr<SearchResults>>;
        fn synthon_space_fingerprint_search_lines(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            generator: &FpGenerator,
            params: &SynthonSpaceSearchParamsData,
            start_line: u64,
            finish_line: u64,
        ) -> Result<UniquePtr<SearchResults>>;
        unsafe fn synthon_space_fingerprint_search_incremental<'a>(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            generator: &FpGenerator,
            callback: &mut SearchCallback<'a>,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<()>;

        fn synthon_space_rascal_search(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            rascal_options: &RascalOptionsData,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<UniquePtr<SearchResults>>;
        fn synthon_space_rascal_search_lines(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            rascal_options: &RascalOptionsData,
            params: &SynthonSpaceSearchParamsData,
            start_line: u64,
            finish_line: u64,
        ) -> Result<UniquePtr<SearchResults>>;
        unsafe fn synthon_space_rascal_search_incremental<'a>(
            space: Pin<&mut SynthonSpace>,
            query: &RWMol,
            rascal_options: &RascalOptionsData,
            callback: &mut SearchCallback<'a>,
            params: &SynthonSpaceSearchParamsData,
        ) -> Result<()>;

        fn substructure_result_get_hit_molecules(
            results: &SearchResults,
        ) -> Result<UniquePtr<MolList>>;
        fn substructure_result_get_max_num_results(results: &SearchResults) -> Result<u64>;
        fn substructure_result_get_timed_out(results: &SearchResults) -> Result<bool>;
        fn substructure_result_get_cancelled(results: &SearchResults) -> Result<bool>;
        fn substructure_result_get_best_hit(results: &SearchResults) -> Result<UniquePtr<RWMol>>;

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
