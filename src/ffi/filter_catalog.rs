// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::pin::Pin;

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::chem::rdfiltercatalog::{CustomFilterMatcher, FilterMatch};
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::mol::ffi::RWMol;

/// The matcher behind a `FilterMatcher` made in Rust. RDKit's catalogs share
/// their matchers between threads, so it is `Send + Sync`.
pub(crate) struct FilterMatcherCallback(pub(crate) Box<dyn CustomFilterMatcher>);

/// Runs a callback, catching a panic at the boundary; the error becomes the
/// exception C++ sees.
fn guarded<T>(call: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    catch_unwind(AssertUnwindSafe(call)).unwrap_or_else(|payload| {
        let message = payload
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "a filter matcher panicked".to_string());
        Err(format!("filter matcher panicked: {message}"))
    })
}

fn filter_matcher_callback_is_valid(callback: &FilterMatcherCallback) -> Result<bool, String> {
    guarded(|| Ok(callback.0.is_valid()))
}

/// Called with the callback's own copy of the molecule.
fn filter_matcher_callback_has_match(
    callback: &FilterMatcherCallback,
    mol: UniquePtr<RWMol>,
) -> Result<bool, String> {
    let mol = Mol::from_ffi(mol);
    guarded(|| Ok(callback.0.has_match(&mol)))
}

/// Called with the callback's own copy of the molecule; moves the matches it
/// returns into `list`.
fn filter_matcher_callback_get_matches(
    callback: &FilterMatcherCallback,
    mol: UniquePtr<RWMol>,
    mut list: Pin<&mut ffi::FilterMatchList>,
) -> Result<bool, String> {
    let mol = Mol::from_ffi(mol);
    guarded(|| {
        let mut matches = Vec::new();
        let matched = callback.0.get_matches(&mol, &mut matches);
        for m in &matches {
            push_match(list.as_mut(), m).map_err(|err| err.to_string())?;
        }
        Ok(matched)
    })
}

pub(crate) fn push_match(
    list: Pin<&mut ffi::FilterMatchList>,
    m: &FilterMatch,
) -> Result<(), cxx::Exception> {
    let pairs: Vec<IndexPair> = m
        .atom_pairs
        .iter()
        .map(|&(first, second)| IndexPair { first, second })
        .collect();
    ffi::filter_match_list_push(list, m.filter_match.raw(), &pairs)
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    extern "Rust" {
        type FilterMatcherCallback;

        fn filter_matcher_callback_is_valid(callback: &FilterMatcherCallback) -> Result<bool>;
        fn filter_matcher_callback_has_match(
            callback: &FilterMatcherCallback,
            mol: UniquePtr<RWMol>,
        ) -> Result<bool>;
        fn filter_matcher_callback_get_matches(
            callback: &FilterMatcherCallback,
            mol: UniquePtr<RWMol>,
            list: Pin<&mut FilterMatchList>,
        ) -> Result<bool>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/filter_catalog.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type FilterMatcherHandle;
        type FilterCatalogEntryHandle;
        type FilterMatchList;
        type FilterMatcherList;
        type FilterCatalogEntryList;
        type FilterCatalogParams;
        type FilterCatalog;

        fn filter_matcher_new(
            name: &str,
            callback: Box<FilterMatcherCallback>,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_matcher_copy(
            matcher: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_matcher_get_name(matcher: &FilterMatcherHandle) -> Result<UniquePtr<CxxString>>;
        fn filter_matcher_is_valid(matcher: &FilterMatcherHandle) -> Result<bool>;
        fn filter_matcher_has_match(matcher: &FilterMatcherHandle, mol: &RWMol) -> Result<bool>;
        fn filter_matcher_get_matches(
            matcher: &FilterMatcherHandle,
            mol: &RWMol,
        ) -> Result<UniquePtr<FilterMatchList>>;

        fn smarts_matcher_new(name: &str) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn smarts_matcher_from_mol(mol: &RWMol) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn smarts_matcher_new_mol(
            name: &str,
            mol: &RWMol,
            min_count: u32,
            max_count: u32,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn smarts_matcher_new_smarts(
            name: &str,
            smarts: &str,
            min_count: u32,
            max_count: u32,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn smarts_matcher_get_pattern(matcher: &FilterMatcherHandle) -> Result<UniquePtr<RWMol>>;
        fn smarts_matcher_set_pattern_mol(
            matcher: Pin<&mut FilterMatcherHandle>,
            mol: &RWMol,
        ) -> Result<()>;
        fn smarts_matcher_set_pattern_smarts(
            matcher: Pin<&mut FilterMatcherHandle>,
            smarts: &str,
        ) -> Result<()>;
        fn smarts_matcher_get_min_count(matcher: &FilterMatcherHandle) -> Result<u32>;
        fn smarts_matcher_set_min_count(
            matcher: Pin<&mut FilterMatcherHandle>,
            count: u32,
        ) -> Result<()>;
        fn smarts_matcher_get_max_count(matcher: &FilterMatcherHandle) -> Result<u32>;
        fn smarts_matcher_set_max_count(
            matcher: Pin<&mut FilterMatcherHandle>,
            count: u32,
        ) -> Result<()>;

        fn exclusion_list_new() -> Result<UniquePtr<FilterMatcherHandle>>;
        fn exclusion_list_add_pattern(
            list: Pin<&mut FilterMatcherHandle>,
            base: &FilterMatcherHandle,
        ) -> Result<()>;
        fn filter_matcher_list_new() -> Result<UniquePtr<FilterMatcherList>>;
        fn filter_matcher_list_push(
            list: Pin<&mut FilterMatcherList>,
            matcher: &FilterMatcherHandle,
        ) -> Result<()>;
        fn exclusion_list_set_exclusion_patterns(
            list: Pin<&mut FilterMatcherHandle>,
            patterns: &FilterMatcherList,
        ) -> Result<()>;

        fn filter_hierarchy_matcher_new() -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_hierarchy_matcher_from_matcher(
            matcher: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_hierarchy_matcher_set_pattern(
            node: Pin<&mut FilterMatcherHandle>,
            matcher: &FilterMatcherHandle,
        ) -> Result<()>;
        fn filter_hierarchy_matcher_add_child(
            node: Pin<&mut FilterMatcherHandle>,
            hierarchy: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;

        fn filter_match_ops_and(
            arg1: &FilterMatcherHandle,
            arg2: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_match_ops_or(
            arg1: &FilterMatcherHandle,
            arg2: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_match_ops_not(
            arg1: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;

        fn filter_match_list_size(list: &FilterMatchList) -> Result<usize>;
        fn filter_match_list_matcher(
            list: &FilterMatchList,
            i: usize,
        ) -> Result<UniquePtr<FilterMatcherHandle>>;
        fn filter_match_list_atom_pairs(list: &FilterMatchList, i: usize)
        -> Result<Vec<IndexPair>>;
        fn filter_match_list_push(
            list: Pin<&mut FilterMatchList>,
            matcher: &FilterMatcherHandle,
            atom_pairs: &[IndexPair],
        ) -> Result<()>;

        fn filter_catalog_entry_new() -> Result<UniquePtr<FilterCatalogEntryHandle>>;
        fn filter_catalog_entry_new_matcher(
            name: &str,
            matcher: &FilterMatcherHandle,
        ) -> Result<UniquePtr<FilterCatalogEntryHandle>>;
        fn filter_catalog_entry_share(
            entry: &FilterCatalogEntryHandle,
        ) -> Result<UniquePtr<FilterCatalogEntryHandle>>;
        fn filter_catalog_entry_same(
            a: &FilterCatalogEntryHandle,
            b: &FilterCatalogEntryHandle,
        ) -> Result<bool>;
        fn filter_catalog_entry_is_valid(entry: &FilterCatalogEntryHandle) -> Result<bool>;
        fn filter_catalog_entry_get_description(
            entry: &FilterCatalogEntryHandle,
        ) -> Result<UniquePtr<CxxString>>;
        fn filter_catalog_entry_set_description(
            entry: Pin<&mut FilterCatalogEntryHandle>,
            description: &str,
        ) -> Result<()>;
        fn filter_catalog_entry_get_filter_matches(
            entry: &FilterCatalogEntryHandle,
            mol: &RWMol,
        ) -> Result<UniquePtr<FilterMatchList>>;
        fn filter_catalog_entry_has_filter_match(
            entry: &FilterCatalogEntryHandle,
            mol: &RWMol,
        ) -> Result<bool>;
        fn filter_catalog_entry_serialize(entry: &FilterCatalogEntryHandle) -> Result<Vec<u8>>;
        fn filter_catalog_entry_get_prop_list(
            entry: &FilterCatalogEntryHandle,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn filter_catalog_entry_set_prop(
            entry: Pin<&mut FilterCatalogEntryHandle>,
            key: &str,
            val: &str,
        ) -> Result<()>;
        fn filter_catalog_entry_get_prop(
            entry: &FilterCatalogEntryHandle,
            key: &str,
        ) -> Result<UniquePtr<CxxString>>;
        fn filter_catalog_entry_clear_prop(
            entry: Pin<&mut FilterCatalogEntryHandle>,
            key: &str,
        ) -> Result<()>;

        fn filter_catalog_entry_list_size(list: &FilterCatalogEntryList) -> Result<usize>;
        fn filter_catalog_entry_list_get(
            list: &FilterCatalogEntryList,
            i: usize,
        ) -> Result<UniquePtr<FilterCatalogEntryHandle>>;
        fn filter_catalog_entry_list_row_ends(list: &FilterCatalogEntryList) -> Result<Vec<usize>>;

        fn filter_catalog_params_new() -> Result<UniquePtr<FilterCatalogParams>>;
        fn filter_catalog_params_from_catalogs(
            catalogs: u32,
        ) -> Result<UniquePtr<FilterCatalogParams>>;
        fn filter_catalog_params_add_catalog(
            params: Pin<&mut FilterCatalogParams>,
            catalogs: u32,
        ) -> Result<bool>;
        // Only the crate's unit tests call this.
        #[allow(dead_code)]
        fn filter_catalogs_values() -> Result<Vec<u32>>;

        fn filter_catalog_new() -> Result<UniquePtr<FilterCatalog>>;
        fn filter_catalog_from_params(
            params: &FilterCatalogParams,
        ) -> Result<UniquePtr<FilterCatalog>>;
        fn filter_catalog_from_catalogs(catalogs: u32) -> Result<UniquePtr<FilterCatalog>>;
        fn filter_catalog_from_pickle(pickle: &[u8]) -> Result<UniquePtr<FilterCatalog>>;
        fn filter_catalog_serialize(catalog: &FilterCatalog) -> Result<Vec<u8>>;
        fn filter_catalog_add_entry(
            catalog: Pin<&mut FilterCatalog>,
            entry: &FilterCatalogEntryHandle,
        ) -> Result<()>;
        fn filter_catalog_remove_entry(catalog: Pin<&mut FilterCatalog>, idx: u32) -> Result<bool>;
        fn filter_catalog_remove_entry_entry(
            catalog: Pin<&mut FilterCatalog>,
            entry: &FilterCatalogEntryHandle,
        ) -> Result<bool>;
        fn filter_catalog_get_num_entries(catalog: &FilterCatalog) -> Result<u32>;
        fn filter_catalog_get_entry(
            catalog: &FilterCatalog,
            idx: u32,
        ) -> Result<UniquePtr<FilterCatalogEntryHandle>>;
        fn filter_catalog_has_match(catalog: &FilterCatalog, mol: &RWMol) -> Result<bool>;
        fn filter_catalog_get_first_match(
            catalog: &FilterCatalog,
            mol: &RWMol,
        ) -> Result<UniquePtr<FilterCatalogEntryHandle>>;
        fn filter_catalog_get_matches(
            catalog: &FilterCatalog,
            mol: &RWMol,
        ) -> Result<UniquePtr<FilterCatalogEntryList>>;
        fn filter_catalog_get_filter_matches(
            catalog: &FilterCatalog,
            mol: &RWMol,
        ) -> Result<UniquePtr<FilterMatchList>>;

        fn filter_catalog_can_serialize() -> Result<bool>;
        fn run_filter_catalog(
            catalog: &FilterCatalog,
            smiles: &[&str],
            num_threads: i32,
        ) -> Result<UniquePtr<FilterCatalogEntryList>>;
        fn get_functional_group_hierarchy() -> Result<UniquePtr<FilterCatalog>>;
        fn get_flattened_functional_group_hierarchy(normalized: bool)
        -> Result<UniquePtr<MolList>>;
    }
}

// SAFETY: a matcher has no thread affinity, and its const methods are what
// RDKit's own threaded RunFilterCatalog runs on shared catalogs: RDKit's
// matchers only match (SubstructMatch locks the recursive queries it writes,
// under RDK_BUILD_THREADSAFE_SSS) and copy (sharing patterns, with atomic
// counts), and Rust's call a Send + Sync callback. The wrappers copy a
// shared pattern only with its recursive queries locked, and change only
// matchers their handle owns alone (filter_catalog.h).
unsafe impl Send for ffi::FilterMatcherHandle {}
unsafe impl Sync for ffi::FilterMatcherHandle {}
// SAFETY: an entry's const methods read its properties and run its matcher,
// as above; a handle changes the entry only when it is its only owner, and
// otherwise copies it first.
unsafe impl Send for ffi::FilterCatalogEntryHandle {}
unsafe impl Sync for ffi::FilterCatalogEntryHandle {}
// SAFETY: parameters are a list of catalog names.
unsafe impl Send for ffi::FilterCatalogParams {}
unsafe impl Sync for ffi::FilterCatalogParams {}
// SAFETY: a catalog owns its parameters and shares its entries, which no one
// changes in place; its const methods read them and run their matchers, as
// RDKit's RunFilterCatalog does from several threads.
unsafe impl Send for ffi::FilterCatalog {}
unsafe impl Sync for ffi::FilterCatalog {}
