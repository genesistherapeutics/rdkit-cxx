// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/frag_catalog.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type FragCatParams;
        type FragCatalog;

        fn frag_cat_params_new(
            l_len: u32,
            u_len: u32,
            fgroup_filename: &[u8],
            tol: f64,
        ) -> Result<UniquePtr<FragCatParams>>;
        fn frag_cat_params_get_type_string(params: &FragCatParams) -> Result<UniquePtr<CxxString>>;
        fn frag_cat_params_get_upper_frag_length(params: &FragCatParams) -> Result<u32>;
        fn frag_cat_params_get_lower_frag_length(params: &FragCatParams) -> Result<u32>;
        fn frag_cat_params_get_tolerance(params: &FragCatParams) -> Result<f64>;
        fn frag_cat_params_get_num_func_groups(params: &FragCatParams) -> Result<u32>;
        fn frag_cat_params_get_func_group(
            params: &FragCatParams,
            fid: u32,
        ) -> Result<UniquePtr<RWMol>>;
        fn frag_cat_params_serialize(params: &FragCatParams) -> Result<Vec<u8>>;

        fn frag_catalog_new() -> Result<UniquePtr<FragCatalog>>;
        fn frag_catalog_from_params(params: &FragCatParams) -> Result<UniquePtr<FragCatalog>>;
        fn frag_catalog_from_pickle(pickle: &[u8]) -> Result<UniquePtr<FragCatalog>>;
        fn frag_catalog_get_num_entries(catalog: &FragCatalog) -> Result<u32>;
        fn frag_catalog_get_fp_length(catalog: &FragCatalog) -> Result<u32>;
        fn frag_catalog_get_catalog_params(
            catalog: &FragCatalog,
        ) -> Result<UniquePtr<FragCatParams>>;
        fn frag_catalog_serialize(catalog: &FragCatalog) -> Result<Vec<u8>>;
        fn frag_catalog_get_bit_description(
            catalog: &FragCatalog,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn frag_catalog_get_bit_order(catalog: &FragCatalog, idx: u32) -> Result<u32>;
        fn frag_catalog_get_bit_func_group_ids(catalog: &FragCatalog, idx: u32)
        -> Result<Vec<i32>>;
        fn frag_catalog_get_bit_entry_id(catalog: &FragCatalog, idx: u32) -> Result<i64>;
        fn frag_catalog_get_bit_discrims(catalog: &FragCatalog, idx: u32) -> Result<Vec<f64>>;
        fn frag_catalog_get_entry_bit_id(catalog: &FragCatalog, idx: u32) -> Result<i32>;
        fn frag_catalog_get_entry_description(
            catalog: &FragCatalog,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn frag_catalog_get_entry_order(catalog: &FragCatalog, idx: u32) -> Result<u32>;
        fn frag_catalog_get_entry_func_group_ids(
            catalog: &FragCatalog,
            idx: u32,
        ) -> Result<Vec<i32>>;
        fn frag_catalog_get_entry_down_ids(catalog: &FragCatalog, idx: u32) -> Result<Vec<u32>>;

        fn frag_cat_generator_add_frags_from_mol(
            mol: &RWMol,
            fcat: Pin<&mut FragCatalog>,
        ) -> Result<u32>;
        fn frag_fp_generator_get_fp_for_mol(
            mol: &RWMol,
            fcat: &FragCatalog,
        ) -> Result<UniquePtr<BitVect>>;
    }
}

// SAFETY: parameters own their functional groups (molecules read from the
// file or a pickle, or copies), and a catalog its parameters, entries, and
// their fragments, with no thread affinity. Neither is Sync: like molecules
// (CONVENTIONS.md S1), and an entry's discriminators, which the bit lookups
// and the generators read, are cached in its properties through const.
unsafe impl Send for ffi::FragCatParams {}
unsafe impl Send for ffi::FragCatalog {}
