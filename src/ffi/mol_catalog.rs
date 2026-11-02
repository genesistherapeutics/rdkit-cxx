// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_catalog.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolCatalog;
        type MolCatalogEntry;

        fn create_mol_catalog() -> Result<UniquePtr<MolCatalog>>;
        fn mol_catalog_from_pickle(pickle: &[u8]) -> Result<UniquePtr<MolCatalog>>;
        fn mol_catalog_get_num_entries(catalog: &MolCatalog) -> Result<u32>;
        fn mol_catalog_get_fp_length(catalog: &MolCatalog) -> Result<u32>;
        fn mol_catalog_serialize(catalog: &MolCatalog) -> Result<Vec<u8>>;
        fn mol_catalog_get_bit_description(
            catalog: &MolCatalog,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_catalog_get_bit_entry_id(catalog: &MolCatalog, idx: u32) -> Result<i64>;
        fn mol_catalog_get_entry_bit_id(catalog: &MolCatalog, idx: u32) -> Result<i32>;
        fn mol_catalog_get_entry_description(
            catalog: &MolCatalog,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_catalog_get_entry_down_ids(catalog: &MolCatalog, idx: u32) -> Result<Vec<u32>>;
        fn mol_catalog_add_entry(
            catalog: Pin<&mut MolCatalog>,
            entry: &MolCatalogEntry,
        ) -> Result<u32>;
        fn mol_catalog_add_edge(catalog: Pin<&mut MolCatalog>, id1: u32, id2: u32) -> Result<()>;

        fn mol_catalog_entry_new() -> Result<UniquePtr<MolCatalogEntry>>;
        fn mol_catalog_entry_from_pickle(pickle: &[u8]) -> Result<UniquePtr<MolCatalogEntry>>;
        fn mol_catalog_entry_get_description(
            entry: &MolCatalogEntry,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_catalog_entry_set_description(
            entry: Pin<&mut MolCatalogEntry>,
            val: &str,
        ) -> Result<()>;
        fn mol_catalog_entry_get_mol(entry: &MolCatalogEntry) -> Result<UniquePtr<RWMol>>;
        fn mol_catalog_entry_set_mol(entry: Pin<&mut MolCatalogEntry>, mol: &RWMol) -> Result<()>;
        fn mol_catalog_entry_get_order(entry: &MolCatalogEntry) -> Result<u32>;
        fn mol_catalog_entry_set_order(entry: Pin<&mut MolCatalogEntry>, order: u32) -> Result<()>;
        fn mol_catalog_entry_serialize(entry: &MolCatalogEntry) -> Result<Vec<u8>>;
    }
}

// SAFETY: a catalog owns its entries, and an entry its molecule (copies the
// wrappers make, or molecules read from a pickle), with no thread affinity.
// Like molecules, they are not Sync (CONVENTIONS.md S1).
unsafe impl Send for ffi::MolCatalog {}
unsafe impl Send for ffi::MolCatalogEntry {}
