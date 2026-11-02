// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/substruct_library.h");
        include!("rdkit-cxx/ffi/substruct.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolBundle = crate::ffi::mol_bundle::ffi::MolBundle;
        type TautomerQuery = crate::ffi::tautomer_query::ffi::TautomerQuery;
        type ExtendedQueryMol = crate::ffi::generalized_substruct::ffi::ExtendedQueryMol;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type SubstructMatchParametersData =
            crate::ffi::substruct::ffi::SubstructMatchParametersData;
        type MolHolderHandle;
        type FpHolderHandle;
        type KeyHolderHandle;
        type SubstructLibrary;

        fn mol_holder_new(kind: u8) -> Result<UniquePtr<MolHolderHandle>>;
        fn mol_holder_kind(holder: &MolHolderHandle) -> Result<u8>;
        fn mol_holder_len(holder: &MolHolderHandle) -> Result<u32>;
        fn mol_holder_add_mol(holder: Pin<&mut MolHolderHandle>, mol: &RWMol) -> Result<u32>;
        fn mol_holder_add_binary(holder: Pin<&mut MolHolderHandle>, pickle: &[u8]) -> Result<u32>;
        fn mol_holder_add_smiles(holder: Pin<&mut MolHolderHandle>, smiles: &str) -> Result<u32>;
        fn mol_holder_get_mol(holder: &MolHolderHandle, idx: u32) -> Result<UniquePtr<RWMol>>;

        fn fp_holder_new(kind: u8, num_bits: u32) -> Result<UniquePtr<FpHolderHandle>>;
        fn fp_holder_kind(holder: &FpHolderHandle) -> Result<u8>;
        fn fp_holder_num_bits(holder: &FpHolderHandle) -> Result<u32>;
        fn fp_holder_len(holder: &FpHolderHandle) -> Result<u32>;
        fn fp_holder_add_mol(holder: Pin<&mut FpHolderHandle>, mol: &RWMol) -> Result<u32>;
        fn fp_holder_add_fingerprint(holder: Pin<&mut FpHolderHandle>, v: &BitVect) -> Result<u32>;
        fn fp_holder_get_fingerprint(
            holder: &FpHolderHandle,
            idx: u32,
        ) -> Result<UniquePtr<BitVect>>;
        fn fp_holder_passes_filter(
            holder: &FpHolderHandle,
            idx: u32,
            query: &BitVect,
        ) -> Result<bool>;
        fn fp_holder_make_fingerprint(
            holder: &FpHolderHandle,
            mol: &RWMol,
        ) -> Result<UniquePtr<BitVect>>;

        fn key_holder_new(propname: &str) -> Result<UniquePtr<KeyHolderHandle>>;
        fn key_holder_len(holder: &KeyHolderHandle) -> Result<u32>;
        fn key_holder_add_mol(holder: Pin<&mut KeyHolderHandle>, mol: &RWMol) -> Result<u32>;
        fn key_holder_add_key(holder: Pin<&mut KeyHolderHandle>, key: &str) -> Result<u32>;
        fn key_holder_get_key(holder: &KeyHolderHandle, idx: u32) -> Result<UniquePtr<CxxString>>;
        fn key_holder_get_keys(
            holder: &KeyHolderHandle,
            indices: &[u32],
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn key_holder_get_prop_name(holder: &KeyHolderHandle) -> Result<UniquePtr<CxxString>>;

        fn substruct_library_view(mols: &MolHolderHandle) -> Result<UniquePtr<SubstructLibrary>>;
        fn substruct_library_view_set_fp_holder(
            view: Pin<&mut SubstructLibrary>,
            fps: &FpHolderHandle,
        ) -> Result<()>;
        fn substruct_library_view_set_key_holder(
            view: Pin<&mut SubstructLibrary>,
            keys: &KeyHolderHandle,
        ) -> Result<()>;
        fn substruct_library_view_set_search_order(
            view: Pin<&mut SubstructLibrary>,
            order: &[u32],
        ) -> Result<()>;
        fn substruct_library_add_mol(view: Pin<&mut SubstructLibrary>, mol: &RWMol) -> Result<u32>;
        fn substruct_library_get_matches(
            view: &SubstructLibrary,
            query: &RWMol,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
            max_results: i32,
        ) -> Result<Vec<u32>>;
        fn substruct_library_get_matches_mol_bundle(
            view: &SubstructLibrary,
            query: &MolBundle,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
            max_results: i32,
        ) -> Result<Vec<u32>>;
        fn substruct_library_count_matches(
            view: &SubstructLibrary,
            query: &RWMol,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<u32>;
        fn substruct_library_count_matches_mol_bundle(
            view: &SubstructLibrary,
            query: &MolBundle,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<u32>;
        fn substruct_library_has_match(
            view: &SubstructLibrary,
            query: &RWMol,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn substruct_library_has_match_mol_bundle(
            view: &SubstructLibrary,
            query: &MolBundle,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn substruct_library_get_matches_tautomer_query(
            view: &SubstructLibrary,
            query: &TautomerQuery,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
            max_results: i32,
        ) -> Result<Vec<u32>>;
        fn substruct_library_count_matches_tautomer_query(
            view: &SubstructLibrary,
            query: &TautomerQuery,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<u32>;
        fn substruct_library_has_match_tautomer_query(
            view: &SubstructLibrary,
            query: &TautomerQuery,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn substruct_library_get_matches_extended_query_mol(
            view: &SubstructLibrary,
            query: &ExtendedQueryMol,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
            max_results: i32,
        ) -> Result<Vec<u32>>;
        fn substruct_library_count_matches_extended_query_mol(
            view: &SubstructLibrary,
            query: &ExtendedQueryMol,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<u32>;
        fn substruct_library_has_match_extended_query_mol(
            view: &SubstructLibrary,
            query: &ExtendedQueryMol,
            ranged: bool,
            start_idx: u32,
            end_idx: u32,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn add_patterns(view: &SubstructLibrary, patterns: Pin<&mut FpHolderHandle>) -> Result<()>;

        fn substruct_library_can_serialize() -> Result<bool>;
        fn substruct_library_serialize(view: &SubstructLibrary) -> Result<Vec<u8>>;
        fn substruct_library_from_pickle(data: &[u8]) -> Result<UniquePtr<SubstructLibrary>>;
        fn substruct_library_mol_holder(
            lib: &SubstructLibrary,
        ) -> Result<UniquePtr<MolHolderHandle>>;
        fn substruct_library_fp_holder(lib: &SubstructLibrary)
        -> Result<UniquePtr<FpHolderHandle>>;
        fn substruct_library_key_holder(
            lib: &SubstructLibrary,
        ) -> Result<UniquePtr<KeyHolderHandle>>;
        fn substruct_library_search_order(lib: &SubstructLibrary) -> Result<Vec<u32>>;
    }
}

// SAFETY: each holder is owned by its one Rust value: its shared_ptr is
// copied only into the views a call builds and drops, or held by a library
// RDKit read from a pickle until Rust takes the holders out. A holder owns
// what it holds (molecules, strings, fingerprints, keys), so moving it moves
// everything RDKit reads. It is not Sync: a search writes the ring info of a
// MolHolder's molecules through const (MolOps::symmetrizeSSSR in
// SubSearcher), as S1 says of molecules.
unsafe impl Send for ffi::MolHolderHandle {}
// SAFETY: as for MolHolderHandle; reading a fingerprint holder writes
// nothing, but it is not Sync either, since nothing needs it to be.
unsafe impl Send for ffi::FpHolderHandle {}
// SAFETY: as for FpHolderHandle.
unsafe impl Send for ffi::KeyHolderHandle {}
