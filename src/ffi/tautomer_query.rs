// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/tautomer_query.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");
        include!("rdkit-cxx/ffi/substruct.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type SubstructMatchParametersData =
            crate::ffi::substruct::ffi::SubstructMatchParametersData;
        type TautomerQuery;

        fn tautomer_query_clone(query: &TautomerQuery) -> Result<UniquePtr<TautomerQuery>>;
        fn tautomer_query_from_mol(
            mol: &RWMol,
            tautomer_transform_file: &[u8],
        ) -> Result<UniquePtr<TautomerQuery>>;
        fn tautomer_query_from_pickle(pickle: &[u8]) -> Result<UniquePtr<TautomerQuery>>;
        fn tautomer_query_to_binary(query: &TautomerQuery) -> Result<Vec<u8>>;

        fn tautomer_query_is_substruct_of(
            query: &TautomerQuery,
            target: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn tautomer_query_get_substruct_match(
            query: &TautomerQuery,
            target: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn tautomer_query_get_substruct_matches(
            query: &TautomerQuery,
            target: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
        fn tautomer_query_get_substruct_matches_with_tautomers(
            query: &TautomerQuery,
            target: &RWMol,
            params: &SubstructMatchParametersData,
            matches: &mut Vec<IndexList>,
        ) -> Result<UniquePtr<MolList>>;

        fn tautomer_query_pattern_fingerprint_template(
            query: &TautomerQuery,
            fingerprint_size: u32,
        ) -> Result<UniquePtr<BitVect>>;
        fn tautomer_query_get_template_molecule(query: &TautomerQuery) -> Result<UniquePtr<RWMol>>;
        fn tautomer_query_get_modified_atoms(query: &TautomerQuery) -> Result<Vec<u32>>;
        fn tautomer_query_get_modified_bonds(query: &TautomerQuery) -> Result<Vec<u32>>;
        fn tautomer_query_get_tautomers(query: &TautomerQuery) -> Result<UniquePtr<MolList>>;

        fn pattern_fingerprint_tautomer_target(
            target: &RWMol,
            fingerprint_size: u32,
        ) -> Result<UniquePtr<BitVect>>;
        fn tautomer_query_can_serialize() -> Result<bool>;
    }
}

// SAFETY: a TautomerQuery owns its template and tautomers; the shared_ptrs to
// the tautomers are held by no other object. It is not Sync: the pattern
// fingerprint finds the template's rings through a const reference.
unsafe impl Send for ffi::TautomerQuery {}
