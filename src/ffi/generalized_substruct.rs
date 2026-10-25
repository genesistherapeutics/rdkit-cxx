// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/generalized_substruct.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");
        include!("rdkit-cxx/ffi/substruct.rs.h");
        include!("rdkit-cxx/ffi/adjust_query.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type SubstructMatchParametersData =
            crate::ffi::substruct::ffi::SubstructMatchParametersData;
        type AdjustQueryParametersData = crate::ffi::adjust_query::ffi::AdjustQueryParametersData;
        type ExtendedQueryMol;

        fn extended_query_mol_from_binary(pickle: &[u8]) -> Result<UniquePtr<ExtendedQueryMol>>;
        fn extended_query_mol_from_json(text: &str) -> Result<UniquePtr<ExtendedQueryMol>>;
        fn extended_query_mol_clone(
            query: &ExtendedQueryMol,
        ) -> Result<UniquePtr<ExtendedQueryMol>>;
        fn extended_query_mol_to_binary(query: &ExtendedQueryMol) -> Result<Vec<u8>>;
        fn extended_query_mol_to_json(query: &ExtendedQueryMol) -> Result<UniquePtr<CxxString>>;
        fn extended_query_mol_pattern_fingerprint_query(
            query: &ExtendedQueryMol,
            fingerprint_size: u32,
        ) -> Result<UniquePtr<BitVect>>;

        fn create_extended_query_mol(
            mol: &RWMol,
            do_enumeration: bool,
            do_tautomers: bool,
            adjust_query_properties: bool,
            adjust_query_parameters: &AdjustQueryParametersData,
        ) -> Result<UniquePtr<ExtendedQueryMol>>;

        fn mol_has_substruct_match_extended(
            mol: &RWMol,
            query: &ExtendedQueryMol,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn mol_get_substruct_match_extended(
            mol: &RWMol,
            query: &ExtendedQueryMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn mol_get_substruct_matches_extended(
            mol: &RWMol,
            query: &ExtendedQueryMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
        fn pattern_fingerprint_target(
            target: &RWMol,
            fingerprint_size: u32,
        ) -> Result<UniquePtr<BitVect>>;
    }
}

// SAFETY: an ExtendedQueryMol owns its molecules: those it reads, makes, or
// copies (extended_query_mol_clone copies a bundle's) are held by no other
// object. It is not Sync: its pattern fingerprint finds their rings through
// a const reference.
unsafe impl Send for ffi::ExtendedQueryMol {}
