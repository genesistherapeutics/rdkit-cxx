// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `SubstructMatchParameters`, without its callback setters.
    struct SubstructMatchParametersData {
        use_chirality: bool,
        use_enhanced_stereo: bool,
        aromatic_matches_conjugated: bool,
        aromatic_matches_single_or_double: bool,
        use_generic_matchers: bool,
        use_query_query_matches: bool,
        recursion_possible: bool,
        uniquify: bool,
        max_matches: u32,
        max_recursive_matches: u32,
        num_threads: i32,
        atom_properties: Vec<String>,
        bond_properties: Vec<String>,
        specified_stereo_query_matches_unspecified: bool,
        extra_atom_check_overrides_default_check: bool,
        extra_bond_check_overrides_default_check: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/substruct.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexList = crate::ffi::geometry::ffi::IndexList;

        fn substruct_match_parameters_default() -> Result<SubstructMatchParametersData>;
        fn mol_has_substruct_match(
            mol: &RWMol,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn mol_get_substruct_match(
            mol: &RWMol,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn mol_get_substruct_matches(
            mol: &RWMol,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
    }
}
