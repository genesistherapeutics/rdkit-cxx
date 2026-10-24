// Copyright (c) 2026 Steven Kearnes. All rights reserved.
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
        /// Whether `atom_coords_matcher` is the extra atom check.
        has_atom_coords_matcher: bool,
        atom_coords_matcher: AtomCoordsMatcherData,
    }

    /// `AtomCoordsMatchFunctor`.
    struct AtomCoordsMatcherData {
        ref_conf_id: i32,
        query_conf_id: i32,
        tol2: f64,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/substruct.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolBundle = crate::ffi::mol_bundle::ffi::MolBundle;
        type ResonanceMolSupplier = crate::ffi::resonance::ffi::ResonanceMolSupplier;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type Atom = crate::ffi::atom::ffi::Atom;

        fn substruct_match_parameters_default() -> Result<SubstructMatchParametersData>;
        fn atom_coords_matcher_default() -> Result<AtomCoordsMatcherData>;
        fn atom_coords_matcher_new(
            ref_conf_id: i32,
            query_conf_id: i32,
            tol: f64,
        ) -> Result<AtomCoordsMatcherData>;
        fn atom_coords_matcher_call(
            matcher: &AtomCoordsMatcherData,
            query_atom: &Atom,
            target_atom: &Atom,
        ) -> Result<bool>;
        fn mol_has_substruct_match(
            target: &RWMol,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn mol_get_substruct_match(
            target: &RWMol,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn mol_get_substruct_matches(
            target: &RWMol,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
        fn mol_has_substruct_match_mol_bundle(
            target: &RWMol,
            query: &MolBundle,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn mol_get_substruct_match_mol_bundle(
            target: &RWMol,
            query: &MolBundle,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn mol_get_substruct_matches_mol_bundle(
            target: &RWMol,
            query: &MolBundle,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
        fn mol_bundle_has_substruct_match(
            target: &MolBundle,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn mol_bundle_get_substruct_match(
            target: &MolBundle,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn mol_bundle_get_substruct_matches(
            target: &MolBundle,
            query: &RWMol,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
        fn mol_bundle_has_substruct_match_mol_bundle(
            target: &MolBundle,
            query: &MolBundle,
            params: &SubstructMatchParametersData,
        ) -> Result<bool>;
        fn mol_bundle_get_substruct_match_mol_bundle(
            target: &MolBundle,
            query: &MolBundle,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<u32>>;
        fn mol_bundle_get_substruct_matches_mol_bundle(
            target: &MolBundle,
            query: &MolBundle,
            params: &SubstructMatchParametersData,
        ) -> Result<Vec<IndexList>>;
        fn resonance_mol_supplier_get_substruct_match(
            suppl: Pin<&mut ResonanceMolSupplier>,
            query: &RWMol,
            use_chirality: bool,
            use_query_query_matches: bool,
        ) -> Result<Vec<u32>>;
        #[allow(clippy::too_many_arguments)]
        fn resonance_mol_supplier_get_substruct_matches(
            suppl: Pin<&mut ResonanceMolSupplier>,
            query: &RWMol,
            uniquify: bool,
            use_chirality: bool,
            use_query_query_matches: bool,
            max_matches: u32,
            num_threads: i32,
        ) -> Result<Vec<IndexList>>;
    }
}
