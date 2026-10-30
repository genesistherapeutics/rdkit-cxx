// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolzipLabel`; cxx checks the values.
    #[repr(u32)]
    enum MolzipLabel {
        AtomMapNumber,
        Isotope,
        FragmentOnBonds,
        AtomType,
    }

    /// `RDKit::MolzipParams`.
    struct MolzipParamsData {
        label: MolzipLabel,
        atom_symbols: Vec<String>,
        enforce_valence_rules: bool,
        generate_coordinates: bool,
        align_coordinates: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/transforms.h");
        include!("rdkit-cxx/ffi/bond.rs.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BondType = crate::ffi::bond::ffi::BondType;
        type MolList = crate::ffi::molops::ffi::MolList;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type MolzipLabel;

        fn delete_substructs(
            mol: &RWMol,
            query: &RWMol,
            only_frags: bool,
            use_chirality: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn replace_substructs(
            mol: &RWMol,
            query: &RWMol,
            replacement: &RWMol,
            replace_all: bool,
            replacement_connection_point: u32,
            use_chirality: bool,
        ) -> Result<UniquePtr<MolList>>;
        fn replace_sidechains(
            mol: &RWMol,
            core_query: &RWMol,
            use_chirality: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn replace_core_with_matches(
            mol: &RWMol,
            core: &RWMol,
            matches: &[IndexPair],
            replace_dummies: bool,
            label_by_index: bool,
            require_dummy_match: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn replace_core(
            mol: &RWMol,
            core_query: &RWMol,
            replace_dummies: bool,
            label_by_index: bool,
            require_dummy_match: bool,
            use_chirality: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn murcko_decompose(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn combine_mols(mol1: &RWMol, mol2: &RWMol, offset: &Point3D) -> Result<UniquePtr<RWMol>>;
        fn get_most_substituted_core_match(
            mol: &RWMol,
            core: &RWMol,
            matches: &[IndexList],
        ) -> Result<Vec<u32>>;
        fn sort_matches_by_degree_of_core_substitution(
            mol: &RWMol,
            core: &RWMol,
            matches: &[IndexList],
        ) -> Result<Vec<IndexList>>;

        fn add_recursive_query(
            mol: Pin<&mut RWMol>,
            query: &RWMol,
            atom_idx: u32,
            preserve_existing_query: bool,
        ) -> Result<()>;
        fn add_recursive_queries(
            mol: Pin<&mut RWMol>,
            queries: Pin<&mut BorrowedMols>,
            prop_name: &str,
        ) -> Result<()>;
        fn parse_mol_query_def_file(
            path: &[u8],
            standardize: bool,
            delimiter: &str,
            comment: &str,
            name_column: u32,
            smarts_column: u32,
        ) -> Result<UniquePtr<MolList>>;
        fn parse_mol_query_def_text(
            text: &str,
            standardize: bool,
            delimiter: &str,
            comment: &str,
            name_column: u32,
            smarts_column: u32,
        ) -> Result<UniquePtr<MolList>>;

        #[allow(clippy::too_many_arguments)]
        fn fragment_on_bonds(
            mol: &RWMol,
            bond_indices: &[u32],
            add_dummies: bool,
            has_dummy_labels: bool,
            dummy_labels: &[IndexPair],
            has_bond_types: bool,
            bond_types: &[BondType],
            cuts_per_atom: &mut Vec<u32>,
        ) -> Result<UniquePtr<RWMol>>;
        #[allow(clippy::too_many_arguments)]
        fn fragment_on_some_bonds(
            mol: &RWMol,
            bond_indices: &[u32],
            num_to_break: u32,
            add_dummies: bool,
            has_dummy_labels: bool,
            dummy_labels: &[IndexPair],
            has_bond_types: bool,
            bond_types: &[BondType],
            return_cuts_per_atom: bool,
            cuts_per_atom: &mut Vec<IndexList>,
        ) -> Result<UniquePtr<MolList>>;
        fn fragment_on_brics_bonds(mol: &RWMol) -> Result<UniquePtr<RWMol>>;

        fn molzip_params_default() -> Result<MolzipParamsData>;
        fn molzip(a: &RWMol, b: &RWMol, params: &MolzipParamsData) -> Result<UniquePtr<RWMol>>;
        fn molzip_single(a: &RWMol, params: &MolzipParamsData) -> Result<UniquePtr<RWMol>>;
        fn molzip_fragments(
            mols: Pin<&mut BorrowedMols>,
            params: &MolzipParamsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn molzip_row(
            row: Pin<&mut BorrowedMols>,
            params: &MolzipParamsData,
        ) -> Result<UniquePtr<RWMol>>;
    }
}
