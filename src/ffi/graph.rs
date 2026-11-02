// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::SubsetMethod`; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum SubsetMethod {
        BONDS_BETWEEN_ATOMS = 0,
        BONDS = 1,
    }

    /// `RDKit::SubsetOptions`.
    struct SubsetOptionsData {
        sanitize: bool,
        clear_computed_props: bool,
        copy_as_query: bool,
        copy_coordinates: bool,
        conformer_idx: i32,
        method: SubsetMethod,
    }

    /// `RDKit::SubsetInfo`'s mappings, as (from, to) pairs.
    struct SubsetInfoData {
        atom_mapping: Vec<IndexPair>,
        bond_mapping: Vec<IndexPair>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/graph.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type SubsetMethod;

        fn get_distance_matrix(mol: &RWMol, use_bo: bool, use_atom_wts: bool) -> Result<Vec<f64>>;
        fn get_3d_distance_matrix(
            mol: &RWMol,
            conf_id: i32,
            use_atom_wts: bool,
        ) -> Result<Vec<f64>>;
        fn get_adjacency_matrix(mol: &RWMol, use_bo: bool, empty_val: i32) -> Result<Vec<f64>>;

        fn find_all_paths_of_length_n(
            mol: &RWMol,
            length: u32,
            use_bonds: bool,
            use_hs: bool,
            rooted_at_atom: i32,
            only_shortest_paths: bool,
        ) -> Result<Vec<IndexList>>;
        fn find_all_subgraphs_of_length_n(
            mol: &RWMol,
            length: u32,
            use_hs: bool,
            rooted_at_atom: i32,
        ) -> Result<Vec<IndexList>>;
        fn find_all_subgraphs_of_length_m_to_n(
            mol: &RWMol,
            min: u32,
            max: u32,
            use_hs: bool,
            rooted_at_atom: i32,
            counts: &mut Vec<u32>,
        ) -> Result<Vec<IndexList>>;
        fn find_unique_subgraphs_of_length_n(
            mol: &RWMol,
            length: u32,
            use_hs: bool,
            use_bo: bool,
            rooted_at_atom: i32,
        ) -> Result<Vec<IndexList>>;
        fn find_atom_environment_of_radius_n(
            mol: &RWMol,
            radius: u32,
            rooted_at_atom: u32,
            use_hs: bool,
            enforce_size: bool,
            atom_map: &mut Vec<IndexPair>,
        ) -> Result<Vec<u32>>;
        fn path_to_submol(
            mol: &RWMol,
            path: &[u32],
            use_query: bool,
            atom_map: &mut Vec<IndexPair>,
        ) -> Result<UniquePtr<RWMol>>;

        fn split_mol_by_pdb_residues(
            mol: &RWMol,
            has_white_list: bool,
            white_list: &[String],
            negate_list: bool,
        ) -> Result<UniquePtr<MolList>>;
        fn split_mol_by_pdb_chain_id(
            mol: &RWMol,
            has_white_list: bool,
            white_list: &[String],
            negate_list: bool,
        ) -> Result<UniquePtr<MolList>>;

        fn renumber_atoms(mol: &RWMol, new_order: &[u32]) -> Result<UniquePtr<RWMol>>;

        fn subset_options_default() -> Result<SubsetOptionsData>;
        fn copy_mol_subset(
            mol: &RWMol,
            path: &[u32],
            options: &SubsetOptionsData,
            info: &mut SubsetInfoData,
        ) -> Result<UniquePtr<RWMol>>;
        fn copy_mol_subset_atoms_bonds(
            mol: &RWMol,
            atom_indices: &[u32],
            bond_indices: &[u32],
            options: &SubsetOptionsData,
            info: &mut SubsetInfoData,
        ) -> Result<UniquePtr<RWMol>>;
    }
}
