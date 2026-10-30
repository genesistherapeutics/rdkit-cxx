// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolOps::SanitizeFlags`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum SanitizeFlags {
        SANITIZE_NONE = 0x0,
        SANITIZE_CLEANUP = 0x1,
        SANITIZE_PROPERTIES = 0x2,
        SANITIZE_SYMMRINGS = 0x4,
        SANITIZE_KEKULIZE = 0x8,
        SANITIZE_FINDRADICALS = 0x10,
        SANITIZE_SETAROMATICITY = 0x20,
        SANITIZE_SETCONJUGATION = 0x40,
        SANITIZE_SETHYBRIDIZATION = 0x80,
        SANITIZE_CLEANUPCHIRALITY = 0x100,
        SANITIZE_ADJUSTHS = 0x200,
        SANITIZE_CLEANUP_ORGANOMETALLICS = 0x400,
        SANITIZE_CLEANUPATROPISOMERS = 0x800,
        SANITIZE_ALL = 0xFFFFFFF,
    }

    /// `AddHs`' arguments, its `onlyOnAtoms` `None` when `has_only_on_atoms`
    /// is false.
    struct AddHsParametersData {
        explicit_only: bool,
        add_coords: bool,
        add_residue_info: bool,
        has_only_on_atoms: bool,
        only_on_atoms: Vec<u32>,
    }

    /// `RemoveHsParameters`, plus `RemoveHs`'s `sanitize`.
    struct RemoveHsParametersData {
        remove_degree_zero: bool,
        remove_higher_degrees: bool,
        remove_only_h_neighbors: bool,
        remove_isotopes: bool,
        remove_and_track_isotopes: bool,
        remove_dummy_neighbors: bool,
        remove_defining_bond_stereo: bool,
        remove_with_wedged_bond: bool,
        remove_with_query: bool,
        remove_mapped: bool,
        remove_in_s_groups: bool,
        show_warnings: bool,
        remove_nonimplicit: bool,
        update_explicit_count: bool,
        remove_hydrides: bool,
        remove_nontetrahedral_neighbors: bool,
        sanitize: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/molops.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type MolList;
        type BorrowedMols;
        type SanitizeFlags;

        fn sanitize_mol(mol: Pin<&mut RWMol>, sanitize_ops: u32, catch_errors: bool)
        -> Result<u32>;
        fn kekulize(mol: Pin<&mut RWMol>, clear_aromatic_flags: bool) -> Result<()>;
        fn add_hs(mol: &RWMol, params: &AddHsParametersData) -> Result<UniquePtr<RWMol>>;
        fn get_mol_frags(mol: &RWMol) -> Result<Vec<IndexList>>;
        fn get_mol_frags_as_mols(mol: &RWMol, sanitize_frags: bool) -> Result<UniquePtr<MolList>>;

        fn remove_hs_parameters_default() -> Result<RemoveHsParametersData>;
        fn remove_hs(mol: &RWMol, params: &RemoveHsParametersData) -> Result<UniquePtr<RWMol>>;
        fn remove_all_hs(mol: &RWMol, sanitize: bool) -> Result<UniquePtr<RWMol>>;
        fn merge_query_hs(
            mol: &RWMol,
            merge_unmapped_only: bool,
            merge_isotopes: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn has_query_hs(mol: &RWMol) -> Result<Vec<bool>>;

        fn mol_list_size(list: &MolList) -> Result<usize>;
        fn mol_list_name(list: &MolList, i: usize) -> Result<&CxxString>;
        fn mol_list_take(list: Pin<&mut MolList>, i: usize) -> Result<UniquePtr<RWMol>>;

        fn borrowed_mols_new() -> Result<UniquePtr<BorrowedMols>>;
        fn borrowed_mols_push(
            list: Pin<&mut BorrowedMols>,
            name: &str,
            mol: Pin<&mut RWMol>,
        ) -> Result<()>;
        fn borrowed_mols_push_shared(
            list: Pin<&mut BorrowedMols>,
            name: &str,
            mol: &RWMol,
        ) -> Result<()>;
    }
}
