// Copyright (c) 2026 Steven Kearnes All rights reserved.
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

    /// `AddHsParameters`, plus `AddHs`'s `onlyOnAtoms` (`None` when
    /// `has_only_on_atoms` is false).
    struct AddHsParametersData {
        explicit_only: bool,
        add_coords: bool,
        add_residue_info: bool,
        skip_queries: bool,
        has_only_on_atoms: bool,
        only_on_atoms: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/molops.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type MolList;
        type SanitizeFlags;

        fn sanitize_mol(mol: Pin<&mut RWMol>, sanitize_ops: u32, catch_errors: bool)
        -> Result<u32>;
        fn kekulize(
            mol: Pin<&mut RWMol>,
            clear_aromatic_flags: bool,
            canonical: bool,
        ) -> Result<()>;
        fn add_hs_parameters_default() -> Result<AddHsParametersData>;
        fn add_hs(mol: &RWMol, params: &AddHsParametersData) -> Result<UniquePtr<RWMol>>;
        fn get_mol_frags(mol: &RWMol) -> Result<Vec<IndexList>>;
        fn get_mol_frags_as_mols(mol: &RWMol, sanitize_frags: bool) -> Result<UniquePtr<MolList>>;

        fn mol_list_size(list: &MolList) -> Result<usize>;
        fn mol_list_take(list: Pin<&mut MolList>, i: usize) -> Result<UniquePtr<RWMol>>;
    }
}
