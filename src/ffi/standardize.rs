// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// An entry of one of `CleanupParameters`' data lists, its strings as
    /// bytes: a normalization or fragment (name, SMARTS), an acid-base pair
    /// (name, acid, base), or a tautomer transform (name, SMARTS, bonds,
    /// charges). Unused fields are empty.
    #[derive(Clone, Debug, Default)]
    struct CleanupDataEntry {
        name: Vec<u8>,
        first: Vec<u8>,
        second: Vec<u8>,
        third: Vec<u8>,
    }

    /// `RDKit::MolStandardize::CleanupParameters` (A7), its strings as bytes.
    #[derive(Clone, Debug, Default)]
    struct CleanupParametersData {
        rdbase: Vec<u8>,
        normalizations: Vec<u8>,
        acidbase_file: Vec<u8>,
        fragment_file: Vec<u8>,
        tautomer_transforms: Vec<u8>,
        max_restarts: i32,
        prefer_organic: bool,
        do_canonical: bool,
        max_tautomers: i32,
        max_transforms: i32,
        tautomer_remove_sp3_stereo: bool,
        tautomer_remove_bond_stereo: bool,
        tautomer_remove_isotopic_hs: bool,
        tautomer_reassign_stereo: bool,
        largest_fragment_chooser_use_atom_count: bool,
        largest_fragment_chooser_count_heavy_atoms_only: bool,
        normalization_data: Vec<CleanupDataEntry>,
        fragment_data: Vec<CleanupDataEntry>,
        acidbase_data: Vec<CleanupDataEntry>,
        tautomer_transform_data: Vec<CleanupDataEntry>,
    }

    /// Which module function a call stands for.
    #[repr(u8)]
    enum StandardizeOp {
        Cleanup,
        Normalize,
        Reionize,
        RemoveFragments,
        CanonicalTautomer,
        TautomerParent,
        FragmentParent,
        StereoParent,
        IsotopeParent,
        ChargeParent,
        SuperParent,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;

        fn cleanup_parameters_default() -> Result<CleanupParametersData>;
        fn update_params_from_json(
            data: &CleanupParametersData,
            json: &str,
        ) -> Result<CleanupParametersData>;

        fn standardize_copy(
            op: StandardizeOp,
            mol: &RWMol,
            data: &CleanupParametersData,
            skip_standardize: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn standardize_in_place(
            op: StandardizeOp,
            mol: Pin<&mut RWMol>,
            data: &CleanupParametersData,
            skip_standardize: bool,
        ) -> Result<()>;
        fn standardize_in_place_mols(
            op: StandardizeOp,
            mols: Pin<&mut BorrowedMols>,
            num_threads: i32,
            data: &CleanupParametersData,
            skip_standardize: bool,
        ) -> Result<()>;

        fn standardize_smiles(smiles: &str) -> Result<UniquePtr<CxxString>>;
    }
}
