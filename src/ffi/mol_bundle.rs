// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_bundle.h");

        type MolBundle;
        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;

        fn mol_bundle_new() -> Result<UniquePtr<MolBundle>>;
        fn mol_bundle_add_borrowed(bundle: Pin<&mut MolBundle>, mol: &RWMol) -> Result<()>;
        fn fixed_mol_size_mol_bundle_check(first: &RWMol, nmol: &RWMol) -> Result<()>;
        fn mol_bundle_from_pkl_string(pkl: &[u8]) -> Result<UniquePtr<MolList>>;
        fn mol_bundle_to_binary(bundle: &MolBundle) -> Result<Vec<u8>>;
        fn mol_bundle_can_serialize() -> Result<bool>;
    }
}
