// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/deprotect.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type DeprotectData;
        type DeprotectList;

        fn deprotect_data_new(
            deprotection_class: &str,
            reaction_smarts: &str,
            abbreviation: &str,
            full_name: &str,
        ) -> Result<UniquePtr<DeprotectData>>;
        fn deprotect_data_clone(data: &DeprotectData) -> Result<UniquePtr<DeprotectData>>;
        fn deprotect_data_deprotection_class(data: &DeprotectData) -> Result<UniquePtr<CxxString>>;
        fn deprotect_data_full_name(data: &DeprotectData) -> Result<UniquePtr<CxxString>>;
        fn deprotect_data_abbreviation(data: &DeprotectData) -> Result<UniquePtr<CxxString>>;
        fn deprotect_data_reaction_smarts(data: &DeprotectData) -> Result<UniquePtr<CxxString>>;
        fn deprotect_data_example(data: &DeprotectData) -> Result<UniquePtr<CxxString>>;
        fn deprotect_data_is_valid(data: &DeprotectData) -> Result<bool>;

        fn num_default_deprotections() -> Result<usize>;
        fn default_deprotection(i: usize) -> Result<UniquePtr<DeprotectData>>;

        fn deprotect_list_new() -> Result<UniquePtr<DeprotectList>>;
        fn deprotect_list_push(list: Pin<&mut DeprotectList>, data: &DeprotectData) -> Result<()>;

        fn deprotect_with(mol: &RWMol, list: &DeprotectList) -> Result<UniquePtr<RWMol>>;
        fn deprotect_default(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn deprotect_in_place_with(mol: Pin<&mut RWMol>, list: &DeprotectList) -> Result<bool>;
        fn deprotect_in_place_default(mol: Pin<&mut RWMol>) -> Result<bool>;
    }
}

// SAFETY: a DeprotectData owns its strings and shares its reaction through a
// shared_ptr (atomic counts). Copies on other threads only run the reaction,
// which RDKit's const runReactant reads once it is initialized.
unsafe impl Send for ffi::DeprotectData {}
