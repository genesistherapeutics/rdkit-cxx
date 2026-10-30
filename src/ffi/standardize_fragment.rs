// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize_fragment.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type FragmentRemover;

        fn fragment_remover_default() -> Result<UniquePtr<FragmentRemover>>;
        fn fragment_remover_new(
            fragment_filename: &[u8],
            leave_last: bool,
            skip_if_all_match: bool,
        ) -> Result<UniquePtr<FragmentRemover>>;
        fn fragment_remover_from_data(
            fragment_data: &str,
            leave_last: bool,
            skip_if_all_match: bool,
        ) -> Result<UniquePtr<FragmentRemover>>;
        fn fragment_remover_remove(
            remover: &FragmentRemover,
            mol: &RWMol,
        ) -> Result<UniquePtr<RWMol>>;
        fn fragment_remover_remove_in_place(
            remover: &FragmentRemover,
            mol: Pin<&mut RWMol>,
        ) -> Result<()>;

        fn largest_fragment_chooser_choose(
            mol: &RWMol,
            prefer_organic: bool,
            use_atom_count: bool,
            count_heavy_atoms_only: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn largest_fragment_chooser_choose_in_place(
            mol: Pin<&mut RWMol>,
            prefer_organic: bool,
            use_atom_count: bool,
            count_heavy_atoms_only: bool,
        ) -> Result<()>;
    }
}

// SAFETY: a FragmentRemover owns its catalog, a deep copy of RDKit's
// definitions, and nothing refers to it.
unsafe impl Send for ffi::FragmentRemover {}
