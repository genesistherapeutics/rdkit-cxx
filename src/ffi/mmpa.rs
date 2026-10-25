// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mmpa.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;

        fn fragment_mol(
            mol: &RWMol,
            min_cuts: u32,
            max_cuts: u32,
            max_cut_bonds: u32,
            pattern: &str,
        ) -> Result<UniquePtr<MolList>>;
        fn fragment_mol_bonds(
            mol: &RWMol,
            bonds_to_cut: &[u32],
            min_cuts: u32,
            max_cuts: u32,
        ) -> Result<UniquePtr<MolList>>;
    }
}
