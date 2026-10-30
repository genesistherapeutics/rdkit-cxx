// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/cip_labeler.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn assign_cip_labels(
            mol: Pin<&mut RWMol>,
            has_atoms: bool,
            atoms: &[u32],
            has_bonds: bool,
            bonds: &[u32],
            max_recursive_iterations: u32,
        ) -> Result<()>;
    }
}
