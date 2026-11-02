// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/partial_charges.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn compute_gasteiger_charges(
            mol: Pin<&mut RWMol>,
            n_iter: i32,
            throw_on_param_failure: bool,
        ) -> Result<()>;
    }
}
