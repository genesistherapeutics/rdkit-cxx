// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/sln_parse.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn mol_from_sln(sln: &str, sanitize: bool, debug_parser: bool) -> Result<UniquePtr<RWMol>>;
        fn mol_from_query_sln(
            sln: &str,
            merge_hs: bool,
            debug_parser: bool,
        ) -> Result<UniquePtr<RWMol>>;
    }
}
