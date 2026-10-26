// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `InchiToMol`'s molecule and `ExtraInchiReturnValues`.
    struct InchiToMolData {
        /// Null when InChI cannot read the string.
        mol: UniquePtr<RWMol>,
        retcode: i32,
        message: UniquePtr<CxxString>,
        log: UniquePtr<CxxString>,
    }

    /// `MolToInchi`'s InChI and `ExtraInchiReturnValues`.
    struct MolToInchiData {
        inchi: UniquePtr<CxxString>,
        retcode: i32,
        message: UniquePtr<CxxString>,
        log: UniquePtr<CxxString>,
        aux_info: UniquePtr<CxxString>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/inchi.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn inchi_to_mol(inchi: &str, sanitize: bool, remove_hs: bool) -> Result<InchiToMolData>;
        fn mol_to_inchi(mol: &RWMol, options: &str) -> Result<MolToInchiData>;
        fn mol_to_inchi_key(mol: &RWMol, options: &str) -> Result<UniquePtr<CxxString>>;
        fn inchi_to_inchi_key(inchi: &str) -> Result<UniquePtr<CxxString>>;
        fn get_inchi_version() -> Result<UniquePtr<CxxString>>;
    }
}
