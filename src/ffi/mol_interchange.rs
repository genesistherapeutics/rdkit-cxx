// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `MolInterchange::JSONParseParameters`.
    struct JSONParseParametersData {
        set_aromatic_bonds: bool,
        strict_valence_check: bool,
        parse_properties: bool,
        parse_conformers: bool,
        use_h_counts: bool,
    }

    /// `MolInterchange::JSONWriteParameters`.
    struct JSONWriteParametersData {
        use_rdkit_extensions: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_interchange.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;

        fn json_parse_parameters_default() -> Result<JSONParseParametersData>;
        fn json_write_parameters_default() -> Result<JSONWriteParametersData>;
        fn json_to_mols(
            json_block: &str,
            params: &JSONParseParametersData,
        ) -> Result<UniquePtr<MolList>>;
        fn mol_to_json(
            mol: &RWMol,
            params: &JSONWriteParametersData,
        ) -> Result<UniquePtr<CxxString>>;
        fn mols_to_json(
            mols: Pin<&mut BorrowedMols>,
            params: &JSONWriteParametersData,
        ) -> Result<UniquePtr<CxxString>>;
    }
}
