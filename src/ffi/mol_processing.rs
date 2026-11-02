// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `GeneralMolSupplier::SupplierOptions`, as Python exposes it.
    struct SupplierOptionsData {
        num_threads: i32,
        sanitize: bool,
        remove_hs: bool,
        strict_parsing: bool,
        delimiter: String,
        smiles_column: i32,
        name_column: i32,
        title_line: bool,
        name_record: String,
        conf_id_2d: i32,
        conf_id_3d: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_processing.h");

        type FpGenerator = crate::ffi::fingerprint_generator::ffi::FpGenerator;
        type FingerprintList = crate::ffi::fingerprint_generator::ffi::FingerprintList;

        fn supplier_options_default() -> Result<SupplierOptionsData>;
        fn get_fingerprints_for_mols_in_file(
            path: &[u8],
            options: &SupplierOptionsData,
        ) -> Result<UniquePtr<FingerprintList>>;
        fn get_fingerprints_for_mols_in_file_with_generator(
            path: &[u8],
            generator: &FpGenerator,
            options: &SupplierOptionsData,
        ) -> Result<UniquePtr<FingerprintList>>;
    }
}
