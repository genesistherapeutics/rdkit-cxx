// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `MolMetadataToPNGString`'s and `MolMetadataToPNGFile`'s flags.
    struct PNGMetadataParamsData {
        include_pkl: bool,
        include_smiles: bool,
        include_mol: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/ffi/molfiles.rs.h");
        include!("rdkit-cxx/png.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type SmilesParserParamsData = crate::ffi::molfiles::ffi::SmilesParserParamsData;
        type PNGMetadata;

        fn mol_from_png_string(
            png: &[u8],
            params: &SmilesParserParamsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_png_file(
            path: &[u8],
            params: &SmilesParserParamsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn mols_from_png_string(
            png: &[u8],
            tag: &str,
            params: &SmilesParserParamsData,
        ) -> Result<UniquePtr<MolList>>;
        fn mols_from_png_file(
            path: &[u8],
            tag: &str,
            params: &SmilesParserParamsData,
        ) -> Result<UniquePtr<MolList>>;
        fn mol_metadata_to_png_string(
            mol: &RWMol,
            png: &[u8],
            params: &PNGMetadataParamsData,
        ) -> Result<Vec<u8>>;
        fn mol_metadata_to_png_file(
            mol: &RWMol,
            path: &[u8],
            params: &PNGMetadataParamsData,
        ) -> Result<Vec<u8>>;
        fn add_metadata_to_png_string(
            keys: &Vec<String>,
            values: &Vec<String>,
            png: &[u8],
        ) -> Result<Vec<u8>>;
        fn add_metadata_to_png_file(
            keys: &Vec<String>,
            values: &Vec<String>,
            path: &[u8],
        ) -> Result<Vec<u8>>;
        fn metadata_from_png_string(png: &[u8]) -> Result<UniquePtr<PNGMetadata>>;
        fn metadata_from_png_file(path: &[u8]) -> Result<UniquePtr<PNGMetadata>>;
        fn png_metadata_size(metadata: &PNGMetadata) -> Result<usize>;
        fn png_metadata_key(metadata: &PNGMetadata, i: usize) -> Result<&CxxString>;
        fn png_metadata_value(metadata: &PNGMetadata, i: usize) -> Result<Vec<u8>>;
    }
}
