// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `SmilesParserParams`, with `replacements` as parallel vectors.
    struct SmilesParserParamsData {
        debug_parse: i32,
        parse_name: bool,
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
        sanitize: bool,
        remove_hs: bool,
        replacement_keys: Vec<String>,
        replacement_values: Vec<String>,
    }

    /// `SmartsParserParams`, with `replacements` as parallel vectors.
    struct SmartsParserParamsData {
        debug_parse: i32,
        parse_name: bool,
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
        merge_hs: bool,
        replacement_keys: Vec<String>,
        replacement_values: Vec<String>,
    }

    /// `SmilesWriteParams`.
    struct SmilesWriteParamsData {
        do_isomeric_smiles: bool,
        do_kekule: bool,
        canonical: bool,
        clean_stereo: bool,
        all_bonds_explicit: bool,
        all_hs_explicit: bool,
        do_random: bool,
        rooted_at_atom: i32,
        include_dative_bonds: bool,
        ignore_atom_map_numbers: bool,
    }

    /// `MolWriterParams`, plus the writers' `confId`.
    struct MolWriterParamsData {
        include_stereo: bool,
        kekulize: bool,
        force_v3000: bool,
        precision: u32,
        conf_id: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/molfiles.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn smiles_parser_params_default() -> Result<SmilesParserParamsData>;
        fn smarts_parser_params_default() -> Result<SmartsParserParamsData>;
        fn smiles_write_params_default() -> Result<SmilesWriteParamsData>;
        fn mol_writer_params_default() -> Result<MolWriterParamsData>;

        fn mol_from_smiles(
            smiles: &str,
            params: &SmilesParserParamsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_to_smiles(
            mol: &RWMol,
            params: &SmilesWriteParamsData,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_from_smarts(
            smarts: &str,
            params: &SmartsParserParamsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_to_smarts(
            mol: &RWMol,
            params: &SmilesWriteParamsData,
        ) -> Result<UniquePtr<CxxString>>;

        fn mol_from_mol_block(
            mol_block: &str,
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_mol_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_to_mol_block(
            mol: &RWMol,
            params: &MolWriterParamsData,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_mol_file(mol: &RWMol, path: &[u8], params: &MolWriterParamsData) -> Result<()>;
    }
}
