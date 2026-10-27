// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
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

    /// `RDKit::SmilesWrite::CXSmilesFields`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum CXSmilesFields {
        CX_NONE = 0,
        CX_ATOM_LABELS = 0x1,
        CX_MOLFILE_VALUES = 0x2,
        CX_COORDS = 0x4,
        CX_RADICALS = 0x8,
        CX_ATOM_PROPS = 0x10,
        CX_LINKNODES = 0x20,
        CX_ENHANCEDSTEREO = 0x40,
        CX_SGROUPS = 0x80,
        CX_POLYMER = 0x100,
        CX_BOND_CFG = 0x200,
        CX_BOND_ATROPISOMER = 0x400,
        CX_COORDINATE_BONDS = 0x800,
        CX_ZERO_BONDS = 0x2000,
        CX_ALL = 0x7fff_ffff,
        CX_ALL_BUT_COORDS = 0x7fff_fffb,
    }

    /// `RDKit::RestoreBondDirOption`; cxx checks the values.
    #[repr(u32)]
    enum RestoreBondDirOption {
        RestoreBondDirOptionTrue = 0,
        RestoreBondDirOptionClear = 1,
    }

    /// `RDKit::v2::FileParsers::SCSRTemplateNames`; cxx checks the values.
    #[repr(i32)]
    enum SCSRTemplateNames {
        AsEntered,
        UseFirstName,
        UseSecondName,
    }

    /// `RDKit::v2::FileParsers::SCSRBaseHbondOptions`; cxx checks the values.
    #[repr(i32)]
    enum SCSRBaseHbondOptions {
        Ignore,
        UseSapAll,
        UseSapOne,
        Auto,
    }

    /// `RDKit::v2::CDXMLParser::CDXMLFormat`; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum CDXMLFormat {
        CDXML = 0,
        CDX = 1,
        Auto = 2,
    }

    /// `RDKit::CDXFormat`, ChemDraw's writer's format; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum CDXFormat {
        CDX = 1,
        CDXML = 2,
    }

    /// `MolFromSCSRParams`.
    struct MolFromSCSRParamsData {
        include_leaving_groups: bool,
        scsr_template_names: SCSRTemplateNames,
        scsr_base_hbond_options: SCSRBaseHbondOptions,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/molfiles.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;
        type Bond = crate::ffi::bond::ffi::Bond;
        type CXSmilesFields;
        type RestoreBondDirOption;
        type SCSRTemplateNames;
        type SCSRBaseHbondOptions;
        type CDXMLFormat;
        type CDXFormat;

        fn smiles_parser_params_default() -> Result<SmilesParserParamsData>;
        fn smarts_parser_params_default() -> Result<SmartsParserParamsData>;
        fn smiles_write_params_default() -> Result<SmilesWriteParamsData>;
        fn mol_writer_params_default() -> Result<MolWriterParamsData>;
        fn mol_from_scsr_params_default() -> Result<MolFromSCSRParamsData>;

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
        fn atom_from_smiles(smiles: &str) -> Result<UniquePtr<Atom>>;
        fn atom_from_smarts(smarts: &str) -> Result<UniquePtr<Atom>>;
        fn bond_from_smiles(smiles: &str) -> Result<UniquePtr<Bond>>;
        fn bond_from_smarts(smarts: &str) -> Result<UniquePtr<Bond>>;
        fn mol_to_cx_smiles(
            mol: &RWMol,
            params: &SmilesWriteParamsData,
            flags: u32,
            restore_bond_dirs: RestoreBondDirOption,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_cx_smarts(mol: &RWMol, isomeric_smiles: bool) -> Result<UniquePtr<CxxString>>;
        #[allow(clippy::too_many_arguments)]
        fn mol_fragment_to_smiles(
            mol: &RWMol,
            params: &SmilesWriteParamsData,
            atoms_to_use: &[u32],
            has_bonds: bool,
            bonds_to_use: &[u32],
            has_atom_symbols: bool,
            atom_symbols: &Vec<String>,
            has_bond_symbols: bool,
            bond_symbols: &Vec<String>,
            cx: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_fragment_to_smarts(
            mol: &RWMol,
            atoms_to_use: &[u32],
            has_bonds: bool,
            bonds_to_use: &[u32],
            isomeric_smarts: bool,
            cx: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_random_smiles_vect(
            mol: &RWMol,
            num_smiles: u32,
            random_seed: u32,
            isomeric_smiles: bool,
            kekule_smiles: bool,
            all_bonds_explicit: bool,
            all_hs_explicit: bool,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn canonical_rank_atoms(
            mol: &RWMol,
            break_ties: bool,
            include_chirality: bool,
            include_isotopes: bool,
            include_atom_maps: bool,
            include_chiral_presence: bool,
        ) -> Result<Vec<u32>>;
        #[allow(clippy::too_many_arguments)]
        fn canonical_rank_atoms_in_fragment(
            mol: &RWMol,
            atoms_to_use: &[u32],
            has_bonds: bool,
            bonds_to_use: &[u32],
            has_atom_symbols: bool,
            atom_symbols: &Vec<String>,
            break_ties: bool,
            include_chirality: bool,
            include_isotopes: bool,
            include_atom_maps: bool,
            include_chiral_presence: bool,
        ) -> Result<Vec<i32>>;
        fn canonicalize_enhanced_stereo(mol: Pin<&mut RWMol>) -> Result<()>;
        /// `kind` is 0 for int, 1 for double, 2 for bool, 3 for string.
        fn create_property_list(
            mol: Pin<&mut RWMol>,
            kind: u8,
            prop_name: &str,
            missing_value_marker: &str,
            line_size: u32,
        ) -> Result<()>;

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
        fn mol_to_v2k_mol_block(
            mol: &RWMol,
            params: &MolWriterParamsData,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_v3k_mol_block(
            mol: &RWMol,
            params: &MolWriterParamsData,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_v3k_mol_file(
            mol: &RWMol,
            path: &[u8],
            params: &MolWriterParamsData,
        ) -> Result<()>;
        fn mol_from_scsr_block(
            mol_block: &str,
            sanitize: bool,
            remove_hs: bool,
            params: &MolFromSCSRParamsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_scsr_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            params: &MolFromSCSRParamsData,
        ) -> Result<UniquePtr<RWMol>>;

        fn mol_from_tpl_block(
            tpl_block: &str,
            sanitize: bool,
            skip_first_conf: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_tpl_file(
            path: &[u8],
            sanitize: bool,
            skip_first_conf: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_to_tpl_block(
            mol: &RWMol,
            partial_charge_prop: &str,
            write_first_conf_twice: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_tpl_file(
            mol: &RWMol,
            path: &[u8],
            partial_charge_prop: &str,
            write_first_conf_twice: bool,
        ) -> Result<()>;
        fn mol_from_mol2_block(
            mol2_block: &str,
            sanitize: bool,
            remove_hs: bool,
            cleanup_substructures: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_mol2_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            cleanup_substructures: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_mrv_block(
            mrv_block: &str,
            sanitize: bool,
            remove_hs: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_mrv_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_to_mrv_block(
            mol: &RWMol,
            include_stereo: bool,
            conf_id: i32,
            kekulize: bool,
            pretty_print: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_mrv_file(
            mol: &RWMol,
            path: &[u8],
            include_stereo: bool,
            conf_id: i32,
            kekulize: bool,
            pretty_print: bool,
        ) -> Result<()>;
        fn mol_from_xyz_block(xyz_block: &str) -> Result<UniquePtr<RWMol>>;
        fn mol_from_xyz_file(path: &[u8]) -> Result<UniquePtr<RWMol>>;
        fn mol_to_xyz_block(
            mol: &RWMol,
            conf_id: i32,
            precision: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_xyz_file(mol: &RWMol, path: &[u8], conf_id: i32, precision: u32) -> Result<()>;
        fn mol_from_pdb_block(
            mol_block: &str,
            sanitize: bool,
            remove_hs: bool,
            flavor: u32,
            proximity_bonding: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_pdb_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            flavor: u32,
            proximity_bonding: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_to_pdb_block(mol: &RWMol, conf_id: i32, flavor: u32)
        -> Result<UniquePtr<CxxString>>;
        fn mol_to_pdb_file(mol: &RWMol, path: &[u8], conf_id: i32, flavor: u32) -> Result<()>;
        fn mol_to_cml_block(
            mol: &RWMol,
            conf_id: i32,
            kekulize: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_to_cml_file(mol: &RWMol, path: &[u8], conf_id: i32, kekulize: bool) -> Result<()>;
        fn mol_from_rdkit_svg(
            svg: &str,
            sanitize: bool,
            remove_hs: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn mol_from_sequence(text: &str, sanitize: bool, flavor: i32) -> Result<UniquePtr<RWMol>>;
        fn mol_from_fasta(text: &str, sanitize: bool, flavor: i32) -> Result<UniquePtr<RWMol>>;
        fn mol_from_helm(text: &str, sanitize: bool) -> Result<UniquePtr<RWMol>>;
        fn mol_to_sequence(mol: &RWMol) -> Result<UniquePtr<CxxString>>;
        fn mol_to_fasta(mol: &RWMol) -> Result<UniquePtr<CxxString>>;
        fn mol_to_helm(mol: &RWMol) -> Result<UniquePtr<CxxString>>;
        fn has_chem_draw_cdx_support() -> Result<bool>;
        fn mol_to_chem_draw_block(mol: &RWMol, format: CDXFormat) -> Result<Vec<u8>>;
    }
}
