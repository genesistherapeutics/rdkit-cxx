// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::Abbreviations::AbbreviationDefinition`'s fields Python
    /// exposes; a null `mol` for none.
    struct AbbreviationDefinitionData {
        label: UniquePtr<CxxString>,
        display_label: UniquePtr<CxxString>,
        display_label_w: UniquePtr<CxxString>,
        includes_x_bonds: bool,
        mol: UniquePtr<RWMol>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/abbreviations.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type AbbreviationList;

        fn abbreviation_list_new() -> Result<UniquePtr<AbbreviationList>>;
        fn abbreviation_list_push(
            list: Pin<&mut AbbreviationList>,
            label: &str,
            display_label: &str,
            display_label_w: &str,
            includes_x_bonds: bool,
        ) -> Result<()>;
        fn abbreviation_list_push_with_mol(
            list: Pin<&mut AbbreviationList>,
            label: &str,
            display_label: &str,
            display_label_w: &str,
            includes_x_bonds: bool,
            mol: &RWMol,
        ) -> Result<()>;

        fn abbreviation_definition_default() -> Result<AbbreviationDefinitionData>;
        fn get_default_abbreviations() -> Result<Vec<AbbreviationDefinitionData>>;
        fn get_default_linkers() -> Result<Vec<AbbreviationDefinitionData>>;
        fn parse_abbreviations(
            text: &str,
            remove_extra_dummies: bool,
            allow_connection_to_dummies: bool,
        ) -> Result<Vec<AbbreviationDefinitionData>>;

        fn condense_mol_abbreviations(
            mol: &RWMol,
            abbrevs: &AbbreviationList,
            max_coverage: f64,
            sanitize: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn label_mol_abbreviations(
            mol: &RWMol,
            abbrevs: &AbbreviationList,
            max_coverage: f64,
        ) -> Result<UniquePtr<RWMol>>;
        fn condense_abbreviation_substance_groups(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
    }
}
