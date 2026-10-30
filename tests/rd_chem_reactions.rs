// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdChemReactions: the expected values come from RDKit's Python bindings.

mod common;

use std::collections::BTreeMap;
use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_reactions::*;
use rdkit_cxx::chem::{
    AtomMonomerInfo, AtomMonomerType, ChiralType, Mol, PropValue, PropertyPickleOptions,
    SmilesWriteParams, mol_from_smarts, mol_from_smiles, mol_to_smiles,
};

const AMIDE: &str = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]";

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn product_smiles(products: &[Vec<Mol>]) -> Vec<Vec<String>> {
    products
        .iter()
        .map(|set| set.iter().map(|mol| mol_to_smiles(mol).unwrap()).collect())
        .collect()
}

fn amide() -> ChemicalReaction {
    reaction_from_smarts(AMIDE).unwrap()
}

#[test]
fn runs_a_reaction() {
    let mut rxn = amide();
    assert_eq!(rxn.get_num_reactant_templates(), 2);
    assert_eq!(rxn.get_num_product_templates(), 1);
    assert_eq!(rxn.get_num_agent_templates(), 0);
    assert!(!rxn.is_initialized());
    let reactants = [smiles("C(=O)O"), smiles("CNC")];
    let products = rxn.run_reactants(&reactants).unwrap();
    assert_eq!(product_smiles(&products), [["CN(C)C=O"]]);
    // Running initializes the reaction, as in Python.
    assert!(rxn.is_initialized());
    assert_eq!(mol_to_smiles(&reactants[1]).unwrap(), "CNC");
}

#[test]
fn run_reactants_stops_at_max_products() {
    let mut rxn = reaction_from_smarts("[O:1]>>[N:1]").unwrap();
    let reactants = [smiles("OCO")];
    assert_eq!(
        product_smiles(&rxn.run_reactants(&reactants).unwrap()),
        [["NCO"], ["NCO"]]
    );
    assert_eq!(
        product_smiles(&rxn.run_reactants_with_max_products(&reactants, 1).unwrap()),
        [["NCO"]]
    );
}

#[test]
fn run_reactants_finds_every_match() {
    let mut rxn = reaction_from_smarts(
        "[C:1]=[C:2].[C:3]=[*:4][C:5]=[C:6]>>[C:1]1[C:2][C:3][*:4]=[C:5][C:6]1",
    )
    .unwrap();
    let products = rxn
        .run_reactants(&[smiles("C=C"), smiles("C=CC=C")])
        .unwrap();
    assert_eq!(products.len(), 4);
    for set in product_smiles(&products) {
        assert_eq!(set, ["C1=CCCCC1"]);
    }
}

#[test]
fn run_reactants_checks_the_reactant_count() {
    let err = amide().run_reactants(&[smiles("C")]).unwrap_err();
    assert!(
        matches!(&err, Error::Value { message } if message.contains("does not match")),
        "{err:?}"
    );
}

#[test]
fn an_invalid_reaction_does_not_run() {
    let mut rxn = reaction_from_smarts("C>>").unwrap();
    assert_eq!(rxn.validate_with_silent(true).unwrap(), (1, 1));
    let err = rxn.run_reactants(&[smiles("C")]).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    assert!(!rxn.is_initialized());
}

#[test]
fn validate_counts_problems() {
    assert_eq!(amide().validate().unwrap(), (0, 0));
    // A product map number repeated, and not in the reactants, trips RDKit's
    // invariant.
    let mut rxn = reaction_from_smarts("[C:1]>>[C:1].[N:2].[N:2]").unwrap();
    assert!(matches!(
        rxn.validate_with_silent(true),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn run_reactant_uses_one_template() {
    let mut rxn = amide();
    let products = rxn.run_reactant(&smiles("CC(=O)O"), 0).unwrap();
    assert_eq!(product_smiles(&products), [["CC(N)=O"]]);
    assert_eq!(
        rxn.run_reactant(&smiles("CC(=O)O"), 5).unwrap_err(),
        Error::Index { index: 5, len: 2 }
    );
}

#[test]
fn run_reactant_in_place() {
    let mut rxn = reaction_from_smarts("[C:1]-[Cl:2]>>[C:1]-[Br:2]").unwrap();
    let mut mol = smiles("CCCl");
    assert!(rxn.run_reactant_in_place(&mut mol).unwrap());
    assert_eq!(mol_to_smiles(&mol).unwrap(), "CCBr");
    let mut unmatched = smiles("CCO");
    assert!(!rxn.run_reactant_in_place(&mut unmatched).unwrap());
    // Only one-reactant, one-product reactions run in place.
    let err = amide().run_reactant_in_place(&mut mol).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}

#[test]
fn template_views() {
    let rxn = amide();
    let template = rxn.get_reactant_template(0).unwrap();
    assert_eq!(template.get_num_atoms(), 3);
    assert_eq!(rxn.get_product_template(0).unwrap().get_num_atoms(), 3);
    assert_eq!(
        rxn.get_reactant_template(2).unwrap_err(),
        Error::Index { index: 2, len: 2 }
    );
    assert_eq!(
        rxn.get_agent_template(0).unwrap_err(),
        Error::Index { index: 0, len: 0 }
    );
    let sizes: Vec<u32> = rxn.get_reactants().map(|t| t.get_num_atoms()).collect();
    assert_eq!(sizes, [3, 1]);
    assert_eq!(rxn.get_reactants().rev().len(), 2);
    assert_eq!(rxn.get_products().len(), 1);
    assert_eq!(rxn.get_agents().len(), 0);
    let copy: Mol = (*template).clone();
    assert_eq!(copy.get_num_atoms(), 3);
}

#[test]
fn added_templates_are_copies() {
    let mut rxn = ChemicalReaction::default();
    let mut reactant = mol_from_smarts("[C:1]=[O:2]").unwrap();
    assert_eq!(rxn.add_reactant_template(&reactant), 1);
    assert_eq!(
        rxn.add_product_template(&mol_from_smarts("[C:1][O:2]").unwrap()),
        1
    );
    assert_eq!(rxn.add_agent_template(&smiles("O")), 1);
    reactant.set_prop("changed", "yes");
    assert!(!rxn.get_reactant_template(0).unwrap().has_prop("changed"));
    let products = rxn.run_reactants(&[smiles("CC=O")]).unwrap();
    assert_eq!(product_smiles(&products), [["CCO"]]);
}

#[test]
fn reacting_atoms_and_template_matches() {
    let mut rxn = amide();
    let acid = smiles("CC(=O)O");
    assert!(matches!(
        rxn.is_molecule_reactant(&acid),
        Err(Error::Value { .. })
    ));
    assert!(matches!(rxn.get_reacting_atoms(), Err(Error::Value { .. })));
    rxn.initialize().unwrap();
    assert!(rxn.is_molecule_reactant(&acid).unwrap());
    assert!(!rxn.is_molecule_product(&acid).unwrap());
    assert!(!rxn.is_molecule_reactant(&smiles("CC")).unwrap());
    assert!(!rxn.is_molecule_agent(&acid).unwrap());
    assert_eq!(rxn.get_reacting_atoms().unwrap(), [vec![0, 2], vec![0]]);
    assert_eq!(
        rxn.get_reacting_atoms_with_mapped_atoms_only(true).unwrap(),
        [vec![0], vec![0]]
    );
}

#[test]
fn remove_unmapped_templates() {
    let mut rxn = reaction_from_smarts("[C:1]=[O:2].CCO>>[C:1][O:2]").unwrap();
    let removed = rxn.remove_unmapped_reactant_templates();
    assert_eq!(removed.len(), 1);
    assert_eq!(mol_to_smiles(&removed[0]).unwrap(), "CCO");
    assert_eq!(rxn.get_num_reactant_templates(), 1);
    assert_eq!(rxn.get_num_agent_templates(), 1);
    assert_eq!(
        reaction_to_smarts(&rxn).unwrap(),
        "[C:1]=[O:2]>CCO>[C:1][O:2]"
    );
    assert!(rxn.remove_unmapped_product_templates().is_empty());
    let agents = rxn.remove_agent_templates();
    assert_eq!(agents.len(), 1);
    assert_eq!(rxn.get_num_agent_templates(), 0);

    let mut rxn = reaction_from_smarts("[C:1]=[O:2].CCO>>[C:1][O:2]").unwrap();
    let params = RemoveUnmappedReactantTemplatesParams {
        move_to_agent_templates: false,
        ..Default::default()
    };
    assert_eq!(
        rxn.remove_unmapped_reactant_templates_with_params(&params)
            .len(),
        1
    );
    assert_eq!(rxn.get_num_agent_templates(), 0);
}

#[test]
fn param_defaults_match_python() {
    let params = RemoveUnmappedReactantTemplatesParams::default();
    assert!((params.threshold_unmapped_atoms - 0.2).abs() < 1e-12);
    assert!(params.move_to_agent_templates);
    let params = AddRecursiveQueriesToReactionParams::default();
    assert!(params.queries.is_empty());
    assert_eq!(params.prop_name, "molFileValue");
    let params = PreprocessReactionParams::default();
    assert!(params.queries.is_empty());
    assert_eq!(params.prop_name, "molFileValue");
    assert_eq!(
        ReactionFromSmartsParams::default(),
        ReactionFromSmartsParams {
            replacements: BTreeMap::new(),
            use_smiles: false
        }
    );
    assert_eq!(
        ReactionFromRxnBlockParams::default(),
        ReactionFromRxnBlockParams {
            sanitize: false,
            remove_hs: false,
            strict_parsing: true
        }
    );
    assert_eq!(
        ReactionFromMrvBlockParams::default(),
        ReactionFromMrvBlockParams {
            sanitize: false,
            remove_hs: false
        }
    );
    assert_eq!(
        ReactionToRxnBlockParams::default(),
        ReactionToRxnBlockParams {
            separate_agents: false,
            force_v3000: false
        }
    );
    let png = ReactionMetadataToPNGStringParams::default();
    assert!(png.include_pkl && png.include_smiles && !png.include_smarts && !png.include_rxn);
    let png = ReactionMetadataToPNGFileParams::default();
    assert!(png.include_pkl && png.include_smiles && !png.include_smarts && !png.include_mol);
    let cx = ReactionToCXSmilesParams::default();
    assert_eq!(cx.params, SmilesWriteParams::default());
    assert_eq!(cx.flags, 0x7fffffff);
    let depict = Compute2DCoordsForReactionParams::default();
    assert!((depict.spacing - 1.0).abs() < 1e-12);
    assert!(depict.update_props && depict.canon_orient && !depict.permute_deg4_nodes);
    assert_eq!(
        (
            depict.n_flips_per_sample,
            depict.n_sample,
            depict.sample_seed
        ),
        (0, 0, 0)
    );
    assert!((depict.bond_length + 1.0).abs() < 1e-12);
    let sanitize = SanitizeRxnParams::default();
    assert_eq!(sanitize.sanitize_ops, SanitizeFlags::SANITIZE_ALL);
    assert_eq!(sanitize.params, get_default_adjust_params());
    assert!(!sanitize.catch_errors);
}

#[test]
fn writers() {
    let rxn = amide();
    assert_eq!(reaction_to_smarts(&rxn).unwrap(), AMIDE);
    assert_eq!(
        reaction_to_smiles(&rxn).unwrap(),
        "O[C:1]=[O:2].[NH3:3]>>[C:1](=[O:2])[N:3]"
    );
    let non_canonical = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    assert_eq!(
        reaction_to_smiles_with_params(&rxn, &non_canonical).unwrap(),
        "[C:1](=[O:2])O.[NH3:3]>>[C:1](=[O:2])[N:3]"
    );
    assert_eq!(reaction_to_cx_smarts(&rxn).unwrap(), AMIDE);
    assert_eq!(
        reaction_to_cx_smiles(&rxn).unwrap(),
        "O[C:1]=[O:2].[NH3:3]>>[C:1](=[O:2])[N:3]"
    );
    let labeled = reaction_from_smarts("[CH3:1][OH:2]>>[CH2:1]=[O:2] |$R1;$|").unwrap();
    assert_eq!(
        reaction_to_cx_smiles(&labeled).unwrap(),
        "[CH3:1][OH:2]>>[CH2:1]=[O:2] |$R1;;;$|"
    );
    assert_eq!(
        reaction_to_smarts(&labeled).unwrap(),
        "[C&H3:1][O&H1:2]>>[C&H2:1]=[O:2]"
    );
    let none = ReactionToCXSmilesParams {
        flags: 0,
        ..Default::default()
    };
    assert_eq!(
        reaction_to_cx_smiles_with_params(&labeled, &none).unwrap(),
        "[CH3:1][OH:2]>>[CH2:1]=[O:2]"
    );
    assert_eq!(
        reaction_to_cx_smarts_with_params(&labeled, &none).unwrap(),
        "[C&H3:1][O&H1:2]>>[C&H2:1]=[O:2]"
    );
}

#[test]
fn smarts_and_smiles_parsers() {
    assert_eq!(
        reaction_to_smiles(&reaction_from_smiles("CC(=O)O.N>>CC(N)=O").unwrap()).unwrap(),
        "CC(=O)O.N>>CC(N)=O"
    );
    let replacements = BTreeMap::from([("{acid}".to_owned(), "[C:1](=[O:2])O".to_owned())]);
    let params = ReactionFromSmartsParams {
        replacements: replacements.clone(),
        use_smiles: false,
    };
    let rxn =
        reaction_from_smarts_with_params("{acid}.[N:3]>>[C:1](=[O:2])[N:3]", &params).unwrap();
    assert_eq!(reaction_to_smarts(&rxn).unwrap(), AMIDE);
    let rxn = reaction_from_smiles_with_replacements("{acid}.N>>CC(N)=O", &replacements).unwrap();
    assert_eq!(rxn.get_num_reactant_templates(), 2);
    for bad in ["CC", "C>>C>>C", "[C:1>>C"] {
        assert!(
            matches!(reaction_from_smarts(bad), Err(Error::Parse { .. })),
            "{bad}"
        );
    }
}

#[test]
fn rxn_blocks() {
    let rxn = amide();
    let expected = "[#6:1](=[#8:2])[#8].[#7:3]>>[#6:1](=[#8:2])[#7:3]";
    let block = reaction_to_rxn_block(&rxn).unwrap();
    assert!(block.starts_with("$RXN\n"));
    let parsed = reaction_from_rxn_block(&block).unwrap();
    assert_eq!(reaction_to_smarts(&parsed).unwrap(), expected);
    let v3000 = reaction_to_v3k_rxn_block(&rxn).unwrap();
    assert!(v3000.starts_with("$RXN V3000"));
    assert_eq!(
        reaction_to_smarts(&reaction_from_rxn_block(&v3000).unwrap()).unwrap(),
        expected
    );
    let params = ReactionToRxnBlockParams {
        separate_agents: true,
        force_v3000: true,
    };
    assert!(
        reaction_to_rxn_block_with_params(&rxn, &params)
            .unwrap()
            .contains("M  V30 BEGIN AGENT")
    );
    assert!(
        reaction_to_v3k_rxn_block_with_separate_agents(&rxn, true)
            .unwrap()
            .contains("M  V30 BEGIN AGENT")
    );
    assert!(matches!(
        reaction_from_rxn_block("junk"),
        Err(Error::Parse { .. })
    ));

    let path = common::temp_path("amide.rxn");
    std::fs::write(&path, &block).unwrap();
    assert_eq!(
        reaction_to_smarts(&reaction_from_rxn_file(&path).unwrap()).unwrap(),
        expected
    );
    let strict = ReactionFromRxnFileParams {
        sanitize: true,
        ..Default::default()
    };
    assert!(reaction_from_rxn_file_with_params(&path, &strict).is_ok());
    std::fs::write(&path, "").unwrap();
    assert!(matches!(
        reaction_from_rxn_file(&path),
        Err(Error::Parse { .. })
    ));
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        reaction_from_rxn_file(&path),
        Err(Error::Os { .. })
    ));
}

#[test]
fn rxn_block_with_a_missing_agent() {
    // RDKit's V2000 parser adds a null agent for an agent block cut short,
    // which Python then dereferences (a segfault in ReactionToSmarts).
    let block = reaction_to_rxn_block(&reaction_from_smarts("[C:1]>>[C:1]O").unwrap()).unwrap();
    let mut lines: Vec<&str> = block.lines().collect();
    lines[4] = "  1  1  1";
    let crafted = format!("{}\n$MOL\n", lines.join("\n"));
    let err = reaction_from_rxn_block(&crafted).unwrap_err();
    assert!(
        matches!(&err, Error::Parse { message } if message.contains("null agent")),
        "{err:?}"
    );
}

#[test]
fn marvin() {
    let rxn = reaction_from_rxn_block(&reaction_to_rxn_block(&amide()).unwrap()).unwrap();
    let expected = "[#6:1](=[#8:2])[#8].[#7:3]>>[#6:1](=[#8:2])[#7:3]";
    let mrv = reaction_to_mrv_block(&rxn).unwrap();
    assert!(mrv_block_is_reaction(&mrv).unwrap());
    assert_eq!(
        reaction_to_smarts(&reaction_from_mrv_block(&mrv).unwrap()).unwrap(),
        expected
    );
    assert!(
        reaction_to_mrv_block_with_pretty_print(&rxn, true)
            .unwrap()
            .contains('\t')
    );
    // The Marvin writer refuses SMARTS's single-or-aromatic bond queries.
    assert!(matches!(
        reaction_to_mrv_block(&amide()),
        Err(Error::Runtime { .. })
    ));
    assert!(!mrv_block_is_reaction("<cml><MDocument/></cml>").unwrap());
    assert!(matches!(
        mrv_block_is_reaction("junk"),
        Err(Error::Runtime { .. })
    ));

    let path = common::temp_path("amide.mrv");
    reaction_to_mrv_file(&rxn, &path).unwrap();
    assert!(mrv_file_is_reaction(&path).unwrap());
    let params = ReactionFromMrvFileParams::default();
    assert_eq!(
        reaction_to_smarts(&reaction_from_mrv_file_with_params(&path, &params).unwrap()).unwrap(),
        expected
    );
    reaction_to_mrv_file_with_pretty_print(&rxn, &path, true).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        reaction_from_mrv_file(&path),
        Err(Error::Os { .. })
    ));
    assert!(matches!(mrv_file_is_reaction(&path), Err(Error::Os { .. })));

    let alex = reaction_from_mrv_file(data("AlexRxn.mrv")).unwrap();
    assert_eq!(
        reaction_to_smarts(&alex).unwrap(),
        "[#6]1=[#6](-[#6]=[#6]-[#6]=[#6]-1)-[#6]>>[#6]1=[#6](-[#6]=[#6]-[#6]=[#6]-1)-[#6]"
    );
}

#[test]
fn marvin_bonds_need_two_atoms() {
    // RDKit's Marvin parser copies a bond's two atom refs before it checks
    // that there are two: Python segfaults on these.
    let rxn = reaction_from_rxn_block(&reaction_to_rxn_block(&amide()).unwrap()).unwrap();
    let mrv = reaction_to_mrv_block(&rxn).unwrap();
    let start = mrv.find("atomRefs2=\"").unwrap() + "atomRefs2=\"".len();
    let end = start + mrv[start..].find('"').unwrap();
    let first = mrv[start..end].split(' ').next().unwrap().to_owned();
    for refs in [first.as_str(), ""] {
        let crafted = format!("{}{refs}{}", &mrv[..start], &mrv[end..]);
        let err = reaction_from_mrv_block(&crafted).unwrap_err();
        assert!(
            matches!(&err, Error::Parse { message } if message.contains("atomRefs2")),
            "{err:?}"
        );
        let path = common::temp_path("one-ref.mrv");
        std::fs::write(&path, &crafted).unwrap();
        assert!(matches!(
            reaction_from_mrv_file(&path),
            Err(Error::Parse { .. })
        ));
        std::fs::remove_file(&path).unwrap();
    }
    let missing = mrv.replacen("atomRefs2=", "atomRefs3=", 1);
    assert!(matches!(
        reaction_from_mrv_block(&missing),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        reaction_from_mrv_block("<cml"),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn cdxml() {
    let rxns = reactions_from_cdxml_file(data("rxn1.cdxml")).unwrap();
    assert_eq!(rxns.len(), 1);
    assert_eq!(
        reaction_to_smarts(&rxns[0]).unwrap(),
        "[#6&D2:2]1:[#6&D2:3]:[#6&D2:4]:[#6&D3:1](:[#6&D2:5]:[#6&D2:6]:1)-[#17&D1].\
         [#6&D3](-[#5&D2]-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1)\
         (-[#8&D1])-[#8&D1]>>[#6&D2:1]1:[#6&D2:5]:[#6&D3:6](:[#6&D2:2]:[#6&D2:3]:[#6&D2:4]:1)\
         -[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1"
    );
    let text = std::fs::read_to_string(data("rxn1.cdxml")).unwrap();
    assert_eq!(reactions_from_cdxml_block(&text).unwrap().len(), 1);
    let params = ReactionsFromCDXMLBlockParams {
        sanitize: true,
        remove_hs: true,
    };
    assert_eq!(
        reactions_from_cdxml_block_with_params(&text, &params)
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        reactions_from_cdxml_file_with_params(data("rxn1.cdxml"), &params)
            .unwrap()
            .len(),
        1
    );
    assert!(matches!(
        reactions_from_cdxml_file(common::temp_path("missing.cdxml")),
        Err(Error::Os { .. })
    ));
    // The ChemDraw library appends to a step's atom map until memory runs
    // out once a token fails to read as an int.
    let auto = r#"ReactionStepAtomMapAuto="300 354 302"#;
    assert!(text.contains(auto));
    for token in ["x", "1.5", "2147483648", "4294967295", "5x", "&#52;"] {
        let bad = text.replacen(
            auto,
            &format!(r#"ReactionStepAtomMapAuto="300 354 {token}"#),
            1,
        );
        assert!(
            matches!(reactions_from_cdxml_block(&bad), Err(Error::Parse { .. })),
            "{token}"
        );
    }
    for token in ["2147483647", "-1", "+7"] {
        let fine = text.replacen(
            auto,
            &format!(r#"ReactionStepAtomMapAuto="300 354 {token}"#),
            1,
        );
        assert!(reactions_from_cdxml_block(&fine).is_ok(), "{token}");
    }
}

#[test]
fn molecules_with_roles() {
    let rxn = amide();
    let mol = reaction_to_molecule(&rxn).unwrap();
    assert_eq!(
        mol_to_smiles(&mol).unwrap(),
        "O[C:1]=[O:2].[C:1](=[O:2])[N:3].[NH3:3]"
    );
    assert_eq!(
        reaction_to_smarts(&reaction_from_molecule(&mol).unwrap()).unwrap(),
        AMIDE
    );
    assert_eq!(
        reaction_from_molecule(&smiles("CC"))
            .unwrap()
            .get_num_reactant_templates(),
        0
    );
}

#[test]
fn pickles() {
    let mut rxn = amide();
    rxn.initialize().unwrap();
    rxn.set_prop("a", "b");
    rxn.set_int_prop("n", 3);
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { ChemicalReaction::from_bin_str(&rxn.to_binary()) }.unwrap();
    assert_eq!(reaction_to_smarts(&copy).unwrap(), AMIDE);
    assert!(copy.is_initialized());
    // SAFETY: a pickle RDKit wrote.
    let with_props = unsafe {
        ChemicalReaction::from_bin_str(
            &rxn.to_binary_with_property_flags(PropertyPickleOptions::AllProps),
        )
    }
    .unwrap();
    assert_eq!(
        with_props.get_props_as_dict().unwrap(),
        BTreeMap::from([
            ("a".to_owned(), PropValue::String("b".to_owned())),
            ("n".to_owned(), PropValue::Int(3)),
        ])
    );
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { ChemicalReaction::from_bin_str(b"junk") },
        Err(Error::Runtime { .. })
    ));
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { ChemicalReaction::from_bin_str(&rxn.to_binary()[..20]) }.is_err());
}

#[test]
fn properties() {
    let mut rxn = amide();
    rxn.set_prop("a", "b");
    rxn.set_int_prop("n", 3);
    rxn.set_double_prop("x", 1.5);
    rxn.set_bool_prop("t", true);
    rxn.set_unsigned_prop("u", 7);
    rxn.set_prop_with_computed("c", "d", true);
    assert!(rxn.has_prop("a"));
    assert_eq!(rxn.get_prop("a").unwrap(), "b");
    assert_eq!(rxn.get_int_prop("n").unwrap(), 3);
    assert!((rxn.get_double_prop("x").unwrap() - 1.5).abs() < 1e-12);
    assert!(rxn.get_bool_prop("t").unwrap());
    assert_eq!(rxn.get_unsigned_prop("u").unwrap(), 7);
    assert_eq!(rxn.get_prop_names().unwrap(), ["a", "n", "x", "t", "u"]);
    assert_eq!(
        rxn.get_prop("missing").unwrap_err(),
        Error::Key {
            key: "missing".to_owned()
        }
    );
    rxn.clear_computed_props();
    assert!(!rxn.has_prop("c"));
    rxn.clear_prop("a");
    assert!(!rxn.has_prop("a"));
    assert_eq!(rxn.get_props_as_dict().unwrap().len(), 4);
}

#[test]
fn clone_and_debug() {
    let rxn = amide();
    let copy = rxn.clone();
    assert_eq!(reaction_to_smarts(&copy).unwrap(), AMIDE);
    assert_eq!(
        format!("{rxn:?}"),
        "ChemicalReaction { num_reactant_templates: 2, num_product_templates: 1, \
         num_agent_templates: 0 }"
    );
    assert!(format!("{:?}", rxn.get_reactant_template(0).unwrap()).starts_with("TemplateRef"));
}

#[test]
fn substruct_params() {
    let mut rxn = amide();
    let mut params = rxn.get_substruct_params();
    assert!(!params.use_chirality);
    assert_eq!(params.max_matches, 1000);
    params.use_chirality = true;
    rxn.set_substruct_params(&params);
    assert!(rxn.get_substruct_params().use_chirality);
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { ChemicalReaction::from_bin_str(&rxn.to_binary()) }.unwrap();
    assert!(copy.get_substruct_params().use_chirality);
}

#[test]
fn png_metadata() {
    let rxn = amide();
    let png = std::fs::read(data("reaction1.no_metadata.png")).unwrap();
    let with_rxn = reaction_metadata_to_png_string(&rxn, &png).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        reaction_to_smarts(&unsafe { reaction_from_png_string(&with_rxn) }.unwrap()).unwrap(),
        AMIDE
    );
    let params = ReactionMetadataToPNGStringParams {
        include_pkl: false,
        include_smiles: false,
        include_smarts: true,
        include_rxn: true,
    };
    let smarts_only = reaction_metadata_to_png_string_with_params(&rxn, &png, &params).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        reaction_to_smarts(&unsafe { reaction_from_png_string(&smarts_only) }.unwrap()).unwrap(),
        AMIDE
    );
    let path = common::temp_path("reaction.png");
    std::fs::write(&path, &with_rxn).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        reaction_to_smarts(&unsafe { reaction_from_png_file(&path) }.unwrap()).unwrap(),
        AMIDE
    );
    std::fs::write(&path, &png).unwrap();
    assert_eq!(
        reaction_metadata_to_png_file(&rxn, &path).unwrap(),
        with_rxn
    );
    let file_params = ReactionMetadataToPNGFileParams {
        include_mol: true,
        ..Default::default()
    };
    assert!(reaction_metadata_to_png_file_with_params(&rxn, &path, &file_params).is_ok());
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { reaction_from_png_file(&path) },
        Err(Error::Os { .. })
    ));
    assert!(matches!(
        reaction_metadata_to_png_file(&rxn, &path),
        Err(Error::Os { .. })
    ));
    // No reaction in the metadata.
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { reaction_from_png_string(&png) },
        Err(Error::Parse { .. })
    ));
}

fn png_chunk(tag: &[u8], data: &[u8]) -> Vec<u8> {
    let mut chunk = (data.len() as u32).to_be_bytes().to_vec();
    chunk.extend_from_slice(tag);
    chunk.extend_from_slice(data);
    chunk.extend_from_slice(&[0, 0, 0, 0]);
    chunk
}

const PNG_SIGNATURE: &[u8] = b"\x89PNG\r\n\x1a\n";

#[test]
fn png_checks() {
    let rxn = amide();
    // RDKit reads a short signature through bytes it never read.
    for short in [&b""[..], b"\x89PNG"] {
        assert!(matches!(
            // SAFETY: too short to be PNG data, which the wrapper refuses before RDKit reads it.
            unsafe { reaction_from_png_string(short) },
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            reaction_metadata_to_png_string(&rxn, short),
            Err(Error::Parse { .. })
        ));
    }
    // A text chunk whose key fills it: RDKit takes front() of an empty
    // string.
    let mut empty_text = PNG_SIGNATURE.to_vec();
    empty_text.extend(png_chunk(b"tEXt", b"rxnSmiles\0"));
    empty_text.extend(png_chunk(b"IEND", b""));
    // SAFETY: a pickle RDKit wrote.
    let err = unsafe { reaction_from_png_string(&empty_text) }.unwrap_err();
    assert!(
        matches!(&err, Error::Parse { message } if message.contains("empty text chunk")),
        "{err:?}"
    );
    // A PNG cut off before its IEND chunk: RDKit's writer reads chunk
    // headers past the end without checking.
    let png = std::fs::read(data("reaction1.no_metadata.png")).unwrap();
    for cut in [png.len() - 12, png.len() - 4, png.len() - 1] {
        let err = reaction_metadata_to_png_string(&rxn, &png[..cut]).unwrap_err();
        assert!(
            matches!(&err, Error::Parse { message } if message.contains("IEND")),
            "{err:?}"
        );
    }
    // A chunk whose length runs past the end.
    let mut long = PNG_SIGNATURE.to_vec();
    long.extend(&png_chunk(b"tEXt", b"rxnSmiles\0C>>C")[..10]);
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { reaction_from_png_string(&long) }.is_err());
}

#[test]
fn fingerprints() {
    let rxn = amide();
    let fp = create_difference_fingerprint_for_reaction(&rxn).unwrap();
    assert_eq!(fp.get_length(), 2048);
    assert_eq!(
        fp.get_nonzero_elements(),
        BTreeMap::from([(367, 10), (399, -10), (497, 10), (1731, -10)])
    );
    let fp = create_structural_fingerprint_for_reaction(&rxn).unwrap();
    assert_eq!(fp.get_num_bits(), 4096);
    assert_eq!(fp.get_num_on_bits(), 28);
    assert_eq!(
        &fp.get_on_bits()[..10],
        [71, 287, 292, 293, 319, 1077, 1260, 1270, 1287, 1301]
    );

    let mut params = ReactionFingerprintParams::new(true, 0.2, 1, 1, 64, FingerprintType::MorganFp);
    let fp =
        create_difference_fingerprint_for_reaction_with_reaction_finger_print_params(&rxn, &params)
            .unwrap();
    assert_eq!(
        fp.get_nonzero_elements(),
        BTreeMap::from([(10, -1), (18, -1), (21, 2), (34, -1), (40, 1), (63, -1)])
    );
    params.fp_type = FingerprintType::RdkitFp;
    let fp =
        create_structural_fingerprint_for_reaction_with_reaction_finger_print_params(&rxn, &params)
            .unwrap();
    assert_eq!(fp.get_num_bits(), 64);
    assert_eq!(fp.get_on_bits(), [2, 5, 8, 11, 15, 28, 30, 31, 43, 51]);
}

#[test]
fn fingerprint_params_defaults() {
    let params = ReactionFingerprintParams::default();
    assert!(!params.include_agents);
    assert!((params.bit_ratio_agents - 0.2).abs() < 1e-12);
    assert_eq!(params.non_agent_weight, 10);
    assert_eq!(params.agent_weight, 1);
    assert_eq!(params.fp_size, 2048);
    assert_eq!(params.fp_type, FingerprintType::AtomPairFp);
}

#[test]
fn fingerprint_checks() {
    let rxn = amide();
    let difference = |params: &ReactionFingerprintParams| {
        create_difference_fingerprint_for_reaction_with_reaction_finger_print_params(&rxn, params)
    };
    let structural = |params: &ReactionFingerprintParams| {
        create_structural_fingerprint_for_reaction_with_reaction_finger_print_params(&rxn, params)
    };
    // Python's weights overflow RDKit's int counts: a weight of 2^31 gives
    // every bit -2^31.
    let mut params =
        ReactionFingerprintParams::new(true, 0.2, 1 << 31, 1, 2048, FingerprintType::AtomPairFp);
    assert!(matches!(difference(&params), Err(Error::Value { .. })));
    params.non_agent_weight = 1 << 30;
    assert_eq!(
        difference(&params).unwrap().get_nonzero_elements(),
        BTreeMap::from([
            (367, 1 << 30),
            (399, -(1 << 30)),
            (497, 1 << 30),
            (1731, -(1 << 30))
        ])
    );
    // A Morgan count of 2 times 2^30: Python's gives -2^31.
    let mut morgan =
        ReactionFingerprintParams::new(true, 0.2, 1 << 30, 1, 64, FingerprintType::MorganFp);
    let err = difference(&morgan).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    morgan.non_agent_weight = 1 << 29;
    assert_eq!(
        difference(&morgan).unwrap().get_nonzero_elements()[&21],
        1 << 30
    );
    params.non_agent_weight = 1;
    params.agent_weight = i32::MIN;
    assert!(difference(&params).is_ok());
    params.fp_size = 0;
    assert!(matches!(difference(&params), Err(Error::Runtime { .. })));
    params.fp_size = 2048;
    params.fp_type = FingerprintType::PatternFp;
    assert!(matches!(difference(&params), Err(Error::Runtime { .. })));

    // RDKit casts ceil(fpSize * bitRatioAgents) to int.
    let mut params = ReactionFingerprintParams::new(
        true,
        f64::NEG_INFINITY,
        1,
        1,
        64,
        FingerprintType::PatternFp,
    );
    assert!(matches!(structural(&params), Err(Error::Value { .. })));
    params.bit_ratio_agents = -1e300;
    assert!(matches!(structural(&params), Err(Error::Value { .. })));
    params.bit_ratio_agents = f64::NAN;
    assert_eq!(structural(&params).unwrap().get_num_bits(), 64);
    params.bit_ratio_agents = 0.5;
    params.fp_size = u32::MAX;
    assert!(matches!(structural(&params), Err(Error::Value { .. })));
    params.fp_size = 0;
    assert!(matches!(structural(&params), Err(Error::Runtime { .. })));
}

#[test]
fn sanitization() {
    let mut rxn = reaction_from_rxn_file(data("boronic1.rxn")).unwrap();
    let expected = "[#35]-[*:1].[*:2]-[#5](-[#8])-[#8]>>[*:1]-[*:2]";
    assert_eq!(reaction_to_smarts(&rxn).unwrap(), expected);
    assert_eq!(
        sanitize_rxn(&mut rxn).unwrap(),
        SanitizeFlags::SANITIZE_NONE
    );
    assert_eq!(reaction_to_smarts(&rxn).unwrap(), expected);
    let params = SanitizeRxnParams {
        sanitize_ops: SanitizeFlags::SANITIZE_ATOM_MAPS | SanitizeFlags::SANITIZE_MERGEHS,
        params: match_only_at_rgroups_adjust_params(),
        catch_errors: true,
    };
    assert_eq!(
        sanitize_rxn_with_params(&mut rxn, &params).unwrap(),
        SanitizeFlags::SANITIZE_NONE
    );
    sanitize_rxn_as_mols(&mut rxn).unwrap();
    sanitize_rxn_as_mols_with_sanitize_ops(&mut rxn, rdkit_cxx::chem::SanitizeFlags::SANITIZE_NONE)
        .unwrap();

    assert!(get_default_adjust_params().aromatize_if_possible);
    assert!(match_only_at_rgroups_adjust_params().adjust_degree);
    assert!(get_chem_draw_rxn_adjust_params().adjust_degree);
    assert_eq!(SanitizeFlags::SANITIZE_ALL.bits(), u32::MAX);
    assert_eq!(SanitizeFlags::SANITIZE_MERGEHS.bits(), 8);
    assert_eq!(
        format!("{:?}", SanitizeFlags::SANITIZE_ATOM_MAPS),
        "SanitizeFlags::SANITIZE_ATOM_MAPS"
    );
}

#[test]
fn sanitize_rxn_checks_atom_map_room() {
    // fixAtomMaps maps the R group to its label plus the largest atom map,
    // summed in int: Python's gives INT_MIN.
    let mut rxn = ChemicalReaction::default();
    let mut reactant = mol_from_smarts("[*][C:1]").unwrap();
    reactant
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_unsigned_prop("_MolFileRLabel", i32::MAX as u32);
    rxn.add_reactant_template(&reactant);
    rxn.add_product_template(&mol_from_smarts("[C:1]").unwrap());
    let params = SanitizeRxnParams {
        sanitize_ops: SanitizeFlags::SANITIZE_ATOM_MAPS,
        ..Default::default()
    };
    let err = sanitize_rxn_with_params(&mut rxn, &params).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    let caught = SanitizeRxnParams {
        catch_errors: true,
        ..params
    };
    assert_eq!(
        sanitize_rxn_with_params(&mut rxn, &caught).unwrap(),
        SanitizeFlags::SANITIZE_ATOM_MAPS
    );
}

#[test]
fn preprocessing() {
    let mut rxn = reaction_from_rxn_file(data("boronic1.rxn")).unwrap();
    rxn.initialize().unwrap();
    let labels = vec![
        vec![(0, "halogen.bromine.aromatic".to_owned())],
        vec![(1, "boronicacid".to_owned())],
    ];
    let result = preprocess_reaction(&mut rxn).unwrap();
    assert_eq!(
        (
            result.n_warn,
            result.n_error,
            result.n_reacts,
            result.n_prods
        ),
        (0, 0, 2, 1)
    );
    assert_eq!(result.reactant_labels, labels);
    let template = rxn.get_reactant_template(0).unwrap();
    assert!(!smiles("CCBr").has_substruct_match(&template).unwrap());
    assert!(smiles("c1ccccc1Br").has_substruct_match(&template).unwrap());

    let mut uninitialized = reaction_from_rxn_file(data("boronic1.rxn")).unwrap();
    let result = preprocess_reaction(&mut uninitialized).unwrap();
    assert_eq!(
        (
            result.n_warn,
            result.n_error,
            result.n_reacts,
            result.n_prods
        ),
        (0, 0, 2, 1)
    );
    assert_eq!(result.reactant_labels, labels);

    let mut invalid = reaction_from_smarts("C>>").unwrap();
    let params = PreprocessReactionParams {
        queries: BTreeMap::from([("x".to_owned(), mol_from_smarts("C").unwrap())]),
        ..Default::default()
    };
    let result = preprocess_reaction_with_params(&mut invalid, &params).unwrap();
    assert_eq!(
        (
            result.n_warn,
            result.n_error,
            result.n_reacts,
            result.n_prods
        ),
        (1, 1, 1, 0)
    );
    assert!(result.reactant_labels.is_empty());
}

#[test]
fn recursive_queries() {
    let mut rxn = reaction_from_rxn_file(data("boronic1.rxn")).unwrap();
    assert!(matches!(
        rxn.add_recursive_queries_to_reaction(),
        Err(Error::Value { .. })
    ));
    rxn.initialize().unwrap();
    // The labels name queries the map lacks.
    assert!(matches!(
        rxn.add_recursive_queries_to_reaction(),
        Err(Error::Key { .. })
    ));
    let params = AddRecursiveQueriesToReactionParams {
        queries: BTreeMap::from([
            (
                "halogen.bromine.aromatic".to_owned(),
                mol_from_smarts("Brc").unwrap(),
            ),
            ("boronicacid".to_owned(), mol_from_smarts("B(O)O").unwrap()),
        ]),
        ..Default::default()
    };
    let mut copy = rxn.clone();
    copy.add_recursive_queries_to_reaction_with_params(&params)
        .unwrap();
    assert_eq!(
        rxn.add_recursive_queries_to_reaction_get_labels_with_params(&params)
            .unwrap(),
        [
            vec![(0, "halogen.bromine.aromatic".to_owned())],
            vec![(1, "boronicacid".to_owned())],
        ]
    );
    let template = rxn.get_reactant_template(0).unwrap();
    assert!(!smiles("CCBr").has_substruct_match(&template).unwrap());
    let mut unlabeled = amide();
    unlabeled.initialize().unwrap();
    assert_eq!(
        unlabeled
            .add_recursive_queries_to_reaction_get_labels()
            .unwrap(),
        [vec![], vec![]]
    );
}

#[test]
fn side_chains_and_atom_maps() {
    let mut rxn = amide();
    let products = rxn
        .run_reactants(&[smiles("CC(=O)O"), smiles("NC")])
        .unwrap();
    assert_eq!(
        mol_to_smiles(&reduce_product_to_side_chains(&products[0][0]).unwrap()).unwrap(),
        "C[*:1].C[*:3]"
    );
    assert_eq!(
        mol_to_smiles(
            &reduce_product_to_side_chains_with_add_dummy_atoms(&products[0][0], false).unwrap()
        )
        .unwrap(),
        "C.C"
    );

    assert!(has_reaction_atom_mapping(&rxn));
    let mut unmapped = amide();
    remove_mapping_numbers_from_reactions(&mut unmapped);
    assert_eq!(reaction_to_smarts(&unmapped).unwrap(), "C(=O)O.N>>C(=O)N");
    assert!(!has_reaction_atom_mapping(&unmapped));
    update_products_stereochemistry(&mut unmapped).unwrap();

    let query = reaction_from_smarts("[C:1]=[O:2]>>[C:1]").unwrap();
    assert!(has_reaction_substruct_match(&rxn, &query).unwrap());
    assert!(has_reaction_substruct_match_with_include_agents(&rxn, &query, true).unwrap());
    assert!(has_reactant_template_substruct_match(&rxn, &query).unwrap());
    assert!(has_product_template_substruct_match(&rxn, &query).unwrap());
    assert!(has_agent_template_substruct_match(&rxn, &query).unwrap());
    assert!(!has_reactant_template_substruct_match(&query, &rxn).unwrap());

    assert!(is_reaction_template_molecule_agent(&smiles("CCO"), 0.2));
    assert!(!is_reaction_template_molecule_agent(
        &mol_from_smarts("[C:1]O").unwrap(),
        0.2
    ));
    assert!(is_reaction_template_molecule_agent(
        &smiles("CCO"),
        f64::NAN
    ));
}

#[test]
fn depiction() {
    let mut rxn = amide();
    compute_2d_coords_for_reaction(&mut rxn).unwrap();
    let counts: Vec<u32> = rxn
        .get_reactants()
        .map(|t| t.get_num_conformers())
        .collect();
    assert_eq!(counts, [1, 1]);

    let params = Compute2DCoordsForReactionParams {
        bond_length: 2.0,
        ..Default::default()
    };
    let mut rxn = amide();
    compute_2d_coords_for_reaction_with_params(&mut rxn, &params).unwrap();
    let template = rxn.get_reactant_template(0).unwrap();
    let conf = template.get_conformer().unwrap();
    let (a, b) = (
        conf.get_atom_position(0).unwrap(),
        conf.get_atom_position(1).unwrap(),
    );
    let length = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt();
    assert!((length - 2.0).abs() < 1e-6, "{length}");
    // The bond length is restored after the call.
    let mut again = amide();
    compute_2d_coords_for_reaction(&mut again).unwrap();
    let template = again.get_reactant_template(0).unwrap();
    let conf = template.get_conformer().unwrap();
    let (a, b) = (
        conf.get_atom_position(0).unwrap(),
        conf.get_atom_position(1).unwrap(),
    );
    let length = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt();
    assert!((length - 1.5).abs() < 1e-6, "{length}");
}

fn square_planar_refused<T>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { message }) if message.contains("square-planar"))
}

/// RDKit's depictor reads the first neighbor of a square-planar atom
/// unchecked: Python segfaults in Compute2DCoordsForReaction and in
/// ReactionToRxnBlock, whose V2000 writer depicts a template without a
/// conformer. The paths that depict a template check it first.
#[test]
fn depiction_refuses_what_the_depictor_would_mishandle() {
    let mut platinum = smiles("[Pt]");
    platinum
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiSquareplanar);
    let mut rxn = reaction_from_smarts("[C:1]>>[C:1]").unwrap();
    rxn.add_product_template(&platinum);
    assert!(square_planar_refused(reaction_to_rxn_block(&rxn)));
    let png = std::fs::read(data("reaction1.no_metadata.png")).unwrap();
    let with_rxn = ReactionMetadataToPNGStringParams {
        include_rxn: true,
        ..Default::default()
    };
    assert!(square_planar_refused(
        reaction_metadata_to_png_string_with_params(&rxn, &png, &with_rxn)
    ));
    assert!(square_planar_refused(compute_2d_coords_for_reaction(
        &mut rxn
    )));
    // The writers that do not depict write it.
    assert!(
        reaction_to_v3k_rxn_block(&rxn)
            .unwrap()
            .starts_with("$RXN V3000")
    );
    reaction_metadata_to_png_string(&rxn, &png).unwrap();

    // So is a bond length past the depictor's coordinate bound.
    let params = Compute2DCoordsForReactionParams {
        bond_length: 1e6,
        ..Default::default()
    };
    assert!(matches!(
        compute_2d_coords_for_reaction_with_params(&mut amide(), &params),
        Err(Error::Value { .. })
    ));
}

#[test]
fn reactions_move_between_threads() {
    let mut rxn = amide();
    let products = std::thread::spawn(move || {
        let products = rxn
            .run_reactants(&[smiles("C(=O)O"), smiles("CNC")])
            .unwrap();
        product_smiles(&products)
    })
    .join()
    .unwrap();
    assert_eq!(products, [["CN(C)C=O"]]);
}

/// A V3000 rxn block whose product has a stereo-group line with `pad` extra
/// spaces.
fn stereo_group_rxn_block(pad: usize) -> String {
    let mut rxn = reaction_from_smarts("[C:1]>>[C:1]").unwrap();
    rxn.add_product_template(&smiles("C[C@H](F)Cl |&1:1|"));
    let block = reaction_to_v3k_rxn_block(&rxn).unwrap();
    assert!(block.contains("MDLV30/STERAC1 ATOMS=(1 2)"));
    block.replace(
        "MDLV30/STERAC1 ATOMS=(1 2)",
        &format!("MDLV30/STERAC1 ATOMS=(1 2{})", " ".repeat(pad)),
    )
}

fn png_with_text(key: &str, value: &[u8]) -> Vec<u8> {
    let mut text = key.as_bytes().to_vec();
    text.push(0);
    text.extend_from_slice(value);
    let mut png = PNG_SIGNATURE.to_vec();
    png.extend(png_chunk(b"tEXt", &text));
    png.extend(png_chunk(b"IEND", b""));
    png
}

/// RDKit reads an rxn block's templates with its mol-block reader, whose
/// regex recurses once per character of a stereo-group line (U12).
#[test]
fn long_stereo_group_lines_are_refused() {
    let short = stereo_group_rxn_block(3000);
    let long = stereo_group_rxn_block(100_000);
    let rxn = reaction_from_rxn_block(&short).unwrap();
    assert_eq!(rxn.get_num_product_templates(), 2);
    assert!(matches!(
        reaction_from_rxn_block(&long),
        Err(Error::Value { .. })
    ));
    let path = common::temp_path("stereo-group.rxn");
    std::fs::write(&path, &long).unwrap();
    assert!(matches!(
        reaction_from_rxn_file(&path),
        Err(Error::Value { .. })
    ));
    std::fs::write(&path, &short).unwrap();
    reaction_from_rxn_file(&path).unwrap();
    // An rxn block in PNG metadata is read the same way.
    // SAFETY: a pickle RDKit wrote.
    unsafe { reaction_from_png_string(&png_with_text("ReactionRxn", short.as_bytes())) }.unwrap();
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { reaction_from_png_string(&png_with_text("ReactionRxn", long.as_bytes())) },
        Err(Error::Value { .. })
    ));
}

/// A reaction whose agent is methane with plain monomer info named "MARK",
/// its pickle, and where the info's type byte is in it.
fn marked_reaction_pickle() -> (ChemicalReaction, Vec<u8>, usize) {
    let mut methane = smiles("C");
    methane
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_monomer_info(&AtomMonomerInfo::new_with_name(AtomMonomerType::Other, "MARK").unwrap());
    let mut rxn = reaction_from_smarts("[C:1]>>[C:1]").unwrap();
    rxn.add_agent_template(&methane);
    let pkl = rxn.to_binary();
    // The info's name (a u32 length, then its bytes), then its type, in 32
    // bits low byte first.
    let name = pkl
        .windows(12)
        .position(|w| w == b"\x04\0\0\0MARK\x02\0\0\0")
        .unwrap();
    (rxn, pkl, name + 8)
}

/// A reaction pickle's templates are molecule pickles: RDKit 2025.09 refuses
/// a monomer type its enum lacks in them.
#[test]
fn reaction_pickles_cannot_give_plain_info_the_pdb_type() {
    let (_, pkl, at) = marked_reaction_pickle();
    // SAFETY: a pickle RDKit wrote.
    unsafe { ChemicalReaction::from_bin_str(&pkl) }.unwrap();
    // SAFETY: as above.
    unsafe { reaction_from_png_string(&png_with_text("rdkitReactionPKL", &pkl)) }.unwrap();
    let mut changed = pkl.clone();
    changed[at] = 3;
    // SAFETY: a pickle RDKit wrote but for a monomer type, which RDKit 2025.09
    // refuses.
    assert!(unsafe { ChemicalReaction::from_bin_str(&changed) }.is_err());
    // SAFETY: as above.
    assert!(
        unsafe { reaction_from_png_string(&png_with_text("rdkitReactionPKL", &changed)) }.is_err()
    );
}
