// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Inputs RDKit 2024.03's parsers and writers misread (CONVENTIONS.md U1):
// each is refused with `Error::Value`, and the valid inputs next to them
// still read and write as 2024.03's Python does.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_reactions::*;
use rdkit_cxx::chem::rdinchi::{mol_to_inchi, mol_to_inchi_key};
use rdkit_cxx::chem::*;

fn is_value_error<T>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

// A refusal of a CXSMILES field by the wrapper, as opposed to RDKit's own
// errors, which may be `Error::Value` too.
fn refused_for_cx<T>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { message }) if message.contains("CXSMILES"))
}

// RDKit's own `Error::Value`, not the wrapper's refusal of a field.
fn rdkit_value_error<T>(result: &Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { message }) if !message.contains("CXSMILES"))
}

fn tmp_path(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("soundness_parsers_{}_{name}", std::process::id()))
}

// CXSMILES extensions.

#[test]
fn cx_fields_that_run_past_the_end_are_refused() {
    for smiles in [
        // parse_data_sgroup_attr reads past the end of the extension text.
        "C |SgD:0:x|",
        "C |SgD:0:x",
        "C |SgD:",
        "C |SgD:0",
        "CC |SgD:0,1:name:value",
        // parse_atom_props dereferences end() + 1.
        "C |atomProp:0.p.v",
        "CO |atomProp:0.p.v:1.q.w",
    ] {
        assert!(is_value_error(mol_from_smiles(smiles)), "{smiles}");
        let lenient = SmilesParserParams {
            strict_cxsmiles: false,
            ..Default::default()
        };
        assert!(
            is_value_error(mol_from_smiles_with_params(smiles, &lenient)),
            "{smiles}"
        );
        assert!(is_value_error(mol_from_smarts(smiles)), "{smiles}");
    }
    // parse_polymer_sgroup steps to end() + 2, which MolFromSmiles copies
    // the text up to when strictCXSMILES is off; strictly, the parse fails.
    for smiles in ["*CC* |Sg:n", "*CC* |Sg:n|"] {
        let lenient = SmilesParserParams {
            strict_cxsmiles: false,
            ..Default::default()
        };
        assert!(
            is_value_error(mol_from_smiles_with_params(smiles, &lenient)),
            "{smiles}"
        );
        assert!(
            matches!(mol_from_smiles(smiles), Err(Error::Parse { .. })),
            "{smiles}"
        );
    }
    // Without CXSMILES the text is a name.
    let no_cx = SmilesParserParams {
        allow_cxsmiles: false,
        ..Default::default()
    };
    let mol = mol_from_smiles_with_params("C |atomProp:0.p.v", &no_cx).unwrap();
    assert_eq!(mol.get_prop("_Name").unwrap(), "|atomProp:0.p.v");
}

#[test]
fn valid_cx_fields_still_parse() {
    // Expected values from RDKit 2024.03.6's Python.
    for (smiles, cx) in [
        ("C |SgD:0:x:y::::|", "C |SgD:0:x:y::::|"),
        ("C |SgD:0:x:y:|", "C |SgD:0:x:y::::|"),
        (
            "CC |SgD:0,1:PieceName:Fragment 1:like:unit:t:(1.,1.)|",
            "CC |SgD:0,1:PieceName:Fragment 1:like:unit:t:|",
        ),
        ("C |atomProp:0.p.v|", "C |atomProp:0.p.v|"),
        ("CO |atomProp:0.p.v:1.q.w|", "CO |atomProp:0.p.v:1.q.w|"),
        ("*CC* |Sg:n:1,2::ht|", "*CC* |Sg:n:1,2::ht:::|"),
        ("C[C@H](O)CC |&1:1|", "CC[C@H](C)O |&1:2|"),
        ("C |SgD:0:x:y::::| name", "C |SgD:0:x:y::::|"),
    ] {
        let mol = mol_from_smiles(smiles).unwrap();
        assert_eq!(mol_to_cx_smiles(&mol).unwrap(), cx, "{smiles}");
    }
    let mol = mol_from_smiles("CO |atomProp:0.p.v:1.q.w|").unwrap();
    assert_eq!(
        mol.get_atom_with_idx(1).unwrap().get_prop("q").unwrap(),
        "w"
    );
    let query = mol_from_smarts("[#6] |atomProp:0.p.v|").unwrap();
    assert_eq!(
        query.get_atom_with_idx(0).unwrap().get_prop("p").unwrap(),
        "v"
    );
}

#[test]
fn reaction_cx_fields_are_checked() {
    for smarts in ["[C:1]>>[C:1]O |atomProp:0.p.v", "[C:1]>>[C:1]O |SgD:0:x|"] {
        assert!(is_value_error(reaction_from_smarts(smarts)), "{smarts}");
    }
    // A reaction's extension that does not parse is RDKit's error, whatever
    // strictCXSMILES says, so a polymer S-group copies nothing.
    let polymer = reaction_from_smarts("[C:1]>>[C:1]O|Sg:n");
    assert!(polymer.is_err() && !refused_for_cx(polymer));
    let rxn = reaction_from_smarts("[C:1]>>[C:1]O |atomProp:0.p.v|").unwrap();
    assert_eq!(rxn.get_num_reactant_templates(), 1);
    assert_eq!(
        rxn.get_reactant_template(0)
            .unwrap()
            .get_atom_with_idx(0)
            .unwrap()
            .get_prop("p")
            .unwrap(),
        "v"
    );
}

#[test]
fn suppliers_parse_smiles_without_cx_extensions() {
    // RDKit's SMILES and TDT suppliers parse with allowCXSMILES off
    // (SmilesMolSupplierParams and TDTMolSupplierParams in MolSupplier.h),
    // so a CXSMILES field is part of the SMILES, which does not parse, as
    // from Python.
    let tabbed = SmilesMolSupplierParams {
        delimiter: "\t".into(),
        title_line: false,
        ..Default::default()
    };
    for text in ["C |atomProp:0.p.v\tfirst\n", "C |atomProp:0.p.v|\tfirst\n"] {
        let mut suppl = smiles_mol_supplier_from_text_with_params(text, &tabbed).unwrap();
        assert!(matches!(suppl.get(0), Err(Error::Parse { .. })), "{text:?}");
    }
    let mut tdt = TDTMolSupplier::default();
    tdt.set_data("$SMI<C |atomProp:0.p.v>\nCAS<1>\n|\n")
        .unwrap();
    assert!(matches!(tdt.get(0), Err(Error::Parse { .. })));
}

// V3000 S-group values.

const SUP_BLOCK: &str = "
  RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1.5 0 0 0
M  V30 3 O 3 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(2 2 3) XBONDS=(1 1) LABEL=EtOH CSTATE=(4 1 1.5 0 0) SAP=(3 2 1 1)
M  V30 END SGROUP
M  V30 END CTAB
M  END
";

const SGROUP_LINE: &str =
    "M  V30 1 SUP 1 ATOMS=(2 2 3) XBONDS=(1 1) LABEL=EtOH CSTATE=(4 1 1.5 0 0) SAP=(3 2 1 1)";

fn with_sgroup_line(line: &str) -> String {
    SUP_BLOCK.replace(SGROUP_LINE, line)
}

fn sd_record(block: &str, i: usize) -> String {
    format!("mol{i}{block}$$$$\n")
}

#[test]
fn malformed_sgroup_values_are_refused() {
    let mol = mol_from_mol_block(SUP_BLOCK).unwrap();
    assert_eq!(mol.get_num_atoms(), 3);
    let bad_lines = [
        "M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE=(",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE=(4",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE=(4 ",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE=",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) SAP=(",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) SAP=(3 2 1",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) SAP=(3 x 1 1)",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) SAP=(99999999999 2 1 1)",
        "M  V30 1 SUP 1 ATOMS=(2",
        "M  V30 1 SUP 1 ATOMS=(2 ",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) PARENT=",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) COMPNO= ",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) BRKXYZ=(9",
        // Joined across a continuation.
        "M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE=(4-\nM  V30 ",
    ];
    for line in bad_lines {
        let block = with_sgroup_line(line);
        assert!(is_value_error(mol_from_mol_block(&block)), "{line}");
        let mut suppl = SDMolSupplier::default();
        let read = suppl
            .set_data(&sd_record(&block, 0))
            .and_then(|()| suppl.get(0));
        assert!(is_value_error(read), "{line}");
    }
    // Values that give RDKit's reads something to read still parse.
    for line in [
        "M  V30 1 SUP 1 ATOMS=(2 2 3) XBONDS=(1 1) CSTATE=(4 1 1.5 0 0)",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) XBONDS=(1 1) CSTATE=(4 1 1.5 0 0) SAP=(3 2 1 1)",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) XBONDS=(1 1) SAP=(3 2 AIDX 1)",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) COMPNO=1",
        "M  V30 1 SUP 1 ATOMS=(0)",
        "M  V30 1 SUP 1 ATOMS=(2 2 3) LABEL=\"CSTATE=(\"",
    ] {
        let block = with_sgroup_line(line);
        assert!(mol_from_mol_block(&block).is_ok(), "{line}");
    }
}

#[test]
fn sgroup_values_are_checked_across_chunks() {
    // An SD file of 64 KiB chunks, where values straddle chunk boundaries:
    // all records read, and a bad value at a boundary is refused.
    let record = sd_record(SUP_BLOCK, 0);
    let count = 3 * 65536 / record.len() + 3;
    let mut text = String::new();
    for i in 0..count {
        text.push_str(&sd_record(SUP_BLOCK, i));
    }
    let path = tmp_path("sgroups.sdf");
    std::fs::write(&path, &text).unwrap();
    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert_eq!(suppl.len().unwrap() as usize, count);
    for i in 0..count as u32 {
        assert_eq!(suppl.get(i).unwrap().get_num_atoms(), 3, "record {i}");
    }
    let mut forward = ForwardSDMolSupplier::from_fileobj(text.as_bytes()).unwrap();
    assert_eq!(forward.mols().filter(|m| m.is_ok()).count(), count);
    // Every cut through "CSTATE=(4 1": the line is cut off by the end of
    // the input there. RDKit misreads it once it has read the label.
    let at = text.find("CSTATE=(4 1").unwrap();
    for end in at..at + 11 {
        let mut suppl = SDMolSupplier::default();
        let read = suppl.set_data(&text[..end]).and_then(|()| suppl.get(0));
        assert_eq!(is_value_error(read), end > at + 5, "{end}");
    }
    // A bad CSTATE that crosses the first 64 KiB boundary, at each of its
    // bytes: refused, with the records before it read.
    let bad_record = sd_record(
        &with_sgroup_line("M  V30 1 SUP 1 ATOMS=(2 2 3) CSTATE=(4"),
        0,
    );
    let token = bad_record.find("CSTATE=(4").unwrap();
    for split in 0..10 {
        let mut prefix = String::new();
        while prefix.len() + record.len() + token + split < 65536 {
            prefix.push_str(&record);
        }
        let pad = 65536 - (prefix.len() + token + split);
        let mut bad = prefix.replacen("mol0", &format!("mol0{}", "x".repeat(pad)), 1);
        let good = bad.matches("$$$$").count();
        bad.push_str(&bad_record);
        bad.push_str(&record);
        assert_eq!(bad.find("CSTATE=(4\n").unwrap() + split, 65536);
        std::fs::write(&path, &bad).unwrap();
        let mut suppl = SDMolSupplier::new(&path).unwrap();
        assert!(suppl.get(0).is_ok(), "{split}");
        let read = (0..good as u32 + 2).try_for_each(|i| suppl.get(i).map(drop));
        assert!(is_value_error(read), "{split}");
        let mut forward = ForwardSDMolSupplier::from_fileobj(bad.as_bytes()).unwrap();
        let read: Result<Vec<_>, _> = forward.mols().collect();
        assert!(is_value_error(read), "{split}");
        let mut forward = ForwardSDMolSupplier::from_fileobj(bad.as_bytes()).unwrap();
        assert_eq!(
            forward.mols().take_while(|m| m.is_ok()).count(),
            good,
            "{split}"
        );
    }
    std::fs::remove_file(&path).unwrap();
}

// SMARTS in mol blocks.

#[test]
fn cx_fields_in_mol_block_smarts_are_refused() {
    let v2000 = "
  RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  MRV SMA   1 [#6] |atomProp:0.p.v
M  END
";
    assert!(is_value_error(mol_from_mol_block(v2000)));
    let fine = v2000.replace(" |atomProp:0.p.v", "");
    assert!(mol_from_mol_block(&fine).is_ok());
    // Data S-group fields, joined across SCD and SED lines as RDKit joins
    // them, and V3000 FIELDDATA.
    let split = v2000.replace(
        "M  MRV SMA   1 [#6] |atomProp:0.p.v\n",
        "M  SCD   1 [#6] |atomPr\nM  SED   1 op:0.p.v\n",
    );
    assert!(is_value_error(mol_from_mol_block(&split)));
    let v3000 = with_sgroup_line(
        "M  V30 1 DAT 1 ATOMS=(1 1) FIELDNAME=SMARTSQ FIELDDATA=\"[#6] |SgD:0:x\"",
    );
    assert!(is_value_error(mol_from_mol_block(&v3000)));
}

// V3000 rxn blocks with agents.

const RXN_NO_AGENTS: &str = "$RXN V3000\n\n      RDKit\n\nM  V30 COUNTS 1 1\nM  V30 BEGIN REACTANT\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 1\nM  V30 END ATOM\nM  V30 END CTAB\nM  V30 END REACTANT\nM  V30 BEGIN PRODUCT\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 1\nM  V30 2 O 0.000000 0.000000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 6 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  V30 END PRODUCT\nM  END\n";

#[test]
fn v3000_rxn_agents_are_refused() {
    let rxn = reaction_from_rxn_block(RXN_NO_AGENTS).unwrap();
    assert_eq!(rxn.get_num_reactant_templates(), 1);
    assert_eq!(rxn.get_num_product_templates(), 1);
    assert_eq!(reaction_to_smarts(&rxn).unwrap(), "[#6:1]>>[#6:1][#8]");
    let agent = "M  V30 BEGIN AGENT\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 O 0.000000 0.000000 0.000000 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  V30 END AGENT\n";
    let with_agent = RXN_NO_AGENTS
        .replace("M  V30 COUNTS 1 1\n", "M  V30 COUNTS 1 1 1\n")
        .replace(
            "M  V30 END PRODUCT\n",
            &format!("M  V30 END PRODUCT\n{agent}"),
        );
    assert!(is_value_error(reaction_from_rxn_block(&with_agent)));
    let path = tmp_path("agents.rxn");
    std::fs::write(&path, &with_agent).unwrap();
    assert!(is_value_error(reaction_from_rxn_file(&path)));
    std::fs::remove_file(&path).unwrap();
    // A zero agent count is no agents.
    let zero = RXN_NO_AGENTS.replace("M  V30 COUNTS 1 1\n", "M  V30 COUNTS 1 1 0\n");
    assert!(reaction_from_rxn_block(&zero).is_ok());
}

// CML.

#[test]
fn cml_tetrahedral_atoms_with_fewer_than_four_neighbors_are_refused() {
    // Three neighbors (an implicit H), and two with a tag set directly.
    let three = mol_from_smiles("F[C@H](Cl)Br").unwrap();
    let mut two = mol_from_smiles("C(F)Cl").unwrap();
    two.get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiTetrahedralCcw);
    for mol in [&three, &two] {
        assert!(is_value_error(mol_to_cml_block(mol)));
        let path = tmp_path("parity.cml");
        assert!(is_value_error(mol_to_cml_file(mol, &path)));
    }
    let cml = mol_to_cml_block(&mol_from_smiles("C[C@](F)(Cl)Br").unwrap()).unwrap();
    assert!(cml.contains("atomParity"));
    assert!(cml.contains("atomRefs4=\"a0 a2 a3 a4\""));
}

// InChI.

#[test]
fn inchi_neighbor_overflows_are_refused() {
    // 21 bonds whose lower-index atom is the dummy: past MAXVAL.
    let crowded = mol_from_smiles(&format!("*{}", "(C)".repeat(21))).unwrap();
    assert!(is_value_error(mol_to_inchi(&crowded)));
    assert!(is_value_error(mol_to_inchi_key(&crowded)));
    // 20 are InChI's to refuse (it does not take '*').
    let full = mol_from_smiles(&format!("*{}", "(C)".repeat(20))).unwrap();
    assert_eq!(mol_to_inchi(&full).unwrap().inchi, "");
    // A chiral tag on an atom of more than four neighbors.
    let octahedral = mol_from_smiles("F[S@OH1](F)(F)(F)(F)F").unwrap();
    assert!(is_value_error(mol_to_inchi(&octahedral)));
    let mut parity = RWMol::from(mol_from_smiles("FS(F)(F)(F)(F)F").unwrap());
    assert!(mol_to_inchi(&parity).is_ok());
    parity
        .get_atom_with_idx_mut(1)
        .unwrap()
        .set_prop("molParity", "1");
    assert!(is_value_error(mol_to_inchi(&parity)));
    // A CW center of two neighbor atoms and two implicit hydrogens.
    let mut chain = mol_from_smiles("CCO").unwrap();
    chain
        .get_atom_with_idx_mut(1)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiTetrahedralCw);
    assert!(is_value_error(mol_to_inchi(&chain)));
    // And a valid one, as from Python.
    let ok = mol_to_inchi(&mol_from_smiles("C[C@](F)(Cl)Br").unwrap()).unwrap();
    assert_eq!(ok.inchi, "InChI=1S/C2H3BrClF/c1-2(3,4)5/h1H3/t2-/m1/s1");
}

// RDKit's text-format readers and the options that take SMARTS parse them
// with CXSMILES on, so text holding a field RDKit could read past the end
// of is refused (or, for a filter matcher, gives no pattern).
#[test]
fn readers_of_smarts_text_refuse_cx_fields() {
    use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
    use rdkit_cxx::chem::rd_deprotect::DeprotectData;
    use rdkit_cxx::chem::rd_fmcs::{MCSParameters, find_mcs_with_params};
    use rdkit_cxx::chem::rd_mol_chemical_features::build_feature_factory_from_string;
    use rdkit_cxx::chem::rd_rascal_mces::{RascalOptions, find_mces_with_opts};
    use rdkit_cxx::chem::rdfiltercatalog::SmartsMatcher;
    use rdkit_cxx::chem::scaffolds::rd_scaffold_network::ScaffoldNetworkParams;

    let params = CleanupParameters::default();
    assert!(is_value_error(normalizer_from_data(
        "n\t[C:1]>>[N:1] |SgD:0:x",
        &params
    )));
    normalizer_from_data("n\t[C:1]>>[N:1]", &params).unwrap();
    let with_data = CleanupParameters {
        normalization_data: vec![("n".into(), "[C:1]>>[N:1] |atomProp:0.p.v".into())],
        ..Default::default()
    };
    assert!(is_value_error(normalizer_from_params(&with_data)));
    // The fragment and acid-base readers erase the spaces from their
    // SMARTS, which leaves no extension: RDKit fails to parse these itself.
    let removed = fragment_remover_from_data("x\tC |Sg:n");
    assert!(rdkit_value_error(&removed));
    fragment_remover_from_data("x\tC").unwrap();
    let reionizer = reionizer_from_data("a\tC(=O)[OH] |SgD:\tC(=O)[O-]");
    assert!(rdkit_value_error(&reionizer));
    assert!(is_value_error(DeprotectData::new(
        "c",
        "[C:1]>>[N:1] |SgD:0:x",
        "a",
        "f"
    )));
    assert!(is_value_error(ScaffoldNetworkParams::new(&[
        "[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2] |atomProp:0.p.v"
    ])));
    assert!(is_value_error(build_feature_factory_from_string(
        "DefineFeature D [N,O;!H0] |SgD:0:x\n  Family Donor\n  Weights 1\nEndFeature\n"
    )));
    let mols = [
        mol_from_smiles("CCO").unwrap(),
        mol_from_smiles("CCN").unwrap(),
    ];
    let mcs = MCSParameters {
        initial_seed: "CC |SgD:0:x".into(),
        ..Default::default()
    };
    assert!(is_value_error(find_mcs_with_params(&mols, &mcs)));
    let rascal = RascalOptions {
        equivalent_atoms: "[N,O]\t|SgD:0:x".into(),
        ..Default::default()
    };
    assert!(is_value_error(find_mces_with_opts(
        &mols[0], &mols[1], &rascal
    )));
    // The same texts without the fields read.
    reionizer_from_data("a\tC(=O)[OH]\tC(=O)[O-]").unwrap();
    DeprotectData::new("c", "[C:1]>>[N:1]", "a", "f").unwrap();
    ScaffoldNetworkParams::new(&["[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]"]).unwrap();
    build_feature_factory_from_string(
        "DefineFeature D [N,O;!H0]\n  Family Donor\n  Weights 1\nEndFeature\n",
    )
    .unwrap();
    let seeded = MCSParameters {
        initial_seed: "CC".into(),
        ..Default::default()
    };
    find_mcs_with_params(&mols, &seeded).unwrap();
    let equivalent = RascalOptions {
        equivalent_atoms: "[N,O]".into(),
        ..Default::default()
    };
    find_mces_with_opts(&mols[0], &mols[1], &equivalent).unwrap();
    assert!(!SmartsMatcher::from_smarts("m", "C |SgD:0:x").is_valid());
    assert!(SmartsMatcher::from_smarts("m", "C").is_valid());
}

const SGROUP_BLOCK: &str = "
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 NSG 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 O 2.598076 -0.000000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
LINES
M  V30 END SGROUP
M  V30 END CTAB
M  END
";

fn sgroup_block(num: &str, lines: &str) -> String {
    SGROUP_BLOCK.replace("NSG", num).replace("LINES\n", lines)
}

// ParseV3000SGroupsBlock reads each S-group line's ids into variables it
// never initializes, which a line without a number, a type, and a number
// leaves unread, as does a block with fewer lines than COUNTS names.
#[test]
fn sgroup_lines_without_ids_are_refused() {
    let ok = sgroup_block("1", "M  V30 1 SUP 0 ATOMS=(1 3) LABEL=OH\n");
    let mol = mol_from_mol_block(&ok).unwrap();
    assert_eq!(mol.get_num_atoms(), 3);
    let with_default = sgroup_block(
        "1",
        "M  V30 DEFAULT CLASS=AA\nM  V30 1 SUP 0 ATOMS=(1 3) LABEL=OH\n",
    );
    mol_from_mol_block(&with_default).unwrap();
    // A third token that is not a number reads as 0.
    mol_from_mol_block(&sgroup_block("1", "M  V30 1 SUP x ATOMS=(1 3)\n")).unwrap();
    for (num, lines) in [
        ("1", "M  V30 1\n"),
        ("1", "M  V30 1 SUP\n"),
        ("1", "M  V30 SUP 1 0\n"),
        ("1", "M  V30 DEFAULT\n"),
        ("2", "M  V30 1 SUP 0 ATOMS=(1 3) LABEL=OH\n"),
    ] {
        assert!(
            is_value_error(mol_from_mol_block(&sgroup_block(num, lines))),
            "{num} {lines}"
        );
    }
    // A cut line is refused too, from an SD file of it.
    let cut = sgroup_block("1", "M  V30 1 SUP 0 ATOMS=(1 3) LABEL=OH\n");
    let cut = &cut[..cut.find("M  V30 1 SUP").unwrap() + 10];
    assert!(mol_from_mol_block(cut).is_err());
}

// Mol2 bond lines.

fn mol2_block(bonds: &str, bond: &str) -> String {
    format!(
        "@<TRIPOS>MOLECULE\nethane\n 2 {bonds} 0 0 0\nSMALL\nNO_CHARGES\n\n\
         @<TRIPOS>ATOM\n\
         \x20     1 C1          0.0000    0.0000    0.0000 C.3     1  MOL       0.0000\n\
         \x20     2 C2          1.5400    0.0000    0.0000 C.3     1  MOL       0.0000\n\
         @<TRIPOS>BOND\n     {bond}\n"
    )
}

// ParseMol2FileBondLine refuses a bond only when both atoms are outside the
// molecule (RDKit #9417), and reads its atom table with each.
#[test]
fn mol2_bond_lines_name_atoms_of_the_molecule() {
    let mol = mol_from_mol2_block(&mol2_block("1", "1 1 2 1")).unwrap();
    assert_eq!(mol.get_num_bonds(), 1);
    let path = tmp_path("ethane.mol2");
    std::fs::write(&path, mol2_block("1", "1 1 2 1")).unwrap();
    assert_eq!(mol_from_mol2_file(&path).unwrap().get_num_bonds(), 1);
    for bond in ["1 1 99 1", "1 0 1 1", "1 1 3 1", "1 -1 2 1"] {
        let block = mol2_block("1", bond);
        assert!(is_value_error(mol_from_mol2_block(&block)), "{bond}");
        std::fs::write(&path, &block).unwrap();
        assert!(is_value_error(mol_from_mol2_file(&path)), "{bond}");
    }
    // Without a BOND header RDKit reads the bond lines from the start.
    let headless = format!(
        "# 1 99 1\n{}",
        mol2_block("1", "").replace("@<TRIPOS>BOND\n", "")
    );
    assert!(is_value_error(mol_from_mol2_block(&headless)));
    // A molecule without bonds reads none.
    mol_from_mol2_block(&mol2_block("0", "1 1 99 1")).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(mol_from_mol2_file(&path), Err(Error::Os { .. })));
}

// Stereo groups that list an atom twice, which RDKit 2024.03 keeps and a
// removal leaves dangling: with its hydrogens removed, Python crashes on
// "[H][C@](C)(F)Cl |&1:0,0,1|".
#[test]
fn stereo_groups_listing_an_atom_twice_are_refused() {
    for text in [
        "C[C@H](F)Cl |&1:1,1|",
        "[H][C@](C)(F)Cl |&1:0,0,1|",
        "C[C@H](F)C[C@H](F)Cl |o1:1,o1:4,o1:1|",
        "C[C@H](F)Cl |a:1,a:1|",
    ] {
        assert!(is_value_error(mol_from_smiles(text)), "{text}");
    }
    assert!(is_value_error(mol_from_smarts("[#6][C@H](F)Cl |&1:1,1|")));
    assert!(is_value_error(reaction_from_smarts(
        "[C:1][C@H](F)Cl>>[C:1] |&1:1,1|"
    )));
    // Groups of different types or numbers merge into none.
    for text in ["C[C@H](F)Cl |&1:1,&2:1|", "C[C@H](F)Cl |a:1,o1:1|"] {
        assert!(mol_from_smiles(text).is_ok(), "{text}");
    }
    let mol = mol_from_smiles("C[C@H](F)Cl |&1:1|").unwrap();
    let block = mol_to_v3k_mol_block(&mol).unwrap();
    assert!(block.contains("ATOMS=(1 2)"));
    assert!(mol_from_mol_block(&block).is_ok());
    let twice = block.replace("ATOMS=(1 2)", "ATOMS=(3 1 1 2)");
    assert!(is_value_error(mol_from_mol_block(&twice)));
    // An SD file of it, read through a checked stream.
    let path = tmp_path("stereo_twice.sdf");
    std::fs::write(&path, format!("{twice}$$$$\n")).unwrap();
    let mut supplier = SDMolSupplier::new(&path).unwrap();
    assert!(is_value_error(supplier.get(0)));
    std::fs::remove_file(&path).unwrap();
    // Marvin groups atoms by the first atom with their id.
    let mrv = mol_to_mrv_block(&mol).unwrap();
    assert!(mol_from_mrv_block(&mrv).is_ok());
    let same_id = mrv.replace(
        "</atomArray>",
        "<atom id=\"a2\" elementType=\"Br\" x2=\"3\" y2=\"3\" mrvStereoGroup=\"and1\"/></atomArray>",
    );
    assert_ne!(same_id, mrv);
    assert!(matches!(
        mol_from_mrv_block(&same_id),
        Err(Error::Parse { .. })
    ));
}

// A V3000 enhanced-stereo collection shorter than its count: RDKit 2024.03
// reads the first atom unset (an uninitialized unsigned), and repeats the
// last for the rest.
#[test]
fn v3000_stereo_collections_shorter_than_their_count_are_refused() {
    let block =
        mol_to_v3k_mol_block(&mol_from_smiles("C[C@H](F)[C@@H](Cl)O |&1:1,3|").unwrap()).unwrap();
    assert!(block.contains("ATOMS=(2 2 4)"));
    for atoms in ["(1 )", "(2 )", "(2 2)", "(3 2 4)"] {
        let bad = block.replace("ATOMS=(2 2 4)", &format!("ATOMS={atoms}"));
        assert!(is_value_error(mol_from_mol_block(&bad)), "{atoms}");
    }
    let fine = block.replace("ATOMS=(2 2 4)", "ATOMS=(2 1 2)");
    mol_from_mol_block(&fine).unwrap();
    // A line RDKit does not read as a stereo group is not checked: outside
    // a collection, or not of a type RDKit reads.
    let outside = block.replace(
        "M  V30 END CTAB",
        "M  V30 BEGIN OBJ3D\nM  V30 MDLV30/STERAC2 ATOMS=(1 )\nM  V30 END OBJ3D\nM  V30 END CTAB",
    );
    assert_ne!(outside, block);
    assert!(mol_from_mol_block(&outside).is_ok());
    // A line cut off by the end of the input is read as a whole one: the
    // line below is not short of atoms (RDKit reads ")" as 0, and refuses
    // that atom itself), but the line cut after its first ")" is.
    let text = sd_record(&block.replace("ATOMS=(2 2 4)", "ATOMS=(2 2) 4)"), 0);
    let at = text.find("ATOMS=(2 2) 4)").unwrap() + "ATOMS=(2 2)".len();
    let mut suppl = SDMolSupplier::default();
    suppl.set_data(&text).unwrap();
    let whole = suppl.get(0);
    assert!(
        matches!(&whole, Err(e) if !e.to_string().contains("enhanced-stereo")),
        "{whole:?}"
    );
    let mut suppl = SDMolSupplier::default();
    let read = suppl.set_data(&text[..at]).and_then(|()| suppl.get(0));
    assert!(is_value_error(read));
}

// RDKit parses an "M  MRV SMA" line's rest as SMARTS, strictly: a CXSMILES
// field is refused only where RDKit would read past the end of the SMARTS,
// and other "M  " lines are not SMARTS.
#[test]
fn mol_block_smarts_are_refused_where_rdkit_reads_past_them() {
    let v2000 = |line: &str| {
        format!(
            "\n  RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    \
             0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n{line}\nM  END\n"
        )
    };
    for sma in [
        "[#6] |atomProp:0.p.v",
        "[#6] |atomProp:0.p.v\r",
        "[#6] |SgD:0:x|",
        "[#6]\t|SgD:0:x",
    ] {
        let block = v2000(&format!("M  MRV SMA   1 {sma}"));
        assert!(is_value_error(mol_from_mol_block(&block)), "{sma:?}");
        // Cut off by the end of the input, too.
        assert!(block.ends_with("\nM  END\n"));
        let cut = &block[..block.len() - 8];
        assert!(is_value_error(mol_from_mol_block(cut)), "{sma:?}");
    }
    for sma in ["[#6] |atomProp:0.p.v|", "[#6] |SgD:0:x::::|"] {
        let mol = mol_from_mol_block(&v2000(&format!("M  MRV SMA   1 {sma}"))).unwrap();
        assert_eq!(mol_to_smarts(&mol).unwrap(), "[#6&$([#6])]", "{sma}");
    }
    // A polymer S-group does not read past the end of strict SMARTS, which
    // RDKit then fails to parse.
    assert!(matches!(
        mol_from_mol_block(&v2000("M  MRV SMA   1 [#6] |Sg:n")),
        Err(Error::Parse { .. })
    ));
    // Only "M  MRV SMA" lines are SMARTS, whose extension starts after a
    // space.
    for line in [
        "M  MRV XYZ   1 [#6] |atomProp:0.p.v",
        "M  MRV SMA   1 [#6]|atomProp:0.p.v",
    ] {
        assert!(!refused_for_cx(mol_from_mol_block(&v2000(line))), "{line}");
    }
    // An SD file cut off inside the line is read as if the line were whole:
    // refused where the cut leaves a field running to the end.
    let text = sd_record(&v2000("M  MRV SMA   1 [#6] |atomProp:0.p.v|"), 0);
    let at = text.find("|atomProp:0.p.v|").unwrap();
    for (end, refused) in [(at + 15, true), (at + 16, false), (text.len(), false)] {
        let mut suppl = SDMolSupplier::default();
        let read = suppl.set_data(&text[..end]).and_then(|()| suppl.get(0));
        assert_eq!(refused_for_cx(read), refused, "{end}");
    }
    for end in at..text.len() {
        let mut suppl = SDMolSupplier::default();
        let _ = suppl.set_data(&text[..end]).and_then(|()| suppl.get(0));
    }
}

// A CXSMILES polymer S-group reads past the end only when strictCXSMILES is
// off.
#[test]
fn polymer_sgroups_are_refused_only_without_strict_cxsmiles() {
    let lax = SmilesParserParams {
        strict_cxsmiles: false,
        ..Default::default()
    };
    assert!(refused_for_cx(mol_from_smiles_with_params(
        "*CC* |Sg:n",
        &lax
    )));
    assert!(matches!(
        mol_from_smiles("*CC* |Sg:n"),
        Err(Error::Parse { .. })
    ));
    let lax = SmartsParserParams {
        strict_cxsmiles: false,
        ..Default::default()
    };
    assert!(refused_for_cx(mol_from_smarts_with_params(
        "*CC* |Sg:n",
        &lax
    )));
    assert!(matches!(
        mol_from_smarts("*CC* |Sg:n"),
        Err(Error::Parse { .. })
    ));
}

// parseQueryDefFile parses its SMARTS column with CXSMILES on.
#[test]
fn query_definitions_are_checked() {
    let fine = "// comment\nAcid\t[CX3](=O)[OX2H1]\nAmine\t[NX3;H2] |atomProp:0.p.v|\n";
    let defs = parse_mol_query_def_file_text(fine).unwrap();
    assert_eq!(defs.len(), 2);
    let bad = format!("{fine}Bad\t[#6] |atomProp:0.p.v\n");
    assert!(refused_for_cx(parse_mol_query_def_file_text(&bad)));
    let path = tmp_path("query_defs.txt");
    std::fs::write(&path, fine).unwrap();
    assert_eq!(parse_mol_query_def_file(&path).unwrap().len(), 2);
    std::fs::write(&path, &bad).unwrap();
    assert!(refused_for_cx(parse_mol_query_def_file(&path)));
    // A commented line, or one of a single column, is not parsed.
    std::fs::write(&path, format!("{fine}//Bad\t[#6] |atomProp:0.p.v\nBad\n")).unwrap();
    assert_eq!(parse_mol_query_def_file(&path).unwrap().len(), 2);
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        parse_mol_query_def_file(&path),
        Err(Error::Os { .. })
    ));
    // addRecursiveQueries parses a value that names no query as SMARTS.
    let mut mol = mol_from_smiles("CC").unwrap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("query", "[#6] |atomProp:0.p.v");
    assert!(refused_for_cx(mol_add_recursive_queries(
        &mut mol,
        &Default::default(),
        "query"
    )));
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("query", "[#6]");
    mol_add_recursive_queries(&mut mol, &Default::default(), "query").unwrap();
}

// The standardizer's files: the normalizations are read once and checked,
// and the fragment, acid-base, and tautomer readers erase spaces from their
// SMARTS, which leaves them no CXSMILES extension.
#[test]
fn standardizer_files_are_read_once_and_checked() {
    use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;

    let path = tmp_path("normalizations.txt");
    std::fs::write(&path, "// comment\nn\t[C:1]>>[N:1]\n").unwrap();
    Normalizer::new(&path, 200).unwrap();
    let params = CleanupParameters {
        normalizations_file: path.clone(),
        ..Default::default()
    };
    normalizer_from_params(&params).unwrap();
    // As the transformations' text (Python's 2024.03 gives "NNO").
    let mut mol = mol_from_smiles("CCO").unwrap();
    let from_data = normalizer_from_data("n\t[C:1]>>[N:1]", &CleanupParameters::default())
        .unwrap()
        .normalize(&mol)
        .unwrap();
    assert_eq!(mol_to_smiles(&from_data).unwrap(), "NNO");
    for normalized in [
        normalize_with_params(&mol, &params).unwrap(),
        Normalizer::new(&path, 200)
            .unwrap()
            .normalize(&mol)
            .unwrap(),
    ] {
        assert_eq!(mol_to_smiles(&normalized).unwrap(), "NNO");
    }
    cleanup_in_place_with_params(&mut mol, &params).unwrap();
    std::fs::write(&path, "n\t[C:1]>>[N:1] |atomProp:0.p.v\n").unwrap();
    assert!(refused_for_cx(Normalizer::new(&path, 200)));
    assert!(refused_for_cx(normalizer_from_params(&params)));
    assert!(refused_for_cx(cleanup_with_params(&mol, &params)));
    // A field RDKit reads within the text is fine.
    std::fs::write(&path, "n\t[C:1]>>[N:1] |SgD:0:x::::|\n").unwrap();
    normalizer_from_params(&params).unwrap();
    // A file of no transformations cannot be handed to RDKit as data.
    std::fs::write(&path, "// none\n").unwrap();
    Normalizer::new(&path, 200).unwrap();
    assert!(is_value_error(normalizer_from_params(&params)));
    // Ops that make no normalizer do not read it.
    let fragments = tmp_path("fragments.txt");
    std::fs::write(&fragments, "x\t[Cl] |atomProp:0.p.v\n").unwrap();
    let params = CleanupParameters {
        normalizations_file: tmp_path("missing.txt"),
        fragment_file: fragments.clone(),
        ..Default::default()
    };
    // RDKit fails to parse "[Cl]|atomProp:0.p.v" itself.
    let salt = mol_from_smiles("CC.Cl").unwrap();
    let removed = remove_fragments_with_params(&salt, &params);
    assert!(rdkit_value_error(&removed));
    assert!(matches!(
        cleanup_with_params(&salt, &params),
        Err(Error::Os { .. })
    ));
    let removed = fragment_remover_from_data("x\t[Cl] |atomProp:0.p.v");
    assert!(rdkit_value_error(&removed));
    std::fs::remove_file(&fragments).unwrap();
    std::fs::remove_file(&path).unwrap();
}
