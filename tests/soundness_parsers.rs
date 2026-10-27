// Copyright (c) 2026 Steven Kearnes
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
        // parse_polymer_sgroup steps to end() + 2.
        "*CC* |Sg:n",
        "*CC* |Sg:n|",
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
    for smarts in [
        "[C:1]>>[C:1]O |atomProp:0.p.v",
        "[C:1]>>[C:1]O |SgD:0:x|",
        "[C:1]>>[C:1]O|Sg:n",
    ] {
        assert!(is_value_error(reaction_from_smarts(smarts)), "{smarts}");
    }
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
