// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_abbreviations::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdAbbreviations.

fn smiles(smi: &str) -> Mol {
    mol_from_smiles(smi).unwrap()
}

fn cx(mol: &Mol) -> String {
    mol_to_cx_smiles(mol).unwrap()
}

// A definition from a SMARTS as written, with the _numDummies property the
// parsers set.
fn definition(label: &str, smarts: &str) -> AbbreviationDefinition {
    let mut mol = mol_from_smarts(smarts).unwrap();
    mol.set_unsigned_prop("_numDummies", 1);
    AbbreviationDefinition {
        label: label.to_owned(),
        mol: Some(mol),
        ..Default::default()
    }
}

const BENZOATE: &str = "CCCOC(=O)c1ccccc1";
const ALL: CondenseMolAbbreviationsParams = CondenseMolAbbreviationsParams {
    max_coverage: 1.0,
    sanitize: true,
};

#[test]
fn default_definition() {
    let def = AbbreviationDefinition::default();
    assert_eq!(def.label, "");
    assert_eq!(def.display_label, "");
    assert_eq!(def.display_label_w, "");
    assert!(def.mol.is_none());
}

#[test]
fn param_defaults() {
    assert_eq!(
        ParseAbbreviationsParams::default(),
        ParseAbbreviationsParams {
            remove_extra_dummies: false,
            allow_connection_to_dummies: false,
        }
    );
    assert_eq!(
        CondenseMolAbbreviationsParams::default(),
        CondenseMolAbbreviationsParams {
            max_coverage: 0.4,
            sanitize: true,
        }
    );
}

#[test]
fn defaults() {
    let abbrevs = get_default_abbreviations();
    assert_eq!(abbrevs.len(), 37);
    let labels: Vec<&str> = abbrevs.iter().take(5).map(|a| a.label.as_str()).collect();
    assert_eq!(labels, ["CO2Et", "COOEt", "OiBu", "nDec", "nNon"]);
    let first = &abbrevs[0];
    assert_eq!(first.display_label, "CO<sub>2</sub>Et");
    assert_eq!(first.display_label_w, "EtO<sub>2</sub>C");
    let mol = first.mol.as_ref().unwrap();
    assert_eq!(
        mol_to_smarts(mol).unwrap(),
        "[!#0][C&D3&R0](=[O&D1&R0])[O&D2&R0][C&D2&R0][C&D1&R0]"
    );
    assert_eq!(mol.get_unsigned_prop("_numDummies").unwrap(), 1);

    let linkers = get_default_linkers();
    let labels: Vec<&str> = linkers.iter().map(|a| a.label.as_str()).collect();
    assert_eq!(
        labels,
        ["PEG6", "PEG5", "PEG4", "PEG3", "Dec", "Non", "Oct", "Hept"]
    );
    let mol = linkers[0].mol.as_ref().unwrap();
    assert!(
        mol_to_smarts(mol)
            .unwrap()
            .starts_with("*[O&D2&R0][C&D2&R0]")
    );
    assert_eq!(mol.get_unsigned_prop("_numDummies").unwrap(), 1);
}

#[test]
fn parse() {
    // The unparsable SMARTS and the lines without one are skipped.
    let defs = parse_abbreviations("Me C\nbad [C\n   \nonly\nEt CC Et EtW\n").unwrap();
    let got: Vec<(&str, &str, &str, String)> = defs
        .iter()
        .map(|d| {
            (
                d.label.as_str(),
                d.display_label.as_str(),
                d.display_label_w.as_str(),
                mol_to_smarts(d.mol.as_ref().unwrap()).unwrap(),
            )
        })
        .collect();
    assert_eq!(
        got,
        [
            ("Me", "", "", "[!#0][C&D1&R0]".to_owned()),
            ("Et", "Et", "EtW", "[!#0][C&D2&R0][C&D1&R0]".to_owned()),
        ]
    );

    let linkers = parse_linkers("PEG2 *OCCOCC*\n").unwrap();
    assert_eq!(linkers.len(), 1);
    let mol = linkers[0].mol.as_ref().unwrap();
    assert_eq!(
        mol_to_smarts(mol).unwrap(),
        "*[O&D2&R0][C&D2&R0][C&D2&R0][O&D2&R0][C&D2&R0][C&D2&R0]"
    );
    assert_eq!(mol.get_unsigned_prop("_numDummies").unwrap(), 1);

    let params = ParseAbbreviationsParams {
        remove_extra_dummies: true,
        allow_connection_to_dummies: true,
    };
    let defs = parse_abbreviations_with_params("Et CC", &params).unwrap();
    assert_eq!(
        mol_to_smarts(defs[0].mol.as_ref().unwrap()).unwrap(),
        "*[C&D2&R0][C&D1&R0]"
    );
}

// RDKit reads a lone label's SMARTS past the end of the line's fields, where
// Python's binding then reads the label again.
#[test]
fn parse_skips_lines_without_smarts() {
    assert!(parse_abbreviations("CC").unwrap().is_empty());
    assert_eq!(parse_abbreviations("CC\n\t \nEt CC").unwrap().len(), 1);
}

#[test]
fn condense() {
    let abbrevs = get_default_abbreviations();
    // CF3 covers 40% of the atoms.
    let mol = smiles("c1ccccc1C(F)(F)F");
    assert_eq!(
        cx(&condense_mol_abbreviations(&mol, &abbrevs).unwrap()),
        "FC(F)(F)c1ccccc1"
    );

    let mol = smiles(BENZOATE);
    let condensed = condense_mol_abbreviations_with_params(&mol, &abbrevs, &ALL).unwrap();
    assert_eq!(cx(&condensed), "*OC(=O)c1ccccc1 |$nPr;;;;;;;;;$|");
    let dummy = condensed.get_atom_with_idx(0).unwrap();
    assert_eq!(dummy.get_atomic_num(), 0);
    assert_eq!(dummy.get_prop("atomLabel").unwrap(), "nPr");
    assert_eq!(dummy.get_prop("_displayLabel").unwrap(), "nPr");
    assert_eq!(
        cx(&condense_mol_abbreviations(&mol, &abbrevs).unwrap()),
        "*OC(=O)c1ccccc1 |$nPr;;;;;;;;;$|"
    );
    let unsanitized = CondenseMolAbbreviationsParams {
        sanitize: false,
        ..ALL
    };
    assert_eq!(
        cx(&condense_mol_abbreviations_with_params(&mol, &abbrevs, &unsanitized).unwrap()),
        "*OC(=O)c1ccccc1 |$nPr;;;;;;;;;$|"
    );
    // The input is unchanged.
    assert_eq!(cx(&mol), "CCCOC(=O)c1ccccc1");

    let peg = smiles("COCCOCCOCCOCCOCCOCCOC");
    assert_eq!(
        cx(&condense_mol_abbreviations_with_params(&peg, &get_default_linkers(), &ALL).unwrap()),
        "C*OC |$;PEG6;;$|"
    );
    assert_eq!(
        cx(&condense_mol_abbreviations(&smiles(""), &abbrevs).unwrap()),
        ""
    );
}

#[test]
fn label_and_condense_groups() {
    let mol = smiles(BENZOATE);
    let labelled =
        label_mol_abbreviations_with_max_coverage(&mol, &get_default_abbreviations(), 1.0).unwrap();
    let groups = get_mol_substance_groups(&labelled);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].get_prop("TYPE").unwrap(), "SUP");
    assert_eq!(groups[0].get_prop("LABEL").unwrap(), "nPr");
    assert_eq!(groups[0].get_atoms(), [2, 1, 0]);
    assert_eq!(groups[0].get_bonds(), [2]);
    // The propyl group covers 30% of the atoms.
    let labelled_default = label_mol_abbreviations(&mol, &get_default_abbreviations()).unwrap();
    assert_eq!(get_mol_substance_groups(&labelled_default).len(), 1);

    let condensed = condense_abbreviation_substance_groups(&labelled).unwrap();
    assert_eq!(cx(&condensed), "*OC(=O)c1ccccc1 |$nPr;;;;;;;;;$|");
    // Without SUP groups, RDKit warns and leaves the molecule as it is.
    assert_eq!(
        cx(&condense_abbreviation_substance_groups(&smiles("CC")).unwrap()),
        "CC"
    );
}

#[test]
fn definitions_without_mols_or_dummy_counts() {
    let mol = smiles("c1ccccc1C(F)(F)F");
    let err = condense_mol_abbreviations(&mol, &[AbbreviationDefinition::default()]).unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
    let err = label_mol_abbreviations(&mol, &[AbbreviationDefinition::default()]).unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");

    let mut def = definition("X", "*C(F)(F)F");
    let mol_without_count = mol_from_smarts("*C(F)(F)F").unwrap();
    def.mol = Some(mol_without_count);
    let err = condense_mol_abbreviations(&mol, std::slice::from_ref(&def)).unwrap_err();
    assert!(matches!(err, Error::Key { .. }), "{err:?}");
    let any = CondenseMolAbbreviationsParams {
        max_coverage: 0.0,
        sanitize: true,
    };
    assert_eq!(
        cx(&condense_mol_abbreviations_with_params(&mol, &[def], &any).unwrap()),
        "*c1ccccc1 |$X;;;;;;$|"
    );
}

// The first abbreviation bonds C3 to its connecting atom, C1, since C2 has a
// bond its query atom lacks; the second then removes C3, and RDKit marks the
// new bond in a bitset sized for the molecule's original bonds (Python raises
// an invariant only after that).
#[test]
fn refuses_abbreviations_bonded_into_later_matches() {
    let mol = smiles("NCCCCO");
    let abbrevs = [definition("X", "[NX3]CC"), definition("Y", "[OX2]CC")];
    let err = condense_mol_abbreviations_with_params(&mol, &abbrevs, &ALL).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    // Either alone is fine.
    assert_eq!(
        cx(&condense_mol_abbreviations_with_params(&mol, &abbrevs[..1], &ALL).unwrap()),
        "N*CCO |$;X;;;$|"
    );
}

// Group A attaches C3 to its connecting atom, C1; group B then removes C3.
#[test]
fn refuses_sup_groups_bonded_into_later_groups() {
    let mut mol = smiles("CCCCCC");
    let bond =
        |mol: &Mol, a: u32, b: u32| mol.get_bond_between_atoms(a, b).unwrap().unwrap().get_idx();
    let (b01, b23, b45) = (bond(&mol, 0, 1), bond(&mol, 2, 3), bond(&mol, 4, 5));
    {
        let mut a = create_mol_substance_group(&mut mol, "SUP");
        a.set_prop("LABEL", "A");
        a.set_atoms(&[1, 2]).unwrap();
        a.set_bonds(&[b01, b23]).unwrap();
    }
    {
        let mut b = create_mol_substance_group(&mut mol, "SUP");
        b.set_prop("LABEL", "B");
        b.set_atoms(&[4, 3]).unwrap();
        b.set_bonds(&[b45]).unwrap();
    }
    let err = condense_abbreviation_substance_groups(&mol).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}

// Removing atoms rewrites every bond's _MolFileBondEndPts, which RDKit reads
// without checking it.
#[test]
fn refuses_bad_end_points() {
    let abbrevs = get_default_abbreviations();
    for (value, ok) in [
        ("", false),
        ("()", false),
        ("( )", false),
        ("(1 -1)", false),
        ("(1 1e10)", false),
        ("(1 nan)", false),
        ("(1 x)", false),
        ("(1 7)", true),
        ("not a list", true),
    ] {
        let mut mol = smiles(BENZOATE);
        mol.get_bond_with_idx_mut(7)
            .unwrap()
            .set_prop("_MolFileBondEndPts", value);
        let res = condense_mol_abbreviations_with_params(&mol, &abbrevs, &ALL);
        assert_eq!(res.is_ok(), ok, "{value:?}: {res:?}");
        if let Err(err) = res {
            assert!(matches!(err, Error::Value { .. }), "{err:?}");
        }
        let mut labelled = label_mol_abbreviations_with_max_coverage(&mol, &abbrevs, 1.0).unwrap();
        labelled
            .get_bond_with_idx_mut(7)
            .unwrap()
            .set_prop("_MolFileBondEndPts", value);
        assert_eq!(
            condense_abbreviation_substance_groups(&labelled).is_ok(),
            ok,
            "{value:?}"
        );
    }
}

#[test]
fn threads() {
    let abbrevs = get_default_abbreviations();
    std::thread::scope(|s| {
        for _ in 0..4 {
            // Definitions hold molecules, which are not Sync (S1).
            let abbrevs = abbrevs.clone();
            s.spawn(move || {
                let mol = smiles(BENZOATE);
                for _ in 0..20 {
                    let res = condense_mol_abbreviations_with_params(&mol, &abbrevs, &ALL).unwrap();
                    assert_eq!(cx(&res), "*OC(=O)c1ccccc1 |$nPr;;;;;;;;;$|");
                }
            });
        }
    });
}

// Removing atom 0 drops the group that holds it, and RDKit then follows the
// other groups' PARENT chains, which here form a cycle it would follow
// forever.
#[test]
fn refuses_parent_cycles() {
    let abbrevs = get_default_abbreviations();
    for (parents, ok) in [([2, 1], false), ([1, 1], false), ([3, 3], true)] {
        let mut mol = smiles(BENZOATE);
        for (i, (atom, parent)) in [(8, Some(parents[0])), (9, Some(parents[1])), (0, None)]
            .into_iter()
            .enumerate()
        {
            let mut group = create_mol_substance_group(&mut mol, "DAT");
            group.add_atom_with_idx(atom).unwrap();
            group.set_unsigned_prop("index", i as u32 + 1);
            if let Some(parent) = parent {
                group.set_unsigned_prop("PARENT", parent);
            }
        }
        let res = condense_mol_abbreviations_with_params(&mol, &abbrevs, &ALL);
        assert_eq!(res.is_ok(), ok, "{parents:?}: {res:?}");
        if let Err(err) = res {
            assert!(matches!(err, Error::Value { .. }), "{err:?}");
        }
    }
}
