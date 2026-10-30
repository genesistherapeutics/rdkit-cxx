// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! `rdMolInterchange`. Expected values come from RDKit's Python wrappers.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

const HEADER: &str = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":["#;

/// A document of molecules (JSON objects).
fn doc(mols: &[String]) -> String {
    format!("{HEADER}{}]}}", mols.join(","))
}

/// But-2-ene, with `extra` members.
fn butene(extra: &str) -> String {
    butene_with_bonds(
        r#"[{"atoms":[0,1]},{"bo":2,"atoms":[1,2]},{"atoms":[2,3]}]"#,
        extra,
    )
}

fn butene_with_bonds(bonds: &str, extra: &str) -> String {
    format!(
        r#"{{"atoms":[{{"impHs":3}},{{"impHs":1}},{{"impHs":1}},{{"impHs":3}}],"bonds":{bonds}{extra}}}"#
    )
}

/// An rdkitQueries extension.
fn queries(atoms: &str, bonds: &str) -> String {
    format!(
        r#","extensions":[{{"name":"rdkitQueries","formatVersion":10,"atomQueries":[{atoms}],"bondQueries":[{bonds}]}}]"#
    )
}

fn atom_query(query: &str) -> String {
    doc(&[butene(&queries(query, ""))])
}

fn bond_query(query: &str) -> String {
    doc(&[butene(&queries("", query))])
}

fn smarts(json: &str) -> Vec<String> {
    json_to_mols(json)
        .unwrap()
        .iter()
        .map(|mol| mol_to_smarts(mol).unwrap())
        .collect()
}

fn is_value(result: Result<Vec<Mol>, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

const STYRENE: &str = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"atoms":[{"impHs":3},{"impHs":1},{"impHs":1},{},{"impHs":1},{"impHs":1},{"impHs":1},{"impHs":1},{"impHs":1}],"bonds":[{"atoms":[0,1]},{"bo":2,"atoms":[1,2],"stereo":"trans","stereoAtoms":[0,3]},{"atoms":[2,3]},{"bo":2,"atoms":[3,4]},{"atoms":[4,5]},{"bo":2,"atoms":[5,6]},{"atoms":[6,7]},{"bo":2,"atoms":[7,8]},{"atoms":[8,3]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6","aromaticAtoms":[3,4,5,6,7,8],"aromaticBonds":[3,4,5,6,7,8],"cipRanks":[0,1,4,6,5,3,2,3,5],"atomRings":[[3,8,7,6,5,4]]}]}]}"#;

#[test]
fn defaults_match_python() {
    assert_eq!(
        JSONParseParameters::default(),
        JSONParseParameters {
            set_aromatic_bonds: true,
            strict_valence_check: false,
            parse_properties: true,
            parse_conformers: true,
            use_h_counts: true,
        }
    );
    assert_eq!(
        JSONWriteParameters::default(),
        JSONWriteParameters {
            use_rdkit_extensions: true
        }
    );
}

#[test]
fn writes_json() {
    let styrene = mol_from_smiles("C/C=C/c1ccccc1").unwrap();
    assert_eq!(mol_to_json(&styrene).unwrap(), STYRENE);
    let commonchem = mol_to_json_with_params(
        &styrene,
        &JSONWriteParameters {
            use_rdkit_extensions: false,
        },
    )
    .unwrap();
    assert_eq!(
        commonchem,
        STYRENE.replace(
            r#"{"rdkitjson":{"version":12}"#,
            r#"{"commonchem":{"version":10}"#
        )
    );
    let two = [mol_from_smiles("C").unwrap(), mol_from_smiles("O").unwrap()];
    assert_eq!(
        mols_to_json(&two).unwrap(),
        doc(&[
            r#"{"atoms":[{"impHs":4}],"bonds":[],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6"}]}"#.to_string(),
            r#"{"atoms":[{"z":8,"impHs":2}],"bonds":[],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6"}]}"#.to_string(),
        ])
    );
    assert_eq!(mols_to_json(&[]).unwrap(), doc(&[]));
    // RDKit reads atom 0 of every molecule, as Python raises.
    assert!(matches!(
        mol_to_json(&Mol::default()),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        mols_to_json(&[styrene, Mol::default()]),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn reads_json() {
    let mols = json_to_mols(STYRENE).unwrap();
    assert_eq!(mols.len(), 1);
    assert_eq!(mol_to_smiles(&mols[0]).unwrap(), "C/C=C/c1ccccc1");
    assert_eq!(json_to_mols(&doc(&[])).unwrap().len(), 0);

    // Names, stereo groups, conformers, properties, CIP labels, and partial
    // charges, as Python's MolToJSON writes them.
    let everything = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"name":"everything","atoms":[{"impHs":3},{"impHs":1,"stereo":"ccw"},{"z":9},{"impHs":2},{"impHs":1,"stereo":"cw"},{"z":17},{"z":8,"impHs":1}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]},{"atoms":[1,3]},{"atoms":[3,4]},{"atoms":[4,5]},{"atoms":[4,6]}],"stereoGroups":[{"type":"or","id":0,"atoms":[1,4]}],"conformers":[{"dim":2,"coords":[[-2.598076211353316E0,-8.571428571428568E-1],[-1.299038105676658E0,-1.0714285714285698E-1],[-1.299038105676658E0,1.3928571428571432E0],[0E0,-8.571428571428574E-1],[1.299038105676658E0,-1.0714285714285743E-1],[1.2990381056766576E0,1.3928571428571432E0],[2.598076211353316E0,-8.571428571428578E-1]]}],"properties":{"count":7,"weight":1.5E0,"note":"hi"},"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6","cipRanks":[0,2,5,1,3,6,4],"cipCodes":[[1,"S"],[4,"R"]]},{"name":"partialCharges","generator":"RDKit","formatVersion":10,"generatorVersion":"2026.03.6","values":[-3.3468141675467966E-2,1.0097975448960565E-1,-2.4764454643317263E-1,1.838815200751798E-2,1.3027253014536672E-1,-9.194848658283385E-2,-3.776002921559416E-1]}]}]}"#;
    let mol = json_to_mols(everything).unwrap().remove(0);
    assert_eq!(
        mol_to_cx_smiles(&mol).unwrap(),
        "C[C@H](F)C[C@H](O)Cl |(-2.59808,-0.857143,;-1.29904,-0.107143,;-1.29904,1.39286,;0,-0.857143,;1.29904,-0.107143,;2.59808,-0.857143,;1.29904,1.39286,),wD:1.0,wU:4.4,o1:1,4|"
    );
    assert_eq!(mol.get_num_conformers(), 1);
    assert_eq!(mol.get_prop("_Name").unwrap(), "everything");
    assert_eq!(mol.get_int_prop("count").unwrap(), 7);
    assert_eq!(mol.get_double_prop("weight").unwrap(), 1.5);
    assert_eq!(mol.get_prop("note").unwrap(), "hi");
    let charges: Vec<f64> = mol
        .get_atoms()
        .map(|atom| atom.get_double_prop("_GasteigerCharge").unwrap())
        .collect();
    assert!((charges[6] + 0.3776002921559416).abs() < 1e-12);
    let codes: Vec<String> = mol
        .get_atoms()
        .filter_map(|atom| atom.get_prop("_CIPCode").ok())
        .collect();
    assert_eq!(codes, ["S", "R"]);

    let params = JSONParseParameters {
        parse_conformers: false,
        parse_properties: false,
        ..Default::default()
    };
    let bare = json_to_mols_with_params(everything, &params)
        .unwrap()
        .remove(0);
    assert_eq!(bare.get_num_conformers(), 0);
    assert!(!bare.has_prop("count"));
    assert!(bare.has_prop("_Name"));
}

#[test]
fn parse_parameters() {
    let benzene = mol_to_json(&mol_from_smiles("c1ccccc1").unwrap()).unwrap();
    let kekule = json_to_mols_with_params(
        &benzene,
        &JSONParseParameters {
            set_aromatic_bonds: false,
            ..Default::default()
        },
    )
    .unwrap()
    .remove(0);
    let types: Vec<BondType> = kekule.get_bonds().map(|b| b.get_bond_type()).collect();
    assert_eq!(
        types,
        [
            BondType::Single,
            BondType::Double,
            BondType::Single,
            BondType::Double,
            BondType::Single,
            BondType::Double
        ]
    );
    assert!(kekule.get_bonds().all(|b| b.get_is_aromatic()));
    let aromatic = json_to_mols(&benzene).unwrap().remove(0);
    assert!(
        aromatic
            .get_bonds()
            .all(|b| b.get_bond_type() == BondType::Aromatic)
    );

    // A carbon whose JSON gives it no hydrogens.
    let carbon = doc(&[r#"{"atoms":[{}],"bonds":[]}"#.to_string()]);
    let hs = |use_h_counts| {
        let params = JSONParseParameters {
            use_h_counts,
            ..Default::default()
        };
        let mol = json_to_mols_with_params(&carbon, &params)
            .unwrap()
            .remove(0);
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_total_num_hs()
            .unwrap()
    };
    assert_eq!(hs(true), 0);
    assert_eq!(hs(false), 4);

    let methane5 = doc(&[r#"{"atoms":[{"impHs":5}],"bonds":[]}"#.to_string()]);
    assert_eq!(smiles_of(&methane5), ["[CH5]"]);
    let strict = JSONParseParameters {
        strict_valence_check: true,
        ..Default::default()
    };
    assert!(matches!(
        json_to_mols_with_params(&methane5, &strict),
        Err(Error::AtomValence { atom_idx: 0, .. })
    ));
}

fn smiles_of(json: &str) -> Vec<String> {
    json_to_mols(json)
        .unwrap()
        .iter()
        .map(|mol| mol_to_smiles(mol).unwrap())
        .collect()
}

#[test]
fn substance_groups() {
    let json = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"name":"","atoms":[{"impHs":3},{"impHs":2},{"z":8,"impHs":1}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]}],"substanceGroups":[{"properties":{"TYPE":"SUP","index":1,"LABEL":"CO"},"atoms":[1,2],"bonds":[0],"brackets":[[[5E-1,0E0,0E0],[5E-1,1E0,0E0],[0E0,0E0,0E0]]],"cstates":[{"bond":0,"vector":[5E-1,0E0,0E0]}],"attachPoints":[{"aIdx":1,"lvIdx":0,"id":"1"}]}],"conformers":[{"dim":2,"coords":[[0E0,0E0],[1E0,0E0],[2E0,0E0]]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6"}]}]}"#;
    let mol = json_to_mols(json).unwrap().remove(0);
    let groups = get_mol_substance_groups(&mol);
    assert_eq!(groups.len(), 1);
    let group = &groups[0];
    assert_eq!(group.get_prop("TYPE").unwrap(), "SUP");
    assert_eq!(group.get_atoms(), [1, 2]);
    assert_eq!(group.get_bonds(), [0]);
    let cstates = group.get_c_states();
    assert_eq!(cstates.len(), 1);
    assert_eq!(cstates[0].bond_idx, 0);
    let attach = group.get_attach_points().unwrap();
    assert_eq!(
        (attach[0].a_idx, attach[0].lv_idx, attach[0].id.as_str()),
        (1, 0, "1")
    );
    assert_eq!(group.get_brackets().len(), 1);
    // The writer gives back the document.
    assert_eq!(mol_to_json(&mol).unwrap(), json);
}

#[test]
fn attachment_points_without_leaving_atoms() {
    // RDKit's writer leaves out an lvIdx of -1, and RDKit 2026.03's reader
    // then leaves it uninitialized; the wrapper gives it -1.
    let json = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"name":"","atoms":[{"impHs":3},{"impHs":2},{"z":8,"impHs":1}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]}],"substanceGroups":[{"properties":{"TYPE":"SUP","index":1,"LABEL":"CO"},"atoms":[1,2],"attachPoints":[{"aIdx":1,"id":"1"}]}]}]}"#;
    let mol = json_to_mols(json).unwrap().remove(0);
    let groups = get_mol_substance_groups(&mol);
    let attach = groups[0].get_attach_points().unwrap();
    assert_eq!((attach[0].a_idx, attach[0].lv_idx), (1, -1));
    // And so does a round trip.
    let json = mol_to_json(&mol).unwrap();
    assert!(
        json.contains(r#""attachPoints":[{"aIdx":1,"id":"1"}]"#),
        "{json}"
    );
    let mol = json_to_mols(&json).unwrap().remove(0);
    let groups = get_mol_substance_groups(&mol);
    assert_eq!(groups[0].get_attach_points().unwrap()[0].lv_idx, -1);
}

#[test]
fn queries_round_trip() {
    // Python's MolToSmarts(JSONToMols(MolToJSON(MolFromSmarts(s)))[0]).
    let cases = [
        ("[C;R2,#7]-[!#1]", "[C;R2,#7]-[!#1]"),
        ("[$(C=O)]-[N;H1,H2]", "[$(C=O)]-[N;H1,H2]"),
        ("[D{2-3}]", "[D{2-3}]"),
        ("[r{5-}]", "[r{5-}]"),
        ("[x{-2}]", "[x{-2}]"),
        ("C~C", "C~C"),
        ("[#6]-,=[#7]", "[#6]-,=[#7]"),
        ("[!@;-]", "[!*@]"),
        ("*-[!C]", "*-[!C]"),
        ("[c,n;H0]:[c]", "[c,n;H0]:c"),
        ("[CX3](=O)[OX2H1]", "[C&X3](=O)[O&X2&H1]"),
        ("[#6;!$([#6]~[#7])]", "[#6&!$([#6]~[#7])]"),
        ("[+,-]", "[+,-]"),
        ("[2H]", "[#1&2*]"),
        ("[13C@@H](F)Cl", "F[C@@&13*&H1]Cl"),
        ("[k5]", "[k5]"),
        ("[R]@[R]", "[R]@[R]"),
        ("[z2]", "[z2]"),
        ("[Z1]", "[Z1]"),
        ("[^3]", "[^3]"),
        ("[v4]", "[v4]"),
        ("[h1]", "[h1]"),
        ("[X{1-2}]", "[X{1-2}]"),
        ("[$([#6]);$(C=O),$(C#N)]", "[$([#6]);$(C=O),$(C#N)]"),
        ("C=,#C", "C=,#C"),
        ("[N;!H0]", "[N&!H0]"),
        ("[#7&+]", "[#7&+]"),
        ("[a]", "a"),
        ("[A]", "A"),
    ];
    for (query, expected) in cases {
        let json = mol_to_json(&mol_from_smarts(query).unwrap()).unwrap();
        assert_eq!(smarts(&json), [expected], "{query}");
    }

    // Python's MolToJSON of a mol block's atom list, A, and Q atoms, ring
    // bond, substitution, and unsaturation counts, and query bonds.
    let mol_block = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"name":"","atoms":[{"z":0},{"z":0},{"z":0},{},{},{"z":7}],"bonds":[{"atoms":[0,1]},{"bo":0,"atoms":[1,2]},{"bo":0,"atoms":[2,3]},{"bo":0,"atoms":[3,4]},{"bo":0,"atoms":[4,5]}],"conformers":[{"dim":2,"coords":[[0E0,0E0],[1E0,0E0],[2E0,0E0],[3E0,0E0],[4E0,0E0],[5E0,0E0]]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6"},{"name":"rdkitQueries","formatVersion":10,"toolkitVersion":"2026.03.6","atomQueries":[{"descr":"AtomOr","tag":31,"children":[{"descr":"AtomOr","tag":31,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomAtomicNum","tag":33,"val":7}]},{"descr":"AtomAtomicNum","tag":33,"val":8}]},{"descr":"AtomAtomicNum","type":"A","negated":true,"tag":33,"val":1},{"descr":"AtomOr","type":"Q","negated":true,"tag":31,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomAtomicNum","tag":33,"val":1}]},{"descr":"AtomAnd","tag":30,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomRingBondCount","tag":41,"val":2}]},{"descr":"AtomAnd","tag":30,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomExplicitDegree","tag":33,"val":2}]},{"descr":"AtomAnd","tag":30,"children":[{"descr":"AtomAtomicNum","tag":33,"val":7},{"descr":"AtomUnsaturated","tag":33,"val":1}]}],"bondQueries":[{},{},{"descr":"SingleOrDoubleBond","type":"BondOrder","tag":33,"val":1},{"descr":"SingleOrAromaticBond","type":"BondOrder","tag":33,"val":1},{"descr":"BondNull","tag":40}]}]}]}"#;
    assert_eq!(
        smarts(mol_block),
        ["[#6,#7,#8]-,=[!#1][!#6&!#1]~[#6&x2]~[#6&D2]~[#7&$(*=,:,#*)]"]
    );
    // CXSMARTS ring-bond counts, one of them RDKit's "as many as now".
    let cx = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"atoms":[{},{},{"z":7}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6"},{"name":"rdkitQueries","formatVersion":10,"toolkitVersion":"2026.03.6","atomQueries":[{"descr":"AtomAnd","tag":30,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomRingBondCount","tag":41,"val":2}]},{"descr":"AtomAnd","tag":30,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomRingBondCount","tag":41,"val":-559038737}]},{"descr":"AtomAtomicNum","tag":33,"val":7}],"bondQueries":[{"descr":"BondOrder","type":"BondOrder","tag":33,"val":1},{"descr":"BondOrder","type":"BondOrder","tag":33,"val":1}]}]}]}"#;
    assert_eq!(smarts(cx), ["[#6&x2]-[#6&x-559038737]-[#7]"]);
    // An XOR query, which ExpandQuery makes.
    let xor = r#"{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"atoms":[{"z":0},{"impHs":2},{"z":8,"impHs":1}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2026.03.6"},{"name":"rdkitQueries","formatVersion":10,"toolkitVersion":"2026.03.6","atomQueries":[{"descr":"AtomXor","tag":32,"children":[{"descr":"AtomAtomicNum","tag":33,"val":6},{"descr":"AtomIsAromatic","tag":33,"val":1}]},{},{}]}]}]}"#;
    assert_eq!(smarts(xor), ["*-[#6]-[#8]"]);

    // Queries RDKit makes that its SMARTS round trips above do not show.
    let accepted = [
        (
            r#"{"descr":"less_AtomRingBondCount","tag":37,"val":2}"#,
            "[x{2-}]-[#6]=[#6]-[#6]",
        ),
        (
            r#"{"descr":"range_AtomTotalDegree","tag":38,"lower":1,"upper":2,"ends":0}"#,
            "[X{1-2}]-[#6]=[#6]-[#6]",
        ),
        (
            r#"{"descr":"AtomInNRings","tag":41,"val":2}"#,
            "[R2]-[#6]=[#6]-[#6]",
        ),
        (
            r#"{"descr":"AtomAtomicNum","tag":33,"val":1,"negated":true,"type":"A"}"#,
            "[!#1]-[#6]=[#6]-[#6]",
        ),
    ];
    for (query, expected) in accepted {
        assert_eq!(smarts(&atom_query(query)), [expected], "{query}");
    }
    // RDKit ignores queries of a later format.
    let later = doc(&[butene(
        r#","extensions":[{"name":"rdkitQueries","formatVersion":11,"atomQueries":[{"tag":39,"descr":"AtomAtomicNum","set":[6]}]}]"#,
    )]);
    assert_eq!(smarts(&later), ["[#6]-[#6]=[#6]-[#6]"]);
}

#[test]
fn malformed_json() {
    let is_parse = |json: &str| matches!(json_to_mols(json), Err(Error::Parse { .. }));
    assert!(is_parse("not json"));
    assert!(is_parse("[]"));
    assert!(is_parse(r#"{"foo": 1}"#));
    assert!(is_parse(r#"{"rdkitjson":{"version":13},"molecules":[]}"#));
    assert!(is_parse(r#"{"rdkitjson":{"version":"12"},"molecules":[]}"#));
    // Without defaults every atom needs a stereo value.
    assert!(is_parse(
        r#"{"rdkitjson":{"version":12},"molecules":[{"atoms":[{"z":6}],"bonds":[]}]}"#
    ));
    // A value of the wrong type, which boost::json refuses.
    assert!(is_parse(&doc(&[butene_with_bonds(
        r#"[{"atoms":[0,"1"]}]"#,
        ""
    )])));
    // RDKit's own checks (invariants, and indices past the end) are runtime
    // errors, as Python raises.
    let rep = |members: &str| {
        doc(&[butene(&format!(
            r#","extensions":[{{"name":"rdkitRepresentation","formatVersion":2{members}}}]"#
        ))])
    };
    for members in [
        r#","aromaticAtoms":[9]"#,
        r#","aromaticAtoms":[-1]"#,
        r#","aromaticBonds":[9]"#,
        r#","atomRings":[[0,1,9]]"#,
        r#","cipRanks":[1,2,3,4,5]"#,
        r#","cipCodes":[[9,"R"]]"#,
    ] {
        assert!(
            matches!(json_to_mols(&rep(members)), Err(Error::Runtime { .. })),
            "{members}"
        );
    }
    // Out-of-range stereo groups and substance-group atoms RDKit checks.
    assert!(
        json_to_mols(&doc(&[butene(
            r#","stereoGroups":[{"type":"or","id":1,"atoms":[9]}]"#
        )]))
        .is_err()
    );
    assert!(
        json_to_mols(&doc(&[butene(
            r#","substanceGroups":[{"properties":{"TYPE":"SRU"},"atoms":[9]}]"#
        )]))
        .is_err()
    );
}

#[test]
fn stereo_atoms_need_two_atoms() {
    // readBondStereo reads stereoAtoms[0] and [1] without checking.
    for atoms in ["[]", "[0]"] {
        let bonds = format!(
            r#"[{{"atoms":[0,1]}},{{"bo":2,"atoms":[1,2],"stereo":"trans","stereoAtoms":{atoms}}},{{"atoms":[2,3]}}]"#
        );
        assert!(is_value(json_to_mols(&doc(&[butene_with_bonds(
            &bonds, ""
        )]))));
    }
    // RDKit does not read the stereo atoms of a bond without stereo.
    let bonds = r#"[{"atoms":[0,1]},{"bo":2,"atoms":[1,2],"stereoAtoms":[]},{"atoms":[2,3]}]"#;
    assert_eq!(smiles_of(&doc(&[butene_with_bonds(bonds, "")])), ["CC=CC"]);
    let bonds = r#"[{"atoms":[0,1]},{"bo":2,"atoms":[1,2],"stereo":"trans","stereoAtoms":[0,3]},{"atoms":[2,3]}]"#;
    assert_eq!(
        smiles_of(&doc(&[butene_with_bonds(bonds, "")])),
        ["C/C=C/C"]
    );
}

#[test]
fn empty_atom_rings_are_refused() {
    // RDKit reads the ring's element size - 1.
    let json = doc(&[butene(
        r#","extensions":[{"name":"rdkitRepresentation","formatVersion":2,"atomRings":[[]]}]"#,
    )]);
    assert!(is_value(json_to_mols(&json)));
}

#[test]
fn substance_group_checks() {
    let sgroup = |members: &str| {
        doc(&[butene(&format!(
            r#","substanceGroups":[{{"properties":{{"TYPE":"SUP"}},"atoms":[1]{members}}}]"#
        ))])
    };
    // RDKit writes each bracket point into an array of three.
    assert!(is_value(json_to_mols(&sgroup(
        r#","brackets":[[[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]]"#
    ))));
    assert_eq!(
        json_to_mols(&sgroup(
            r#","brackets":[[[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]]"#
        ))
        .unwrap()
        .len(),
        1
    );
    // RDKit does not check these indices; the crate trusts them.
    for members in [
        r#","cstates":[{"bond":3}]"#,
        r#","cstates":[{"bond":-1}]"#,
        r#","attachPoints":[{"aIdx":4}]"#,
        r#","attachPoints":[{"aIdx":1,"lvIdx":4}]"#,
        r#","attachPoints":[{"aIdx":1,"lvIdx":-2}]"#,
    ] {
        assert!(is_value(json_to_mols(&sgroup(members))), "{members}");
    }
    let mol = json_to_mols(&sgroup(
        r#","cstates":[{"bond":2}],"attachPoints":[{"aIdx":3,"lvIdx":-1}]"#,
    ))
    .unwrap()
    .remove(0);
    let group = &get_mol_substance_groups(&mol)[0];
    assert_eq!(group.get_c_states()[0].bond_idx, 2);
    assert_eq!(group.get_attach_points().unwrap()[0].lv_idx, -1);
}

#[test]
fn bond_queries_on_stereo_group_bonds_are_refused() {
    // replaceBond leaves the group pointing at the bond it deletes.
    let bond_order = r#"{"descr":"BondOrder","type":"BondOrder","tag":33,"val":1}"#;
    let groups = r#","stereoGroups":[{"type":"or","id":1,"bonds":[1]}]"#;
    let json = doc(&[butene(&format!(
        "{groups}{}",
        queries("", &format!("{bond_order},{bond_order}"))
    ))]);
    assert!(is_value(json_to_mols(&json)));
    // A query on another bond, or on a group's atoms, is fine: replaceAtom
    // updates the groups.
    let json = doc(&[butene(&format!("{groups}{}", queries("", bond_order)))]);
    let mol = json_to_mols(&json).unwrap().remove(0);
    assert_eq!(mol.get_stereo_groups()[0].get_bonds(), [1]);
    let carbon = r#"{"descr":"AtomAtomicNum","tag":33,"val":6}"#;
    let json = doc(&[butene(&format!(
        r#","stereoGroups":[{{"type":"or","id":1,"atoms":[1,2]}}]{}"#,
        queries(&format!("{carbon},{carbon},{carbon}"), "")
    ))]);
    let mol = json_to_mols(&json).unwrap().remove(0);
    assert_eq!(mol.get_stereo_groups()[0].get_atoms(), [1, 2]);
    assert_eq!(mol_to_smarts(&mol).unwrap(), "[#6]-[#6]=[#6]-[#6]");
}

#[test]
fn queries_rdkit_does_not_make_are_refused() {
    let carbon = r#"{"descr":"AtomAtomicNum","tag":33,"val":6}"#;
    let atom_queries = [
        // AND, OR, and XOR with fewer than two children, which RDKit's
        // writers read.
        r#"{"descr":"AtomAnd","tag":30}"#.to_string(),
        format!(r#"{{"descr":"AtomAnd","tag":30,"children":[{carbon}]}}"#),
        format!(r#"{{"descr":"AtomOr","tag":31,"children":[{carbon}]}}"#),
        format!(r#"{{"descr":"AtomXor","tag":32,"children":[{carbon}]}}"#),
        // Children of a query that is not AND, OR, or XOR.
        format!(r#"{{"descr":"AtomAtomicNum","tag":33,"val":6,"children":[{carbon},{carbon}]}}"#),
        // Descriptions with classes RDKit never gives them.
        r#"{"descr":"AtomAtomicNum","tag":39,"set":[6]}"#.to_string(),
        r#"{"descr":"AtomRingSize","tag":40}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":42,"subquery":{"atoms":[{}],"bonds":[]}}"#.to_string(),
        r#"{"descr":"RecursiveStructure","tag":33,"val":6}"#.to_string(),
        r#"{"descr":"AtomAnd","tag":31,"children":[{"descr":"AtomNull","tag":40},{"descr":"AtomNull","tag":40}]}"#.to_string(),
        r#"{"descr":"AtomNull","tag":33,"val":0}"#.to_string(),
        r#"{"descr":"less_AtomTotalDegree","tag":38,"lower":1,"upper":2}"#.to_string(),
        r#"{"descr":"range_AtomNull","tag":38,"lower":1,"upper":2}"#.to_string(),
        r#"{"descr":"HasProp","tag":40}"#.to_string(),
        r#"{"descr":"Foo","tag":33,"val":1}"#.to_string(),
        // A recursive query without atoms, whose match reads its first.
        r#"{"descr":"RecursiveStructure","tag":42,"subquery":{"atoms":[],"bonds":[]}}"#.to_string(),
        // A recursive query's molecule is checked as any other.
        r#"{"descr":"RecursiveStructure","tag":42,"subquery":{"atoms":[{},{},{}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]},{"atoms":[2,0]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"atomRings":[[]]}]}}"#.to_string(),
        // Values beyond 2^30 (U7), as rdqueries refuses them.
        r#"{"descr":"AtomType","tag":33,"val":-2147483648}"#.to_string(),
        r#"{"descr":"AtomType","tag":33,"val":2147483647}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33,"val":6,"tol":-2147483648}"#.to_string(),
        r#"{"descr":"range_AtomTotalDegree","tag":38,"lower":-2147483648,"upper":2}"#.to_string(),
    ];
    for query in &atom_queries {
        assert!(is_value(json_to_mols(&atom_query(query))), "{query}");
        // Inside an AND query too.
        let nested = format!(r#"{{"descr":"AtomAnd","tag":30,"children":[{carbon},{query}]}}"#);
        assert!(is_value(json_to_mols(&atom_query(&nested))), "{nested}");
    }
    let single = r#"{"descr":"BondOrder","type":"BondOrder","tag":33,"val":1}"#;
    let bond_queries = [
        format!(r#"{{"descr":"BondOr","tag":31,"children":[{single}]}}"#),
        format!(r#"{{"descr":"BondXor","tag":32,"children":[{single},{single}]}}"#),
        // SmartsWrite.cpp casts AND and OR queries' other children to
        // EqualityQuery.
        format!(
            r#"{{"descr":"BondOr","tag":31,"children":[{single},{{"descr":"BondNull","tag":40}}]}}"#
        ),
        r#"{"descr":"BondOrder","tag":40}"#.to_string(),
        r#"{"descr":"BondRingSize","tag":39,"set":[5]}"#.to_string(),
        r#"{"descr":"BondOrder","tag":33,"val":2147483647}"#.to_string(),
    ];
    for query in &bond_queries {
        assert!(is_value(json_to_mols(&bond_query(query))), "{query}");
    }
    // A null bond query alone is fine, as in SMARTS's "~".
    assert_eq!(
        smarts(&bond_query(r#"{"descr":"BondNull","tag":40}"#)),
        ["[#6]~[#6]=[#6]-[#6]"]
    );
}
