// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! `json_to_mols` against RDKit 2024.09's reader (MolInterchange/Parser.cpp),
//! which reads with rapidjson and checks few of the types it reads: a value
//! of another type, a missing key it reads without checking, or a short
//! array is read as if it were what RDKit wants, so the wrapper refuses it
//! first (`Error::Value`). Where RDKit 2024.09 checks a value itself it
//! still reports it (`Error::Parse`, `Error::Runtime`). Valid documents come
//! from RDKit 2024.09.6's Python `MolToJSON`.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

const HEADER: &str = r#"{"rdkitjson":{"version":11},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":["#;

/// A document of molecules (JSON values).
fn doc(mols: &[String]) -> String {
    format!("{HEADER}{}]}}", mols.join(","))
}

fn mol_doc(mol: &str) -> String {
    doc(&[mol.to_string()])
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

fn butene_doc(extra: &str) -> String {
    doc(&[butene(extra)])
}

/// But-2-ene with a trans double bond of `stereo_atoms`.
fn trans(stereo_atoms: &str) -> String {
    doc(&[butene_with_bonds(
        &format!(
            r#"[{{"atoms":[0,1]}},{{"bo":2,"atoms":[1,2],"stereo":"trans","stereoAtoms":{stereo_atoms}}},{{"atoms":[2,3]}}]"#
        ),
        "",
    )])
}

/// A substance group of atom 1 with `members`.
fn sgroup(members: &str) -> String {
    butene_doc(&format!(
        r#","substanceGroups":[{{"properties":{{"TYPE":"SUP"}},"atoms":[1]{members}}}]"#
    ))
}

fn conformer(conf: &str) -> String {
    butene_doc(&format!(r#","conformers":[{conf}]"#))
}

/// An rdkitRepresentation extension with `members`.
fn representation(members: &str) -> String {
    butene_doc(&format!(
        r#","extensions":[{{"name":"rdkitRepresentation","formatVersion":2{members}}}]"#
    ))
}

/// An rdkitQueries extension.
fn queries(atoms: &str, bonds: &str) -> String {
    format!(
        r#","extensions":[{{"name":"rdkitQueries","formatVersion":10,"atomQueries":[{atoms}],"bondQueries":[{bonds}]}}]"#
    )
}

fn atom_query(query: &str) -> String {
    butene_doc(&queries(query, ""))
}

fn bond_query(query: &str) -> String {
    butene_doc(&queries("", query))
}

const CARBON: &str = r#"{"descr":"AtomAtomicNum","tag":33,"val":6}"#;
const SINGLE: &str = r#"{"descr":"BondOrder","type":"BondOrder","tag":33,"val":1}"#;

fn is_value(json: &str) -> bool {
    matches!(json_to_mols(json), Err(Error::Value { .. }))
}

fn is_parse(json: &str) -> bool {
    matches!(json_to_mols(json), Err(Error::Parse { .. }))
}

fn assert_value(json: &str) {
    assert!(is_value(json), "{json}");
}

fn assert_parse(json: &str) {
    assert!(is_parse(json), "{json}");
}

fn smarts(json: &str) -> Vec<String> {
    json_to_mols(json)
        .unwrap()
        .iter()
        .map(|mol| mol_to_smarts(mol).unwrap())
        .collect()
}

const BUTENE_SMARTS: &str = "[#6H3]-[#6H]=[#6H]-[#6H3]";

#[test]
fn rdkit_documents_still_parse() {
    // Python's MolToJSON of C[C@H](F)C/C=C/[C@@H](O)Cl |o1:1,6| with 2D
    // coordinates.
    let stereo = r#"{"rdkitjson":{"version":11},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"atoms":[{"impHs":3},{"impHs":1,"stereo":"ccw"},{"z":9},{"impHs":2},{"impHs":1},{"impHs":1},{"impHs":1,"stereo":"cw"},{"z":8,"impHs":1},{"z":17}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]},{"atoms":[1,3]},{"atoms":[3,4]},{"bo":2,"atoms":[4,5],"stereo":"trans","stereoAtoms":[3,6]},{"atoms":[5,6]},{"atoms":[6,7]},{"atoms":[6,8]}],"stereoGroups":[{"type":"or","id":0,"atoms":[1,6]}],"conformers":[{"dim":2,"coords":[[-3.8971,-0.75],[-2.5981,0.0],[-2.5981,1.5],[-1.299,-0.75],[0.0,-0.0],[1.299,-0.75],[2.5981,0.0],[3.8971,-0.75],[2.5981,1.5]]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2024.09.6","cipRanks":[0,4,7,1,2,3,5,6,8],"cipCodes":[[1,"S"],[6,"S"]]}]}]}"#;
    let mol = json_to_mols(stereo).unwrap().remove(0);
    assert_eq!(
        mol_to_smarts(&mol).unwrap(),
        "[#6H3]-[#6@H](-[#9])-[#6H2]/[#6H]=[#6H]/[#6@@H](-[#8H])-[#17]"
    );
    assert_eq!(mol.get_num_conformers(), 1);
    assert_eq!(mol.get_stereo_groups()[0].get_atoms(), [1, 6]);

    // Python's MolToJSON of MolFromSmarts("[$(C=O);R1,r{5-6}]-[N;H1,H2]"):
    // a recursive query, AND and OR, ring and range queries.
    let query = r#"{"rdkitjson":{"version":11},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"atoms":[{"z":0},{"z":7}],"bonds":[{"atoms":[0,1]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2024.09.6"},{"name":"rdkitQueries","formatVersion":10,"toolkitVersion":"2024.09.6","atomQueries":[{"descr":"AtomAnd","tag":30,"children":[{"descr":"RecursiveStructure","tag":42,"subquery":{"atoms":[{},{"z":8}],"bonds":[{"bo":2,"atoms":[0,1]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2024.09.6"},{"name":"rdkitQueries","formatVersion":10,"toolkitVersion":"2024.09.6","atomQueries":[{"descr":"AtomType","tag":33,"val":6},{"descr":"AtomType","tag":33,"val":8}],"bondQueries":[{"descr":"BondOrder","type":"BondOrder","tag":33,"val":2}]}]}},{"descr":"AtomOr","tag":31,"children":[{"descr":"AtomInNRings","tag":41,"val":1},{"descr":"range_AtomMinRingSize","tag":38,"lower":5,"upper":6,"ends":0}]}]},{"descr":"AtomAnd","tag":30,"children":[{"descr":"AtomType","tag":33,"val":7},{"descr":"AtomOr","tag":31,"children":[{"descr":"AtomHCount","tag":33,"val":1},{"descr":"AtomHCount","tag":33,"val":2}]}]}],"bondQueries":[{"descr":"BondOrder","type":"BondOrder","tag":33,"val":1}]}]}]}"#;
    assert_eq!(smarts(query), ["[$(C=O);R1,r{5-6}]-[N;H1,H2]"]);

    // Python's MolToJSON of a V3000 mol block's SUP group.
    let sup = r#"{"rdkitjson":{"version":11},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},"molecules":[{"name":"","atoms":[{"impHs":3},{"impHs":2},{"z":8,"impHs":1}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]}],"substanceGroups":[{"properties":{"TYPE":"SUP","index":1,"LABEL":"CO","DATAFIELDS":"[]"},"atoms":[1,2],"bonds":[0],"brackets":[[[0.5,0.0,0.0],[0.5,1.0,0.0],[0.0,0.0,0.0]]],"cstates":[{"bond":0,"vector":[0.5,0.0,0.0]}],"attachPoints":[{"aIdx":1,"lvIdx":0,"id":"1"}]}],"conformers":[{"dim":2,"coords":[[0.0,0.0],[1.0,0.0],[2.0,0.0]]}],"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2024.09.6","cipRanks":[0,1,2]}]}]}"#;
    let mol = json_to_mols(sup).unwrap().remove(0);
    assert_eq!(get_mol_substance_groups(&mol).len(), 1);
    assert_eq!(mol_to_json(&mol).unwrap(), sup);

    // A commonchem document without defaults, whose atoms RDKit reads from
    // an empty DefaultValueCache, and defaults that are not an object, which
    // RDKit's JSON pointer does not look into.
    let methane = r#"[{"atoms":[{"z":6,"impHs":4,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"}],"bonds":[]}]"#;
    for defaults in ["", r#","defaults":[]"#] {
        let json = format!(r#"{{"commonchem":{{"version":10}}{defaults},"molecules":{methane}}}"#);
        assert_eq!(smarts(&json), ["[#6H4]"], "{json}");
    }
}

#[test]
fn headers_are_ints() {
    // DocToMols reads the version with GetInt.
    for header in [
        r#"{"commonchem":{"version":"10"}"#,
        r#"{"commonchem":{"version":10.0}"#,
        r#"{"commonchem":{"version":null}"#,
        r#"{"rdkitjson":{"version":"11"}"#,
        r#"{"rdkitjson":{"version":11.0}"#,
        r#"{"rdkitjson":{"version":4294967307}"#,
    ] {
        assert_value(&format!(r#"{header},"molecules":[]}}"#));
    }
    // RDKit checks the version's value, and the defaults' types.
    assert_parse(r#"{"rdkitjson":{"version":12},"molecules":[]}"#);
    assert_parse(r#"{"commonchem":{"version":10},"defaults":{"atom":[]},"molecules":[]}"#);
}

#[test]
fn molecules_atoms_and_bonds() {
    // processMol, readAtom, and readBond read their members without checking
    // that they have one.
    for mol in ["1", "[]", "null", r#""mol""#] {
        assert_value(&mol_doc(mol));
    }
    assert_value(&butene_doc(r#","name":5"#));
    assert_value(&mol_doc(r#"{"atoms":{},"bonds":[]}"#));
    assert_value(&mol_doc(r#"{"atoms":[],"bonds":0}"#));
    assert_value(&mol_doc(r#"{"atoms":[1],"bonds":[]}"#));
    assert_value(&mol_doc(r#"{"atoms":[null],"bonds":[]}"#));
    assert_value(&doc(&[butene_with_bonds("[[0,1]]", "")]));
    // readBond reads atoms[0] and atoms[1] with GetInt.
    for bond in [
        r#"{"bo":1}"#,
        r#"{"atoms":{"0":0}}"#,
        r#"{"atoms":[]}"#,
        r#"{"atoms":[0]}"#,
        r#"{"atoms":[0,"1"]}"#,
        r#"{"atoms":[0,1.0]}"#,
        r#"{"atoms":[0,2147483648]}"#,
        r#"{"atoms":[-2147483649,1]}"#,
    ] {
        assert_value(&doc(&[butene_with_bonds(&format!("[{bond}]"), "")]));
    }
    // RWMol::addBond checks the atoms.
    for bond in [r#"{"atoms":[0,9]}"#, r#"{"atoms":[0,-1]}"#] {
        assert!(matches!(
            json_to_mols(&doc(&[butene_with_bonds(&format!("[{bond}]"), "")])),
            Err(Error::Runtime { .. })
        ));
    }
    // RDKit checks the missing bonds, and a bond's stereo, before it reads
    // the next bond.
    assert_parse(&mol_doc(r#"{"atoms":[]}"#));
    assert_parse(&doc(&[butene_with_bonds(
        r#"[{"atoms":[0,1],"stereo":1},{"atoms":[1,"x"]}]"#,
        "",
    )]));
}

#[test]
fn stereo_atoms() {
    // readBondStereo reads stereoAtoms[0] and [1] with GetInt.
    for atoms in ["0", "{}", "[]", "[0]", r#"[0,"3"]"#, "[0.0,3]", "[0,null]"] {
        assert_value(&trans(atoms));
    }
    assert_eq!(smarts(&trans("[0,3]")), ["[#6H3]/[#6H]=[#6H]/[#6H3]"]);
    // Bond::setStereoAtoms checks the atoms.
    assert!(matches!(
        json_to_mols(&trans("[0,9]")),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn stereo_groups() {
    // readStereoGroups reads the type with GetString and the id and atoms
    // with GetUint.
    for groups in [
        "{}",
        "1",
        "[1]",
        "[[]]",
        r#"[{"type":1,"atoms":[1]}]"#,
        r#"[{"type":"or","id":-1,"atoms":[1]}]"#,
        r#"[{"type":"or","id":"1","atoms":[1]}]"#,
        r#"[{"type":"and","id":1.0,"atoms":[1]}]"#,
        r#"[{"type":"or","atoms":1}]"#,
        r#"[{"type":"or","atoms":[-1]}]"#,
        r#"[{"type":"or","atoms":["1"]}]"#,
        r#"[{"type":"or","atoms":[1.0]}]"#,
        r#"[{"type":"or","atoms":[4294967296]}]"#,
    ] {
        assert_value(&butene_doc(&format!(r#","stereoGroups":{groups}"#)));
    }
    // RDKit does not read an absolute group's id.
    let mol = json_to_mols(&butene_doc(
        r#","stereoGroups":[{"type":"abs","id":-1,"atoms":[1]}]"#,
    ))
    .unwrap()
    .remove(0);
    assert_eq!(mol.get_stereo_groups()[0].get_atoms(), [1]);
    // RDKit 2024.09 needs a group's atoms, and checks its type and atoms.
    assert_parse(&butene_doc(
        r#","stereoGroups":[{"type":"or","id":1,"bonds":[1]}]"#,
    ));
    assert_parse(&butene_doc(
        r#","stereoGroups":[{"type":"xor","atoms":["x"]}]"#,
    ));
    assert!(
        json_to_mols(&butene_doc(
            r#","stereoGroups":[{"type":"or","atoms":[9]}]"#
        ))
        .is_err()
    );
}

#[test]
fn substance_groups() {
    // readSubstanceGroups reads the TYPE with GetString, indices with
    // GetUint, coordinates with GetDouble, and an attachment point's id with
    // GetString, from members it reads without checking.
    for groups in [
        "{}",
        "[1]",
        r#"[{"properties":[]}]"#,
        r#"[{"properties":"SUP"}]"#,
        r#"[{"properties":{"TYPE":1}}]"#,
    ] {
        assert_value(&butene_doc(&format!(r#","substanceGroups":{groups}"#)));
    }
    for members in [
        r#","atoms":1"#,
        r#","atoms":[-1]"#,
        r#","atoms":["1"]"#,
        r#","bonds":{}"#,
        r#","bonds":[1.0]"#,
        r#","parentAtoms":[-1]"#,
        r#","brackets":{}"#,
        r#","brackets":[1]"#,
        r#","brackets":[[1]]"#,
        r#","brackets":[[[0.5,"0.5",0.5]]]"#,
        r#","brackets":[[[0.5,null,0.5]]]"#,
        // RDKit writes each point into an array of three.
        r#","brackets":[[[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]]"#,
        r#","cstates":{}"#,
        r#","cstates":[1]"#,
        r#","cstates":[{}]"#,
        r#","cstates":[{"bond":"0"}]"#,
        r#","cstates":[{"bond":-1}]"#,
        r#","cstates":[{"bond":3}]"#,
        r#","cstates":[{"bond":0,"vector":1}]"#,
        r#","cstates":[{"bond":0,"vector":[0,"0",0]}]"#,
        r#","attachPoints":{}"#,
        r#","attachPoints":[1]"#,
        r#","attachPoints":[{}]"#,
        r#","attachPoints":[{"aIdx":"1"}]"#,
        r#","attachPoints":[{"aIdx":4}]"#,
        r#","attachPoints":[{"aIdx":1,"lvIdx":"0"}]"#,
        // RDKit reads lvIdx with GetUint; its writer leaves out -1.
        r#","attachPoints":[{"aIdx":1,"lvIdx":-1}]"#,
        r#","attachPoints":[{"aIdx":1,"lvIdx":-2}]"#,
        r#","attachPoints":[{"aIdx":1,"lvIdx":4}]"#,
        r#","attachPoints":[{"aIdx":1,"lvIdx":2147483648}]"#,
        r#","attachPoints":[{"aIdx":1,"id":1}]"#,
    ] {
        assert_value(&sgroup(members));
    }
    // RDKit checks the TYPE's presence, a point's size before it writes it,
    // and the atoms and bonds.
    assert_parse(&butene_doc(r#","substanceGroups":[{"properties":{}}]"#));
    assert_parse(&sgroup(
        r#","brackets":[[[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5]]]"#,
    ));
    assert_parse(&sgroup(r#","cstates":[{"bond":0,"vector":[0,0]}]"#));
    assert!(json_to_mols(&sgroup(r#","bonds":[5]"#)).is_err());

    let mol = json_to_mols(&sgroup(
        r#","bonds":[0],"parentAtoms":[2],"brackets":[[[0.5,0,0],[1,1,0],[0,0,0]]],"cstates":[{"bond":2,"vector":[0,1,0]}],"attachPoints":[{"aIdx":3},{"aIdx":1,"lvIdx":4294967295,"id":"2"}]"#,
    ))
    .unwrap()
    .remove(0);
    let group = &get_mol_substance_groups(&mol)[0];
    assert_eq!(group.get_c_states()[0].bond_idx, 2);
    let attach = group.get_attach_points().unwrap();
    assert_eq!((attach[0].a_idx, attach[0].lv_idx), (3, -1));
    // 4294967295 is -1 in AttachPoint's int.
    assert_eq!((attach[1].a_idx, attach[1].lv_idx), (1, -1));
}

#[test]
fn conformers() {
    let coords = "[[0,0],[1,0],[2,0],[3,0]]";
    // readConformer reads dim with GetInt and coordinates with GetFloat, a
    // cast to float that is undefined beyond float's range.
    for conf in [
        "1".to_string(),
        format!(r#"{{"dim":"2","coords":{coords}}}"#),
        format!(r#"{{"dim":2.0,"coords":{coords}}}"#),
        r#"{"dim":2,"coords":{}}"#.to_string(),
        r#"{"dim":2,"coords":[1,[0,0],[0,0],[0,0]]}"#.to_string(),
        r#"{"dim":2,"coords":[[0,"0"],[1,0],[2,0],[3,0]]}"#.to_string(),
        r#"{"dim":2,"coords":[[0,1e39],[1,0],[2,0],[3,0]]}"#.to_string(),
        r#"{"dim":3,"coords":[[0,0,-1e300],[1,0,0],[2,0,0],[3,0,0]]}"#.to_string(),
        r#"{"dim":2,"coords":[[0,1e400],[1,0],[2,0],[3,0]]}"#.to_string(),
    ] {
        assert_value(&conformer(&conf));
    }
    assert_value(&butene_doc(r#","conformers":{}"#));
    let mol = json_to_mols(&conformer(
        r#"{"dim":2,"coords":[[0,3.3e38],[1,-3.3e38],[2,0],[3,0]]}"#,
    ))
    .unwrap()
    .remove(0);
    assert_eq!(mol.get_num_conformers(), 1);
    // RDKit checks the counts.
    assert_parse(&conformer(r#"{"dim":2,"coords":[[0,0],[1,0],[2,0]]}"#));
    assert_parse(&conformer(r#"{"dim":2,"coords":[[0,0,0],[1,"x"]]}"#));
    assert!(
        json_to_mols(&conformer(
            r#"{"dim":2,"coords":[[0,0],[1,0],[2,0],[3,0],[4,0]]}"#
        ))
        .is_err()
    );
    // Without parse_conformers RDKit reads none.
    let params = JSONParseParameters {
        parse_conformers: false,
        ..Default::default()
    };
    let mol = json_to_mols_with_params(&butene_doc(r#","conformers":{}"#), &params)
        .unwrap()
        .remove(0);
    assert_eq!(mol.get_num_conformers(), 0);
}

#[test]
fn properties() {
    // parseProperties reads them with GetObject, and checks each value's
    // type.
    assert_value(&butene_doc(r#","properties":[]"#));
    let params = JSONParseParameters {
        parse_properties: false,
        ..Default::default()
    };
    assert_eq!(
        json_to_mols_with_params(&butene_doc(r#","properties":[]"#), &params)
            .unwrap()
            .len(),
        1
    );
    let mol = json_to_mols(&butene_doc(
        r#","properties":{"a":1,"b":2.5,"c":"x","d":[1],"e":null,"f":4294967296}"#,
    ))
    .unwrap()
    .remove(0);
    assert_eq!(mol.get_int_prop("a").unwrap(), 1);
    assert!(!mol.has_prop("d"));
    assert!(!mol.has_prop("f"));
}

#[test]
fn extensions() {
    // processMol reads the name with GetString, and each extension its
    // formatVersion with GetInt.
    for extensions in [
        "{}",
        "[1]",
        r#"[{"name":1}]"#,
        r#"[{"name":"rdkitRepresentation","formatVersion":"2"}]"#,
        r#"[{"name":"rdkitRepresentation","formatVersion":2.0}]"#,
        r#"[{"name":"partialCharges","formatVersion":"10"}]"#,
        r#"[{"name":"rdkitQueries","formatVersion":10.0}]"#,
    ] {
        assert_value(&butene_doc(&format!(r#","extensions":{extensions}"#)));
    }
    assert_parse(&butene_doc(r#","extensions":[{"formatVersion":1}]"#));
    // RDKit ignores other extensions, and extensions of later versions.
    for extensions in [
        r#"[{"name":"other","formatVersion":"x"}]"#,
        r#"[{"name":"rdkitRepresentation","formatVersion":3,"atomRings":[[]]}]"#,
    ] {
        assert_eq!(
            smarts(&butene_doc(&format!(r#","extensions":{extensions}"#))),
            [BUTENE_SMARTS]
        );
    }
    // RDKit compares the name as a C string.
    assert_value(&butene_doc(
        r#","extensions":[{"name":"rdkitRepresentation\u0000x","formatVersion":2,"cipCodes":[[1]]}]"#,
    ));
    let mol = json_to_mols(&butene_doc(
        r#","extensions":[{"name":"rdkitRepresentation\u0000x","formatVersion":2,"cipCodes":[[1,"R"]]}]"#,
    ))
    .unwrap()
    .remove(0);
    assert_eq!(
        mol.get_atom_with_idx(1)
            .unwrap()
            .get_prop("_CIPCode")
            .unwrap(),
        "R"
    );
}

#[test]
fn representations() {
    // readRDKitRepresentation reads its lists with GetArray, a CIP code's
    // val[0] and val[1] with GetInt and GetString, and a ring's atoms with
    // GetInt.
    for members in [
        r#","aromaticAtoms":{}"#,
        r#","aromaticBonds":1"#,
        r#","cipRanks":"0""#,
        r#","cipCodes":{}"#,
        r#","cipCodes":[[]]"#,
        r#","cipCodes":[[1]]"#,
        r#","cipCodes":[["1","R"]]"#,
        r#","cipCodes":[[1.0,"R"]]"#,
        r#","cipCodes":[[1,2]]"#,
        r#","atomRings":{}"#,
        r#","atomRings":[[]]"#,
        r#","atomRings":[[0,"1",2]]"#,
        r#","atomRings":[[0,1,2.0]]"#,
    ] {
        assert_value(&representation(members));
    }
    // RDKit checks these types.
    for members in [
        r#","aromaticAtoms":["1"]"#,
        r#","cipCodes":[1]"#,
        r#","atomRings":[1]"#,
    ] {
        assert_parse(&representation(members));
    }
    assert_eq!(
        smarts(&representation(
            r#","aromaticAtoms":[],"cipRanks":[0,1,2,3],"cipCodes":[[1,"R"]]"#
        )),
        [BUTENE_SMARTS]
    );
}

#[test]
fn partial_charges() {
    let charges = |values: &str| {
        butene_doc(&format!(
            r#","extensions":[{{"name":"partialCharges","formatVersion":10,"values":{values}}}]"#
        ))
    };
    // readPartialCharges reads the values with GetArray, and checks their
    // count and types.
    assert_value(&charges("{}"));
    assert_parse(&charges("[0,0.1,0.1,0.1]"));
    assert_parse(&charges("[0.1]"));
    let mol = json_to_mols(&charges("[0.0,0.1,0.1,-0.2]"))
        .unwrap()
        .remove(0);
    assert_eq!(
        mol.get_atom_with_idx(3)
            .unwrap()
            .get_double_prop("_GasteigerCharge")
            .unwrap(),
        -0.2
    );
}

#[test]
fn query_values() {
    // readQueries reads its lists with GetArray, and checks each query is
    // an object with a tag; readQuery and readBaseQuery read the tag and the
    // values with GetInt, finishQuery the description and type with
    // GetString and the negation with GetBool.
    assert_value(&butene_doc(
        r#","extensions":[{"name":"rdkitQueries","formatVersion":10,"atomQueries":{}}]"#,
    ));
    assert_value(&butene_doc(
        r#","extensions":[{"name":"rdkitQueries","formatVersion":10,"bondQueries":1}]"#,
    ));
    let subquery = r#"{"atoms":[{}],"bonds":[]}"#;
    let atom_queries = [
        r#"{"descr":"AtomAtomicNum","tag":"33","val":6}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33.0,"val":6}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33,"val":"6"}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33,"val":6.0}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33,"val":6,"tol":"0"}"#.to_string(),
        r#"{"descr":"range_AtomTotalDegree","tag":38,"upper":2}"#.to_string(),
        r#"{"descr":"range_AtomTotalDegree","tag":38,"lower":1,"upper":"2"}"#.to_string(),
        r#"{"descr":"range_AtomTotalDegree","tag":38,"lower":1,"upper":2,"ends":"0"}"#.to_string(),
        r#"{"descr":"AtomInNRings","tag":41}"#.to_string(),
        r#"{"descr":"AtomInNRings","tag":41,"val":1,"tol":null}"#.to_string(),
        // finishQuery reads the description without checking it is there
        // for a recursive or ring query.
        format!(r#"{{"tag":42,"subquery":{subquery}}}"#),
        r#"{"tag":41,"val":1}"#.to_string(),
        r#"{"descr":5,"tag":33,"val":6}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33,"val":6,"type":1}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":33,"val":6,"negated":1}"#.to_string(),
        r#"{"descr":"AtomAnd","tag":30,"children":{}}"#.to_string(),
        r#"{"descr":"AtomAnd","tag":30,"children":[1,2]}"#.to_string(),
        format!(r#"{{"descr":"AtomAnd","tag":30,"children":[{CARBON},null]}}"#),
        r#"{"descr":"RecursiveStructure","tag":42,"subquery":[]}"#.to_string(),
        r#"{"descr":"RecursiveStructure","tag":42,"subquery":{"atoms":[1],"bonds":[]}}"#
            .to_string(),
        r#"{"descr":"AtomAtomicNum","tag":39,"set":["6"]}"#.to_string(),
        r#"{"descr":"AtomAtomicNum","tag":39,"set":6}"#.to_string(),
        // RDKit 2024.09's finalizeQueryFromDescription casts an AtomRingSize
        // query to EqualityQuery, which a range query is not.
        r#"{"descr":"range_AtomRingSize","tag":38,"lower":5,"upper":6}"#.to_string(),
    ];
    for query in &atom_queries {
        assert_value(&atom_query(query));
        let nested = format!(r#"{{"descr":"AtomAnd","tag":30,"children":[{CARBON},{query}]}}"#);
        assert_value(&atom_query(&nested));
    }
    let bond_queries = [
        r#"{"descr":"BondOrder","tag":"33","val":1}"#.to_string(),
        r#"{"descr":"BondOrder","tag":33}"#.to_string(),
        r#"{"descr":"BondOrder","tag":33,"val":[1]}"#.to_string(),
        r#"{"descr":[],"tag":33,"val":1}"#.to_string(),
        r#"{"descr":"BondOrder","tag":33,"val":1,"negated":"no"}"#.to_string(),
        format!(r#"{{"descr":"BondOr","tag":31,"children":[{SINGLE},1]}}"#),
    ];
    for query in &bond_queries {
        assert_value(&bond_query(query));
    }

    // RDKit checks a query is an object, that there are no more queries
    // than atoms, a child's tag, the description of a base query, and the
    // tag.
    assert_parse(&atom_query("1"));
    assert_parse(&atom_query(&[CARBON; 5].join(",")));
    assert_parse(&atom_query(&format!(
        r#"{{"descr":"AtomAnd","tag":30,"children":[{CARBON},{{}}]}}"#
    )));
    assert_parse(&atom_query(r#"{"tag":33,"val":"x"}"#));
    assert_parse(&atom_query(
        r#"{"descr":"AtomAtomicNum","tag":99,"val":"x"}"#,
    ));
    assert_parse(&bond_query(r#"{"descr":"BondOrder","tag":42}"#));

    // Queries RDKit makes.
    let accepted = [
        (
            r#"{"descr":"less_AtomRingSize","tag":37,"val":5}"#.to_string(),
            "*{5-}-[#6H]=[#6H]-[#6H3]",
        ),
        (
            r#"{"descr":"AtomInNRings","tag":41,"val":1,"tol":0}"#.to_string(),
            "[R1]-[#6H]=[#6H]-[#6H3]",
        ),
        (
            format!(
                r#"{{"descr":"RecursiveStructure","tag":42,"subquery":{subquery},"negated":false,"type":"x"}}"#
            ),
            "[$([#6])]-[#6H]=[#6H]-[#6H3]",
        ),
    ];
    for (query, expected) in &accepted {
        assert_eq!(smarts(&atom_query(query)), [*expected], "{query}");
    }
    assert_eq!(
        smarts(&bond_query(&format!(
            r#"{{"descr":"BondOr","tag":31,"children":[{SINGLE},{{"descr":"BondOrder","tag":33,"val":2}}]}}"#
        ))),
        ["[#6H3]-,=[#6H]=[#6H]-[#6H3]"]
    );
}
