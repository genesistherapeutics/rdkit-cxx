// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_generalized_substruct::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::all_probe_bits_match;

// Expected values come from Python's rdGeneralizedSubstruct.

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn xqm_type(query: &ExtendedQueryMol) -> String {
    let json = query.to_json().unwrap();
    let start = json.find("\"xqm_type\": \"").unwrap() + 13;
    json[start..start + 1].to_string()
}

fn num_matches(mol: &str, query: &ExtendedQueryMol) -> usize {
    mol_get_substruct_matches(&smiles(mol), query)
        .unwrap()
        .len()
}

#[test]
fn link_node_queries() {
    let query = create_extended_query_mol(&smiles("COCc1n[nH]c(F)c1 |LN:1:1.3|")).unwrap();
    // A bundle of tautomer queries, one per link-node repeat.
    assert_eq!(xqm_type(&query), "4");
    let mol = smiles("COOCc1n[nH]c(F)c1");
    let expected = [0, 1, 2, 3, 4, 5, 9, 6, 7, 8];
    assert_eq!(
        mol_get_substruct_matches(&mol, &query).unwrap(),
        [expected.to_vec()]
    );
    assert_eq!(mol_get_substruct_match(&mol, &query).unwrap(), expected);
    assert!(mol_has_substruct_match(&mol, &query).unwrap());
    assert!(!mol_has_substruct_match(&smiles("CCCC"), &query).unwrap());
    assert!(
        mol_get_substruct_match(&smiles("CCCC"), &query)
            .unwrap()
            .is_empty()
    );

    let fp = query.pattern_fingerprint_query().unwrap();
    assert_eq!(fp.get_num_bits(), 2048);
    assert_eq!(fp.get_num_on_bits(), 135);
    assert_eq!(
        query
            .pattern_fingerprint_query_with_fingerprint_size(1024)
            .unwrap()
            .get_num_on_bits(),
        132
    );
    let target = pattern_fingerprint_target(&mol).unwrap();
    assert_eq!(target.get_num_on_bits(), 187);
    assert_eq!(
        pattern_fingerprint_target_with_fingerprint_size(&mol, 1024)
            .unwrap()
            .get_num_on_bits(),
        180
    );
    assert!(all_probe_bits_match(&fp, &target).unwrap());
    assert!(matches!(
        query.pattern_fingerprint_query_with_fingerprint_size(1 << 31),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        pattern_fingerprint_target_with_fingerprint_size(&mol, u32::MAX),
        Err(Error::Value { .. })
    ));
}

#[test]
fn match_parameters() {
    let query =
        create_extended_query_mol(&smiles("COCc1n[nH]c(F)c1[C@H](F)Cl |LN:1:1.3|")).unwrap();
    let chiral = SubstructMatchParameters {
        use_chirality: true,
        ..Default::default()
    };
    let same = smiles("COOCc1n[nH]c(F)c1[C@H](F)Cl");
    let inverted = smiles("COOCc1n[nH]c(F)c1[C@@H](F)Cl");
    assert_eq!(mol_get_substruct_matches(&same, &query).unwrap().len(), 1);
    assert_eq!(
        mol_get_substruct_matches(&inverted, &query).unwrap().len(),
        1
    );
    assert_eq!(
        mol_get_substruct_matches_with_params(&same, &query, &chiral)
            .unwrap()
            .len(),
        1
    );
    assert!(
        mol_get_substruct_matches_with_params(&inverted, &query, &chiral)
            .unwrap()
            .is_empty()
    );
    assert!(!mol_has_substruct_match_with_params(&inverted, &query, &chiral).unwrap());
    assert!(
        mol_get_substruct_match_with_params(&inverted, &query, &chiral)
            .unwrap()
            .is_empty()
    );

    let bundle = create_extended_query_mol(&smiles("COC |LN:1:1.3|")).unwrap();
    assert_eq!(xqm_type(&bundle), "2");
    let all = SubstructMatchParameters {
        uniquify: false,
        ..Default::default()
    };
    assert_eq!(
        mol_get_substruct_matches_with_params(&smiles("COOCC"), &bundle, &all).unwrap(),
        [vec![0, 1, 2, 3], vec![3, 2, 1, 0]]
    );
    assert_eq!(
        mol_get_substruct_matches_with_params(&smiles("COCC"), &bundle, &all).unwrap(),
        [vec![0, 1, 2], vec![2, 1, 0]]
    );
}

#[test]
fn creation_steps() {
    let mol = smiles("COCC1OC(N)=N1 |LN:1:1.3|");
    let step = |do_enumeration, do_tautomers| {
        create_extended_query_mol_with_params(
            &mol,
            &CreateExtendedQueryMolParams {
                do_enumeration,
                do_tautomers,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let queries = [
        step(true, true),
        step(false, true),
        step(true, false),
        step(false, false),
    ];
    let types: Vec<_> = queries.iter().map(xqm_type).collect();
    assert_eq!(types, ["4", "3", "2", "1"]);
    for (target, expected) in [
        ("COCC1OC(N)=N1", [1, 1, 1, 1]),
        ("COCC1OC(=N)N1", [1, 1, 0, 0]),
        ("COOCC1OC(N)=N1", [1, 0, 1, 0]),
    ] {
        let found: Vec<_> = queries.iter().map(|q| num_matches(target, q)).collect();
        assert_eq!(found, expected, "{target}");
    }

    let mol = smiles("COCc1n[nH]cc1 |LN:1:1.3|");
    let adjusted = |params: AdjustQueryParameters| {
        create_extended_query_mol_with_params(
            &mol,
            &CreateExtendedQueryMolParams {
                adjust_query_properties: true,
                adjust_query_parameters: params,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let queries = [
        create_extended_query_mol(&mol).unwrap(),
        adjusted(AdjustQueryParameters::default()),
        adjusted(AdjustQueryParameters::no_adjustments()),
    ];
    for (target, expected) in [
        ("COCc1n[nH]cc1", [1, 1, 1]),
        ("COCc1n[nH]c(F)c1", [1, 0, 1]),
    ] {
        let found: Vec<_> = queries.iter().map(|q| num_matches(target, q)).collect();
        assert_eq!(found, expected, "{target}");
    }
}

#[test]
fn defaults() {
    let params = CreateExtendedQueryMolParams::default();
    assert!(params.do_enumeration);
    assert!(params.do_tautomers);
    assert!(!params.adjust_query_properties);
    assert_eq!(
        params.adjust_query_parameters,
        AdjustQueryParameters::default()
    );
}

#[test]
fn pickles_and_json() {
    let mol = smiles("COCC1OC(N)=N1 |LN:1:1.3|");
    let target = smiles("COCC1OC(=N)N1");
    for (do_enumeration, do_tautomers) in
        [(true, true), (false, true), (true, false), (false, false)]
    {
        let query = create_extended_query_mol_with_params(
            &mol,
            &CreateExtendedQueryMolParams {
                do_enumeration,
                do_tautomers,
                ..Default::default()
            },
        )
        .unwrap();
        let expected = mol_has_substruct_match(&target, &query).unwrap();
        let pickle = query.to_binary().unwrap();
        let json = query.to_json().unwrap();
        let copies = [
            // SAFETY: a pickle RDKit wrote.
            unsafe { ExtendedQueryMol::from_data(&pickle) }.unwrap(),
            // SAFETY: as above.
            unsafe { ExtendedQueryMol::new_with_is_json(&json, true) }.unwrap(),
            query.clone(),
        ];
        for copy in &copies {
            assert_eq!(xqm_type(copy), xqm_type(&query));
            assert_eq!(mol_has_substruct_match(&target, copy).unwrap(), expected);
        }
        let mut replaced = create_extended_query_mol(&smiles("CC")).unwrap();
        // SAFETY: a pickle RDKit wrote.
        unsafe { replaced.init_from_binary(&pickle) }.unwrap();
        assert_eq!(xqm_type(&replaced), xqm_type(&query));
        let mut replaced = create_extended_query_mol(&smiles("CC")).unwrap();
        // SAFETY: a pickle RDKit wrote.
        unsafe { replaced.init_from_json(&json) }.unwrap();
        assert_eq!(xqm_type(&replaced), xqm_type(&query));
    }

    let plain = create_extended_query_mol(&smiles("CC")).unwrap();
    assert_eq!(xqm_type(&plain), "1");
    // A failed read leaves the query as it was.
    let mut query = plain.clone();
    // SAFETY: not a pickle; RDKit refuses it before reading a molecule.
    assert!(unsafe { query.init_from_binary(b"junk") }.is_err());
    // SAFETY: not JSON; RDKit refuses it before reading a molecule.
    assert!(unsafe { query.init_from_json("{") }.is_err());
    assert_eq!(xqm_type(&query), "1");
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { ExtendedQueryMol::new_with_is_json("{\"xqm_type\": \"9\"}", true) },
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { ExtendedQueryMol::new_with_is_json("not json", true) },
        Err(Error::Value { .. })
    ));
    // The string form of the binary constructor reads its bytes.
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { ExtendedQueryMol::new("not a pickle") }.is_err());
}

#[test]
fn json_is_checked() {
    let json = |mol: &str| format!("{{\"xqm_type\": \"1\", \"mol\": {{{mol}}}}}");
    // SAFETY: a pickle RDKit wrote.
    let read = |text: &str| unsafe { ExtendedQueryMol::new_with_is_json(text, true) };
    // RDKit holds, then dereferences, the null molecule SmilesToMol and
    // SmartsToMol return for text that does not parse.
    assert!(matches!(
        read(&json("\"smiles\": \"C1CC\"")),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        read(&json("\"smarts\": \"[C\"")),
        Err(Error::Parse { .. })
    ));
    let ok = read(&json("\"smarts\": \"[#6]~[#7]\"")).unwrap();
    assert!(mol_has_substruct_match(&smiles("CN"), &ok).unwrap());
    assert!(read(&json("\"smiles\": \"CN\"")).is_ok());

    // RDKit's base64 decoder writes before its buffer for text under two
    // characters, and past it for some unpadded lengths.
    for bad in ["A", "AAA", "AAAAAAA"] {
        assert!(
            matches!(
                read(&json(&format!("\"pkl\": \"{bad}\""))),
                Err(Error::Value { .. })
            ),
            "{bad}"
        );
    }
    let pkl = rdkit_cxx::chem::Mol::to_binary(&smiles("CC"));
    let encoded = base64(&pkl);
    let from_pkl = read(&json(&format!("\"pkl\": \"{encoded}\""))).unwrap();
    assert!(mol_has_substruct_match(&smiles("CCC"), &from_pkl).unwrap());
    // A pickle older than format 4.0, which RDKit would give stereochemistry
    // as it reads it, before its relevant cycles can be counted.
    for version in common::OLD_PICKLE_VERSIONS {
        let old = base64(&common::with_pickle_version(&pkl, 0, version));
        assert!(
            matches!(
                read(&json(&format!("\"pkl\": \"{old}\""))),
                Err(Error::Value { message }) if message.contains("relevant cycles")
            ),
            "{version:?}"
        );
    }

    // An empty bundle has no pattern fingerprint (RDKit returns null).
    assert!(matches!(
        read("{\"xqm_type\": \"2\", \"mols\": []}"),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        read("{\"xqm_type\": \"4\", \"tautomerQueries\": []}"),
        Err(Error::Value { .. })
    ));
    // The matcher reads each tautomer at the template's modified atoms.
    let tautomer_query = |tautomer: &str, atoms: &str| {
        format!(
            "{{\"xqm_type\": \"3\", \"tautomers\": [{{\"smiles\": \"{tautomer}\"}}], \
             \"template\": {{\"smiles\": \"CC=O\"}}, \"modifiedAtoms\": [{atoms}], \
             \"modifiedBonds\": []}}"
        )
    };
    assert!(read(&tautomer_query("CC=O", "\"0\", \"2\"")).is_ok());
    for (tautomer, atoms) in [("C", "\"0\""), ("CC=O", "\"3\""), ("CC=O", "\"-1\"")] {
        assert!(
            matches!(
                read(&tautomer_query(tautomer, atoms)),
                Err(Error::Value { .. })
            ),
            "{tautomer} {atoms}"
        );
    }
}

/// Standard base64, padded.
fn base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    for chunk in bytes.chunks(3) {
        let n = chunk.len();
        let b = [
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ];
        let v = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        for i in 0..4 {
            if i <= n {
                out.push(ALPHABET[((v >> (18 - 6 * i)) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

#[test]
fn link_nodes_are_checked() {
    let enumerate = |text: &str, nodes: &str| {
        let mut mol = smiles(text);
        mol.set_prop("_molLinkNodes", nodes);
        create_extended_query_mol(&mol)
    };
    assert!(enumerate("C1CCCC1", "1 3 2 2 1 2 3").is_ok());
    // RDKit's link-node enumeration takes the first product of a reaction
    // that matches nothing: past an empty list (a segfault in Python).
    for (text, nodes) in [
        // A double bond.
        ("C1=CCC1", "1 3 2 1 2 1 4"),
        // Isotope labels 1000 * (index + 1) wrap past the 65th atom.
        (
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            "1 2 2 66 65 66 67",
        ),
        // Both bonds to one atom.
        ("C1CCCC1", "1 3 2 2 1 2 1"),
        // Two nodes on one atom.
        ("C1CCCC1", "2 3 2 2 1 2 3|2 4 2 2 3 2 1"),
    ] {
        assert!(
            matches!(enumerate(text, nodes), Err(Error::Value { .. })),
            "{nodes}"
        );
    }
    // The 65th atom still fits.
    assert!(
        enumerate(
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            "1 2 2 64 63 64 65",
        )
        .is_ok()
    );
    // Without enumeration the link nodes are not read.
    let mut mol = smiles("C1=CCC1");
    mol.set_prop("_molLinkNodes", "1 3 2 1 2 1 4");
    assert!(
        create_extended_query_mol_with_params(
            &mol,
            &CreateExtendedQueryMolParams {
                do_enumeration: false,
                ..Default::default()
            },
        )
        .is_ok()
    );
}

#[test]
fn queries_move_between_threads() {
    let query = create_extended_query_mol(&smiles("COC |LN:1:1.3|")).unwrap();
    let copy = query.clone();
    let handle =
        std::thread::spawn(move || mol_has_substruct_match(&smiles("COOC"), &copy).unwrap());
    assert!(handle.join().unwrap());
    assert!(mol_has_substruct_match(&smiles("COOOC"), &query).unwrap());
}

// RDKit writes a tautomer bundle's size in 16 bits, so a bundle of more
// queries, which JSON can hold, would read back as only some of them.
#[test]
fn large_tautomer_bundles_are_not_pickled() {
    let entry = r#"{"tautomers": [{"smiles": "C"}], "template": {"smiles": "C"}, "modifiedAtoms": [], "modifiedBonds": []}"#;
    let entries = vec![entry; 65536].join(", ");
    let json = format!(r#"{{"xqm_type": "4", "tautomerQueries": [{entries}]}}"#);
    // SAFETY: JSON with SMILES and no molecule pickles.
    let query = unsafe { ExtendedQueryMol::new_with_is_json(&json, true) }.unwrap();
    assert!(matches!(query.to_binary(), Err(Error::Value { .. })));
    assert!(query.to_json().is_ok());
    let small = format!(r#"{{"xqm_type": "4", "tautomerQueries": [{entry}]}}"#);
    // SAFETY: as above.
    let small = unsafe { ExtendedQueryMol::new_with_is_json(&small, true) }.unwrap();
    // SAFETY: as above.
    unsafe { ExtendedQueryMol::from_data(&small.to_binary().unwrap()) }.unwrap();
}

// Enumerating takes RDKit's enumerator a recursion level per variation
// point, as rdMolEnumerator's own checks say (CONVENTIONS.md U12).
#[test]
fn enumeration_refuses_too_many_variation_points() {
    let mut mol = smiles(&"C".repeat(4002));
    for idx in 0..4001 {
        let mut bond = mol.get_bond_with_idx_mut(idx).unwrap();
        bond.set_prop("_MolFileBondEndPts", "(1 1)");
        bond.set_prop("_MolFileBondAttach", "ANY");
    }
    match create_extended_query_mol(&mol).unwrap_err() {
        Error::Value { message } => assert!(message.contains("recursion deeper"), "{message}"),
        err => panic!("{err:?}"),
    }
}
