// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_tautomer_query::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::all_probe_bits_match;

// Expected values come from Python's rdTautomerQuery.

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn cyclohexanone() -> TautomerQuery {
    TautomerQuery::new(&smiles("O=C1CCCCC1")).unwrap()
}

fn canonical(mols: &[Mol]) -> Vec<String> {
    mols.iter().map(|m| mol_to_smiles(m).unwrap()).collect()
}

#[test]
fn query_from_a_molecule() {
    let query = cyclohexanone();
    assert_eq!(query.get_modified_atoms(), [0, 2, 6]);
    assert_eq!(query.get_modified_bonds(), [0, 1, 6]);
    assert_eq!(
        canonical(&query.get_tautomers()),
        ["O=C1CCCCC1", "OC1=CCCCC1"]
    );
    assert_eq!(
        mol_to_smarts(&query.get_template_molecule()).unwrap(),
        "[#8]-,=,:[#6]1-,=,:[#6]-[#6]-[#6]-[#6]-[#6]-,=,:1"
    );

    let target = smiles("OC1=CCCC(CC)C1");
    assert!(query.is_substruct_of(&target).unwrap());
    assert_eq!(
        query.get_substruct_match(&target).unwrap(),
        [0, 1, 2, 3, 4, 5, 8]
    );
    assert_eq!(
        query.get_substruct_matches(&target).unwrap(),
        [vec![0, 1, 2, 3, 4, 5, 8]]
    );
    let with_tautomers = query.get_substruct_matches_with_tautomers(&target).unwrap();
    assert_eq!(with_tautomers.len(), 1);
    assert_eq!(with_tautomers[0].0, [0, 1, 2, 3, 4, 5, 8]);
    assert_eq!(mol_to_smiles(&with_tautomers[0].1).unwrap(), "OC1=CCCCC1");

    let hexane = smiles("CCCCCC");
    assert!(!query.is_substruct_of(&hexane).unwrap());
    assert!(query.get_substruct_match(&hexane).unwrap().is_empty());

    let acetaldehyde = TautomerQuery::new(&smiles("CC=O")).unwrap();
    assert!(acetaldehyde.is_substruct_of(&smiles("OC(C)=O")).unwrap());
}

#[test]
fn match_parameters() {
    let query = cyclohexanone();
    let target = smiles("O=C1CCCC(CC)C1");
    let params = SubstructMatchParameters {
        uniquify: false,
        ..Default::default()
    };
    assert_eq!(
        query
            .get_substruct_matches_with_params(&target, &params)
            .unwrap(),
        [vec![0, 1, 2, 3, 4, 5, 8], vec![0, 1, 8, 5, 4, 3, 2]]
    );
    let with_tautomers = query
        .get_substruct_matches_with_tautomers_with_params(&target, &params)
        .unwrap();
    assert_eq!(with_tautomers.len(), 2);
    assert_eq!(with_tautomers[1].0, [0, 1, 8, 5, 4, 3, 2]);
    for (_, tautomer) in &with_tautomers {
        assert_eq!(mol_to_smiles(tautomer).unwrap(), "O=C1CCCCC1");
    }
    assert_eq!(
        query.get_substruct_matches(&target).unwrap(),
        [vec![0, 1, 2, 3, 4, 5, 8]]
    );
    assert_eq!(
        query
            .get_substruct_match_with_params(&target, &params)
            .unwrap()
            .len(),
        7
    );
    assert!(query.is_substruct_of_with_params(&target, &params).unwrap());
}

#[test]
fn pattern_fingerprints() {
    let query = cyclohexanone();
    let target = smiles("OC1=CCCC(CC)C1");
    let template_fp = query.pattern_fingerprint_template().unwrap();
    let target_fp = pattern_fingerprint_tautomer_target(&target).unwrap();
    assert_eq!(template_fp.get_num_bits(), 2048);
    assert_eq!(template_fp.get_num_on_bits(), 75);
    assert_eq!(target_fp.get_num_on_bits(), 132);
    assert_eq!(&template_fp.get_on_bits()[..4], [56, 173, 217, 296]);
    assert!(all_probe_bits_match(&template_fp, &target_fp).unwrap());

    let small = query
        .pattern_fingerprint_template_with_fingerprint_size(1024)
        .unwrap();
    assert_eq!(small.get_num_bits(), 1024);
    assert_eq!(small.get_num_on_bits(), 75);
    let small_target =
        pattern_fingerprint_tautomer_target_with_fingerprint_size(&target, 1024).unwrap();
    assert_eq!(small_target.get_num_on_bits(), 129);

    // ExplicitBitVects hold at most i32::MAX bits (U3); RDKit refuses 0.
    for size in [u32::MAX, 1 << 31] {
        assert!(matches!(
            query.pattern_fingerprint_template_with_fingerprint_size(size),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            pattern_fingerprint_tautomer_target_with_fingerprint_size(&target, size),
            Err(Error::Value { .. })
        ));
    }
    assert!(matches!(
        query.pattern_fingerprint_template_with_fingerprint_size(0),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn pickles() {
    assert!(tautomer_query_can_serialize());
    let query = cyclohexanone();
    let pickle = query.to_binary().unwrap();
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { TautomerQuery::from_pickle(&pickle) }.unwrap();
    assert_eq!(copy.get_modified_atoms(), [0, 2, 6]);
    assert_eq!(copy.get_tautomers().len(), 2);
    assert!(copy.is_substruct_of(&smiles("OC1=CCCC(CC)C1")).unwrap());
    assert_eq!(copy.get_modified_bonds(), [0, 1, 6]);
    // SAFETY: not a pickle; RDKit refuses it before reading a molecule.
    assert!(unsafe { TautomerQuery::from_pickle(b"") }.is_err());

    let clone = query.clone();
    drop(query);
    assert_eq!(clone.get_modified_atoms(), [0, 2, 6]);
}

#[test]
fn pickle_indices_are_checked() {
    // The archive ends with the modified atoms and bonds: a count, an item
    // version, and the indices.
    let pickle = cyclohexanone().to_binary().unwrap();
    let tail = b" 3 0 0 2 6 3 0 0 1 6\n";
    assert!(pickle.ends_with(tail));
    let head = &pickle[..pickle.len() - tail.len()];
    // RDKit's matcher indexes the template's match by each modified atom,
    // and each tautomer by each modified bond, without checking either.
    for bad in [&b" 3 0 0 2 99 3 0 0 1 6\n"[..], b" 3 0 0 2 6 3 0 0 1 60\n"] {
        let crafted = [head, bad].concat();
        assert!(matches!(
            // SAFETY: a pickle RDKit wrote but for a modified atom or bond index, which the wrapper
            // checks.
            unsafe { TautomerQuery::from_pickle(&crafted) },
            Err(Error::Value { .. })
        ));
    }
    let fine = [head, &b" 1 0 0 1 0 0\n"[..]].concat();
    assert_eq!(
        // SAFETY: as above, with indices the wrapper accepts.
        unsafe { TautomerQuery::from_pickle(&fine) }
            .unwrap()
            .get_modified_atoms(),
        [0]
    );
}

/// `CCCC(=O)C` with the first atom a query, as Python builds it with
/// `RWMol.InsertMol`.
fn with_first_atom(smarts: &str) -> Mol {
    let mut mol = RWMol::from(mol_from_smarts(smarts).unwrap());
    mol.insert_mol(&smiles("CCC(=O)C"));
    mol.add_bond_with_order(0, 1, BondType::Single).unwrap();
    mol.update_property_cache_with_strict(false).unwrap();
    mol.into()
}

#[test]
fn recursive_queries_in_the_template() {
    // RDKit hands the tautomer check the recursive SMARTS's own one-atom
    // matches, which it reads at the modified atoms 2, 4, and 5, past their
    // end. Here those matches pass unchecked, so the recursive query matches
    // as [CH3] does; the expected values are Python's for [CH3].
    let plain = TautomerQuery::new(&with_first_atom("[CH3]")).unwrap();
    let recursive = TautomerQuery::new(&with_first_atom("[$([CH3])]")).unwrap();
    assert_eq!(recursive.get_modified_atoms(), [2, 4, 5]);
    assert_eq!(recursive.get_modified_bonds(), [1, 2, 3]);
    assert_eq!(recursive.get_tautomers().len(), 3);
    let expected: [(&str, Vec<Vec<u32>>); 5] = [
        ("CCCC(=O)C", vec![vec![0, 1, 2, 3, 4, 5]]),
        ("CCCC(O)=C", vec![vec![0, 1, 2, 3, 4, 5]]),
        ("OCCC(=O)C", vec![]),
        ("CCC=C(O)C", vec![vec![0, 1, 2, 3, 4, 5]]),
        (
            "CC(C)CC(=O)C",
            vec![vec![0, 1, 3, 4, 5, 6], vec![2, 1, 3, 4, 5, 6]],
        ),
    ];
    for (target, matches) in expected {
        let target = smiles(target);
        assert_eq!(plain.get_substruct_matches(&target).unwrap(), matches);
        assert_eq!(recursive.get_substruct_matches(&target).unwrap(), matches);
        let with_tautomers = recursive
            .get_substruct_matches_with_tautomers(&target)
            .unwrap();
        assert_eq!(with_tautomers.len(), matches.len());
    }
}

#[test]
fn tautomer_transform_files() {
    let mol = smiles("O=C1CCCCC1");
    let write = |name: &str, text: &str| {
        let path = common::temp_path(name);
        std::fs::write(&path, text).unwrap();
        path
    };
    let keto = write(
        "keto.in",
        "1,3 (thio)keto/enol f\t[CX4!H0]-[C]=[O,S,Se,Te;X1]\n\
         1,3 (thio)keto/enol r\t[O,S,Se,Te;X2!H0]-[C]=[C]\n",
    );
    let query = TautomerQuery::new_with_tautomer_transform_file(&mol, &keto).unwrap();
    assert_eq!(query.get_modified_atoms(), [0, 2, 6]);
    assert_eq!(
        canonical(&query.get_tautomers()),
        ["O=C1CCCCC1", "OC1=CCCCC1"]
    );
    let default = TautomerQuery::new_with_tautomer_transform_file(&mol, "").unwrap();
    assert_eq!(default.get_tautomers().len(), 2);

    // RDKit sets each bond's type and adds each atom's charge from the lists
    // by position, reading past a short list (in Python, garbage bond types
    // and charges).
    let short_bonds = write(
        "short-bonds.in",
        "short bonds\t[CX4!H0]-[C]=[O,S,Se,Te;X1]\t=\n",
    );
    let short_charges = write(
        "short-charges.in",
        "short charges\t[CX4!H0]-[C]=[O,S,Se,Te;X1]\t=-\t+\n",
    );
    for path in [&short_bonds, &short_charges] {
        assert!(matches!(
            TautomerQuery::new_with_tautomer_transform_file(&mol, path),
            Err(Error::Value { .. })
        ));
    }
    let bad_smarts = write("bad-smarts.in", "bad\t[C\n");
    assert!(matches!(
        TautomerQuery::new_with_tautomer_transform_file(&mol, &bad_smarts),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        TautomerQuery::new_with_tautomer_transform_file(&mol, common::temp_path("missing.in")),
        Err(Error::Os { .. })
    ));
    assert!(matches!(
        TautomerQuery::new_with_tautomer_transform_file(&mol, "a\0b"),
        Err(Error::Value { .. })
    ));
    for path in [keto, short_bonds, short_charges, bad_smarts] {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn queries_move_between_threads() {
    let query = cyclohexanone();
    let target = smiles("OC1=CCCC(CC)C1");
    let handle = std::thread::spawn(move || query.is_substruct_of(&target).unwrap());
    assert!(handle.join().unwrap());
}
