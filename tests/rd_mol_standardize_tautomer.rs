// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::{
    Atom, Mol, RWMol, SmilesWriteParams, mol_from_smarts, mol_from_smiles, mol_to_smiles,
    mol_to_smiles_with_params,
};

// Expected values come from Python's rdMolStandardize.

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn smiles(mol: &Mol) -> String {
    mol_to_smiles(mol).unwrap()
}

fn all_smiles(mols: &[Mol]) -> Vec<String> {
    mols.iter().map(smiles).collect()
}

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

struct Case {
    input: &'static str,
    tautomers: &'static [&'static str],
    modified_atoms: &'static [u32],
    modified_bonds: &'static [u32],
    canonical: &'static str,
    kekulized: &'static [&'static str],
    scores: &'static [i32],
    rings: &'static [i32],
    substructs: &'static [i32],
}

const CASES: [Case; 5] = [
    Case {
        input: "Oc1ccccn1",
        tautomers: &["O=C1CC=CC=N1", "O=c1cccc[nH]1", "Oc1ccccn1"],
        modified_atoms: &[0, 2, 6],
        modified_bonds: &[0, 1, 2, 3, 4, 5, 6],
        canonical: "O=c1cccc[nH]1",
        kekulized: &["O=C1CC=CC=N1", "O=C1C=CC=CN1", "OC1=CC=CC=N1"],
        scores: &[6, 102, 100],
        rings: &[0, 100, 100],
        substructs: &[6, 2, 0],
    },
    Case {
        input: "CC(=O)C",
        tautomers: &["C=C(C)O", "CC(C)=O"],
        modified_atoms: &[0, 2, 3],
        modified_bonds: &[0, 1, 2],
        canonical: "CC(C)=O",
        kekulized: &["C=C(C)O", "CC(C)=O"],
        scores: &[1, 5],
        rings: &[0, 0],
        substructs: &[1, 5],
    },
    Case {
        input: "OC(C)=C(C)C",
        tautomers: &["C=C(O)C(C)C", "CC(=O)C(C)C", "CC(C)=C(C)O"],
        modified_atoms: &[0, 2, 3],
        modified_bonds: &[0, 1, 2],
        canonical: "CC(=O)C(C)C",
        kekulized: &["C=C(O)C(C)C", "CC(=O)C(C)C", "CC(C)=C(C)O"],
        scores: &[2, 6, 3],
        rings: &[0, 0, 0],
        substructs: &[2, 6, 3],
    },
    Case {
        input: "c1ccc2[nH]ncc2c1",
        tautomers: &["c1ccc2[nH]ncc2c1", "c1ccc2n[nH]cc2c1"],
        modified_atoms: &[4, 5],
        modified_bonds: &[0, 1, 2, 3, 5, 6, 7, 8, 9],
        canonical: "c1ccc2[nH]ncc2c1",
        kekulized: &["C1=CC=C2NN=CC2=C1", "C1=CC2=CNN=C2C=C1"],
        scores: &[350, 350],
        rings: &[350, 350],
        substructs: &[0, 0],
    },
    Case {
        input: "C[C@H](N)C(=O)O",
        tautomers: &["C=C(N)C(O)O", "CC(=N)C(O)O", "CC(N)=C(O)O", "CC(N)C(=O)O"],
        modified_atoms: &[0, 1, 2, 3, 4, 5],
        modified_bonds: &[0, 1, 2, 3, 4],
        canonical: "CC(N)C(=O)O",
        kekulized: &["C=C(N)C(O)O", "CC(=N)C(O)O", "CC(N)=C(O)O", "CC(N)C(=O)O"],
        scores: &[0, 2, 1, 6],
        rings: &[0, 0, 0, 0],
        substructs: &[0, 2, 1, 6],
    },
];

#[test]
fn enumeration() {
    let enumerator = TautomerEnumerator::default();
    let kekule = SmilesWriteParams {
        do_kekule: true,
        ..Default::default()
    };
    for case in CASES {
        let result = enumerator.enumerate(&mol(case.input)).unwrap();
        assert_eq!(result.smiles().unwrap(), case.tautomers, "{}", case.input);
        assert_eq!(all_smiles(&result.tautomers()), case.tautomers);
        assert_eq!(result.status(), TautomerEnumeratorStatus::Completed);
        assert_eq!(result.modified_atoms(), case.modified_atoms);
        assert_eq!(result.modified_bonds(), case.modified_bonds);
        assert_eq!(result.len(), case.tautomers.len());
        assert!(!result.is_empty());
        let last = case.tautomers.len() as i32 - 1;
        assert_eq!(smiles(&result.get(0).unwrap()), case.tautomers[0]);
        assert_eq!(
            smiles(&result.get(-1).unwrap()),
            case.tautomers[last as usize]
        );
        assert!(matches!(result.get(last + 1), Err(Error::Index { .. })));
        assert!(matches!(result.get(-last - 2), Err(Error::Index { .. })));

        let map = result.smiles_tautomer_map().unwrap();
        assert_eq!(map.len(), case.tautomers.len());
        assert_eq!(map.keys(), case.tautomers);
        for ((key, tautomer), (want, kekulized)) in map
            .items()
            .into_iter()
            .zip(case.tautomers.iter().zip(case.kekulized))
        {
            assert_eq!(key, *want);
            assert_eq!(smiles(&tautomer.tautomer()), *want);
            assert_eq!(
                mol_to_smiles_with_params(&tautomer.kekulized().unwrap(), &kekule).unwrap(),
                *kekulized
            );
        }
        assert_eq!(map.values().len(), case.tautomers.len());
        let first = map.get(case.tautomers[0]).unwrap();
        assert_eq!(smiles(&first.tautomer()), case.tautomers[0]);
        assert!(map.contains_key(case.tautomers[0]));
        assert!(map.get("not a tautomer").is_none());
        assert!(!map.contains_key("not a tautomer"));

        assert_eq!(
            smiles(&enumerator.canonicalize(&mol(case.input)).unwrap()),
            case.canonical
        );
        assert_eq!(
            smiles(
                &enumerator
                    .pick_canonical_tautomer_enumerator_result(&result)
                    .unwrap()
            ),
            case.canonical
        );
        let mols: Vec<Mol> = case.tautomers.iter().map(|s| mol(s)).collect();
        assert_eq!(
            smiles(&enumerator.pick_canonical(&mols).unwrap()),
            case.canonical
        );

        for (i, tautomer) in case.tautomers.iter().enumerate() {
            let m = mol(tautomer);
            assert_eq!(
                TautomerEnumerator::score_tautomer(&m).unwrap(),
                case.scores[i]
            );
            assert_eq!(score_rings(&m).unwrap(), case.rings[i]);
            assert_eq!(score_substructs(&m).unwrap(), case.substructs[i]);
            assert_eq!(score_hetero_hs(&m).unwrap(), 0);
        }
    }
    assert_eq!(score_rings(&mol("c1ccccc1-c1ccccn1")).unwrap(), 350);
    assert_eq!(score_hetero_hs(&mol("SCC[PH2]")).unwrap(), -3);
    assert_eq!(TautomerEnumerator::tautomer_score_version(), "1.0.0");
}

#[test]
fn scoring_functions() {
    let enumerator = TautomerEnumerator::default();
    let m = mol("Oc1ccccn1");
    let aliphatic = |m: &Mol| {
        let aromatic = (0..m.get_num_atoms())
            .filter(|&i| m.get_atom_with_idx(i).unwrap().get_is_aromatic())
            .count();
        m.get_num_atoms() as i32 - aromatic as i32
    };
    assert_eq!(
        smiles(
            &enumerator
                .canonicalize_with_score_func(&m, aliphatic)
                .unwrap()
        ),
        "O=C1CC=CC=N1"
    );
    let result = enumerator.enumerate(&m).unwrap();
    let negated = |m: &Mol| -TautomerEnumerator::score_tautomer(m).unwrap();
    assert_eq!(
        smiles(
            &enumerator
                .pick_canonical_tautomer_enumerator_result_with_score_func(&result, negated)
                .unwrap()
        ),
        "O=C1CC=CC=N1"
    );
    let mols = result.tautomers();
    let length = |m: &Mol| smiles(m).len() as i32;
    assert_eq!(
        smiles(
            &enumerator
                .pick_canonical_with_score_func(&mols, length)
                .unwrap()
        ),
        "O=c1cccc[nH]1"
    );
    // A closure may borrow.
    let offsets = [3, 2, 1];
    let by_offset = |m: &Mol| offsets[m.get_num_atoms() as usize % 3];
    enumerator
        .pick_canonical_with_score_func(&mols, by_offset)
        .unwrap();
}

#[test]
fn score_functions_that_return_i32_min() {
    // RDKit starts its search at i32::MIN and dereferences null when no
    // tautomer scores above it (Python segfaults); i32::MIN counts as
    // i32::MIN + 1, so the first tautomer is picked.
    let enumerator = TautomerEnumerator::default();
    let min = |_: &Mol| i32::MIN;
    let mols = [mol("CC"), mol("CCC")];
    assert_eq!(
        smiles(
            &enumerator
                .pick_canonical_with_score_func(&mols, min)
                .unwrap()
        ),
        "CC"
    );
    assert_eq!(
        smiles(
            &enumerator
                .canonicalize_with_score_func(&mol("Oc1ccccn1"), min)
                .unwrap()
        ),
        "O=C1CC=CC=N1"
    );
    let result = enumerator.enumerate(&mol("Oc1ccccn1")).unwrap();
    assert_eq!(
        smiles(
            &enumerator
                .pick_canonical_tautomer_enumerator_result_with_score_func(&result, min)
                .unwrap()
        ),
        "O=C1CC=CC=N1"
    );
}

#[test]
fn picking_from_nothing_is_refused() {
    // pickCanonical dereferences null when it has no tautomers (Python
    // segfaults for PickCanonical([])).
    let enumerator = TautomerEnumerator::default();
    assert!(is_value(enumerator.pick_canonical(&[])));
    assert!(is_value(
        enumerator.pick_canonical_with_score_func(&[], |_| 0)
    ));
    // One molecule is returned unscored.
    assert_eq!(
        smiles(&enumerator.pick_canonical(&[mol("Oc1ccccn1")]).unwrap()),
        "Oc1ccccn1"
    );
}

#[test]
fn panics_in_score_functions_are_errors() {
    let enumerator = TautomerEnumerator::default();
    let err = enumerator
        .canonicalize_with_score_func(&mol("Oc1ccccn1"), |_| panic!("no score"))
        .unwrap_err();
    assert!(
        matches!(&err, Error::Runtime { message } if message.contains("no score")),
        "{err:?}"
    );
}

#[test]
fn options_and_limits() {
    let mut enumerator = TautomerEnumerator::default();
    assert_eq!(enumerator.get_max_tautomers(), 1000);
    assert_eq!(enumerator.get_max_transforms(), 1000);
    assert!(enumerator.get_remove_sp3_stereo());
    assert!(enumerator.get_remove_bond_stereo());
    assert!(enumerator.get_reassign_stereo());
    assert!(enumerator.get_callback().is_none());
    let original = enumerator.clone();
    enumerator.set_max_tautomers(2);
    enumerator.set_max_transforms(3);
    enumerator.set_remove_sp3_stereo(false);
    enumerator.set_remove_bond_stereo(false);
    enumerator.set_reassign_stereo(false);
    assert_eq!(enumerator.get_max_tautomers(), 2);
    assert_eq!(enumerator.get_max_transforms(), 3);
    assert!(!enumerator.get_remove_sp3_stereo());
    assert!(!enumerator.get_remove_bond_stereo());
    assert!(!enumerator.get_reassign_stereo());
    assert_eq!(original.get_max_tautomers(), 1000);

    let m = mol("Oc1nc(O)ccc1");
    let result = enumerator.enumerate(&m).unwrap();
    assert_eq!(
        result.smiles().unwrap(),
        ["O=C1CC=CC(O)=N1", "Oc1cccc(O)n1"]
    );
    assert_eq!(
        result.status(),
        TautomerEnumeratorStatus::MaxTautomersReached
    );
    let mut transforms = TautomerEnumerator::default();
    transforms.set_max_transforms(1);
    let result = transforms.enumerate(&m).unwrap();
    assert_eq!(result.smiles().unwrap(), ["Oc1cccc(O)n1"]);
    assert_eq!(
        result.status(),
        TautomerEnumeratorStatus::MaxTransformsReached
    );
    assert_eq!(
        original.enumerate(&m).unwrap().smiles().unwrap(),
        [
            "O=C1C=CCC(=O)N1",
            "O=C1C=CCC(O)=N1",
            "O=C1CC=CC(O)=N1",
            "O=c1cccc(O)[nH]1",
            "Oc1cccc(O)n1"
        ]
    );

    let params = CleanupParameters {
        max_tautomers: 3,
        tautomer_remove_sp3_stereo: false,
        ..Default::default()
    };
    let from_params = TautomerEnumerator::new(&params).unwrap();
    assert_eq!(from_params.get_max_tautomers(), 3);
    assert!(!from_params.get_remove_sp3_stereo());
    assert_eq!(
        from_params
            .enumerate(&mol("C[C@H](N)C(=O)O"))
            .unwrap()
            .smiles()
            .unwrap(),
        ["CC(=N)C(O)O", "CC(N)=C(O)O", "C[C@H](N)C(=O)O"]
    );
    assert_eq!(
        from_params.enumerate(&m).unwrap().status(),
        TautomerEnumeratorStatus::MaxTautomersReached
    );

    assert_eq!(
        get_v1_tautomer_enumerator()
            .enumerate(&mol("c1ccc2[nH]ncc2c1"))
            .unwrap()
            .smiles()
            .unwrap(),
        ["c1ccc2[nH]ncc2c1", "c1ccc2n[nH]cc2c1"]
    );
}

#[test]
fn callbacks() {
    let calls = Arc::new(AtomicUsize::new(0));
    let sizes = Arc::new(std::sync::Mutex::new(Vec::new()));
    let (counter, seen) = (calls.clone(), sizes.clone());
    let callback = TautomerEnumeratorCallback::new(move |m, res| {
        assert_eq!(m.get_num_atoms(), 8);
        seen.lock().unwrap().push(res.len());
        counter.fetch_add(1, Ordering::SeqCst) + 1 < 3
    });
    let mut enumerator = TautomerEnumerator::default();
    enumerator.set_callback(Some(&callback));
    assert!(enumerator.get_callback().is_some());
    let m = mol("Oc1nc(O)ccc1");
    let result = enumerator.enumerate(&m).unwrap();
    assert_eq!(
        result.smiles().unwrap(),
        ["O=C1CC=CC(O)=N1", "O=c1cccc(O)[nH]1", "Oc1cccc(O)n1"]
    );
    assert_eq!(result.status(), TautomerEnumeratorStatus::Canceled);
    assert_eq!(calls.load(Ordering::SeqCst), 3);
    assert_eq!(*sizes.lock().unwrap(), [1, 2, 3]);
    // Canonicalize enumerates with the callback too.
    assert_eq!(
        smiles(&enumerator.canonicalize(&m).unwrap()),
        "Oc1cccc(O)n1"
    );
    assert_eq!(calls.load(Ordering::SeqCst), 4);
    // The callback can be called directly, and from clones on other threads.
    assert!(
        !enumerator
            .get_callback()
            .unwrap()
            .call(&m, &result)
            .unwrap()
    );
    std::thread::scope(|scope| {
        for _ in 0..3 {
            let enumerator = enumerator.clone();
            scope.spawn(move || {
                let m = mol("Oc1nc(O)ccc1");
                enumerator.enumerate(&m).unwrap();
            });
        }
    });
    assert_eq!(calls.load(Ordering::SeqCst), 8);
    enumerator.set_callback(None);
    assert!(enumerator.get_callback().is_none());
    assert_eq!(enumerator.enumerate(&m).unwrap().len(), 5);
    drop(enumerator);
    drop(callback);
    assert_eq!(Arc::strong_count(&calls), 1);
}

#[test]
fn callbacks_see_the_result_so_far() {
    // The result lent to a callback shares RDKit's tautomers; its accessors
    // copy them.
    let callback = TautomerEnumeratorCallback::new(|_, res| {
        let tautomers = res.tautomers();
        let map = res.smiles_tautomer_map().unwrap();
        assert_eq!(map.len(), tautomers.len());
        for tautomer in map.values() {
            tautomer.kekulized().unwrap();
        }
        TautomerEnumerator::default()
            .pick_canonical_tautomer_enumerator_result(res)
            .unwrap();
        true
    });
    let mut enumerator = TautomerEnumerator::default();
    enumerator.set_callback(Some(&callback));
    assert_eq!(enumerator.enumerate(&mol("Oc1nc(O)ccc1")).unwrap().len(), 5);
}

#[test]
fn closures_can_match_substructures() {
    // Enumerating and scoring hold control_c_mutex, released while a Rust
    // callback or score runs, so the closures can call wrappers that take it.
    let callback = TautomerEnumeratorCallback::new(|m, _| {
        let nitrogen = mol_from_smarts("[#7]").unwrap();
        m.has_substruct_match(&nitrogen).unwrap()
    });
    let mut enumerator = TautomerEnumerator::default();
    enumerator.set_callback(Some(&callback));
    let m = mol("Oc1ccccn1");
    assert_eq!(enumerator.enumerate(&m).unwrap().len(), 3);
    let hydroxyls = |t: &Mol| {
        let query = mol_from_smarts("[OX2H]").unwrap();
        t.get_substruct_matches(&query).unwrap().len() as i32
    };
    assert_eq!(
        smiles(
            &enumerator
                .canonicalize_with_score_func(&m, hydroxyls)
                .unwrap()
        ),
        "Oc1ccccn1"
    );
}

#[test]
fn panics_in_callbacks_are_errors() {
    let callback = TautomerEnumeratorCallback::new(|_, _| panic!("no callback"));
    let mut enumerator = TautomerEnumerator::default();
    enumerator.set_callback(Some(&callback));
    let err = enumerator.enumerate(&mol("Oc1ccccn1")).unwrap_err();
    assert!(
        matches!(&err, Error::Runtime { message } if message.contains("no callback")),
        "{err:?}"
    );
}

#[test]
fn substruct_terms() {
    let terms = get_default_tautomer_score_substructs();
    assert_eq!(terms.len(), 12);
    assert_eq!(
        terms[0],
        SubstructTerm::new(
            "benzoquinone",
            "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
            25
        )
    );
    assert_eq!(terms[1], SubstructTerm::new("oxim", "[#6]=[N][OH]", 4));
    let custom = [
        SubstructTerm::new("C=O", "[#6]=,:[#8]", 7),
        SubstructTerm::new("bad", "[C", 5),
    ];
    let acetone = mol("CC(=O)C");
    assert_eq!(score_substructs_with_terms(&acetone, &custom).unwrap(), 7);
    assert_eq!(score_substructs(&acetone).unwrap(), 5);
    // RDKit counts at most 1000 matches of a term.
    let chain = mol(&"C".repeat(3000));
    let carbon = [SubstructTerm::new("c", "[#6]", 1)];
    assert_eq!(score_substructs_with_terms(&chain, &carbon).unwrap(), 1000);
}

#[test]
fn large_scores_are_refused() {
    // RDKit adds each term's score times its matches (up to 1000) in int;
    // Python gets -1000 here.
    let chain = mol(&"C".repeat(3000));
    let big = [SubstructTerm::new("c", "[#6]", i32::MAX)];
    assert!(is_value(score_substructs_with_terms(&chain, &big)));
    let big = [
        SubstructTerm::new("c", "[#6]", i32::MAX / 2000 + 1),
        SubstructTerm::new("n", "[#7]", -(i32::MAX / 2000 + 1)),
    ];
    assert!(is_value(score_substructs_with_terms(&chain, &big)));

    // The default score adds up to 382 per atom for hydrogens on P, S, Se,
    // and Te, and 250 per ring, in int. (Rings enough to overflow
    // score_rings, some 8.6 million, take gigabytes to enumerate.)
    let atom = Atom::from_num(6);
    let mut many = RWMol::default();
    for _ in 0..(i32::MAX as u32 - 47000) / 383 + 1 {
        many.add_atom(&atom);
    }
    let enumerator = TautomerEnumerator::default();
    assert!(is_value(TautomerEnumerator::score_tautomer(&many)));
    assert!(is_value(enumerator.canonicalize(&many)));
    let many = Mol::from(many);
    assert!(is_value(enumerator.pick_canonical(&[many, mol("CC")])));
}

#[test]
fn short_transforms_are_refused() {
    // RDKit reads a transform's bonds and charges past their ends when
    // they are shorter than its SMARTS (Python reads garbage).
    let params = CleanupParameters {
        tautomer_transform_data: vec![(
            "x".to_string(),
            "[CX4!H0]-[C]=[O]".to_string(),
            "=".to_string(),
            String::new(),
        )],
        ..Default::default()
    };
    assert!(is_value(TautomerEnumerator::new(&params)));
}
