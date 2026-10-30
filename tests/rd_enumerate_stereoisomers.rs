// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeSet;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_enumerate_stereoisomers::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdEnumerateStereoisomers. The isomers
// come in an order drawn from a random generator, whose distributions differ
// between C++ standard libraries, so they are compared as sets.

fn isomers(smi: &str, options: &StereoEnumerationOptions) -> (u64, Vec<String>) {
    let mol = mol_from_smiles(smi).unwrap();
    let mut enumerator = StereoisomerEnumerator::new_with_params(&mol, options).unwrap();
    let count = enumerator.get_stereoisomer_count();
    let mut smiles: Vec<String> = enumerator
        .by_ref()
        .map(|m| mol_to_smiles(&m.unwrap()).unwrap())
        .collect();
    // Once done, it stays done.
    assert!(enumerator.next().is_none());
    smiles.sort();
    (count, smiles)
}

fn with(f: impl FnOnce(&mut StereoEnumerationOptions)) -> StereoEnumerationOptions {
    let mut options = StereoEnumerationOptions::default();
    f(&mut options);
    options
}

#[test]
fn default_options() {
    assert_eq!(
        StereoEnumerationOptions::default(),
        StereoEnumerationOptions {
            try_embedding: false,
            only_unassigned: true,
            only_stereo_groups: false,
            unique: true,
            max_isomers: 0,
            random_seed: -1,
            verbose: true,
        }
    );
}

#[test]
fn enumerate() {
    let defaults = StereoEnumerationOptions::default();
    assert_eq!(
        isomers("CC(F)Cl", &defaults),
        (2, vec!["C[C@@H](F)Cl".to_owned(), "C[C@H](F)Cl".to_owned()])
    );
    assert_eq!(
        isomers("CC=CC", &defaults),
        (2, vec!["C/C=C/C".to_owned(), "C/C=C\\C".to_owned()])
    );
    assert_eq!(
        isomers("C[C@H](F)Cl", &defaults),
        (1, vec!["C[C@H](F)Cl".to_owned()])
    );
    assert_eq!(
        isomers("C[C@H](F)Cl", &with(|o| o.only_unassigned = false)),
        (2, vec!["C[C@@H](F)Cl".to_owned(), "C[C@H](F)Cl".to_owned()])
    );
    let meso = [
        "C[C@@H](F)[C@@H](C)F",
        "C[C@H](F)[C@@H](C)F",
        "C[C@H](F)[C@H](C)F",
    ];
    assert_eq!(
        isomers("CC(F)C(F)C", &defaults),
        (4, meso.map(String::from).to_vec())
    );
    let (count, all) = isomers("CC(F)C(F)C", &with(|o| o.unique = false));
    assert_eq!(count, 4);
    assert_eq!(all.len(), 4);
    assert_eq!(all.iter().cloned().collect::<BTreeSet<_>>().len(), 3);
    assert_eq!(isomers("CCO", &defaults), (1, vec!["CCO".to_owned()]));
    assert_eq!(isomers("", &defaults), (1, vec![String::new()]));
    assert_eq!(
        isomers("C[C@H](F)C[C@@H](C)F |o1:1,4|", &defaults),
        (2, vec!["C[C@H](F)C[C@@H](C)F".to_owned()])
    );
}

#[test]
fn max_isomers_and_seed() {
    let (count, one) = isomers("CC(F)Cl", &with(|o| o.max_isomers = 1));
    assert_eq!(count, 2);
    assert_eq!(one.len(), 1);
    assert!(["C[C@@H](F)Cl", "C[C@H](F)Cl"].contains(&one[0].as_str()));
    // A seed fixes the order.
    let seeded = with(|o| o.random_seed = 42);
    let order = |options: &StereoEnumerationOptions| -> Vec<String> {
        let mol = mol_from_smiles("CC(F)C(Cl)C(Br)C").unwrap();
        StereoisomerEnumerator::new_with_params(&mol, options)
            .unwrap()
            .map(|m| mol_to_smiles(&m.unwrap()).unwrap())
            .collect()
    };
    let first = order(&seeded);
    assert_eq!(first.len(), 8);
    assert_eq!(order(&seeded), first);
}

#[test]
fn try_embedding() {
    // Only the exo and endo isomers of 2-methylnorbornane exist in 3D.
    let options = with(|o| o.try_embedding = true);
    let (count, embedded) = isomers("CC1CC2CCC1C2", &options);
    assert_eq!(count, 8);
    assert_eq!(
        embedded,
        [
            "C[C@@H]1C[C@@H]2CC[C@H]1C2",
            "C[C@@H]1C[C@H]2CC[C@@H]1C2",
            "C[C@H]1C[C@@H]2CC[C@H]1C2",
            "C[C@H]1C[C@H]2CC[C@@H]1C2",
        ]
    );
    let mol = mol_from_smiles("CC1CC2CCC1C2").unwrap();
    let isomer = StereoisomerEnumerator::new_with_params(&mol, &options)
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(isomer.get_num_atoms(), 8);
    assert_eq!(isomer.get_num_conformers(), 1);
    assert!(isomer.get_conformer().unwrap().is_3d());
    // Without stereo elements RDKit 2025.03's enumerator still embeds the
    // one molecule it returns, as Python's rdEnumerateStereoisomers does.
    let plain = StereoisomerEnumerator::new_with_params(&mol_from_smiles("CCO").unwrap(), &options)
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(plain.get_num_conformers(), 1);
    // verbose defaults to true, as in Python.
    assert_eq!(StereoisomerEnumerator::new(&mol).unwrap().count(), 8);
}

// With its hydrogens the chain has more atoms than RDKit's minimizer can
// size its Hessian for (BFGSOpt sizes it in unsigned int), which the
// embedder's checks refuse before RDKit embeds anything.
#[test]
fn embedding_checks() {
    let smi = format!("FC(Cl){}", "C".repeat(5500));
    let mol = mol_from_smiles(&smi).unwrap();
    let options = with(|o| o.try_embedding = true);
    let mut enumerator = StereoisomerEnumerator::new_with_params(&mol, &options).unwrap();
    assert_eq!(enumerator.get_stereoisomer_count(), 2);
    for _ in 0..2 {
        let err = enumerator.next().unwrap().unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{err:?}");
    }
}

#[test]
fn threads() {
    std::thread::scope(|s| {
        for _ in 0..4 {
            s.spawn(|| {
                let options = with(|o| o.try_embedding = true);
                for _ in 0..5 {
                    assert_eq!(isomers("CC(F)C(Cl)Br", &options).1.len(), 4);
                }
            });
        }
    });
    // An enumerator moves between threads.
    let mol = mol_from_smiles("CC(F)Cl").unwrap();
    let enumerator = StereoisomerEnumerator::new(&mol).unwrap();
    let count = std::thread::spawn(move || enumerator.count())
        .join()
        .unwrap();
    assert_eq!(count, 2);
}
