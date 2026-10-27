// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's property registry. Expected values are from Python's
//! rdMolDescriptors. The registry is process-wide, so each test registers
//! properties under its own names.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_descriptors::*;
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles};

fn apap() -> Mol {
    mol_from_smiles("CC(=O)Nc1ccc(O)cc1").unwrap()
}

const BUILT_IN: [&str; 43] = [
    "exactmw",
    "amw",
    "lipinskiHBA",
    "lipinskiHBD",
    "NumRotatableBonds",
    "NumHBD",
    "NumHBA",
    "NumHeavyAtoms",
    "NumAtoms",
    "NumHeteroatoms",
    "NumAmideBonds",
    "FractionCSP3",
    "NumRings",
    "NumAromaticRings",
    "NumAliphaticRings",
    "NumSaturatedRings",
    "NumHeterocycles",
    "NumAromaticHeterocycles",
    "NumSaturatedHeterocycles",
    "NumAliphaticHeterocycles",
    "NumSpiroAtoms",
    "NumBridgeheadAtoms",
    "NumAtomStereoCenters",
    "NumUnspecifiedAtomStereoCenters",
    "labuteASA",
    "tpsa",
    "CrippenClogP",
    "CrippenMR",
    "chi0v",
    "chi1v",
    "chi2v",
    "chi3v",
    "chi4v",
    "chi0n",
    "chi1n",
    "chi2n",
    "chi3n",
    "chi4n",
    "hallKierAlpha",
    "kappa1",
    "kappa2",
    "kappa3",
    "Phi",
];

const VALUES: [f64; 43] = [
    151.063328528,
    151.165,
    3.0,
    2.0,
    1.0,
    2.0,
    2.0,
    11.0,
    20.0,
    3.0,
    1.0,
    0.125,
    1.0,
    1.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    64.66690538413799,
    49.33,
    1.3505999999999998,
    42.41050000000002,
    6.164862962722325,
    3.2490981480278296,
    // RDKit 2024.03 computes chi2v and chi2n as chi3v and chi3n (#9643).
    1.1876733148984595,
    1.1876733148984595,
    0.7480265991539318,
    6.164862962722325,
    3.2490981480278296,
    1.1876733148984595,
    1.1876733148984595,
    0.7480265991539318,
    -1.5099999999999998,
    7.595374077976819,
    3.0531381005525,
    2.292303972882069,
    2.1081569077654403,
];

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn built_in_properties() {
    let available = Properties::get_available_properties();
    assert_eq!(available[..43], BUILT_IN);
    let all = Properties::default();
    let names = all.get_property_names();
    assert_eq!(names[..43], BUILT_IN);
    let values = all.compute_properties(&apap()).unwrap();
    for (value, expected) in values.iter().zip(VALUES) {
        assert_close(*value, expected);
    }

    let subset = Properties::new(&["exactmw", "tpsa"]).unwrap();
    assert_eq!(subset.get_property_names(), vec!["exactmw", "tpsa"]);
    let mut mol = apap();
    let values = subset
        .compute_properties_with_annotate_mol(&mut mol, false)
        .unwrap();
    assert_close(values[1], 49.33);
    assert!(!mol.has_prop("tpsa"));
    subset
        .compute_properties_with_annotate_mol(&mut mol, true)
        .unwrap();
    assert_close(mol.get_double_prop("tpsa").unwrap(), 49.33);
    let mut other = apap();
    subset.annotate_properties(&mut other).unwrap();
    assert_close(other.get_double_prop("exactmw").unwrap(), 151.063328528);
    assert!(matches!(Properties::new(&["nope"]), Err(Error::Key { .. })));

    let functor = Properties::get_property("NumHBD").unwrap();
    assert_eq!(functor.get_name(), "NumHBD");
    assert_eq!(functor.get_version(), "2.0.1");
    assert_eq!(functor.call(&apap()).unwrap(), 2.0);
    assert!(matches!(
        Properties::get_property("nope"),
        Err(Error::Key { .. })
    ));
}

#[test]
fn range_queries() {
    let mol = apap();
    assert!(
        make_property_range_query("exactmw", 0.0, 500.0)
            .unwrap()
            .match_(&mol)
            .unwrap()
    );
    assert!(
        !make_property_range_query("exactmw", 0.0, 100.0)
            .unwrap()
            .match_(&mol)
            .unwrap()
    );
    assert!(
        !make_property_range_query("exactmw", 151.063328528, 151.063328528)
            .unwrap()
            .match_(&mol)
            .unwrap()
    );
    assert!(matches!(
        make_property_range_query("nope", 0.0, 1.0),
        Err(Error::Key { .. })
    ));
}

#[test]
fn rust_property_functors() {
    let calls = Arc::new(AtomicUsize::new(0));
    let counted = Arc::clone(&calls);
    let heavy = PropertyFunctor::new("rustHeavyAtoms", "1.0.0", move |mol| {
        counted.fetch_add(1, Ordering::SeqCst);
        f64::from(mol.get_num_heavy_atoms())
    });
    assert_eq!(heavy.get_name(), "rustHeavyAtoms");
    assert_eq!(heavy.get_version(), "1.0.0");
    assert_eq!(heavy.call(&apap()).unwrap(), 11.0);
    let idx = Properties::register_property(&heavy);
    assert!(idx >= 43);
    assert_eq!(Properties::register_property(&heavy), idx);
    let available = Properties::get_available_properties();
    assert_eq!(available[idx as usize], "rustHeavyAtoms");
    let properties = Properties::new(&["rustHeavyAtoms", "tpsa"]).unwrap();
    let values = properties.compute_properties(&apap()).unwrap();
    assert_eq!(values[0], 11.0);
    assert_close(values[1], 49.33);
    assert_eq!(
        Properties::get_property("rustHeavyAtoms")
            .unwrap()
            .call(&apap())
            .unwrap(),
        11.0
    );
    assert_eq!(calls.load(Ordering::SeqCst), 3);
    // The query calls RDKit's descriptor function, which a Rust property
    // does not have.
    let query = make_property_range_query("rustHeavyAtoms", 0.0, 100.0).unwrap();
    assert!(matches!(query.match_(&apap()), Err(Error::Runtime { .. })));

    // A replacement keeps the index; earlier Properties keep the old one.
    let replacement = PropertyFunctor::new("rustHeavyAtoms", "2.0.0", |_| -1.0);
    assert_eq!(Properties::register_property(&replacement), idx);
    assert_eq!(properties.compute_properties(&apap()).unwrap()[0], 11.0);
    let updated = Properties::new(&["rustHeavyAtoms"]).unwrap();
    assert_eq!(updated.compute_properties(&apap()).unwrap(), vec![-1.0]);
}

#[test]
fn rust_property_functors_can_match_substructures() {
    // Computing properties holds control_c_mutex, released while a Rust
    // functor runs, so its closure can call wrappers that take it.
    let hydroxyls = PropertyFunctor::new("rustHydroxyls", "1.0.0", |mol| {
        let query = mol_from_smarts("[OX2H]").unwrap();
        mol.get_substruct_matches(&query).unwrap().len() as f64
    });
    Properties::register_property(&hydroxyls);
    let properties = Properties::new(&["NumHBD", "rustHydroxyls", "NumHBA"]).unwrap();
    assert_eq!(
        properties.compute_properties(&apap()).unwrap(),
        vec![2.0, 1.0, 2.0]
    );
}

#[test]
fn panics_are_caught() {
    let panicking = PropertyFunctor::new("rustPanics", "1.0.0", |_| panic!("no value"));
    match panicking.call(&apap()) {
        Err(Error::Runtime { message }) => assert!(message.contains("no value"), "{message}"),
        other => panic!("{other:?}"),
    }
    Properties::register_property(&panicking);
    let properties = Properties::new(&["tpsa", "rustPanics"]).unwrap();
    assert!(matches!(
        properties.compute_properties(&apap()),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn closures_get_a_copy() {
    // The closure sees the molecule's atoms but not the original object.
    let functor = PropertyFunctor::new("rustAtoms", "1.0.0", |mol| f64::from(mol.get_num_atoms()));
    let mut mol = apap();
    mol.set_prop("x", "y");
    assert_eq!(functor.call(&mol).unwrap(), 11.0);
}

#[test]
fn registry_in_threads() {
    // RDKit writes the registry without a lock; the wrappers hold one.
    let threads: Vec<_> = (0..4)
        .map(|i| {
            std::thread::spawn(move || {
                let name = format!("rustThread{i}");
                let functor = PropertyFunctor::new(&name, "1.0.0", move |_| f64::from(i));
                for _ in 0..20 {
                    Properties::register_property(&functor);
                    let properties = Properties::new(&[name.as_str()]).unwrap();
                    assert_eq!(
                        properties.compute_properties(&apap()).unwrap(),
                        vec![f64::from(i)]
                    );
                    assert!(Properties::get_available_properties().len() >= 44);
                }
            })
        })
        .collect();
    for thread in threads {
        thread.join().unwrap();
    }
}
