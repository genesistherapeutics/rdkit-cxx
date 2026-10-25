// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// rdMolStandardize on several threads at once. Most of its calls hold
// control_c_mutex for the whole RDKit call (CONVENTIONS.md U4), and the
// closures RDKit calls during them release it, so a closure may call
// functions that take it.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles};

// Expected values come from Python's rdMolStandardize, besides those of the
// Rust closures.

const BENZOATE: &str = "
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    2.5833    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8333    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5833   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3333   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4167    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9167    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6667    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9167   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4167   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.2990    0.0000 Na  0  0  0  0  0 15  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9  4  1  0
M  CHG  2   3  -1  10   1
M  END
";

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn smiles(mol: &Mol) -> String {
    mol_to_smiles(mol).unwrap()
}

// What one thread does, with the objects all threads share: a pipeline
// (Sync), a validation whose Rust method charge-parents the molecule inside
// a MolVSValidation whose FragmentValidation matches, an enumerator (each
// thread's clone shares the catalog and a callback that normalizes), and a
// score function that scores substructures. The normalizer and enumerator
// of the multi-molecule calls are shared by their own threads.
struct Expected {
    pipeline: String,
    validation: Vec<String>,
    canonical: String,
    scored: String,
    normalized: Vec<String>,
    tautomer_parents: Vec<String>,
}

fn work(
    pipeline: &Pipeline,
    validation: &MolVSValidation,
    enumerator: &TautomerEnumerator,
) -> Expected {
    let pyridone = mol("Oc1ccccn1");
    let mut normalized: Vec<Mol> = (0..8).map(|_| mol("CN(=O)=O")).collect();
    normalize_in_place_mols(&mut normalized, 4).unwrap();
    let mut tautomer_parents: Vec<Mol> = (0..8).map(|_| mol("[Na+].Oc1ccccn1")).collect();
    tautomer_parent_in_place_mols(&mut tautomer_parents, 4).unwrap();
    Expected {
        pipeline: pipeline.run(BENZOATE).unwrap().output_mol_data,
        validation: validation
            .validate_with_report_all_failures(&mol("CC(=O)[O-].[Na+]"), true)
            .unwrap(),
        canonical: smiles(&enumerator.canonicalize(&pyridone).unwrap()),
        scored: smiles(
            &enumerator
                .canonicalize_with_score_func(&pyridone, |m| score_substructs(m).unwrap())
                .unwrap(),
        ),
        normalized: normalized.iter().map(smiles).collect(),
        tautomer_parents: tautomer_parents.iter().map(smiles).collect(),
    }
}

#[test]
fn shared_objects_on_threads() {
    let pipeline = Pipeline::default();
    let parents = Arc::new(AtomicUsize::new(0));
    let counter = parents.clone();
    let nested = ValidationMethod::new(move |mol, _| {
        counter.fetch_add(1, Ordering::Relaxed);
        vec![smiles(&charge_parent(mol).unwrap())]
    });
    let validation = MolVSValidation::new(&[FragmentValidation.into(), nested]).unwrap();
    let normalizations = Arc::new(AtomicUsize::new(0));
    let counter = normalizations.clone();
    let callback = TautomerEnumeratorCallback::new(move |mol, _| {
        counter.fetch_add(1, Ordering::Relaxed);
        Normalizer::default().normalize(mol).is_ok()
    });
    let mut enumerator = TautomerEnumerator::default();
    enumerator.set_callback(Some(&callback));

    let expected = work(&pipeline, &validation, &enumerator);
    assert_eq!(
        expected.validation,
        [
            "INFO: [FragmentValidation] sodium is present",
            "INFO: [FragmentValidation] acetate/acetic acid is present",
            "CC(=O)O",
        ]
    );
    assert_eq!(expected.canonical, "O=c1cccc[nH]1");
    assert_eq!(expected.scored, "O=C1CC=CC=N1");
    assert_eq!(expected.normalized, vec!["C[N+](=O)[O-]"; 8]);
    assert_eq!(expected.tautomer_parents, vec!["O=c1cccc[nH]1.[Na+]"; 8]);
    assert!(parents.load(Ordering::Relaxed) > 0);
    assert!(normalizations.load(Ordering::Relaxed) > 0);

    std::thread::scope(|scope| {
        for _ in 0..4 {
            let enumerator = enumerator.clone();
            let (pipeline, validation, expected) = (&pipeline, &validation, &expected);
            scope.spawn(move || {
                for _ in 0..5 {
                    let got = work(pipeline, validation, &enumerator);
                    assert_eq!(got.pipeline, expected.pipeline);
                    assert_eq!(got.validation, expected.validation);
                    assert_eq!(got.canonical, expected.canonical);
                    assert_eq!(got.scored, expected.scored);
                    assert_eq!(got.normalized, expected.normalized);
                    assert_eq!(got.tautomer_parents, expected.tautomer_parents);
                }
            });
        }
    });
}

#[test]
fn closures_called_directly() {
    // Called from Rust, outside any RDKit call, the closures hold no lock to
    // release.
    let callback = TautomerEnumeratorCallback::new(|mol, _| {
        smiles(&charge_parent(mol).unwrap()) == "O=c1cccc[nH]1"
    });
    let enumerator = TautomerEnumerator::default();
    let result = enumerator.enumerate(&mol("Oc1ccccn1")).unwrap();
    assert!(callback.call(&mol("O=c1cccc[nH]1"), &result).unwrap());
    let method = ValidationMethod::new(|mol, _| vec![standardize_smiles(&smiles(mol)).unwrap()]);
    assert_eq!(
        method.validate(&mol("C[N](=O)=O")).unwrap(),
        ["C[N+](=O)[O-]"]
    );
    // A score function given to pick_canonical, which matches nothing.
    let picked = enumerator
        .pick_canonical_with_score_func(&[mol("Oc1ccccn1"), mol("O=c1cccc[nH]1")], |m| {
            score_tautomer_of(m)
        })
        .unwrap();
    assert_eq!(smiles(&picked), "O=c1cccc[nH]1");
}

fn score_tautomer_of(mol: &Mol) -> i32 {
    TautomerEnumerator::score_tautomer(mol).unwrap()
}
