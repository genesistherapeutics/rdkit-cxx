// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::{
    AsAtom, Atom, Mol, RWMol, SmilesParserParams, mol_from_smarts, mol_from_smiles,
    mol_from_smiles_with_params,
};

// Expected values come from Python's rdMolStandardize.

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn unsanitized(smiles: &str) -> Mol {
    let mut m = mol_from_smiles_with_params(
        smiles,
        &SmilesParserParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    m.update_property_cache_with_strict(false).unwrap();
    m
}

fn check(method: &ValidationMethod, m: &Mol, first: &[&str], all: &[&str]) {
    assert_eq!(method.validate(m).unwrap(), first);
    assert_eq!(
        method.validate_with_report_all_failures(m, true).unwrap(),
        all
    );
}

#[test]
fn rdkit_validation() {
    let validation = RDKitValidation;
    let message = "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is greater than \
                   permitted";
    check(
        &validation.into(),
        &unsanitized("CO(C)C"),
        &[message],
        &[message],
    );
    assert_eq!(
        validation.validate(&unsanitized("CO(C)C")).unwrap(),
        [message]
    );
    let empty = Mol::default();
    assert_eq!(
        validation.validate(&empty).unwrap(),
        ["ERROR: [NoAtomValidation] Molecule has no atoms"]
    );
    assert_eq!(
        NoAtomValidation.validate(&empty).unwrap(),
        ["ERROR: [NoAtomValidation] Molecule has no atoms"]
    );
    assert!(NoAtomValidation.validate(&mol("C")).unwrap().is_empty());
}

#[test]
fn fragment_neutral_and_isotope_validations() {
    let frags = [
        "INFO: [FragmentValidation] chlorine is present",
        "INFO: [FragmentValidation] dichloromethane is present",
    ];
    let m = mol("COC.ClCCl.Cl");
    assert_eq!(FragmentValidation.validate(&m).unwrap(), frags);
    assert_eq!(
        FragmentValidation
            .validate_with_report_all_failures(&m, true)
            .unwrap(),
        frags
    );
    assert_eq!(
        NeutralValidation.validate(&mol("CC(=O)[O-]")).unwrap(),
        ["INFO: [NeutralValidation] Not an overall neutral system (-1)"]
    );
    assert_eq!(
        NeutralValidation
            .validate(&mol("C[N+](C)(C)C.[Na+]"))
            .unwrap(),
        ["INFO: [NeutralValidation] Not an overall neutral system (+2)"]
    );
    let isotopes = mol("[2H]C[13CH3]");
    let validation = IsotopeValidation;
    check(
        &validation.into(),
        &isotopes,
        // RDKit 2024.03 reports every isotope, in sorted order, either way.
        &[
            "INFO: [IsotopeValidation] Molecule contains isotope 13C",
            "INFO: [IsotopeValidation] Molecule contains isotope 2H",
        ],
        &[
            "INFO: [IsotopeValidation] Molecule contains isotope 13C",
            "INFO: [IsotopeValidation] Molecule contains isotope 2H",
        ],
    );
}

#[test]
fn mol_vs_validation() {
    let default = MolVSValidation::default();
    assert_eq!(
        default
            .validate(&mol("COc1cccc(C=N[N-]C(N)=O)c1[O-].O.O.O.O=[U+2]=O"))
            .unwrap(),
        ["INFO: [FragmentValidation] water/hydroxide is present"]
    );
    let some = MolVSValidation::new(&[NoAtomValidation.into(), NeutralValidation.into()]).unwrap();
    assert_eq!(
        some.validate(&mol("[2H]C(=O)[O-]")).unwrap(),
        ["INFO: [NeutralValidation] Not an overall neutral system (-1)"]
    );
    assert!(matches!(
        MolVSValidation::new(&[]),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        validate_smiles("ClCCCl.c1ccccc1O").unwrap(),
        ["INFO: [FragmentValidation] 1,2-dichloroethane is present"]
    );
    assert_eq!(
        validate_smiles("[2H]C(=O)[O-]").unwrap(),
        [
            "INFO: [NeutralValidation] Not an overall neutral system (-1)",
            "INFO: [IsotopeValidation] Molecule contains isotope 2H"
        ]
    );
    assert!(matches!(validate_smiles("C1CC"), Err(Error::Value { .. })));
    assert!(matches!(
        validate_smiles("CO(C)(C)C"),
        Err(Error::AtomValence { .. })
    ));
}

#[test]
fn atom_validations() {
    let (c, n, o, f) = (
        Atom::from_num(6),
        Atom::from_num(7),
        Atom::from_num(8),
        Atom::from_num(9),
    );
    let allowed = AllowedAtomsValidation::new(&[&c, &n, &o]).unwrap();
    assert_eq!(
        allowed.validate(&mol("CC(=O)CF")).unwrap(),
        ["INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list"]
    );
    assert_eq!(
        allowed
            .validate_with_report_all_failures(&mol("CC(=O)CF.ClBr"), true)
            .unwrap(),
        [
            "INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list",
            "INFO: [AllowedAtomsValidation] Atom Cl is not in allowedAtoms list",
            "INFO: [AllowedAtomsValidation] Atom Br is not in allowedAtoms list",
        ]
    );
    let disallowed = DisallowedAtomsValidation::new(&[&n, &f]).unwrap();
    assert_eq!(
        disallowed
            .validate_with_report_all_failures(&mol("CC(=O)CF.CN"), true)
            .unwrap(),
        [
            "INFO: [DisallowedAtomsValidation] Atom F is in disallowedAtoms list",
            "INFO: [DisallowedAtomsValidation] Atom N is in disallowedAtoms list",
        ]
    );
    // A query atom, from a molecule.
    let query = mol_from_smarts("[#6;R]").unwrap();
    let ring_carbon = query.get_atom_with_idx(0).unwrap();
    let ring = AllowedAtomsValidation::new(&[&ring_carbon as &dyn AsAtom]).unwrap();
    assert_eq!(
        ring.validate(&mol("c1ccccc1C")).unwrap(),
        ["INFO: [AllowedAtomsValidation] Atom C is not in allowedAtoms list"]
    );
    assert!(matches!(
        AllowedAtomsValidation::new(&[]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        DisallowedAtomsValidation::new(&[]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn rust_validation_methods() {
    let calls = Arc::new(AtomicUsize::new(0));
    let counter = calls.clone();
    let method = ValidationMethod::new(move |m, report_all_failures| {
        counter.fetch_add(1, Ordering::SeqCst);
        vec![format!(
            "INFO: [V] {} atoms {}",
            m.get_num_atoms(),
            if report_all_failures { "True" } else { "False" }
        )]
    });
    check(
        &method,
        &mol("CC"),
        &["INFO: [V] 2 atoms False"],
        &["INFO: [V] 2 atoms True"],
    );
    let composite = MolVSValidation::new(&[method.clone(), NeutralValidation.into()]).unwrap();
    assert_eq!(
        composite.validate(&mol("C[O-]")).unwrap(),
        ["INFO: [V] 2 atoms False"]
    );
    assert_eq!(
        composite
            .validate_with_report_all_failures(&mol("C[O-]"), true)
            .unwrap(),
        [
            "INFO: [V] 2 atoms True",
            "INFO: [NeutralValidation] Not an overall neutral system (-1)"
        ]
    );
    assert_eq!(calls.load(Ordering::SeqCst), 4);

    // Clones share the closure, from any thread.
    std::thread::scope(|scope| {
        for _ in 0..4 {
            let composite = composite.clone();
            scope.spawn(move || {
                for _ in 0..10 {
                    composite.validate(&mol("CC")).unwrap();
                }
            });
        }
    });
    assert_eq!(calls.load(Ordering::SeqCst), 44);
    drop(composite);
    drop(method);
    // The closure is dropped with the last method that shares it.
    assert_eq!(Arc::strong_count(&calls), 1);
}

#[test]
fn rust_validation_methods_can_match_substructures() {
    // Validating holds control_c_mutex, released while a Rust method runs,
    // so the method can call wrappers that take it.
    let carboxyls = ValidationMethod::new(|m, _| {
        let query = mol_from_smarts("C(=O)[OH]").unwrap();
        let n = m.get_substruct_matches(&query).unwrap().len();
        vec![format!("INFO: [V] {n} carboxyl groups")]
    });
    let composite = MolVSValidation::new(&[FragmentValidation.into(), carboxyls]).unwrap();
    assert_eq!(
        composite
            .validate_with_report_all_failures(&mol("OC(=O)CC(=O)O.Cl"), true)
            .unwrap(),
        [
            "INFO: [FragmentValidation] chlorine is present",
            "INFO: [V] 2 carboxyl groups"
        ]
    );
}

#[test]
fn panics_in_validation_methods_are_errors() {
    let method = ValidationMethod::new(|_, _| panic!("no validation"));
    let err = method.validate(&mol("C")).unwrap_err();
    assert!(
        matches!(&err, Error::Runtime { message } if message.contains("no validation")),
        "{err:?}"
    );
    let composite = MolVSValidation::new(&[method]).unwrap();
    assert!(matches!(
        composite.validate(&mol("C")),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn large_charge_sums_are_refused() {
    // NeutralValidation sums the formal charges in int (U9): 2^24 atoms of
    // charge -128 overflow its negation.
    let mut atom = Atom::from_num(6);
    atom.set_formal_charge(-128);
    let mut big = RWMol::default();
    for _ in 0..1 << 24 {
        big.add_atom(&atom);
    }
    let is_value =
        |result: rdkit_cxx::Result<Vec<String>>| matches!(result, Err(Error::Value { .. }));
    assert!(is_value(NeutralValidation.validate(&big)));
    assert!(is_value(MolVSValidation::default().validate(&big)));
    let composite = MolVSValidation::new(&[NeutralValidation.into()]).unwrap();
    assert!(is_value(composite.validate(&big)));
    // Validations that do not sum charges are not refused.
    assert!(NoAtomValidation.validate(&big).unwrap().is_empty());
}
