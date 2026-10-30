// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/EState unit tests (UnitTestEState.py,
// UnitTestTypes.py, UnitTestFingerprints.py, and UnitTestVSA.py, whose
// reference file is tests/data/EState_VSA.csv), with values from the papers
// they cite, and pins values from Python's modules in RDKit 2026.03.6.

use rdkit_cxx::chem::e_state::e_state_vsa;
use rdkit_cxx::chem::e_state::fingerprinter::fingerprint_mol;
use rdkit_cxx::chem::e_state::*;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

/// As UnitTestEState's `_validate`, to 0.01.
fn validate(data: &[(&str, &[f64])]) {
    for &(smi, ans) in data {
        let m = mol(smi);
        let inds = e_state_indices(&m).unwrap();
        assert_eq!(inds.len(), ans.len(), "{smi}");
        for (ind, a) in inds.iter().zip(ans) {
            assert!(
                (ind - a).abs() < 1e-2,
                "bad EStates for smiles: {smi}: {inds:?}"
            );
        }
        let fold = |f: fn(f64, f64) -> f64, abs: bool| {
            ans.iter()
                .map(|a| if abs { a.abs() } else { *a })
                .reduce(f)
                .unwrap()
        };
        assert!((max_e_state_index(&m).unwrap() - fold(f64::max, false)).abs() < 1e-2);
        assert!((min_e_state_index(&m).unwrap() - fold(f64::min, false)).abs() < 1e-2);
        assert!((max_abs_e_state_index(&m).unwrap() - fold(f64::max, true)).abs() < 1e-2);
        assert!((min_abs_e_state_index(&m).unwrap() - fold(f64::min, true)).abs() < 1e-2);
    }
}

#[test]
fn simple_molecules() {
    validate(&[
        ("CCCC", &[2.18, 1.32, 1.32, 2.18]),
        ("CCCCC", &[2.21, 1.34, 1.39, 1.34, 2.21]),
        ("CCCCCCC", &[2.24, 1.36, 1.42, 1.44, 1.42, 1.36, 2.24]),
        (
            "CCCCCCCCCC",
            &[2.27, 1.37, 1.44, 1.46, 1.47, 1.47, 1.46, 1.44, 1.37, 2.27],
        ),
    ]);
}

#[test]
fn isomers() {
    validate(&[
        ("CCCCCC", &[2.23, 1.36, 1.41, 1.41, 1.36, 2.23]),
        ("CCC(C)CC", &[2.23, 1.33, 0.94, 2.28, 1.33, 2.23]),
        ("CC(C)CCC", &[2.25, 0.90, 2.25, 1.38, 1.33, 2.22]),
        ("CC(C)(C)CC", &[2.24, 0.54, 2.24, 2.24, 1.27, 2.20]),
    ]);
}

#[test]
fn heteroatoms() {
    validate(&[
        (
            "CCCCOCCCC",
            &[2.18, 1.24, 1.21, 0.95, 5.31, 0.95, 1.21, 1.24, 2.18],
        ),
        (
            "CCC(C)OC(C)CC",
            &[2.15, 1.12, 0.43, 2.12, 5.54, 0.43, 2.12, 1.12, 2.15],
        ),
        (
            "CC(C)(C)OC(C)(C)C",
            &[2.07, -0.02, 2.07, 2.07, 5.63, -0.02, 2.07, 2.07, 2.07],
        ),
        ("CC(C)CC", &[2.22, 0.88, 2.22, 1.31, 2.20]),
        ("CC(C)CN", &[2.10, 0.66, 2.10, 0.81, 5.17]),
        ("CC(C)CO", &[1.97, 0.44, 1.97, 0.31, 8.14]),
        ("CC(C)CF", &[1.85, 0.22, 1.85, -0.19, 11.11]),
        ("CC(C)CCl", &[2.09, 0.65, 2.09, 0.78, 5.34]),
        ("CC(C)CBr", &[2.17, 0.80, 2.17, 1.11, 3.31]),
        ("CC(C)CI", &[2.21, 0.87, 2.21, 1.28, 2.38]),
        ("CC(N)C(=O)O", &[1.42, -0.73, 4.84, -0.96, 9.57, 7.86]),
        ("CCOCC", &[1.99, 0.84, 4.83, 0.84, 1.99]),
        // As in RDKit's test, this does not match the paper's values.
        ("CCSCC", &[2.17, 1.26, 1.96, 1.26, 2.17]),
        ("CC(=O)OC", &[1.36, -0.24, 9.59, 4.11, 1.35]),
        ("CC(=S)OC", &[1.73, 0.59, 4.47, 4.48, 1.56]),
    ]);
}

#[test]
fn aromatics() {
    validate(&[
        (
            "Fc1ccc(C)cc1",
            &[12.09, -0.17, 1.45, 1.75, 1.09, 1.93, 1.75, 1.45],
        ),
        (
            "Clc1ccc(C)cc1",
            &[5.61, 0.80, 1.89, 1.99, 1.24, 2.04, 1.99, 1.89],
        ),
        (
            "Brc1ccc(C)cc1",
            &[3.35, 1.14, 2.04, 2.07, 1.30, 2.08, 2.07, 2.04],
        ),
        (
            "Ic1ccc(C)cc1",
            &[2.30, 1.30, 2.10, 2.11, 1.32, 2.09, 2.11, 2.10],
        ),
    ]);
}

#[test]
fn principle_quantum_number() {
    let ranges = [
        (1, 2),
        (3, 10),
        (11, 18),
        (19, 36),
        (37, 54),
        (55, 86),
        (87, 120),
    ];
    for (principal_qn, (nmin, nmax)) in (1..).zip(ranges) {
        for n in nmin..=nmax {
            assert_eq!(get_principle_quantum_number(n), principal_qn);
        }
    }
}

#[test]
fn indices_match_python() {
    // Python's EState.EStateIndices, printed with repr.
    let inds = e_state_indices(&mol("CC(N)C(=O)[O-].[Na+]")).unwrap();
    let expected = [
        1.3564814814814814,
        -0.8425925925925928,
        4.773148148148148,
        -1.212962962962963,
        9.462962962962962,
        9.462962962962962,
        0.0,
    ];
    assert_eq!(inds, expected);
    assert!(max_e_state_index(&mol("")).is_err());
}

#[test]
fn atom_types_simple_molecules() {
    for (smi, ans) in [
        ("CC", vec!["sCH3", "sCH3"]),
        ("CCC", vec!["sCH3", "ssCH2", "sCH3"]),
        ("CCOC", vec!["sCH3", "ssCH2", "ssO", "sCH3"]),
        (
            "c1ccccc1[NH3+]",
            vec!["aaCH", "aaCH", "aaCH", "aaCH", "aaCH", "aasC", "sNH3"],
        ),
        (
            "c1ccccc1N",
            vec!["aaCH", "aaCH", "aaCH", "aaCH", "aaCH", "aasC", "sNH2"],
        ),
        ("C#C", vec!["tCH", "tCH"]),
        ("C=C=C", vec!["dCH2", "ddC", "dCH2"]),
        (
            "c1[nH]cnc1CC(N)C(O)=O",
            vec![
                "aaCH", "aaNH", "aaCH", "aaN", "aasC", "ssCH2", "sssCH", "sNH2", "dssC", "sOH",
                "dO",
            ],
        ),
        (
            "c1nc[n-]c1CC(N)C(O)=O",
            vec![
                "aaCH", "aaN", "aaCH", "aaN", "aasC", "ssCH2", "sssCH", "sNH2", "dssC", "sOH", "dO",
            ],
        ),
        (
            "c1cccc2c1cccc2",
            vec![
                "aaCH", "aaCH", "aaCH", "aaCH", "aaaC", "aaaC", "aaCH", "aaCH", "aaCH", "aaCH",
            ],
        ),
    ] {
        let types = type_atoms(&mol(smi)).unwrap();
        let expected: Vec<Vec<String>> = ans.iter().map(|t| vec![(*t).to_owned()]).collect();
        assert_eq!(types, expected, "bad types for SMILES: {smi}");
    }
}

#[test]
fn atom_type_patterns() {
    let patts = es_patterns();
    assert_eq!(patts.len(), 79);
    assert_eq!(patts[0].0, "sLi");
    assert_eq!(patts[78].0, "ssssPb");
    assert_eq!(build_patts().len(), 79);
    let custom = build_patts_with_raw_v(&[("sCH3", "[CD1H3]-*"), ("bad", "[C")]);
    assert_eq!(custom[0].as_ref().unwrap().0, "sCH3");
    assert!(custom[1].is_none());
    assert_eq!(type_atoms(&mol("")).unwrap(), Vec::<Vec<String>>::new());
}

/// As UnitTestFingerprints' `_validate`, on the nonzero entries, to 0.01.
#[test]
fn fingerprints() {
    for (smi, c, v) in [
        (
            "c1[nH]cnc1CC(N)C(O)=O",
            vec![1, 2, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![
                0.26, 3.12, -0.86, -1.01, 0.67, 5.25, 2.71, 3.84, 8.42, 10.26,
            ],
        ),
        (
            "NCCc1ccc(O)c(O)c1",
            vec![2, 3, 3, 1, 2],
            vec![1.26, 4.71, 0.75, 5.30, 17.97],
        ),
    ] {
        let (counts, vals) = fingerprint_mol(&mol(smi)).unwrap();
        assert_eq!(counts.len(), 79);
        let counts: Vec<i64> = counts.into_iter().filter(|&x| x != 0).collect();
        let vals: Vec<f64> = vals.into_iter().filter(|&x| x != 0.0).collect();
        assert_eq!(counts, c, "bad counts for SMILES: {smi}");
        assert_eq!(vals.len(), v.len(), "bad val len for SMILES: {smi}");
        for (a, b) in vals.iter().zip(&v) {
            assert!((a - b).abs() < 1e-2, "bad val for SMILES: {smi}: {vals:?}");
        }
    }
}

/// UnitTestVSA: every EState_VSA and VSA_EState descriptor of the
/// molecules in EState_VSA.csv, to 1e-4.
#[test]
fn vsa_reference_file() {
    let text = include_str!("data/EState_VSA.csv");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split(',').skip(1).collect();
    let mut rows = 0;
    for line in lines {
        let mut fields = line.split(',').map(str::trim);
        let smiles = fields.next().unwrap();
        let m = mol(smiles);
        let estate_vsa = e_state_vsa::e_state_vsa(&m).unwrap();
        let vsa_estate = e_state_vsa::vsa_e_state(&m).unwrap();
        for (name, value) in header.iter().zip(fields) {
            let expected: f64 = value.parse().unwrap();
            let calc = if let Some(bin) = name.strip_prefix("EState_VSA") {
                estate_vsa[bin.parse::<usize>().unwrap() - 1]
            } else {
                let bin = name.strip_prefix("VSA_EState").unwrap();
                vsa_estate[bin.parse::<usize>().unwrap() - 1]
            };
            assert!(
                (calc - expected).abs() <= 1e-4,
                "{smiles}: {name}: {calc:.4} != {expected:.4}"
            );
        }
        rows += 1;
    }
    assert_eq!(rows, 97);
}

/// Agrees with Python's values to the last bits another platform's build of
/// RDKit (fused multiply-adds) may change in the surface areas.
fn assert_all_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert!(
            (a - e).abs() <= 1e-12 * e.abs().max(1.0),
            "{actual:?} != {expected:?}"
        );
    }
}

#[test]
fn vsa_functions_match_python() {
    // Python's EState_VSA functions, printed with repr.
    let m = mol("c1ccccc1O");
    assert_eq!(e_state_vsa::e_state_vsa1(&m).unwrap(), 0.0);
    assert_eq!(e_state_vsa::vsa_bins().len(), 9);
    assert_eq!(e_state_vsa::estate_bins().len(), 10);
    let all = e_state_vsa::e_state_vsa(&m).unwrap();
    assert_all_close(
        &all,
        &[
            0.0,
            0.0,
            5.749511833283905,
            0.0,
            0.0,
            24.26546827384644,
            6.06636706846161,
            0.0,
            5.106527394840706,
            0.0,
            0.0,
        ],
    );
    let fns = [
        e_state_vsa::e_state_vsa1,
        e_state_vsa::e_state_vsa2,
        e_state_vsa::e_state_vsa3,
        e_state_vsa::e_state_vsa4,
        e_state_vsa::e_state_vsa5,
        e_state_vsa::e_state_vsa6,
        e_state_vsa::e_state_vsa7,
        e_state_vsa::e_state_vsa8,
        e_state_vsa::e_state_vsa9,
        e_state_vsa::e_state_vsa10,
        e_state_vsa::e_state_vsa11,
    ];
    for (f, value) in fns.iter().zip(&all) {
        assert_eq!(f(&m).unwrap(), *value);
    }
    let all = e_state_vsa::vsa_e_state(&m).unwrap();
    assert_all_close(
        &all,
        &[
            0.0,
            0.0,
            8.632222222222222,
            0.0,
            0.3217592592592595,
            8.712685185185185,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    );
    let fns = [
        e_state_vsa::vsa_e_state1,
        e_state_vsa::vsa_e_state2,
        e_state_vsa::vsa_e_state3,
        e_state_vsa::vsa_e_state4,
        e_state_vsa::vsa_e_state5,
        e_state_vsa::vsa_e_state6,
        e_state_vsa::vsa_e_state7,
        e_state_vsa::vsa_e_state8,
        e_state_vsa::vsa_e_state9,
        e_state_vsa::vsa_e_state10,
    ];
    for (f, value) in fns.iter().zip(&all) {
        assert_eq!(f(&m).unwrap(), *value);
    }
    // Unsorted bins bin as Python's bisect does.
    let custom = e_state_vsa::e_state_vsa_with_bins(&m, &[5.0, 1.0]).unwrap();
    assert_all_close(&custom, &[5.749511833283905, 0.0, 35.43836273714876]);
}
