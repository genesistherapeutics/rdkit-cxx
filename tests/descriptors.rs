// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestDescriptors.py and the doctests of
// rdkit.Chem.Descriptors, and checks CalcMolDescriptors against Python's in
// RDKit 2026.03.6 on the molecules of tests/data/calc_mol_descriptors.tsv,
// which Python's CalcMolDescriptors wrote (repr of each value, None where a
// descriptor failed).

use std::collections::BTreeMap;

use rdkit_cxx::chem::descriptors::*;
use rdkit_cxx::chem::rd_mol_descriptors::bcut_2d;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

/// Agrees with Python's value to the last few bits, NaN with NaN: the
/// descriptors computed in Rust add in Python's order, but `Ipc`'s matrix
/// products sum in another order than NumPy's BLAS.
fn same(actual: f64, expected: f64) -> bool {
    (actual.is_nan() && expected.is_nan())
        || actual == expected
        || (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0)
}

/// Agrees with Python's value to the last bits another platform's build of
/// RDKit (fused multiply-adds, its math library) may change.
fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn calc_mol_descriptors_match_python() {
    let text = include_str!("data/calc_mol_descriptors.tsv");
    let mut lines = text.lines();
    let names: Vec<&str> = lines.next().unwrap().split('\t').skip(1).collect();
    let listed: Vec<&str> = desc_list().iter().map(|(name, _)| *name).collect();
    assert_eq!(listed, names, "descList's names and order");
    let mut rows = 0;
    for line in lines {
        let mut fields = line.split('\t');
        let smiles = fields.next().unwrap();
        let m = mol(smiles);
        let calculated = calc_mol_descriptors(&m);
        // In descList's order, as Python's dictionary has them.
        let calculated_names: Vec<&str> = calculated.iter().map(|(name, _)| *name).collect();
        assert_eq!(calculated_names, names);
        for ((name, value), (_, actual)) in names.iter().zip(fields).zip(calculated) {
            let expected = (value != "None").then(|| value.parse::<f64>().unwrap());
            match (actual, expected) {
                (Some(a), Some(e)) => assert!(same(a, e), "{smiles}: {name}: {a} != {e}"),
                _ => assert_eq!(actual, expected, "{smiles}: {name}"),
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 40);
}

#[test]
fn missing_val() {
    let m = mol("");
    let descs: BTreeMap<&str, Option<f64>> = calc_mol_descriptors(&m).into_iter().collect();
    assert_eq!(descs["BalabanJ"], None);
    let filled = calc_mol_descriptors_with_missing_val(&m, -1.0);
    let names: Vec<&str> = desc_list().iter().map(|(name, _)| *name).collect();
    let filled_names: Vec<&str> = filled.iter().map(|(name, _)| *name).collect();
    assert_eq!(filled_names, names);
    let filled: BTreeMap<&str, f64> = filled.into_iter().collect();
    assert_eq!(filled["BalabanJ"], -1.0);
    assert_eq!(filled["MolWt"], 0.0);
}

/// testGithub1287 and testBadAtomHandling: every descriptor runs on these.
#[test]
fn every_descriptor_runs() {
    for smiles in ["CCC", "CC[Pu]", "CC[*]"] {
        let m = mol(smiles);
        for (name, function) in desc_list() {
            assert!(function(&m).is_ok(), "SMILES: {smiles}; Descriptor: {name}");
        }
    }
}

#[test]
fn fp_density_morgan() {
    let close = |a: f64, b: f64| assert!((a - b).abs() < 5e-8, "{a} != {b}");
    close(fp_density_morgan2(&mol("C")).unwrap(), 1.0);
    close(fp_density_morgan2(&mol("CC")).unwrap(), 1.0);
    close(fp_density_morgan2(&mol("CCC")).unwrap(), 4.0 / 3.0);
    close(
        fp_density_morgan2(&mol(&"C".repeat(10))).unwrap(),
        8.0 / 10.0,
    );
    close(
        fp_density_morgan2(&mol(&"C".repeat(100))).unwrap(),
        8.0 / 100.0,
    );
    let m = mol("CCCc1ccccc1");
    let fpd1 = fp_density_morgan1(&m).unwrap();
    let fpd2 = fp_density_morgan2(&m).unwrap();
    let fpd3 = fp_density_morgan3(&m).unwrap();
    close(fpd1, 10.0 / 9.0);
    assert!(fpd1 < fpd2);
    assert!(fpd2 < fpd3);
    // Github #4567.
    assert_eq!(fp_density_morgan1(&mol("[HH]")).unwrap(), 0.0);
}

#[test]
fn vector_descriptors() {
    let m = mol("CCCc1ccccc1");
    let results = bcut_2d(&m).unwrap();
    let names = [
        "BCUT2D_MWHI",
        "BCUT2D_MWLOW",
        "BCUT2D_CHGHI",
        "BCUT2D_CHGLO",
        "BCUT2D_LOGPHI",
        "BCUT2D_LOGPLOW",
        "BCUT2D_MRHI",
        "BCUT2D_MRLOW",
    ];
    let list = desc_list();
    for (i, name) in names.iter().enumerate() {
        let (_, f) = list.iter().find(|(n, _)| n == name).unwrap();
        assert_eq!(f(&m).unwrap(), results[i]);
    }
    // A molecule BCUT2D refuses gives NaN, as Python's wrapper does.
    let (_, f) = list.iter().find(|(n, _)| *n == "BCUT2D_MWHI").unwrap();
    assert!(f(&mol("CC[*]")).unwrap().is_nan());
    assert!(!list.iter().any(|(n, _)| n.starts_with("AUTOCORR2D_")));
}

#[test]
fn get_mol_descriptors() {
    let descs: BTreeMap<&str, Option<f64>> =
        calc_mol_descriptors(&mol("CCCO")).into_iter().collect();
    assert!(descs.contains_key("MolLogP"));
    assert_eq!(descs["NumHDonors"], Some(1.0));
}

#[test]
fn doctests() {
    let close = |a: f64, b: f64| assert!((a - b).abs() < 1e-2, "{a} != {b}");
    assert_eq!(mol_wt(&mol("CC")).unwrap(), 30.07);
    close(mol_wt(&mol("[NH4+].[Cl-]")).unwrap(), 53.49);
    close(heavy_atom_mol_wt(&mol("CC")).unwrap(), 24.02);
    assert_eq!(heavy_atom_mol_wt(&mol("[NH4+].[Cl-]")).unwrap(), 49.46);
    close(exact_mol_wt(&mol("CC")).unwrap(), 30.04);
    close(exact_mol_wt(&mol("[13CH3]C")).unwrap(), 31.05);
    assert_eq!(
        exact_mol_wt_with_only_heavy(&mol("CC"), true).unwrap(),
        rdkit_cxx::chem::rd_mol_descriptors::calc_exact_mol_wt_with_only_heavy(&mol("CC"), true)
            .unwrap()
    );
    assert_eq!(
        mol_wt_with_only_heavy(&mol("CC"), true).unwrap(),
        heavy_atom_mol_wt(&mol("CC")).unwrap()
    );
    assert_eq!(num_valence_electrons(&mol("CC")).unwrap(), 14);
    assert_eq!(num_valence_electrons(&mol("C(=O)O")).unwrap(), 18);
    assert_eq!(num_valence_electrons(&mol("C(=O)[O-]")).unwrap(), 18);
    assert_eq!(num_valence_electrons(&mol("C(=O)")).unwrap(), 12);
    assert_eq!(num_radical_electrons(&mol("CC")), 0);
    assert_eq!(num_radical_electrons(&mol("C[CH3]")), 0);
    assert_eq!(num_radical_electrons(&mol("C[CH2]")), 1);
    assert_eq!(num_radical_electrons(&mol("C[CH]")), 2);
    assert_eq!(num_radical_electrons(&mol("C[C]")), 3);
}

#[test]
fn partial_charges_match_python() {
    // Python's Max/MinPartialCharge and Max/MinAbsPartialCharge, with repr.
    let m = mol("CC(=O)Nc1ccc(O)cc1");
    assert_close(max_partial_charge(&m).unwrap(), 0.22074743631242583);
    assert_close(min_partial_charge(&m).unwrap(), -0.5079643157925783);
    assert_close(max_abs_partial_charge(&m).unwrap(), 0.5079643157925783);
    assert_close(min_abs_partial_charge(&m).unwrap(), 0.22074743631242583);
    // RDKit has no Gasteiger parameters for plutonium, which gets NaN, and
    // Python's min and max pass NaN on.
    let m = mol("CC[Pu]");
    assert!(max_partial_charge(&m).unwrap().is_nan());
    assert!(min_partial_charge(&m).unwrap().is_nan());
}
