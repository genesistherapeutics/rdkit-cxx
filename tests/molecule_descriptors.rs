// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/ML/Descriptors/UnitTestMolDescriptors.py, and pins
// values from Python's MolecularDescriptorCalculator in RDKit 2026.03.6.
// Also tests the ports of rdkit.Chem.ReducedGraphs and
// rdkit.Chem.ChemUtils.DescriptorUtilities.

use rdkit_cxx::chem::chem_utils::descriptor_utilities::VectorDescriptorWrapper;
use rdkit_cxx::chem::descriptors::desc_list;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_mol_descriptors::bcut_2d;
use rdkit_cxx::chem::reduced_graphs::{get_erg_fingerprint, tanimoto_similarity};
use rdkit_cxx::ml::descriptors::molecule_descriptors::MolecularDescriptorCalculator;

/// Agrees with Python's value, NaN with NaN, to the last bits that BCUT2D's
/// eigensolver may change across platforms.
fn same(actual: f64, expected: f64) -> bool {
    (actual.is_nan() && expected.is_nan())
        || (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0)
}

#[test]
fn calc_vals() {
    let calc = MolecularDescriptorCalculator::new(&["MolLogP", "Chi1v"]);
    assert_eq!(calc.get_descriptor_names(), ["MolLogP", "Chi1v"]);
    for (smi, vals) in [
        ("CCOC", [0.6527, 1.40403]),
        ("CC=O", [0.2052, 0.81305]),
        ("CCC(=O)O", [0.481, 1.48839]),
    ] {
        let res = calc.calc_descriptors(&mol_from_smiles(smi).unwrap());
        for (r, v) in res.iter().zip(vals) {
            assert!(
                (r - v).abs() < 1e-4,
                "bad descriptor values for SMILES {smi}"
            );
        }
    }
    assert_eq!(calc.get_descriptor_versions(), ["1.2.0", "1.2.0"]);
    assert_eq!(
        calc.get_descriptor_summaries(),
        ["Wildman-Crippen LogP value", "N/A"]
    );
    assert_eq!(calc.get_descriptor_funcs().len(), 2);
}

#[test]
fn matches_python() {
    let names = [
        "MolLogP",
        "Chi1v",
        "TPSA",
        "NotADescriptor",
        "AUTOCORR2D_5",
        "NumHDonors",
        "fr_ArN",
        "BCUT2D_MWHI",
        "FractionCSP3",
    ];
    let calc = MolecularDescriptorCalculator::new(&names);
    assert_eq!(calc.get_descriptor_names(), names);
    assert_eq!(
        calc.get_descriptor_versions(),
        [
            "1.2.0", "1.2.0", "2.0.0", "N/A", "1.0.0", "1.0.0", "N/A", "1.0.0", "1.0.0"
        ]
    );
    assert_eq!(
        calc.get_descriptor_summaries(),
        [
            "Wildman-Crippen LogP value",
            "N/A",
            "N/A",
            "N/A",
            "N/A",
            "Number of Hydrogen Bond Donors",
            "Number of N functional groups attached to aromatics",
            "N/A",
            "CalcFractionCSP3( (rdkit.Chem.rdchem.Mol)mol) -> float : returns the fraction of C \
             atoms that are SP3 hybridized",
        ]
    );
    let values = calc.calc_descriptors(&mol_from_smiles("CCC(=O)O").unwrap());
    let expected = [
        0.481,
        1.4883911147617317,
        37.3,
        777.0,
        0.0,
        1.0,
        0.0,
        16.365026720551658,
        0.6666666666666666,
    ];
    for (v, e) in values.iter().zip(expected) {
        assert!(same(*v, e), "{v} != {e}");
    }
    // Python reads RDKit's autocorrelation tables past their end for a
    // dummy atom (CONVENTIONS.md U1), and gets 0 for AUTOCORR2D_5; here the
    // function refuses the molecule, so the descriptor is NaN.
    let values = calc.calc_descriptors(&mol_from_smiles("CC[*]").unwrap());
    let expected = [
        0.9739,
        0.7071067811865475,
        0.0,
        777.0,
        f64::NAN,
        0.0,
        0.0,
        f64::NAN,
        1.0,
    ];
    for (v, e) in values.iter().zip(expected) {
        assert!(same(*v, e), "{v} != {e}");
    }
}

/// Github #3511: every descriptor of descList runs.
#[test]
fn every_descriptor() {
    let names: Vec<&str> = desc_list().iter().map(|(name, _)| *name).collect();
    let calc = MolecularDescriptorCalculator::new(&names);
    let values = calc.calc_descriptors(&mol_from_smiles("C").unwrap());
    assert_eq!(values.len(), 217);
    assert!(values.iter().all(|&v| v != -666.0 && v != 777.0));
    // Python's fragment counts, lambdas, have no version.
    let versions = calc.get_descriptor_versions();
    assert_eq!(versions.iter().filter(|v| *v == "N/A").count(), 85);
}

#[test]
fn tanimoto_similarity_matches_python() {
    let a = get_erg_fingerprint(&mol_from_smiles("c1ccccc1CCN").unwrap()).unwrap();
    let b = get_erg_fingerprint(&mol_from_smiles("c1ccccc1CCO").unwrap()).unwrap();
    assert!((tanimoto_similarity(&a, &b).unwrap() - 0.6666666666666667).abs() < 1e-15);
    assert_eq!(tanimoto_similarity(&a, &a).unwrap(), 1.0);
    assert_eq!(tanimoto_similarity(&[0.0, 0.0], &[0.0, 0.0]).unwrap(), 0.0);
    assert!(tanimoto_similarity(&[1.0], &[1.0, 2.0]).is_err());
}

#[test]
fn vector_descriptor_wrapper() {
    let names: Vec<String> = (0..8).map(|i| format!("BCUT{i}")).collect();
    let wrapper = VectorDescriptorWrapper::new(bcut_2d, names.clone(), "1.0.0");
    assert_eq!(wrapper.names(), names);
    assert_eq!(wrapper.version(), "1.0.0");
    let m = mol_from_smiles("CCCc1ccccc1").unwrap();
    let results = bcut_2d(&m).unwrap();
    for (i, expected) in results.iter().enumerate() {
        assert_eq!(wrapper.call_desc(&m, i as u32).unwrap(), *expected);
    }
    assert!(wrapper.call_desc(&m, 8).is_err());
    // The function fails for a dummy atom: NaN, as in Python.
    let dummy = mol_from_smiles("CC[*]").unwrap();
    assert!(wrapper.call_desc(&dummy, 0).unwrap().is_nan());
}
