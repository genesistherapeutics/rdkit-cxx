// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestSpacialScore.py (RDKit 2025.09.6).

use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::spacial_score::{sps, sps_with_normalize};

const TEST_CASES: [(&str, f64, f64); 29] = [
    (r"C/C=C\C1C=CCC(Br)C1C2=CC(C#C)=CC=C2", 491.0, 27.28),
    ("C=CCCBr", 37.0, 7.40),
    ("CCCCC", 42.0, 8.40),
    ("CCCCBr", 42.0, 8.40),
    ("CC(C)CBr", 48.0, 9.60),
    ("C/C=C/CBr", 50.0, 10.00),
    ("BrC(C)(C)(C)", 60.0, 12.00),
    ("CC[C@@H](C)Br", 75.0, 15.00),
    ("c1ccccc1", 48.0, 8.00),
    ("c1cnccc1", 48.0, 8.00),
    ("C1=CC=CCC1", 112.0, 18.67),
    ("C1=CCCCC1", 128.0, 21.33),
    ("C1CCCCC1", 144.0, 24.00),
    ("C1CCCC1C", 153.0, 25.50),
    ("C1CCC1(C)C", 174.0, 29.00),
    (
        "O=C1COCCN1c4ccc(N3C[C@H](CNC(=O)c2ccc(Cl)s2)OC3=O)cc4",
        563.0,
        19.41,
    ),
    (
        "CC(C)CCC[C@@H](C)[C@H]3CC[C@H]4C2CC=C1C[C@@H](O)CC[C@]1(C)[C@H]2CC[C@]34C",
        1303.0,
        46.54,
    ),
    ("CC=CC", 38.0, 9.50),
    ("C/C=C/C", 38.0, 9.50),
    ("C=C1CCCCC1", 158.0, 22.57),
    ("C=C1CCCNC1", 158.0, 22.57),
    ("CC=C1CCCCC1", 167.0, 20.88),
    (r"C/C=C\1CCCCC1", 167.0, 20.88),
    ("CC=C1CCCNC1", 211.0, 26.38),
    ("CC=C(C1CCCNC1)C2CCCNC2", 485.0, 32.33),
    ("C/C=C(C1CCCNC1)/C2CCCNC2", 485.0, 32.33),
    ("CC=C(C1CCCNC1)C2CCCCN2", 511.0, 34.07),
    ("C/C=C(C1CCCNC1)/C2CCCCN2", 511.0, 34.07),
    ("CC1=CCCCC1", 151.0, 21.57),
];

#[test]
fn cases() {
    for (idx, (smiles, expected_sps, expected_nsps)) in TEST_CASES.into_iter().enumerate() {
        let m = mol_from_smiles(smiles).unwrap();
        assert_eq!(
            sps_with_normalize(&m, false).unwrap(),
            expected_sps,
            "entry {idx}"
        );
        let nsps = sps_with_normalize(&m, true).unwrap();
        // assertAlmostEqual(places=2).
        assert!((nsps - expected_nsps).abs() < 5e-3, "entry {idx}: {nsps}");
        assert_eq!(sps(&m).unwrap(), nsps);
    }
}

#[test]
fn no_heavy_atoms() {
    let m = mol_from_smiles("[HH]").unwrap();
    assert_eq!(sps_with_normalize(&m, false).unwrap(), 0.0);
    assert!(sps(&m).is_err());
}
