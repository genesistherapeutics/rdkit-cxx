// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestQED.py, with its reference file
// tests/data/NCI_200_qed.csv (rdkit/Chem/test_data/QED), and the doctests of
// rdkit.Chem.QED in RDKit 2025.09.6.

use rdkit_cxx::chem::qed::*;
use rdkit_cxx::chem::{Mol, add_hs, mol_from_smiles};

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

/// Agrees to Python's assertAlmostEqual's 7 places.
fn almost_equal(a: f64, b: f64, what: &str) {
    assert!((a - b).abs() < 5e-8, "{what}: {a} != {b}");
}

#[test]
fn nci_200() {
    let mut rows = 0;
    for line in include_str!("data/NCI_200_qed.csv")
        .lines()
        .filter(|line| !line.starts_with('#'))
    {
        let (smiles, expected) = line.trim().split_once(',').unwrap();
        let expected: f64 = expected.parse().unwrap();
        let m = mol(smiles);
        almost_equal(qed(&m).unwrap(), expected, smiles);
        // Hydrogens in the graph do not change the result.
        almost_equal(qed(&add_hs(&m).unwrap()).unwrap(), expected, smiles);
        rows += 1;
    }
    assert_eq!(rows, 200);
}

#[test]
fn properties_match_python() {
    let m = mol("N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O");
    for m in [m.clone(), add_hs(&m).unwrap()] {
        let p = properties(&m).unwrap();
        almost_equal(p.mw, 337.456, "MW");
        almost_equal(p.alogp, -0.55833, "ALOGP");
        assert_eq!(p.hba, 6.0);
        assert_eq!(p.hbd, 5.0);
        almost_equal(p.psa, 173.33, "PSA");
        assert_eq!(p.rotb, 7.0);
        assert_eq!(p.arom, 1.0);
        assert_eq!(p.alerts, 3.0);
    }
}

#[test]
fn examples() {
    for (smiles, expected) in [
        // Paroxetine, 0.935 in the paper.
        ("c1cc2OCOc2cc1OCC1CNCCC1c1ccc(F)cc1", 0.934),
        // Leflunomide, 0.929.
        ("C1=NOC(C)=C1C(=O)Nc1ccc(cc1)C(F)(F)F", 0.911),
        // Clomipramine, 0.779.
        ("CN(C)CCCN1c2ccccc2CCc2ccc(Cl)cc21", 0.818),
        // Tegaserod, 0.213.
        ("CCCCCNC(=N)NN=CC1=CNc2ccc(CO)cc21", 0.235),
    ] {
        let value = qed(&mol(smiles)).unwrap();
        assert!((value - expected).abs() < 5e-4, "{smiles}: {value}");
    }
}

#[test]
fn doctests() {
    for (smiles, prefix) in [
        ("N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O", "0.253"),
        ("CNC(=NCCSCc1nc[nH]c1C)NC#N", "0.234"),
        ("CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12", "0.234"),
    ] {
        let value = qed(&mol(smiles)).unwrap().to_string();
        assert!(value.starts_with(prefix), "{smiles}: {value}");
    }
}

/// Agrees with Python's value to the last bits that the platform's `exp`
/// and `log` may change.
fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn weights_match_python() {
    // Python's weights_max, weights_mean, weights_none, and default, printed
    // with repr.
    let m = mol("c1cc2OCOc2cc1OCC1CNCCC1c1ccc(F)cc1");
    assert_close(weights_max(&m).unwrap(), 0.9245739032271649);
    assert_close(weights_mean(&m).unwrap(), 0.9339413767145557);
    assert_close(weights_none(&m).unwrap(), 0.937329699400303);
    assert_close(default(&m).unwrap(), 0.9339413767145557);
    let given = QedParams {
        qed_properties: Some(QEDproperties {
            mw: 300.0,
            alogp: 2.0,
            hba: 3.0,
            hbd: 1.0,
            psa: 60.0,
            rotb: 4.0,
            arom: 2.0,
            alerts: 0.0,
        }),
        ..QedParams::default()
    };
    assert_close(qed_with_params(&m, &given).unwrap(), 0.9391885840701519);
    assert_close(ads(2.0, &ads_parameters()["ALOGP"]), 0.9559756451041136);
}

#[test]
fn constants_match_python() {
    assert_eq!(QedParams::default().w, weight_mean());
    assert_eq!(QedParams::default().qed_properties, None);
    assert_eq!(weight_max().alerts, 1.0);
    assert_eq!(weight_none().mw, 1.0);
    assert_eq!(acceptor_smarts().len(), 11);
    assert_eq!(acceptors().len(), 11);
    assert_eq!(structural_alert_smarts().len(), 116);
    assert_eq!(structural_alerts().len(), 116);
    assert_eq!(
        structural_alert_smarts()[26],
        "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]"
    );
    let params = ads_parameters();
    assert_eq!(params.len(), 8);
    assert_eq!(
        params["MW"],
        ADSparameter {
            a: 2.817065973,
            b: 392.5754953,
            c: 290.7489764,
            d: 2.419764353,
            e: 49.22325677,
            f: 65.37051707,
            dmax: 104.9805561,
        }
    );
    assert_eq!(aliphatic_rings().get_num_atoms(), 1);
}
