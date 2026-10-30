// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestSurf.py and UnitTestCrippen.py, with
// their reference files tests/data/first_200.tpsa.csv (RDKit's
// Data/NCI/first_200.tpsa.csv), PP_descrs_regress.csv, and Crippen.csv, and
// pins values from Python's rdkit.Chem.MolSurf and rdkit.Chem.Crippen in
// RDKit 2026.09.1, each computed on a new molecule, since Python reads the
// results of earlier calls back from the molecule.

use rdkit_cxx::chem::crippen::*;
use rdkit_cxx::chem::mol_surf::*;
use rdkit_cxx::chem::rd_mol_descriptors::{peoe_vsa, slogp_vsa, smr_vsa};
use rdkit_cxx::chem::{BondType, Mol, add_hs, mol_from_smiles};

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

/// Agrees to Python's assertAlmostEqual's 7 places.
fn almost_equal(a: f64, b: f64, what: &str) {
    assert!((a - b).abs() < 5e-8, "{what}: {a} != {b}");
}

fn almost_equal_all(a: &[f64], b: &[f64], what: &str) {
    assert_eq!(a.len(), b.len(), "{what}");
    for (x, y) in a.iter().zip(b) {
        almost_equal(*x, *y, what);
    }
}

/// The SMILES and TPSA of each molecule in first_200.tpsa.csv.
fn nci_200() -> Vec<(&'static str, f64)> {
    include_str!("data/first_200.tpsa.csv")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let (smiles, expected) = line.trim().split_once(',').unwrap();
            (smiles, expected.parse().unwrap())
        })
        .collect()
}

#[test]
fn tpsa_short() {
    let data = nci_200();
    assert_eq!(data.len(), 200);
    for (smiles, expected) in data {
        let m = mol(smiles);
        let calc = tpsa(&m).unwrap();
        assert!(
            (calc - expected).abs() <= 1e-4,
            "bad TPSA for SMILES {smiles}"
        );
        // With hydrogens in the graph (sf.net issue 1969745).
        let calc = tpsa(&add_hs(&m).unwrap()).unwrap();
        assert!(
            (calc - expected).abs() <= 1e-4,
            "bad TPSA with Hs for SMILES {smiles}"
        );
    }
    let m = mol("c1c[nH]cc1");
    assert_eq!(tpsa(&m).unwrap(), tpsa(&add_hs(&m).unwrap()).unwrap());
}

/// UnitTestSurf's TestCase_python: Python's reference implementations agree
/// with rdMolDescriptors' on the NCI molecules.
#[test]
fn reference_implementations_agree() {
    for (smiles, _) in nci_200() {
        let m = mol(smiles);
        almost_equal(labute_asa(&m).unwrap(), py_labute_asa(&m).unwrap(), smiles);
        almost_equal_all(&peoe_vsa(&m).unwrap(), &py_peoe_vsa(&m).unwrap(), smiles);
        almost_equal_all(&slogp_vsa(&m).unwrap(), &py_slogp_vsa(&m).unwrap(), smiles);
        almost_equal_all(&smr_vsa(&m).unwrap(), &py_smr_vsa(&m).unwrap(), smiles);
    }
}

/// UnitTestSurf's descriptor regression: TPSA (column 28) of
/// PP_descrs_regress.csv, to 1e-4, and LabuteASA (column 6), which RDKit
/// runs only in its long tests.
#[test]
fn regression_file() {
    for line in include_str!("data/PP_descrs_regress.csv")
        .lines()
        .filter(|line| !line.starts_with('#'))
    {
        let fields: Vec<&str> = line.split(',').collect();
        let m = mol(fields[0]);
        let col = |i: usize| fields[i].trim().parse::<f64>().unwrap();
        assert!((tpsa(&m).unwrap() - col(28)).abs() <= 1e-4, "{}", fields[0]);
        assert!(
            (labute_asa(&m).unwrap() - col(6)).abs() <= 1e-4,
            "{}",
            fields[0]
        );
    }
}

struct Case {
    smiles: &'static str,
    labute_asa: [f64; 2],
    py_labute_asa: [f64; 2],
    tpsa: [f64; 2],
    py_smr_vsa: &'static [f64],
    py_slogp_vsa: &'static [f64],
    py_peoe_vsa: &'static [f64],
    py_peoe_vsa_bins: &'static [f64],
    crippen: [f64; 4],
}

/// LabuteASA and pyLabuteASA with and without includeHs, TPSA with and
/// without includeSandP, the reference VSA vectors, pyPEOE_VSA_ with bins
/// [0.1, -0.1], and MolLogP and MolMR with and without includeHs.
const CASES: [Case; 4] = [
    Case {
        smiles: "c1ccccc1C(=O)O",
        labute_asa: [52.75210108613335, 51.928733633924615],
        py_labute_asa: [52.75210108613335, 51.928733633924615],
        tpsa: [37.3, 37.3],
        py_smr_vsa: &[
            9.901064578912528,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            35.89528683400505,
            0.0,
            0.0,
            5.969305287951849,
        ],
        py_slogp_vsa: &[
            0.0,
            11.075832682792555,
            0.0,
            0.0,
            10.357988675768818,
            30.33183534230805,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        py_peoe_vsa: &[
            5.106527394840706,
            0.0,
            4.794537184071822,
            0.0,
            0.0,
            18.19910120538483,
            12.13273413692322,
            0.0,
            5.563451491696996,
            0.0,
            0.0,
            0.0,
            0.0,
            5.969305287951849,
        ],
        py_peoe_vsa_bins: &[9.901064578912528, 0.0, 41.864592121956896],
        crippen: [1.3848, 33.40129999999999, 0.4718, 26.311299999999996],
    },
    Case {
        smiles: "CC[N+](C)(C)C.[Cl-]",
        labute_asa: [52.32191789476409, 51.71234295404747],
        py_labute_asa: [52.32191789476409, 51.71234295404747],
        tpsa: [0.0, 0.0],
        py_smr_vsa: &[
            16.8900640188746,
            0.0,
            0.0,
            0.0,
            6.923737199690624,
            27.687772353944144,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        py_slogp_vsa: &[
            12.407033161179976,
            32.17080321163877,
            0.0,
            0.0,
            6.923737199690624,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        py_peoe_vsa: &[
            16.8900640188746,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            6.923737199690624,
            0.0,
            27.687772353944144,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        py_peoe_vsa_bins: &[16.8900640188746, 0.0, 34.61150955363477],
        crippen: [
            -2.2834999999999965,
            28.573399999999985,
            -4.0055,
            13.775400000000001,
        ],
    },
    Case {
        smiles: "OP(=O)(O)OCC",
        labute_asa: [41.355816550600146, 40.322615841872846],
        py_labute_asa: [41.35581655060015, 40.322615841872846],
        tpsa: [66.76, 76.57000000000001],
        py_smr_vsa: &[
            18.8756186521481,
            0.0,
            0.0,
            0.0,
            6.923737199690624,
            6.606881964512918,
            0.0,
            0.0,
            0.0,
            7.822697123132171,
        ],
        py_slogp_vsa: &[
            0.0,
            16.393705170017867,
            9.088795446643147,
            0.0,
            6.923737199690624,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            7.822697123132171,
        ],
        py_peoe_vsa: &[
            9.78682320550495,
            4.5237471617118175,
            4.565048284931329,
            0.0,
            0.0,
            0.0,
            6.923737199690624,
            0.0,
            6.606881964512918,
            0.0,
            0.0,
            0.0,
            0.0,
            7.822697123132171,
        ],
        py_peoe_vsa_bins: &[18.8756186521481, 0.0, 21.353316287335712],
        crippen: [
            0.11559999999999998,
            23.670099999999994,
            0.03599999999999996,
            15.595100000000002,
        ],
    },
    Case {
        smiles: "C[S](=O)C",
        labute_asa: [28.478447345802074, 27.661714937277722],
        py_labute_asa: [28.478447345802074, 27.661714937277722],
        tpsa: [17.07, 36.28],
        py_smr_vsa: &[
            4.208898492164469,
            0.0,
            0.0,
            0.0,
            0.0,
            12.511538367022808,
            0.0,
            0.0,
            0.0,
            10.799569042904308,
        ],
        py_slogp_vsa: &[
            0.0,
            16.720436859187277,
            10.799569042904308,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        py_peoe_vsa: &[
            0.0,
            4.208898492164469,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            23.311107409927118,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        py_peoe_vsa_bins: &[4.208898492164469, 0.0, 23.311107409927118],
        crippen: [-0.005299999999999971, 19.990399999999998, -0.7433, 13.6484],
    },
];

#[test]
fn values_match_python() {
    for case in &CASES {
        let m = mol(case.smiles);
        let what = case.smiles;
        almost_equal(labute_asa(&m).unwrap(), case.labute_asa[0], what);
        almost_equal(
            labute_asa_with_include_hs(&m, false).unwrap(),
            case.labute_asa[1],
            what,
        );
        almost_equal(py_labute_asa(&m).unwrap(), case.py_labute_asa[0], what);
        almost_equal(
            py_labute_asa_with_include_hs(&m, false).unwrap(),
            case.py_labute_asa[1],
            what,
        );
        assert_eq!(tpsa(&m).unwrap(), case.tpsa[0], "{what}");
        assert_eq!(
            tpsa_with_include_s_and_p(&m, true).unwrap(),
            case.tpsa[1],
            "{what}"
        );
        almost_equal_all(&py_smr_vsa(&m).unwrap(), case.py_smr_vsa, what);
        almost_equal_all(&py_slogp_vsa(&m).unwrap(), case.py_slogp_vsa, what);
        almost_equal_all(&py_peoe_vsa(&m).unwrap(), case.py_peoe_vsa, what);
        almost_equal_all(
            &py_peoe_vsa_with_bins(&m, &[0.1, -0.1]).unwrap(),
            case.py_peoe_vsa_bins,
            what,
        );
        almost_equal(mol_log_p(&m).unwrap(), case.crippen[0], what);
        almost_equal(mol_mr(&m).unwrap(), case.crippen[1], what);
        almost_equal(
            mol_log_p_with_include_hs(&m, false).unwrap(),
            case.crippen[2],
            what,
        );
        almost_equal(
            mol_mr_with_include_hs(&m, false).unwrap(),
            case.crippen[3],
            what,
        );
    }
}

#[test]
fn vsa_bin_functions() {
    let m = mol("c1ccccc1C(=O)O");
    let smr = smr_vsa(&m).unwrap();
    let smr_fns = [
        smr_vsa1, smr_vsa2, smr_vsa3, smr_vsa4, smr_vsa5, smr_vsa6, smr_vsa7, smr_vsa8, smr_vsa9,
        smr_vsa10,
    ];
    assert_eq!(smr_fns.len(), mr_bins().len() + 1);
    for (f, v) in smr_fns.iter().zip(&smr) {
        assert_eq!(f(&m).unwrap(), *v);
    }
    let slogp = slogp_vsa(&m).unwrap();
    let slogp_fns = [
        slogp_vsa1,
        slogp_vsa2,
        slogp_vsa3,
        slogp_vsa4,
        slogp_vsa5,
        slogp_vsa6,
        slogp_vsa7,
        slogp_vsa8,
        slogp_vsa9,
        slogp_vsa10,
        slogp_vsa11,
        slogp_vsa12,
    ];
    assert_eq!(slogp_fns.len(), logp_bins().len() + 1);
    for (f, v) in slogp_fns.iter().zip(&slogp) {
        assert_eq!(f(&m).unwrap(), *v);
    }
    let peoe = peoe_vsa(&m).unwrap();
    let peoe_fns = [
        peoe_vsa1, peoe_vsa2, peoe_vsa3, peoe_vsa4, peoe_vsa5, peoe_vsa6, peoe_vsa7, peoe_vsa8,
        peoe_vsa9, peoe_vsa10, peoe_vsa11, peoe_vsa12, peoe_vsa13, peoe_vsa14,
    ];
    assert_eq!(peoe_fns.len(), chg_bins().len() + 1);
    for (f, v) in peoe_fns.iter().zip(&peoe) {
        assert_eq!(f(&m).unwrap(), *v);
    }
    assert_eq!(chg_bins()[0], -0.3);
    assert_eq!(logp_bins()[2], 0.0);
    assert_eq!(mr_bins()[8], 4.0);
    assert_eq!(
        bond_scale_facts(),
        [
            (BondType::Aromatic, 0.1),
            (BondType::Single, 0.0),
            (BondType::Double, 0.2),
            (BondType::Triple, 0.3)
        ]
    );
}

/// UnitTestCrippen's testLogP and testRepeat, against Crippen.csv, to 1e-5.
#[test]
fn crippen_reference_file() {
    for line in include_str!("data/Crippen.csv").lines() {
        let fields: Vec<&str> = line.split(',').collect();
        let m = mol(fields[0]);
        for _ in 0..2 {
            assert!((mol_log_p(&m).unwrap() - fields[1].parse::<f64>().unwrap()).abs() <= 1e-5);
            assert!((mol_mr(&m).unwrap() - fields[2].parse::<f64>().unwrap()).abs() <= 1e-5);
        }
    }
    // Issue 80: another descriptor's matching does not change the result.
    let m = mol("CCOC");
    let reference = mol_log_p(&m).unwrap();
    rdkit_cxx::chem::lipinski::nhoh_count(&m).unwrap();
    assert_eq!(mol_log_p(&m).unwrap(), reference);
}
