// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestGraphDescriptors_2.py, whose values come
// from the papers it cites and from tests/data/PP_descrs_regress.csv, and
// pins values from Python's rdkit.Chem.GraphDescriptors in RDKit 2025.09.6.

use rdkit_cxx::chem::graph_descriptors::*;
use rdkit_cxx::chem::{
    GetDistanceMatrixParams, Mol, SanitizeFlags, SanitizeMolParams, SmilesParserParams,
    get_distance_matrix, get_distance_matrix_with_params, mol_from_smiles,
    mol_from_smiles_with_params, sanitize_mol_with_params,
};

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

/// Agrees with Python's value to the last few bits: the characteristic
/// polynomial's matrix products sum in another order than NumPy's BLAS.
fn assert_close(actual: f64, expected: f64, what: &str) {
    let tol = 1e-12 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: {actual} != {expected}"
    );
}

fn feq(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: {actual} != {expected}"
    );
}

/// The molecule sanitized without the property checks, as RDKit's kappa
/// tests do for SMILES with bogus valences.
fn loose_mol(smiles: &str) -> Mol {
    let mut m = mol_from_smiles_with_params(
        smiles,
        &SmilesParserParams {
            sanitize: false,
            ..SmilesParserParams::default()
        },
    )
    .unwrap();
    m.update_property_cache_with_strict(false).unwrap();
    sanitize_mol_with_params(
        &mut m,
        &SanitizeMolParams {
            sanitize_ops: SanitizeFlags::from_bits_retain(
                SanitizeFlags::SANITIZE_ALL.bits() ^ SanitizeFlags::SANITIZE_PROPERTIES.bits(),
            ),
            ..SanitizeMolParams::default()
        },
    )
    .unwrap();
    m
}

#[test]
fn bertz_ct_values() {
    for (smi, expected) in [
        ("C=CC=C", 21.01955),
        ("O=CC=O", 25.01955),
        ("FCC(=O)CF", 46.7548875),
        ("O=C1C=CC(=O)C=C1", 148.705216),
        ("C12C(F)=C(O)C(F)C1C(F)=C(O)C(F)2", 315.250442),
        ("C12CC=CCC1C(=O)C3CC=CCC3C(=O)2", 321.539522),
    ] {
        feq(super_bertz(smi), expected, 1e-3, smi);
    }
    // Github #8421.
    feq(
        super_bertz("Fc1cccc(C2(c3nnc(Cc4cccc5ccccc45)o3)CCOCC2)c1"),
        1143.0567597744753,
        1e-4,
        "Github #8421",
    );
    // RDKit 2025.09's Python takes a bond's order as its BondType's number:
    // 17 for a dative bond.
    feq(
        super_bertz("C[NH2]->[Cu]"),
        1192.1337526548132,
        1e-9,
        "dative",
    );
}

fn super_bertz(smi: &str) -> f64 {
    bertz_ct(&mol(smi)).unwrap()
}

#[test]
fn ipc_values() {
    for (smi, res1, res2) in [
        ("CCCCC", 1.40564, 11.24511),
        ("CCC(C)C", 1.37878, 9.65148),
        ("CC(C)(C)C", 0.72193, 3.60964),
        ("CN(CC)CCC", 1.67982, 31.91664),
        ("C1CCCCC1", 1.71997, 34.39946),
        ("CC1CCCCC1", 1.68562, 47.19725),
        ("Cc1ccccc1", 1.68562, 47.19725),
        ("CC(C)=C(C)C", 1.36096, 13.60964),
        ("C#N", 1.00000, 2.00000),
        ("OC#N", 0.91830, 2.75489),
    ] {
        let m = mol(smi);
        let forced = IpcParams {
            force_d_mat: true,
            ..IpcParams::default()
        };
        let value = ipc_with_params(&m, &forced).unwrap();
        let avg = ipc_with_params(
            &m,
            &IpcParams {
                avg: true,
                ..forced.clone()
            },
        )
        .unwrap();
        feq(avg, res1, 1e-3, smi);
        let avg_ipc_forced = avg_ipc_with_params(
            &m,
            &AvgIpcParams {
                force_d_mat: true,
                ..AvgIpcParams::default()
            },
        )
        .unwrap();
        assert_eq!(avg, avg_ipc_forced);
        feq(value, res2, 1e-3, smi);
        feq(ipc(&m).unwrap(), res2, 1e-3, smi);
        feq(avg_ipc(&m).unwrap(), res1, 1e-3, smi);
    }
}

#[test]
fn kappas() {
    for (smi, res) in [
        ("C12CC2C3CC13", 2.344),
        ("C1CCC12CC2", 3.061),
        ("C1CCCCC1", 4.167),
        ("CCCCCC", 6.000),
        ("CCC(C)C1CCC(C)CC1", 9.091),
        ("CC(C)CC1CCC(C)CC1", 9.091),
        ("CC(C)C1CCC(C)CCC1", 9.091),
    ] {
        feq(kappa1(&mol(smi)).unwrap(), res, 1e-3, smi);
    }
    for (smi, res) in [
        ("[C+2](C)(C)(C)(C)(C)C", 0.667),
        ("[C+](C)(C)(C)(C)(CC)", 1.240),
        ("C(C)(C)(C)(CCC)", 2.3444),
        ("CC(C)CCCC", 4.167),
        ("CCCCCCC", 6.000),
        ("CCCCCC", 5.000),
        ("C1CCCC1", 1.440),
        ("C1CCCC1C", 1.633),
        ("C1CCCCC1", 2.222),
        ("C1CCCCCC1", 3.061),
        ("CCCCC", 4.00),
        ("CC=CCCC", 4.740),
        ("C1=CN=CN1", 0.884),
        ("c1ccccc1", 1.606),
        ("c1cnccc1", 1.552),
        ("n1ccncc1", 1.500),
        ("CCCCF", 3.930),
        ("CCCCCl", 4.290),
        ("CCCCBr", 4.480),
        ("CCC(C)C1CCC(C)CC1", 4.133),
        ("CC(C)CC1CCC(C)CC1", 4.133),
        ("CC(C)C1CCC(C)CCC1", 4.133),
    ] {
        feq(kappa2(&loose_mol(smi)).unwrap(), res, 1e-3, smi);
    }
    for (smi, res) in [
        ("C[C+](C)(C)(C)C(C)(C)C", 2.000),
        ("CCC(C)C(C)(C)(CC)", 2.380),
        ("CCC(C)CC(C)CC", 4.500),
        ("CC(C)CCC(C)CC", 5.878),
        ("CC(C)CCCC(C)C", 8.000),
        ("CCC(C)C1CCC(C)CC1", 2.500),
        ("CC(C)CC1CCC(C)CC1", 3.265),
        ("CC(C)C1CCC(C)CCC1", 2.844),
    ] {
        feq(kappa3(&loose_mol(smi)).unwrap(), res, 1e-3, smi);
    }
}

#[test]
fn balaban_j_values() {
    for (smi, res) in [
        ("CC", 1.0),
        ("CCC", 1.6330),
        ("CCCC", 1.9747),
        ("CC(C)C", 2.3238),
        ("CCCCC", 2.1906),
        ("CC(C)CC", 2.5396),
        ("CC(C)(C)C", 3.0237),
        ("CCCCCC", 2.3391),
        ("CC(C)CCC", 2.6272),
        ("CCC(C)CC", 2.7542),
        ("CC(C)(C)CC", 3.1685),
        ("CC(C)C(C)C", 2.9935),
        ("C1CCCCC1", 2.0000),
        ("C1C(C)CCCC1", 2.1229),
        ("C1C(CC)CCCC1", 2.1250),
        ("C1C(C)C(C)CCC1", 2.2794),
        ("C1C(C)CC(C)CC1", 2.2307),
        ("C1C(C)CCC(C)C1", 2.1924),
        ("C1C(CCC)CCCC1", 2.0779),
        ("C1C(C(C)C)CCCC1", 2.2284),
        ("C1C(CC)C(C)CCC1", 2.2973),
        ("C1C(CC)CC(C)CC1", 2.2317),
        ("C1C(CC)CCC(C)C1", 2.1804),
        ("C1C(C)C(C)C(C)CC1", 2.4133),
        ("C1C(C)C(C)CC(C)C1", 2.3462),
        ("C1C(C)CC(C)CC1(C)", 2.3409),
        ("c1ccccc1", 3.0000),
        ("c1c(C)cccc1", 3.0215),
        ("c1c(CC)cccc1", 2.8321),
        ("c1c(C)c(C)ccc1", 3.1349),
        ("c1c(C)cc(C)cc1", 3.0777),
        ("c1c(C)ccc(C)c1", 3.0325),
        ("c1c(CCC)cccc1", 2.6149),
        ("c1c(C(C)C)cccc1", 2.8483),
        ("c1c(CC)c(C)ccc1", 3.0065),
        ("c1c(CC)cc(C)cc1", 2.9369),
        ("c1c(CC)ccc(C)c1", 2.8816),
        ("c1c(C)c(C)c(C)cc1", 3.2478),
        ("c1c(C)c(C)cc(C)c1", 3.1717),
        ("c1c(C)cc(C)cc1(C)", 3.1657),
    ] {
        let m = mol(smi);
        let forced = BalabanJParams {
            force_d_mat: true,
            ..BalabanJParams::default()
        };
        feq(balaban_j_with_params(&m, &forced).unwrap(), res, 1e-4, smi);
        feq(balaban_j(&m).unwrap(), res, 1e-4, smi);
    }
}

#[test]
fn order_dependence() {
    // RDKit's testOrderDepend and testIssue125: no call changes another's
    // result, as Python's cached matrices once did.
    for (smi, ct, bal) in [
        ("C=CC=C", 21.01955, 2.73205),
        ("O=CC=O", 25.01955, 2.73205),
        ("FCC(=O)CF", 46.7548875, 2.98816),
        ("O=C1C=CC(=O)C=C1", 148.705216, 2.8265),
        ("C12C(F)=C(O)C(F)C1C(F)=C(O)C(F)2", 315.250442, 2.4509),
        ("C12CC=CCC1C(=O)C3CC=CCC3C(=O)2", 321.539522, 1.95986),
    ] {
        let m = mol(smi);
        for _ in 0..2 {
            feq(balaban_j(&m).unwrap(), bal, 1e-4, smi);
            feq(bertz_ct(&m).unwrap(), ct, 1e-4, smi);
        }
    }
}

#[test]
fn connectivity_indices() {
    type Case = (&'static str, f64);
    let check = |f: fn(&Mol) -> rdkit_cxx::Result<f64>, data: &[Case], what: &str| {
        for &(smi, res) in data {
            feq(f(&mol(smi)).unwrap(), res, 1e-3, &format!("{what} {smi}"));
        }
    };
    let chi0v_data = [
        ("CCCCCC", 4.828),
        ("CCC(C)CC", 4.992),
        ("CC(C)CCC", 4.992),
        ("CC(C)C(C)C", 5.155),
        ("CC(C)(C)CC", 5.207),
        ("CCCCCO", 4.276),
        ("CCC(O)CC", 4.439),
        ("CC(O)(C)CC", 4.654),
        ("c1ccccc1O", 3.834),
        ("CCCl", 2.841),
        ("CCBr", 3.671),
        ("CCI", 4.242),
    ];
    check(chi0v, &chi0v_data, "Chi0v");
    let chi1 = [
        ("CCCCCC", 2.914),
        ("CCC(C)CC", 2.808),
        ("CC(C)CCC", 2.770),
        ("CC(C)C(C)C", 2.643),
        ("CC(C)(C)CC", 2.561),
        ("CCCCCO", 2.523),
        ("CCC(O)CC", 2.489),
        ("CC(O)(C)CC", 2.284),
        ("c1ccccc1O", 2.134),
    ];
    check(chi1v, &chi1, "Chi1v");
    check(chi1n, &chi1, "Chi1n");
    let chi2 = [
        ("CCCCCC", 1.707),
        ("CCC(C)CC", 1.922),
        ("CC(C)CCC", 2.183),
        ("CC(C)C(C)C", 2.488),
        ("CC(C)(C)CC", 2.914),
        ("CCCCCO", 1.431),
        ("CCC(O)CC", 1.470),
        ("CC(O)(C)CC", 2.166),
        ("c1ccccc1O", 1.336),
    ];
    check(chi2v, &chi2, "Chi2v");
    check(chi2n, &chi2, "Chi2n");
    let chi3 = [
        ("CCCCCC", 0.957),
        ("CCC(C)CC", 1.394),
        ("CC(C)CCC", 0.866),
        ("CC(C)C(C)C", 1.333),
        ("CC(C)(C)CC", 1.061),
        ("CCCCCO", 0.762),
        ("CCC(O)CC", 0.943),
        ("CC(O)(C)CC", 0.865),
        ("c1ccccc1O", 0.756),
    ];
    check(chi3v, &chi3, "Chi3v");
    check(chi3n, &chi3, "Chi3n");
    let chi4 = [
        ("CCCCCC", 0.500),
        ("CCC(C)CC", 0.289),
        ("CC(C)CCC", 0.577),
        ("CC(C)C(C)C", 0.000),
        ("CC(C)(C)CC", 0.000),
        ("CCCCCO", 0.362),
        ("CCC(O)CC", 0.289),
        ("CC(O)(C)CC", 0.000),
        ("c1ccccc1O", 0.428),
    ];
    check(chi4v, &chi4, "Chi4v");
    check(chi4n, &chi4, "Chi4n");
    for (smi, res) in [
        ("CCCCCC", 0.250),
        ("CCC(C)CC", 0.000),
        ("CC(C)CCC", 0.000),
        ("CC(C)C(C)C", 0.000),
        ("CC(C)(C)CC", 0.000),
        ("CCCCCO", 0.112),
        ("CCC(O)CC", 0.000),
        ("CC(O)(C)CC", 0.000),
        ("c1ccccc1O", 0.242),
    ] {
        feq(chi_nv(&mol(smi), 5).unwrap(), res, 1e-3, smi);
    }
    let mut chi0n_data = chi0v_data;
    chi0n_data[9].1 = 2.085;
    chi0n_data[10].1 = 2.085;
    chi0n_data[11].1 = 2.085;
    check(chi0n, &chi0n_data, "Chi0n");
    assert_eq!(
        chi_nn(&mol("CCCCCC"), 4).unwrap(),
        chi4n(&mol("CCCCCC")).unwrap()
    );
}

/// UnitTestGraphDescriptors' short regression tests, against
/// PP_descrs_regress.csv: Chi0 (column 2), HallKierAlpha (3), and Kappa3
/// (5), to 1e-4.
#[test]
fn regression_file() {
    let text = include_str!("data/PP_descrs_regress.csv");
    let mut rows = 0;
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split(',').collect();
        let m = mol(fields[0]);
        let col = |i: usize| fields[i].trim().parse::<f64>().unwrap();
        feq(chi0(&m), col(2), 1e-4, fields[0]);
        feq(hall_kier_alpha(&m).unwrap(), col(3), 1e-4, fields[0]);
        feq(kappa3(&m).unwrap(), col(5), 1e-4, fields[0]);
        rows += 1;
    }
    assert_eq!(rows, 98);
}

#[test]
fn values_match_python() {
    // Python's Ipc, AvgIpc, BalabanJ, BertzCT, Chi0, and Chi1, printed with
    // repr; None where Python raises.
    let cases: [(&str, [Option<f64>; 6]); 9] = [
        (
            "CCCC",
            [
                Some(6.854752972273345),
                Some(1.370950594454669),
                Some(1.9747448713915894),
                Some(2.0),
                Some(3.414213562373095),
                Some(1.914213562373095),
            ],
        ),
        (
            "c1ccccc1C(=O)O",
            [
                Some(129.22533870330605),
                Some(1.9003726279897948),
                Some(2.98145461404113),
                Some(203.415952604261),
                Some(6.69023444431199),
                Some(4.304530452640311),
            ],
        ),
        (
            "CC.CCO",
            [
                Some(8.75488750216347),
                Some(1.459147917027245),
                Some(0.0),
                Some(3.6096404744368114),
                Some(4.707106781186548),
                Some(2.414213562373095),
            ],
        ),
        ("C", [Some(0.0); 6]),
        (
            "",
            [Some(0.0), Some(0.0), None, Some(0.0), Some(0.0), Some(0.0)],
        ),
        (
            "[Na+].[Cl-]",
            [
                Some(0.0),
                Some(0.0),
                Some(0.0),
                Some(2.0),
                Some(0.0),
                Some(0.0),
            ],
        ),
        (
            "C1CC1C#N",
            [
                Some(29.69332843436964),
                Some(2.1209520310264027),
                Some(2.2993499215373854),
                Some(66.20736108013043),
                Some(3.6986706127492686),
                Some(2.431851652578137),
            ],
        ),
        (
            "O=C(OC)C1=C(C)NC(C)=C(C(OC)=O)C1C2=CC=CC=C2[N+]([O-])=O",
            [
                Some(291739.89613473933),
                Some(2.5932665724561046),
                Some(2.803490087611988),
                Some(771.1110283918142),
                Some(18.723250160202095),
                Some(11.82780492029403),
            ],
        ),
        (
            "Fc1cccc(C2(c3nnc(Cc4cccc5ccccc45)o3)CCOCC2)c1",
            [
                Some(9626644.371953605),
                Some(3.236081803570815),
                Some(1.5167597530215096),
                Some(1143.0567597744753),
                Some(19.68358750805833),
                Some(14.24395005272655),
            ],
        ),
    ];
    for (smi, expected) in cases {
        let m = mol(smi);
        let actual = [
            ipc(&m).ok(),
            avg_ipc(&m).ok(),
            balaban_j(&m).ok(),
            bertz_ct(&m).ok(),
            Some(chi0(&m)),
            Some(chi1(&m)),
        ];
        for (a, e) in actual.iter().zip(expected) {
            match (a, e) {
                (Some(a), Some(e)) => assert_close(*a, e, smi),
                _ => assert_eq!(a.is_some(), e.is_some(), "{smi}"),
            }
        }
    }
}

#[test]
fn matrix_arguments_match_python() {
    let m = mol("c1ccccc1C(=O)O");
    let ct = |params: BertzCTParams| bertz_ct_with_params(&m, &params).unwrap();
    assert_close(
        ct(BertzCTParams {
            cutoff: 3,
            ..BertzCTParams::default()
        }),
        173.9527932179028,
        "cutoff 3",
    );
    assert_close(
        ct(BertzCTParams {
            cutoff: -2,
            ..BertzCTParams::default()
        }),
        203.415952604261,
        "cutoff -2",
    );
    assert_close(
        ct(BertzCTParams {
            force_d_mat: false,
            ..BertzCTParams::default()
        }),
        200.66106510209752,
        "forceDMat=0",
    );
    let bo = get_distance_matrix_with_params(
        &m,
        &GetDistanceMatrixParams {
            use_bo: true,
            use_atom_wts: false,
        },
    )
    .unwrap();
    assert_close(
        ct(BertzCTParams {
            d_mat: Some(bo),
            force_d_mat: false,
            ..BertzCTParams::default()
        }),
        203.415952604261,
        "dMat",
    );
    let d = get_distance_matrix(&m).unwrap();
    let ipc_d = ipc_with_params(
        &m,
        &IpcParams {
            d_mat: Some(d.clone()),
            ..IpcParams::default()
        },
    )
    .unwrap();
    assert_close(ipc_d, 129.22533870330605, "Ipc dMat");
    let j = balaban_j_with_params(
        &m,
        &BalabanJParams {
            d_mat: Some(d),
            ..BalabanJParams::default()
        },
    )
    .unwrap();
    assert_close(j, 2.2283621149350434, "BalabanJ dMat");
    let bad = BalabanJParams {
        d_mat: Some(vec![vec![0.0]]),
        ..BalabanJParams::default()
    };
    assert!(balaban_j_with_params(&m, &bad).is_err());
}

#[test]
fn hall_kier_alpha_table() {
    let alphas = hall_kier_alphas();
    assert_eq!(alphas.len(), 10);
    assert_eq!(alphas["C"], [Some(-0.22), Some(-0.13), Some(0.0)]);
    assert_eq!(alphas["O"], [None, Some(-0.2), Some(-0.04)]);
}

#[test]
fn defaults_match_python() {
    assert_eq!(
        BertzCTParams::default(),
        BertzCTParams {
            cutoff: 100,
            d_mat: None,
            force_d_mat: true
        }
    );
    assert_eq!(
        IpcParams::default(),
        IpcParams {
            avg: false,
            d_mat: None,
            force_d_mat: false
        }
    );
    assert_eq!(
        BalabanJParams::default(),
        BalabanJParams {
            d_mat: None,
            force_d_mat: false
        }
    );
    assert_eq!(
        AvgIpcParams::default(),
        AvgIpcParams {
            d_mat: None,
            force_d_mat: false
        }
    );
}
