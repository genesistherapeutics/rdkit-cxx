// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.TorsionFingerprints in RDKit 2024.09.6, and
// the tests of UnitTestMol3D.py. The ligands are that test's PDB files with
// their bond orders assigned from its template (AllChem's
// AssignBondOrdersFromTemplate), written as mol files; Python gives the same
// values for both.

use rdkit_cxx::Error;
use rdkit_cxx::chem::torsion_fingerprints::*;
use rdkit_cxx::chem::{Mol, mol_from_mol_file, mol_from_smiles};

fn ligand(name: &str) -> Mol {
    mol_from_mol_file(format!(
        "{}/tests/data/{name}.mol",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

/// Angles and deviations come from trigonometry in the platform's math
/// library, so they are compared to a few ulps.
fn close(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1e-300)
}

fn assert_floats(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!(close(*a, *e), "{actual:?} != {expected:?}");
    }
}

fn assert_torsions(actual: &[Torsion], expected: &[Torsion]) {
    assert_eq!(actual.len(), expected.len());
    for ((aq, ad), (eq, ed)) in actual.iter().zip(expected) {
        assert_eq!(aq, eq);
        assert!(close(*ad, *ed), "{ad} != {ed}");
    }
}

fn assert_angles(actual: &[TorsionAngles], expected: &[TorsionAngles]) {
    assert_eq!(actual.len(), expected.len());
    for ((aa, ad), (ea, ed)) in actual.iter().zip(expected) {
        assert_floats(aa, ea);
        assert_eq!(ad, ed);
    }
}

#[test]
fn ligand_matches_python() {
    let mol = ligand("1DWD_ligand");
    let mol2 = ligand("1PPC_ligand");
    let (tl, tlr) = calculate_torsion_lists(&mol).unwrap();
    assert_eq!((tl.len(), tlr.len()), (11, 4));
    assert_torsions(
        &tl,
        &[
            (vec![[13, 0, 4, 3]], 180.0),
            (vec![[4, 0, 13, 14]], 180.0),
            (vec![[0, 13, 14, 15]], 180.0),
            (vec![[13, 14, 15, 17]], 180.0),
            (vec![[14, 15, 17, 18]], 180.0),
            (vec![[15, 17, 18, 21]], 180.0),
            (vec![[17, 18, 19, 20]], 180.0),
            (vec![[17, 18, 21, 22]], 180.0),
            (vec![[20, 19, 31, 32], [20, 19, 31, 36]], 180.0),
            (vec![[18, 21, 22, 23], [18, 21, 22, 24]], 180.0),
            (vec![[25, 27, 28, 30], [26, 27, 28, 30]], 180.0),
        ],
    );
    assert_torsions(
        &tlr,
        &[
            (
                vec![
                    [3, 11, 12, 6],
                    [11, 12, 6, 5],
                    [12, 6, 5, 4],
                    [6, 5, 4, 3],
                    [5, 4, 3, 11],
                    [4, 3, 11, 12],
                ],
                36.34137323903797,
            ),
            (
                vec![
                    [7, 12, 11, 10],
                    [12, 11, 10, 9],
                    [11, 10, 9, 8],
                    [10, 9, 8, 7],
                    [9, 8, 7, 12],
                    [8, 7, 12, 11],
                ],
                36.34137323903797,
            ),
            (
                vec![
                    [23, 25, 27, 26],
                    [25, 27, 26, 24],
                    [27, 26, 24, 22],
                    [26, 24, 22, 23],
                    [24, 22, 23, 25],
                    [22, 23, 25, 27],
                ],
                36.34137323903797,
            ),
            (
                vec![
                    [32, 33, 34, 35],
                    [33, 34, 35, 36],
                    [34, 35, 36, 31],
                    [35, 36, 31, 32],
                    [36, 31, 32, 33],
                    [31, 32, 33, 34],
                ],
                36.34137323903797,
            ),
        ],
    );
    let spec = CalculateTorsionListsParams {
        max_dev: MaxDev::Spec,
        ..Default::default()
    };
    let (tl_spec, _) = calculate_torsion_lists_with_params(&mol, &spec).unwrap();
    assert_eq!(tl_spec.last().unwrap().1, 90.0);
    assert_torsions(
        &tl_spec,
        &[
            (vec![[13, 0, 4, 3]], 180.0),
            (vec![[4, 0, 13, 14]], 180.0),
            (vec![[0, 13, 14, 15]], 180.0),
            (vec![[13, 14, 15, 17]], 180.0),
            (vec![[14, 15, 17, 18]], 180.0),
            (vec![[15, 17, 18, 21]], 180.0),
            (vec![[17, 18, 19, 20]], 180.0),
            (vec![[17, 18, 21, 22]], 180.0),
            (vec![[20, 19, 31, 32], [20, 19, 31, 36]], 90.0),
            (vec![[18, 21, 22, 23], [18, 21, 22, 24]], 90.0),
            (vec![[25, 27, 28, 30], [26, 27, 28, 30]], 90.0),
        ],
    );
    let r0 = CalculateTorsionListsParams {
        symm_radius: 0,
        ..Default::default()
    };
    let (tl_r0, _) = calculate_torsion_lists_with_params(&mol, &r0).unwrap();
    assert_eq!(tl_r0[0].0.len(), 2);
    assert_torsions(
        &tl_r0,
        &[
            (vec![[13, 0, 4, 3], [13, 0, 4, 5]], 180.0),
            (vec![[4, 0, 13, 14]], 180.0),
            (vec![[0, 13, 14, 15]], 180.0),
            (vec![[13, 14, 15, 17]], 180.0),
            (vec![[16, 15, 17, 18]], 180.0),
            (vec![[15, 17, 18, 21]], 180.0),
            (vec![[17, 18, 19, 20]], 180.0),
            (vec![[17, 18, 21, 22]], 180.0),
            (vec![[20, 19, 31, 32], [20, 19, 31, 36]], 180.0),
            (vec![[18, 21, 22, 23], [18, 21, 22, 24]], 180.0),
            (vec![[25, 27, 28, 30], [26, 27, 28, 30]], 180.0),
        ],
    );

    let weights = calculate_torsion_weights(&mol).unwrap();
    assert_eq!(weights[4], 1.0);
    assert_eq!(weights.len(), tl.len() + tlr.len());
    assert_floats(
        &weights,
        &[
            0.1621349131306841,
            0.3593813663804627,
            0.6345548031926865,
            0.8925186172992409,
            1.0,
            0.8925186172992409,
            0.6345548031926865,
            0.6345548031926865,
            0.3593813663804627,
            0.3593813663804627,
            0.003804056104782509,
            0.05004301611600176,
            0.0020731342544747337,
            0.17480479985664862,
            0.17480479985664862,
        ],
    );
    let params = CalculateTorsionWeightsParams {
        aid1: Some(15),
        aid2: Some(14),
        ..Default::default()
    };
    let weights_15_14 = calculate_torsion_weights_with_params(&mol, &params).unwrap();
    assert_eq!(weights_15_14[3], 1.0);
    assert_floats(
        &weights_15_14,
        &[
            0.3593813663804627,
            0.6345548031926865,
            0.8925186172992409,
            1.0,
            0.8925186172992409,
            0.6345548031926865,
            0.3593813663804627,
            0.3593813663804627,
            0.1621349131306841,
            0.1621349131306841,
            0.000691044751491578,
            0.17480479985664862,
            0.011412168314347527,
            0.05004301611600176,
            0.05004301611600176,
        ],
    );
    for (aid1, aid2) in [(Some(15), Some(3)), (Some(15), None)] {
        let params = CalculateTorsionWeightsParams {
            aid1,
            aid2,
            ..Default::default()
        };
        assert!(matches!(
            calculate_torsion_weights_with_params(&mol, &params),
            Err(Error::Value { .. })
        ));
    }

    let torsions = calculate_torsion_angles(&mol, &tl, &tlr).unwrap();
    assert_eq!(torsions.len(), weights.len());
    assert!((torsions[2].0[0] - 232.5346).abs() < 1e-4);
    assert_angles(
        &torsions,
        &[
            (vec![69.604352958508], 180.0),
            (vec![82.40015795594695], 180.0),
            (vec![232.53464807020788], 180.0),
            (vec![189.19010043401678], 180.0),
            (vec![181.12758376870792], 180.0),
            (vec![239.1918122803744], 180.0),
            (vec![99.5442945702578], 180.0),
            (vec![70.21428692407258], 180.0),
            (vec![179.02322702164676, 357.8171737710665], 180.0),
            (vec![87.78049328843373, 271.075152181605], 180.0),
            (vec![18.032158497437244, 196.20444062588834], 180.0),
            (vec![2.2924952491480357], 36.34137323903797),
            (vec![0.3108357033217675], 36.34137323903797),
            (vec![1.8463860510709473], 36.34137323903797),
            (vec![53.29120673414943], 36.34137323903797),
        ],
    );
    let torsions2 = calculate_torsion_angles(&mol2, &tl, &tlr).unwrap();
    assert_angles(
        &torsions2,
        &[
            (vec![249.63720827416375], 180.0),
            (vec![60.81150721153571], 180.0),
            (vec![248.9338077521336], 180.0),
            (vec![196.21377507169868], 180.0),
            (vec![182.24115150479017], 180.0),
            (vec![251.64771080235596], 180.0),
            (vec![102.85277711172922], 180.0),
            (vec![68.22756801701306], 180.0),
            (vec![170.30431359896414, 0.3252046705816744], 180.0),
            (vec![277.558509585149, 93.19134093741678], 180.0),
            (vec![194.44513617098428, 15.836061852720821], 180.0),
            (vec![1.4625336776928173], 36.34137323903797),
            (vec![1.5075146458903361], 36.34137323903797),
            (vec![1.9232525279250072], 36.34137323903797),
            (vec![60.77005278532204], 36.34137323903797),
        ],
    );
    assert_eq!(calculate_tfd(&torsions, &torsions).unwrap(), 0.0);
    let tfd_weighted = calculate_tfd_with_weights(&torsions, &torsions2, Some(&weights)).unwrap();
    assert!(close(tfd_weighted, 0.06912369904282437));
    assert!(close(
        calculate_tfd(&torsions, &torsions2).unwrap(),
        0.11147525399171747
    ));
    assert!(matches!(
        calculate_tfd(&torsions, &torsions2[1..]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        calculate_tfd_with_weights(&torsions, &torsions2, Some(&weights[1..])),
        Err(Error::Value { .. })
    ));

    assert!(close(
        get_tfd_between_molecules(&mol, &mol2).unwrap(),
        0.06912369904282437
    ));
    assert!(close(
        get_best_tfd_between_molecules(&mol2, &mol).unwrap(),
        0.06912369904282437
    ));
    let params = GetTFDBetweenMoleculesParams {
        use_weights: false,
        max_dev: MaxDev::Spec,
        symm_radius: 0,
        ..Default::default()
    };
    assert!(close(
        get_tfd_between_molecules_with_params(&mol, &mol2, &params).unwrap(),
        0.05032161138750459
    ));
    let benzene = mol_from_smiles("c1ccccc1").unwrap();
    assert!(matches!(
        get_tfd_between_molecules(&mol, &benzene),
        Err(Error::Value { .. })
    ));

    // UnitTestMol3D.testTorsionFingerprintsAtomReordering.
    let reordered = ligand("1DWD_ligand_reordered");
    assert_eq!(get_tfd_between_molecules(&mol, &reordered).unwrap(), 0.0);
    assert!(close(
        get_tfd_between_molecules(&reordered, &mol2).unwrap(),
        0.06912369904282437
    ));

    let mut m = mol.clone();
    let conf = rdkit_cxx::chem::Conformer::from(mol2.get_conformer().unwrap());
    m.add_conformer_with_assign_id(&conf, true).unwrap();
    m.add_conformer_with_assign_id(&conf, true).unwrap();
    let tfd = get_tfd_between_conformers(&m, &[0], &[1, 2]).unwrap();
    assert_floats(&tfd, &[0.06912369904282437, 0.06912369904282437]);
    let params = GetTFDBetweenConformersParams {
        use_weights: false,
        ..Default::default()
    };
    assert_floats(
        &get_tfd_between_conformers_with_params(&m, &[0, 1], &[2], &params).unwrap(),
        &[0.11147525399171747, 0.0],
    );
    assert_floats(
        &get_tfd_matrix(&m).unwrap(),
        &[0.06912369904282437, 0.06912369904282437, 0.0],
    );
    let params = GetTFDMatrixParams {
        max_dev: MaxDev::Spec,
        ..Default::default()
    };
    assert_floats(
        &get_tfd_matrix_with_params(&m, &params).unwrap(),
        &[0.07162550588035656, 0.07162550588035656, 0.0],
    );
    assert_eq!(get_best_tfd_between_molecules(&mol2, &m).unwrap(), 0.0);
}

#[test]
fn torsion_lists_match_python() {
    // Lists and weights, with and without colinear bonds; the first four are
    // UnitTestMol3D's testTorsionFingerprintsColinearBonds, testGithub4720,
    // and testTorsionAngleLeargerThan14.
    {
        let m = mol_from_smiles("CCC#CCC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[]);
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[],
        );
    }
    {
        let m = mol_from_smiles("CCC#CCC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[(vec![[0, 1, 4, 5]], 180.0)]);
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.00010000000000000009],
        );
    }
    {
        let m = mol_from_smiles("C#CCC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[]);
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[],
        );
    }
    {
        let m = mol_from_smiles("C#CCC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[]);
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[],
        );
    }
    {
        let m = mol_from_smiles("S(F)(F)(F)(F)(Cl)c1ccccc1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[(vec![[5, 0, 6, 7], [5, 0, 6, 11]], 180.0)]);
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [6, 11, 10, 9],
                    [11, 10, 9, 8],
                    [10, 9, 8, 7],
                    [9, 8, 7, 6],
                    [8, 7, 6, 11],
                    [7, 6, 11, 10],
                ],
                36.34137323903797,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.3593813663804628, 0.17480479985664862],
        );
    }
    {
        let m = mol_from_smiles("S(F)(F)(F)(F)(Cl)c1ccccc1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[(vec![[5, 0, 6, 7], [5, 0, 6, 11]], 180.0)]);
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [6, 11, 10, 9],
                    [11, 10, 9, 8],
                    [10, 9, 8, 7],
                    [9, 8, 7, 6],
                    [8, 7, 6, 11],
                    [7, 6, 11, 10],
                ],
                36.34137323903797,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.3593813663804628, 0.17480479985664862],
        );
    }
    {
        let m = mol_from_smiles("C1CCCCCCCCCCCCCC1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[]);
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [0, 1, 2, 3],
                    [1, 2, 3, 4],
                    [2, 3, 4, 5],
                    [3, 4, 5, 6],
                    [4, 5, 6, 7],
                    [5, 6, 7, 8],
                    [6, 7, 8, 9],
                    [7, 8, 9, 10],
                    [8, 9, 10, 11],
                    [9, 10, 11, 12],
                    [10, 11, 12, 13],
                    [11, 12, 13, 14],
                    [12, 13, 14, 0],
                    [13, 14, 0, 1],
                    [14, 0, 1, 2],
                ],
                180.0,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.49691917232902644],
        );
    }
    {
        let m = mol_from_smiles("C1CCCCCCCCCCCCCC1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[]);
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [0, 1, 2, 3],
                    [1, 2, 3, 4],
                    [2, 3, 4, 5],
                    [3, 4, 5, 6],
                    [4, 5, 6, 7],
                    [5, 6, 7, 8],
                    [6, 7, 8, 9],
                    [7, 8, 9, 10],
                    [8, 9, 10, 11],
                    [9, 10, 11, 12],
                    [10, 11, 12, 13],
                    [11, 12, 13, 14],
                    [12, 13, 14, 0],
                    [13, 14, 0, 1],
                    [14, 0, 1, 2],
                ],
                180.0,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.49691917232902644],
        );
    }
    {
        let m = mol_from_smiles("CC(C)(C)c1ccc(cc1)C(=O)NCC=C=CC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (
                    vec![
                        [0, 1, 4, 5],
                        [0, 1, 4, 9],
                        [2, 1, 4, 5],
                        [2, 1, 4, 9],
                        [3, 1, 4, 5],
                        [3, 1, 4, 9],
                    ],
                    180.0,
                ),
                (vec![[6, 7, 10, 11], [8, 7, 10, 11]], 180.0),
                (vec![[11, 10, 12, 13]], 180.0),
                (vec![[10, 12, 13, 14]], 180.0),
                (vec![[12, 13, 14, 15]], 180.0),
            ],
        );
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [4, 9, 8, 7],
                    [9, 8, 7, 6],
                    [8, 7, 6, 5],
                    [7, 6, 5, 4],
                    [6, 5, 4, 9],
                    [5, 4, 9, 8],
                ],
                36.34137323903797,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                0.0166810053720006,
                1.0,
                0.7742636826811271,
                0.3593813663804628,
                0.10000000000000006,
                1.0781440991413884,
            ],
        );
    }
    {
        let m = mol_from_smiles("CC(C)(C)c1ccc(cc1)C(=O)NCC=C=CC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (
                    vec![
                        [0, 1, 4, 5],
                        [0, 1, 4, 9],
                        [2, 1, 4, 5],
                        [2, 1, 4, 9],
                        [3, 1, 4, 5],
                        [3, 1, 4, 9],
                    ],
                    180.0,
                ),
                (vec![[6, 7, 10, 11], [8, 7, 10, 11]], 180.0),
                (vec![[11, 10, 12, 13]], 180.0),
                (vec![[10, 12, 13, 14]], 180.0),
                (vec![[12, 13, 14, 15]], 180.0),
                (vec![[13, 14, 16, 17]], 180.0),
            ],
        );
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [4, 9, 8, 7],
                    [9, 8, 7, 6],
                    [8, 7, 6, 5],
                    [7, 6, 5, 4],
                    [6, 5, 4, 9],
                    [5, 4, 9, 8],
                ],
                36.34137323903797,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                0.0166810053720006,
                1.0,
                0.7742636826811271,
                0.3593813663804628,
                0.10000000000000006,
                0.0166810053720006,
                1.0781440991413884,
            ],
        );
    }
    {
        let m = mol_from_smiles("c1ccccc1-c1ccccc1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[(
                vec![[4, 5, 6, 7], [4, 5, 6, 11], [0, 5, 6, 7], [0, 5, 6, 11]],
                180.0,
            )],
        );
        assert_torsions(
            &tlr,
            &[
                (
                    vec![
                        [0, 5, 4, 3],
                        [5, 4, 3, 2],
                        [4, 3, 2, 1],
                        [3, 2, 1, 0],
                        [2, 1, 0, 5],
                        [1, 0, 5, 4],
                    ],
                    36.34137323903797,
                ),
                (
                    vec![
                        [7, 8, 9, 10],
                        [8, 9, 10, 11],
                        [9, 10, 11, 6],
                        [10, 11, 6, 7],
                        [11, 6, 7, 8],
                        [6, 7, 8, 9],
                    ],
                    36.34137323903797,
                ),
            ],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[1.0, 0.30000000000000004, 0.30000000000000004],
        );
    }
    {
        let m = mol_from_smiles("c1ccccc1-c1ccccc1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[(
                vec![[4, 5, 6, 7], [4, 5, 6, 11], [0, 5, 6, 7], [0, 5, 6, 11]],
                180.0,
            )],
        );
        assert_torsions(
            &tlr,
            &[
                (
                    vec![
                        [0, 5, 4, 3],
                        [5, 4, 3, 2],
                        [4, 3, 2, 1],
                        [3, 2, 1, 0],
                        [2, 1, 0, 5],
                        [1, 0, 5, 4],
                    ],
                    36.34137323903797,
                ),
                (
                    vec![
                        [7, 8, 9, 10],
                        [8, 9, 10, 11],
                        [9, 10, 11, 6],
                        [10, 11, 6, 7],
                        [11, 6, 7, 8],
                        [6, 7, 8, 9],
                    ],
                    36.34137323903797,
                ),
            ],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[1.0, 0.30000000000000004, 0.30000000000000004],
        );
    }
    {
        let m = mol_from_smiles("CCCCC(CC)CC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (vec![[0, 1, 2, 3]], 180.0),
                (vec![[1, 2, 3, 4]], 180.0),
                (vec![[2, 3, 4, 5], [2, 3, 4, 7]], 180.0),
                (vec![[3, 4, 5, 6]], 180.0),
                (vec![[3, 4, 7, 8]], 180.0),
            ],
        );
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                0.0166810053720006,
                0.3593813663804628,
                1.0,
                0.3593813663804628,
                0.3593813663804628,
            ],
        );
    }
    {
        let m = mol_from_smiles("CCCCC(CC)CC").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (vec![[0, 1, 2, 3]], 180.0),
                (vec![[1, 2, 3, 4]], 180.0),
                (vec![[2, 3, 4, 5], [2, 3, 4, 7]], 180.0),
                (vec![[3, 4, 5, 6]], 180.0),
                (vec![[3, 4, 7, 8]], 180.0),
            ],
        );
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                0.0166810053720006,
                0.3593813663804628,
                1.0,
                0.3593813663804628,
                0.3593813663804628,
            ],
        );
    }
    {
        let m = mol_from_smiles("C1CC1CC1CC1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (vec![[1, 2, 3, 4], [0, 2, 3, 4]], 180.0),
                (vec![[2, 3, 4, 5], [2, 3, 4, 6]], 180.0),
            ],
        );
        assert_torsions(
            &tlr,
            &[
                (
                    vec![[0, 1, 2, 0], [1, 2, 0, 1], [2, 0, 1, 2]],
                    8.740407828601791,
                ),
                (
                    vec![[5, 4, 6, 5], [4, 6, 5, 4], [6, 5, 4, 6]],
                    8.740407828601791,
                ),
            ],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                1.0,
                0.10000000000000002,
                0.0250215080580009,
                5.390720495706936e-06,
            ],
        );
    }
    {
        let m = mol_from_smiles("C1CC1CC1CC1").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (vec![[1, 2, 3, 4], [0, 2, 3, 4]], 180.0),
                (vec![[2, 3, 4, 5], [2, 3, 4, 6]], 180.0),
            ],
        );
        assert_torsions(
            &tlr,
            &[
                (
                    vec![[0, 1, 2, 0], [1, 2, 0, 1], [2, 0, 1, 2]],
                    8.740407828601791,
                ),
                (
                    vec![[5, 4, 6, 5], [4, 6, 5, 4], [6, 5, 4, 6]],
                    8.740407828601791,
                ),
            ],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                1.0,
                0.10000000000000002,
                0.0250215080580009,
                5.390720495706936e-06,
            ],
        );
    }
    {
        let m = mol_from_smiles("OC(=O)CC#CC#CCN").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(&tl, &[(vec![[2, 1, 3, 4]], 180.0)]);
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.10000000000000002],
        );
    }
    {
        let m = mol_from_smiles("OC(=O)CC#CC#CCN").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[(vec![[2, 1, 3, 4]], 180.0), (vec![[1, 3, 8, 9]], 180.0)],
        );
        assert_torsions(&tlr, &[]);
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[0.10000000000000002, 0.10000000000000002],
        );
    }
    {
        let m = mol_from_smiles("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (vec![[0, 1, 3, 4], [2, 1, 3, 4]], 180.0),
                (vec![[1, 3, 4, 5], [1, 3, 4, 9]], 180.0),
                (vec![[6, 7, 10, 12], [8, 7, 10, 12]], 180.0),
                (vec![[7, 10, 12, 14]], 180.0),
            ],
        );
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [4, 9, 8, 7],
                    [9, 8, 7, 6],
                    [8, 7, 6, 5],
                    [7, 6, 5, 4],
                    [6, 5, 4, 9],
                    [5, 4, 9, 8],
                ],
                36.34137323903797,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: true,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                0.03630780547701015,
                0.22908676527677732,
                0.22908676527677732,
                0.03630780547701015,
                1.078144099141388,
            ],
        );
    }
    {
        let m = mol_from_smiles("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O").unwrap();
        let params = CalculateTorsionListsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        let (tl, tlr) = calculate_torsion_lists_with_params(&m, &params).unwrap();
        assert_torsions(
            &tl,
            &[
                (vec![[0, 1, 3, 4], [2, 1, 3, 4]], 180.0),
                (vec![[1, 3, 4, 5], [1, 3, 4, 9]], 180.0),
                (vec![[6, 7, 10, 12], [8, 7, 10, 12]], 180.0),
                (vec![[7, 10, 12, 14]], 180.0),
            ],
        );
        assert_torsions(
            &tlr,
            &[(
                vec![
                    [4, 9, 8, 7],
                    [9, 8, 7, 6],
                    [8, 7, 6, 5],
                    [7, 6, 5, 4],
                    [6, 5, 4, 9],
                    [5, 4, 9, 8],
                ],
                36.34137323903797,
            )],
        );
        let params = CalculateTorsionWeightsParams {
            ignore_colinear_bonds: false,
            ..Default::default()
        };
        assert_floats(
            &calculate_torsion_weights_with_params(&m, &params).unwrap(),
            &[
                0.03630780547701015,
                0.22908676527677732,
                0.22908676527677732,
                0.03630780547701015,
                1.078144099141388,
            ],
        );
    }
}
