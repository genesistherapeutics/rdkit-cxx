// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's 3D descriptors. Expected values are from Python's
//! rdMolDescriptors, for paracetamol embedded with ETKDG (randomSeed=42).

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_descriptors::*;
use rdkit_cxx::chem::{
    Conformer, Mol, MolFromMolBlockParams, mol_from_mol_block_with_params, mol_from_smiles,
};
use rdkit_cxx::geometry::Point3D;

const APAP_3D: &str = r"
     RDKit          3D

 20 20  0  0  0  0  0  0  0  0999 V2000
    3.7421   -0.2862   -0.1068 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3052   -0.1363   -0.5344 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9839   -0.0125   -1.7458 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3305   -0.1397    0.4734 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0678   -0.0098    0.2545 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6530    1.2425    0.2442 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0225    1.3686    0.0298 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7857    0.2117   -0.1731 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1511    0.3852   -0.3836 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1870   -1.0401   -0.1605 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8111   -1.1585    0.0555 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8086   -1.3020    0.3658 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.4421   -0.1647   -0.9572 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.9396    0.5018    0.6203 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6520   -0.2465    1.4745 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0106    2.1003    0.4071 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5035    2.3374    0.0172 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8160    0.4153    0.3833 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8125   -1.9168   -0.3219 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3830   -2.1498    0.0579 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  1  0
  8 10  2  0
 10 11  1  0
 11  5  2  0
  1 12  1  0
  1 13  1  0
  1 14  1  0
  4 15  1  0
  6 16  1  0
  7 17  1  0
  9 18  1  0
 10 19  1  0
 11 20  1  0
M  END
";

/// Paracetamol with its hydrogens and two conformers: the embedded one
/// (id 0) and a copy stretched twofold along x (id 1).
fn apap() -> Mol {
    let params = MolFromMolBlockParams {
        remove_hs: false,
        ..Default::default()
    };
    let mut mol = mol_from_mol_block_with_params(APAP_3D, &params).unwrap();
    let positions = mol.get_conformer().unwrap().get_positions();
    let mut stretched = Conformer::new(mol.get_num_atoms());
    for (i, p) in positions.iter().enumerate() {
        stretched
            .set_atom_position(
                i as u32,
                Point3D {
                    x: 2.0 * p.x,
                    y: p.y,
                    z: p.z,
                },
            )
            .unwrap();
    }
    stretched.set_3d(true);
    assert_eq!(
        mol.add_conformer_with_assign_id(&stretched, true).unwrap(),
        1
    );
    mol
}

/// The molecule with each atom's `f64` property q set to 0.1 * its index.
fn apap_with_q() -> Mol {
    let mut mol = apap();
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)
            .unwrap()
            .set_double_prop("q", 0.1 * f64::from(idx));
    }
    mol
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn assert_all_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert_close(*a, *e, 1e-9);
    }
}

#[test]
fn moments() {
    let mol = apap();
    type Moment = fn(&Mol) -> rdkit_cxx::Result<f64>;
    let defaults: [(Moment, f64); 9] = [
        (calc_pmi1, 146.23636151183038),
        (calc_pmi2, 949.7245895519983),
        (calc_pmi3, 984.2463819125318),
        (calc_npr1, 0.14857698661555876),
        (calc_npr2, 0.9649256598805548),
        (calc_radius_of_gyration, 2.6230869660176164),
        (calc_inertial_shape_factor, 0.0065983976208440695),
        (calc_eccentricity, 0.9889008438909536),
        (calc_asphericity, 0.6235100671716063),
    ];
    for (f, expected) in defaults {
        assert_close(f(&mol).unwrap(), expected, 1e-7);
    }
    let conf1 = CalcPMI1Params {
        conf_id: 1,
        ..Default::default()
    };
    assert_close(
        calc_pmi1_with_params(&mol, &conf1).unwrap(),
        146.71352597029002,
        1e-7,
    );
    let unweighted = CalcPMI1Params {
        use_atomic_masses: false,
        ..Default::default()
    };
    assert_close(
        calc_pmi1_with_params(&mol, &unweighted).unwrap(),
        33.65336672555779,
        1e-7,
    );
    let conf1 = CalcPMI3Params {
        conf_id: 1,
        ..Default::default()
    };
    assert_close(
        calc_pmi3_with_params(&mol, &conf1).unwrap(),
        3645.5508099716258,
        1e-7,
    );
    let p = CalcNPR1Params {
        conf_id: 1,
        use_atomic_masses: true,
    };
    assert_close(
        calc_npr1_with_params(&mol, &p).unwrap(),
        0.040244542901167774,
        1e-7,
    );
    let p = CalcRadiusOfGyrationParams {
        use_atomic_masses: false,
        ..Default::default()
    };
    assert_close(
        calc_radius_of_gyration_with_params(&mol, &p).unwrap(),
        3.036808849643816,
        1e-7,
    );
    let p = CalcAsphericityParams {
        conf_id: 1,
        ..Default::default()
    };
    assert_close(
        calc_asphericity_with_params(&mol, &p).unwrap(),
        0.8846921596545829,
        1e-7,
    );
    let p = CalcNPR2Params::default();
    assert_close(
        calc_npr2_with_params(&mol, &p).unwrap(),
        0.9649256598805548,
        1e-7,
    );
    let p = CalcInertialShapeFactorParams {
        use_atomic_masses: false,
        ..Default::default()
    };
    assert_close(
        calc_inertial_shape_factor_with_params(&mol, &p).unwrap(),
        0.026680263839480477,
        1e-7,
    );
    let p = CalcEccentricityParams {
        conf_id: 1,
        ..Default::default()
    };
    assert_close(
        calc_eccentricity_with_params(&mol, &p).unwrap(),
        0.9991898602201065,
        1e-7,
    );
    let p = CalcPMI2Params {
        conf_id: 1,
        ..Default::default()
    };
    assert_close(
        calc_pmi2_with_params(&mol, &p).unwrap(),
        3610.708708805409,
        1e-7,
    );
    assert_close(
        calc_spherocity_index(&mol).unwrap(),
        0.12699269845318972,
        1e-7,
    );
    assert_close(
        calc_spherocity_index_with_conf_id(&mol, 1).unwrap(),
        0.03692236795845263,
        1e-7,
    );
    // The caches are gone after each call (U8).
    assert!(!mol.has_prop("_PMI1_mass"));
    assert!(!mol.has_prop("_principalAxes_cov"));
}

#[test]
fn moment_params_defaults() {
    let expected = CalcPMI1Params {
        conf_id: -1,
        use_atomic_masses: true,
    };
    assert_eq!(CalcPMI1Params::default(), expected);
    assert_eq!(CalcAsphericityParams::default().conf_id, -1);
    assert!(CalcEccentricityParams::default().use_atomic_masses);
    assert_eq!(
        CalcWHIMParams::default(),
        CalcWHIMParams {
            conf_id: -1,
            thresh: 0.001,
            custom_atom_property: String::new(),
        }
    );
    assert_eq!(CalcGETAWAYParams::default().precision, 2.0);
    assert_eq!(CalcRDFParams::default().conf_id, -1);
    assert!(CalcMORSEParams::default().custom_atom_property.is_empty());
    assert_eq!(CalcAUTOCORR3DParams::default().conf_id, -1);
}

#[test]
fn missing_conformers() {
    let flat = mol_from_smiles("CCO").unwrap();
    assert!(matches!(calc_pmi1(&flat), Err(Error::Runtime { .. })));
    assert!(matches!(calc_whim(&flat), Err(Error::Runtime { .. })));
    assert!(matches!(
        calc_eem_charges(&flat),
        Err(Error::Runtime { .. })
    ));
    let mol = apap();
    let p = CalcPMI1Params {
        conf_id: 7,
        ..Default::default()
    };
    assert!(matches!(
        calc_pmi1_with_params(&mol, &p),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        calc_pbf_with_conf_id(&mol, 7),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        calc_coulomb_mat_with_conf_id(&mol, 7),
        Err(Error::Value { .. })
    ));
}

#[test]
fn plane_of_best_fit() {
    // Python's CalcPBF(mol, confId=1) after CalcPBF(mol) reads conformer 0's
    // cached plane (0.4411...); each call here uses its own conformer.
    let mol = apap();
    assert_close(calc_pbf(&mol).unwrap(), 0.4396263888229508, 1e-7);
    assert_close(
        calc_pbf_with_conf_id(&mol, 1).unwrap(),
        0.4395703959190972,
        1e-7,
    );
}

#[test]
fn autocorrelations() {
    let mol = apap();
    let values = calc_autocorr_3d(&mol).unwrap();
    assert_eq!(values.len(), 80);
    assert_all_close(&values[..4], &[0.132, 0.364, 0.604, 0.706]);
    assert_all_close(&values[40..42], &[0.085, 0.198]);
    let custom = CalcAUTOCORR3DParams {
        custom_atom_property: "q".to_string(),
        ..Default::default()
    };
    assert_all_close(
        &calc_autocorr_3d_with_params(&apap_with_q(), &custom).unwrap(),
        &[
            0.066, 0.224, 0.445, 0.735, 0.657, 0.574, 0.647, 0.139, 0.289, 0.0,
        ],
    );

    let values = calc_rdf(&mol).unwrap();
    assert_eq!(values.len(), 210);
    assert_all_close(&values[..4], &[4.668, 4.056, 2.616, 11.375]);
    assert_all_close(&values[60..62], &[1.111, 3.742]);
    let custom = CalcRDFParams {
        custom_atom_property: "q".to_string(),
        ..Default::default()
    };
    let values = calc_rdf_with_params(&apap_with_q(), &custom).unwrap();
    assert_eq!(values.len(), 30);
    assert_all_close(&values[..5], &[4.234, 1.148, 1.848, 9.881, 9.044]);

    let values = calc_morse(&mol).unwrap();
    assert_eq!(values.len(), 224);
    assert_all_close(&values[..4], &[190.0, 18.411, -1.025, 0.842]);
    assert_all_close(&values[100..102], &[-2.439, -1.811]);
    let custom = CalcMORSEParams {
        custom_atom_property: "q".to_string(),
        ..Default::default()
    };
    let values = calc_morse_with_params(&apap_with_q(), &custom).unwrap();
    assert_eq!(values.len(), 32);
    assert_all_close(&values[..5], &[168.15, 9.396, -1.57, 2.075, -0.818]);
    assert!(!mol.has_prop("_3DDistanceMatrix_Conf0"));
    assert!(!mol.has_prop("DistanceMatrix"));
}

#[test]
fn whim_and_getaway() {
    let mol = apap();
    let values = calc_whim(&mol).unwrap();
    assert_eq!(values.len(), 114);
    assert_all_close(&values[..4], &[7.54, 1.292, 0.39, 0.818]);
    assert_all_close(&values[100..103], &[0.362, 0.215, 0.218]);
    let loose = CalcWHIMParams {
        thresh: 0.01,
        ..Default::default()
    };
    assert_all_close(
        &calc_whim_with_params(&mol, &loose).unwrap()[5..8],
        &[0.188, 0.23, 0.188],
    );
    let custom = CalcWHIMParams {
        custom_atom_property: "q".to_string(),
        ..Default::default()
    };
    assert_all_close(
        &calc_whim_with_params(&apap_with_q(), &custom).unwrap(),
        &[
            8.249, 1.974, 0.326, 0.782, 0.187, 0.188, 0.188, 0.188, 0.634, 0.829, 0.134, 10.549,
            19.623, 0.188, 0.673, 0.532, 35.488,
        ],
    );

    let values = calc_getaway(&mol).unwrap();
    assert_eq!(values.len(), 273);
    assert_all_close(&values[..4], &[34.054, 0.895, 4.019, 11.411]);
    assert_all_close(&values[200..203], &[0.096, 1.697, 1.452]);
    let precise = CalcGETAWAYParams {
        precision: 3.0,
        ..Default::default()
    };
    assert_all_close(
        &calc_getaway_with_params(&mol, &precise).unwrap()[..4],
        &[38.054, 1.0, 4.019, 11.411],
    );
    let custom = CalcGETAWAYParams {
        custom_atom_property: "q".to_string(),
        ..Default::default()
    };
    let values = calc_getaway_with_params(&apap_with_q(), &custom).unwrap();
    assert_eq!(values.len(), 45);
    assert_all_close(&values[..6], &[34.054, 0.895, 4.019, 11.411, 4.659, 1.022]);
}

#[test]
fn getaway_precision_is_checked() {
    // RDKit converts the precision to unsigned int (undefined for these)
    // and prints that many digits into a 32-byte buffer.
    let mol = apap();
    for precision in [-1.0, 25.0, 1e20, f64::NAN, f64::INFINITY] {
        let p = CalcGETAWAYParams {
            precision,
            ..Default::default()
        };
        assert!(
            matches!(calc_getaway_with_params(&mol, &p), Err(Error::Value { .. })),
            "{precision}"
        );
    }
    let p = CalcGETAWAYParams {
        precision: 24.9,
        ..Default::default()
    };
    assert_eq!(calc_getaway_with_params(&mol, &p).unwrap().len(), 273);
}

#[test]
fn getaway_coordinates_are_checked() {
    // NaN coordinates make NaN leverages, which RDKit sorts.
    let mut mol = apap();
    mol.get_conformer_mut()
        .unwrap()
        .set_atom_position(
            0,
            Point3D {
                x: f64::NAN,
                y: 0.0,
                z: 0.0,
            },
        )
        .unwrap();
    assert!(matches!(calc_getaway(&mol), Err(Error::Value { .. })));
    let p = CalcGETAWAYParams {
        conf_id: 1,
        ..Default::default()
    };
    assert_eq!(calc_getaway_with_params(&mol, &p).unwrap().len(), 273);
}

/// Eigen's SVD and eigensolvers leave their output uninitialized for NaN
/// input, and RDKit reads it, so the descriptors that solve matrices of the
/// coordinates need them finite, and those that divide by distances need
/// the atoms apart.
#[test]
fn solved_descriptors_check_coordinates() {
    let nan = Point3D {
        x: f64::NAN,
        y: 0.0,
        z: 0.0,
    };
    let mut mol = apap();
    mol.get_conformer_mut()
        .unwrap()
        .set_atom_position(0, nan)
        .unwrap();
    assert!(matches!(calc_whim(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_pbf(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_pmi1(&mol), Err(Error::Value { .. })));
    assert!(matches!(
        calc_spherocity_index(&mol),
        Err(Error::Value { .. })
    ));
    assert!(matches!(calc_eem_charges(&mol), Err(Error::Value { .. })));
    // Two atoms at one position, as embedding [Na+].[Cl-] leaves them.
    let mut mol = apap();
    let p = mol.get_conformer().unwrap().get_atom_position(1).unwrap();
    mol.get_conformer_mut()
        .unwrap()
        .set_atom_position(0, p)
        .unwrap();
    for result in [calc_getaway(&mol), calc_eem_charges(&mol)] {
        assert!(
            matches!(&result, Err(Error::Value { message }) if message.contains("atoms 0 and 1")),
            "{result:?}"
        );
    }
    assert_eq!(calc_whim(&mol).unwrap().len(), 114);
    // WHIM weights the covariance matrix it solves by a custom property.
    let mut mol = apap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_double_prop("w", f64::INFINITY);
    let p = CalcWHIMParams {
        custom_atom_property: "w".to_string(),
        ..Default::default()
    };
    assert!(matches!(
        calc_whim_with_params(&mol, &p),
        Err(Error::Value { .. })
    ));
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_double_prop("w", 2.0);
    assert_eq!(calc_whim_with_params(&mol, &p).unwrap().len(), 17);
}

#[test]
fn no_atoms() {
    // AUTOCORR3D, WHIM, and GETAWAY take &v[0] of empty vectors (GETAWAY
    // crashes in Python).
    let mut mol = Mol::default();
    mol.add_conformer(&Conformer::new(0)).unwrap();
    assert!(matches!(calc_getaway(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_whim(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_autocorr_3d(&mol), Err(Error::Value { .. })));
    assert_eq!(calc_rdf(&mol).unwrap().len(), 210);
    assert_eq!(calc_morse(&mol).unwrap().len(), 224);
    assert!(calc_coulomb_mat(&mol).unwrap().is_empty());
    assert_eq!(calc_pbf(&mol).unwrap(), 0.0);
}

#[test]
fn property_tables() {
    // The tables cover atomic numbers 1 to 110.
    let mut mol = apap();
    mol.get_atom_with_idx_mut(0).unwrap().set_atomic_num(0);
    assert!(matches!(calc_rdf(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_morse(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_whim(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_getaway(&mol), Err(Error::Value { .. })));
    assert!(matches!(calc_autocorr_3d(&mol), Err(Error::Value { .. })));
    let custom = CalcRDFParams {
        custom_atom_property: "q".to_string(),
        ..Default::default()
    };
    assert_eq!(calc_rdf_with_params(&mol, &custom).unwrap().len(), 30);
}

#[test]
fn coulomb_matrix() {
    let mol = apap();
    let matrix = calc_coulomb_mat(&mol).unwrap();
    assert_eq!(matrix.len(), 20);
    assert!(matrix.iter().all(|row| row.len() == 20));
    assert_all_close(
        &matrix[0][..3],
        &[36.85810519942594, 23.894074058763216, 19.84131850783672],
    );
    assert_close(matrix[2][1], 38.113847145857925, 1e-9);
    assert_close(
        calc_coulomb_mat_with_conf_id(&mol, 1).unwrap()[0][1],
        12.37410204369376,
        1e-9,
    );
}

#[test]
fn eem_charges() {
    let mol = apap();
    let charges = calc_eem_charges(&mol).unwrap();
    assert_eq!(charges.len(), 20);
    assert_all_close(
        &charges[..4],
        &[
            -0.571887634198,
            0.430354242696,
            -0.584150376929,
            -0.589517105879,
        ],
    );
    assert_close(
        calc_eem_charges_with_conf_id(&mol, 1).unwrap()[0],
        -0.3430800128730414,
        1e-9,
    );
    assert!(!mol.get_atom_with_idx(0).unwrap().has_prop("_EEMCharge"));
}

#[test]
fn eem_tables_are_checked() {
    // EEM reads past its 36-entry tables for iodine (Python gives NaN) and
    // uninitialized memory for a quadruple bond.
    let mut iodo = apap();
    iodo.get_atom_with_idx_mut(0).unwrap().set_atomic_num(53);
    assert!(matches!(calc_eem_charges(&iodo), Err(Error::Value { .. })));
    let mut quadruple = mol_from_smiles("[C]$[C]").unwrap();
    let mut conf = Conformer::new(2);
    conf.set_atom_position(
        1,
        Point3D {
            x: 1.2,
            y: 0.0,
            z: 0.0,
        },
    )
    .unwrap();
    quadruple.add_conformer(&conf).unwrap();
    assert!(matches!(
        calc_eem_charges(&quadruple),
        Err(Error::Value { .. })
    ));
    // Bromine (35) is the last element the tables cover.
    let mut bromo = apap();
    bromo.get_atom_with_idx_mut(0).unwrap().set_atomic_num(35);
    assert_eq!(calc_eem_charges(&bromo).unwrap().len(), 20);
}

const USR: [f64; 12] = [
    2.7999016659236347,
    1.175907585838725,
    -0.47153046570251544,
    2.7938731637763445,
    1.2189913995719395,
    -0.6154407115370574,
    5.07634266039561,
    2.6393739732648003,
    -0.3054338256853875,
    4.8169357564172355,
    2.587446019349795,
    -0.4586430361058503,
];

#[test]
fn usr() {
    let mol = apap();
    let usr = get_usr(&mol).unwrap();
    assert_all_close(&usr, &USR);
    assert_ne!(get_usr_with_conf_id(&mol, 1).unwrap(), usr);
    assert!(matches!(
        get_usr(&mol_from_smiles("CCO").unwrap()),
        Err(Error::Value { .. })
    ));
    let mut two = mol_from_smiles("CO").unwrap();
    two.add_conformer(&Conformer::new(2)).unwrap();
    assert!(matches!(get_usr(&two), Err(Error::Value { .. })));

    let cat = get_usrcat(&mol).unwrap();
    assert_eq!(cat.len(), 60);
    assert_all_close(
        &cat[12..16],
        &[
            2.2927727604779036,
            0.8602982266839411,
            0.826826276799466,
            2.2793797955145485,
        ],
    );
    assert_all_close(&cat[48..50], &[2.802823155653427, 1.3837331566664002]);
    // A skewness of zero is the cube root of a third moment that cancels to
    // rounding error, which the cube root magnifies (3e-6 on aarch64).
    assert_close(cat[50], 0.0, 1e-5);
    assert_close(cat[51], 2.786484859935878, 1e-9);
    let selected = GetUSRCATParams {
        atom_selections: Some(vec![vec![1, 2, 3], vec![5, 6, 7, 8]]),
        ..Default::default()
    };
    let cat = get_usrcat_with_params(&mol, &selected).unwrap();
    assert_eq!(cat.len(), 36);
    assert_all_close(
        &cat[12..15],
        &[2.9224927289376783, 0.5987658959877975, 0.8429171414834553],
    );
    assert_all_close(
        &cat[24..27],
        &[1.7324056476401157, 0.9866472965298341, -0.7509311965521294],
    );
    assert_eq!(GetUSRCATParams::default().atom_selections, None);
    assert_eq!(GetUSRCATParams::default().conf_id, -1);
}

#[test]
fn usrcat_atom_numbers_are_checked() {
    // RDKit reads the coordinates at each selected atom unchecked. Numbers
    // start at 1; 0 wraps to u32::MAX, as in Python.
    let mol = apap();
    let selection = |numbers: Vec<u32>| GetUSRCATParams {
        atom_selections: Some(vec![numbers]),
        ..Default::default()
    };
    assert!(matches!(
        get_usrcat_with_params(&mol, &selection(vec![21])),
        Err(Error::Index { index: 20, len: 20 })
    ));
    assert!(matches!(
        get_usrcat_with_params(&mol, &selection(vec![0])),
        Err(Error::Index {
            index: 4_294_967_295,
            len: 20
        })
    ));
    assert_eq!(
        get_usrcat_with_params(&mol, &selection(vec![20]))
            .unwrap()
            .len(),
        24
    );
    let empty = GetUSRCATParams {
        atom_selections: Some(Vec::new()),
        ..Default::default()
    };
    assert!(matches!(
        get_usrcat_with_params(&mol, &empty),
        Err(Error::Value { .. })
    ));
}

#[test]
fn usr_distributions() {
    let mol = apap();
    let coords = mol.get_conformer().unwrap().get_positions();
    let (distances, outputs) = get_usr_distributions_with_outputs(&coords).unwrap();
    assert_eq!(distances.len(), 4);
    assert!(distances.iter().all(|d| d.len() == 20));
    assert_all_close(
        &distances[0][..3],
        &[3.754537724171246, 2.370247099191348, 2.6426933634504404],
    );
    assert_all_close(&distances[3][..2], &[1.1081256291594377, 2.178510502614114]);
    assert_eq!(outputs.points.len(), 4);
    assert_eq!(
        outputs.points[1],
        Point3D {
            x: -0.0678,
            y: -0.0098,
            z: 0.2545
        }
    );
    assert_eq!(
        outputs.points[2],
        Point3D {
            x: -4.816,
            y: 0.4153,
            z: 0.3833
        }
    );
    assert_eq!(get_usr_distributions(&coords).unwrap(), distances);
    assert!(matches!(
        get_usr_distributions(&[]),
        Err(Error::Value { .. })
    ));

    let origin = Point3D {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    let ones = Point3D {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };
    let from_points = get_usr_distributions_from_points(&coords, &[origin, ones]).unwrap();
    assert_eq!(from_points.len(), 2);
    assert_all_close(
        &from_points[1][..3],
        &[3.2246595308652353, 2.3128138900482247, 3.0874962510098696],
    );
    assert!(matches!(
        get_usr_distributions_from_points(&coords, &[]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        get_usr_distributions_from_points(&[], &[origin]),
        Err(Error::Value { .. })
    ));

    assert_all_close(&get_usr_from_distributions(&distances).unwrap(), &USR);
    assert!(matches!(
        get_usr_from_distributions(&[]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        get_usr_from_distributions(&[vec![1.0], vec![]]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        get_usr_from_distributions(&distances[..3]),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn usr_score() {
    let mut other = USR;
    other[0] += 1.2;
    assert_close(
        get_usr_score(&USR, &other).unwrap(),
        0.9090909090909091,
        1e-12,
    );
    let cat = get_usrcat(&apap()).unwrap();
    assert_eq!(get_usr_score(&cat, &cat).unwrap(), 1.0);
    let shifted: Vec<f64> = cat.iter().map(|c| c + 0.5).collect();
    assert_close(
        get_usr_score_with_weights(&cat, &shifted, &[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap(),
        0.11764705882352941,
        1e-12,
    );
    assert!(matches!(
        get_usr_score(&cat, &USR),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        get_usr_score_with_weights(&cat, &cat, &[1.0]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn double_cubic_lattice_volume() {
    let mol = apap();
    let assert_values = |dclv: DoubleCubicLatticeVolume, expected: [f64; 5]| {
        let actual = [
            dclv.get_surface_area(),
            dclv.get_volume(),
            dclv.get_vdw_volume(),
            dclv.get_compactness(),
            dclv.get_packing_density(),
        ];
        for (a, e) in actual.into_iter().zip(expected) {
            assert_close(a, e, 1e-9);
        }
    };
    assert_values(
        DoubleCubicLatticeVolume::new(&mol).unwrap(),
        [
            298.4404626262069,
            411.9602195996842,
            137.902015045098,
            1.1146453948109476,
            0.33474594993444295,
        ],
    );
    let probe = DoubleCubicLatticeVolumeParams {
        probe_radius: 1.4,
        ..Default::default()
    };
    assert_values(
        DoubleCubicLatticeVolume::new_with_params(&mol, &probe).unwrap(),
        [
            324.85594105067,
            474.59129705476926,
            137.902015045098,
            1.1040621188945368,
            0.29057004606046055,
        ],
    );
    let shallow = DoubleCubicLatticeVolumeParams {
        depth: 2,
        ..Default::default()
    };
    assert_values(
        DoubleCubicLatticeVolume::new_with_params(&mol, &shallow).unwrap(),
        [
            296.9432863397013,
            410.48591111676745,
            137.56011729647295,
            1.1117075345034586,
            0.3351153196031188,
        ],
    );
    // A protein's surface area is its residues', and this has none.
    let protein = DoubleCubicLatticeVolumeParams {
        is_protein: true,
        ..Default::default()
    };
    assert_values(
        DoubleCubicLatticeVolume::new_with_params(&mol, &protein).unwrap(),
        [
            0.0,
            503.6006013573,
            147.5474438172982,
            0.0,
            0.29298504294798217,
        ],
    );
    assert_eq!(
        DoubleCubicLatticeVolumeParams::default(),
        DoubleCubicLatticeVolumeParams {
            is_protein: false,
            include_ligand: true,
            probe_radius: 1.2,
            depth: 4,
            dot_density: 0,
        }
    );
}

#[test]
fn double_cubic_lattice_volume_checks() {
    let mol = apap();
    let value = |result: rdkit_cxx::Result<DoubleCubicLatticeVolume>| {
        matches!(result, Err(Error::Value { .. }))
    };
    // RDKit casts scaled probe radii and coordinates to grid indices.
    for bad in [f64::NAN, f64::INFINITY, 1e7] {
        let probe = DoubleCubicLatticeVolumeParams {
            probe_radius: bad,
            ..Default::default()
        };
        assert!(value(DoubleCubicLatticeVolume::new_with_params(
            &mol, &probe
        )));
    }
    let mut far = apap();
    far.get_conformer_mut()
        .unwrap()
        .set_atom_position(
            2,
            Point3D {
                x: 1e9,
                y: 0.0,
                z: 0.0,
            },
        )
        .unwrap();
    assert!(value(DoubleCubicLatticeVolume::new(&far)));
    // RDKit 2024.03 never sizes the standard dots for a dot density.
    let dots = DoubleCubicLatticeVolumeParams {
        dot_density: 1,
        ..Default::default()
    };
    assert!(value(DoubleCubicLatticeVolume::new_with_params(
        &mol, &dots
    )));
    // Python's RuntimeError.
    let deep = DoubleCubicLatticeVolumeParams {
        depth: 9,
        ..Default::default()
    };
    assert!(matches!(
        DoubleCubicLatticeVolume::new_with_params(&mol, &deep),
        Err(Error::Runtime { .. })
    ));
    assert!(value(DoubleCubicLatticeVolume::new(
        &mol_from_smiles("CCO").unwrap()
    )));
    // As in Python, a molecule without atoms has no surface or volume.
    let mut empty = Mol::default();
    empty.add_conformer(&Conformer::new(0)).unwrap();
    let dclv = DoubleCubicLatticeVolume::new(&empty).unwrap();
    assert_eq!((dclv.get_surface_area(), dclv.get_volume()), (0.0, 0.0));
}

#[test]
fn double_cubic_lattice_volume_in_threads() {
    // Each computation keeps its own state.
    let threads: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                let mol = apap();
                let dclv = DoubleCubicLatticeVolume::new(&mol).unwrap();
                (dclv.get_surface_area(), dclv.get_volume())
            })
        })
        .collect();
    for thread in threads {
        let (area, volume) = thread.join().unwrap();
        assert_close(area, 298.4404626262069, 1e-9);
        assert_close(volume, 411.9602195996842, 1e-9);
    }
}
