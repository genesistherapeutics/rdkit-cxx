// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{assert_close, ethanol, ethanol_two_confs};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_transforms::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values below come from Python's rdMolTransforms on the same
// molecules.

fn assert_point(actual: Point3D, expected: [f64; 3]) {
    assert_close(actual.x, expected[0]);
    assert_close(actual.y, expected[1]);
    assert_close(actual.z, expected[2]);
}

fn assert_rows(actual: &[Vec<f64>], expected: &[&[f64]]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert_eq!(a.len(), e.len());
        for (x, y) in a.iter().zip(e.iter()) {
            assert!((x - y).abs() < 1e-12, "{actual:?}");
        }
    }
}

#[test]
fn params_defaults() {
    assert_eq!(
        ComputeCentroidParams::default(),
        ComputeCentroidParams {
            ignore_hs: true,
            weights: None
        }
    );
    assert_eq!(
        ComputeCanonicalTransformParams::default(),
        ComputeCanonicalTransformParams {
            center: None,
            normalize_covar: false,
            ignore_hs: true,
        }
    );
    assert_eq!(
        ComputePrincipalAxesAndMomentsParams::default(),
        ComputePrincipalAxesAndMomentsParams {
            ignore_hs: true,
            weights: None
        }
    );
    assert_eq!(
        ComputePrincipalAxesAndMomentsFromGyrationMatrixParams::default(),
        ComputePrincipalAxesAndMomentsFromGyrationMatrixParams {
            ignore_hs: true,
            weights: None
        }
    );
    assert_eq!(
        CanonicalizeConformerParams::default(),
        CanonicalizeConformerParams {
            center: None,
            normalize_covar: false,
            ignore_hs: true,
        }
    );
    assert_eq!(
        CanonicalizeMolParams::default(),
        CanonicalizeMolParams {
            normalize_covar: false,
            ignore_hs: true,
        }
    );
}

#[test]
fn centroids_and_axes() {
    let mol = ethanol();
    let conf = mol.get_conformer().unwrap();
    assert_point(
        compute_centroid(&conf).unwrap(),
        [0.2546, 0.016066666666666674, 0.20553333333333335],
    );
    let all = ComputeCentroidParams {
        ignore_hs: false,
        ..Default::default()
    };
    assert_point(
        compute_centroid_with_params(&conf, &all).unwrap(),
        [-1.1111111111060544e-05, 0.0, -1.1111111111097552e-05],
    );
    let weighted = ComputeCentroidParams {
        weights: Some((1..=9).map(f64::from).collect()),
        ..Default::default()
    };
    assert_point(
        compute_centroid_with_params(&conf, &weighted).unwrap(),
        [0.6192166666666666, 0.07048333333333333, 0.3285666666666667],
    );
    let short = ComputeCentroidParams {
        weights: Some(vec![1.0]),
        ..Default::default()
    };
    assert!(matches!(
        compute_centroid_with_params(&conf, &short),
        Err(Error::Runtime { .. })
    ));

    assert_rows(
        &compute_canonical_transform(&conf).unwrap(),
        &[
            &[
                0.9456520915724026,
                0.12673756086503438,
                0.29946571151088786,
                -0.30434945856476975,
            ],
            &[
                -0.31544733007111075,
                0.5811578273693904,
                0.7501656894568678,
                -0.08320843356333164,
            ],
            &[
                -0.07896267254689594,
                -0.8038614124046954,
                0.5895521401801372,
                -0.0881533467552824,
            ],
            &[0.0, 0.0, 0.0, 1.0],
        ],
    );
    let centered = ComputeCanonicalTransformParams {
        center: Some(Point3D::new(1.0, 0.0, 0.0)),
        normalize_covar: true,
        ignore_hs: false,
    };
    let trans = compute_canonical_transform_with_params(&conf, &centered).unwrap();
    assert_rows(
        &trans[..1],
        &[&[
            0.9966927012968062,
            0.0742106466439149,
            0.033112521896824156,
            -0.9966927012968062,
        ]],
    );

    let (axes, moments) = compute_principal_axes_and_moments(&conf).unwrap().unwrap();
    assert_rows(
        &axes,
        &[
            &[0.9456520915724026, 0.3154473300711106, -0.07896267254689619],
            &[0.12673756086503438, -0.581157827369391, -0.8038614124046948],
            &[0.2994657115108878, -0.750165689456867, 0.5895521401801378],
        ],
    );
    for (m, e) in moments
        .iter()
        .zip([0.49606496359837693, 2.725197209734958, 3.2212621733333324])
    {
        assert_close(*m, e);
    }
    let weighted = ComputePrincipalAxesAndMomentsParams {
        ignore_hs: false,
        weights: Some(vec![1.0; 9]),
    };
    let (_, moments) = compute_principal_axes_and_moments_with_params(&conf, &weighted)
        .unwrap()
        .unwrap();
    for (m, e) in moments
        .iter()
        .zip([8.659080778170027, 15.91327223250775, 18.029691104877777])
    {
        assert_close(*m, e);
    }
    let (axes, moments) = compute_principal_axes_and_moments_from_gyration_matrix(&conf)
        .unwrap()
        .unwrap();
    assert!(moments[0].abs() < 1e-12);
    assert_close(moments[1], 0.16535498786612557);
    assert_close(moments[2], 0.9083990699116524);
    assert_rows(
        &axes[..1],
        &[&[
            0.07896267254689594,
            -0.31544733007111075,
            0.9456520915724026,
        ]],
    );
    // Python requires one weight per atom here.
    let short = ComputePrincipalAxesAndMomentsFromGyrationMatrixParams {
        weights: Some(vec![1.0]),
        ..Default::default()
    };
    assert!(matches!(
        compute_principal_axes_and_moments_from_gyration_matrix_with_params(&conf, &short),
        Err(Error::Value { .. })
    ));
}

#[test]
fn internal_coordinates() {
    let mut mol = ethanol();
    {
        let conf = mol.get_conformer().unwrap();
        assert_close(get_bond_length(&conf, 0, 1).unwrap(), 1.5231098876968794);
        assert_close(get_angle_deg(&conf, 0, 1, 2).unwrap(), 106.73896473395976);
        assert_close(get_angle_rad(&conf, 0, 1, 2).unwrap(), 1.862946374777711);
        assert_close(
            get_dihedral_deg(&conf, 3, 0, 1, 2).unwrap(),
            106.12122690450155,
        );
        assert_close(
            get_dihedral_rad(&conf, 3, 0, 1, 2).unwrap(),
            1.852164815739542,
        );
        // Python raises RDKit's RuntimeError for these.
        let index = Error::Index { index: 9, len: 9 };
        assert_eq!(get_bond_length(&conf, 0, 9), Err(index.clone()));
        assert_eq!(get_angle_deg(&conf, 9, 0, 1), Err(index.clone()));
        assert_eq!(get_dihedral_rad(&conf, 0, 1, 2, 9), Err(index));
    }
    let mut conf = mol.get_conformer_mut().unwrap();
    set_bond_length(&mut conf, 0, 1, 1.6).unwrap();
    assert_close(get_bond_length(&conf, 0, 1).unwrap(), 1.6);
    assert_point(
        conf.get_atom_position(2).unwrap(),
        [1.2978130858541683, 0.3252639230417486, 0.7668520258838434],
    );
    set_angle_deg(&mut conf, 0, 1, 2, 120.0).unwrap();
    assert_close(get_angle_deg(&conf, 0, 1, 2).unwrap(), 120.0);
    assert_point(
        conf.get_atom_position(8).unwrap(),
        [2.1777779572953686, 0.8050407671262714, 0.19609624590506405],
    );
    set_dihedral_deg(&mut conf, 3, 0, 1, 2, 60.0).unwrap();
    assert_close(get_dihedral_deg(&conf, 3, 0, 1, 2).unwrap(), 60.0);
    assert_point(
        conf.get_atom_position(2).unwrap(),
        [1.5224955437473264, 0.6817945132589264, -0.17886017889626696],
    );
    set_angle_rad(&mut conf, 0, 1, 2, 2.0).unwrap();
    set_dihedral_rad(&mut conf, 3, 0, 1, 2, 1.0).unwrap();
    assert_close(get_angle_rad(&conf, 0, 1, 2).unwrap(), 2.0);
    assert_close(get_dihedral_rad(&conf, 3, 0, 1, 2).unwrap(), 1.0);
    assert!(matches!(
        set_bond_length(&mut conf, 0, 2, 1.5),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        set_dihedral_deg(&mut conf, 0, 1, 2, 9, 1.0),
        Err(Error::Index { index: 9, len: 9 })
    );

    let mut ring = mol_from_smiles("C1CC1").unwrap();
    let mut conf = Conformer::new(3);
    for i in 0..3 {
        let x = f64::from(i);
        conf.set_atom_position(i, Point3D::new(x, x * x, 0.0))
            .unwrap();
    }
    ring.add_conformer(&conf).unwrap();
    assert!(matches!(
        set_bond_length(&mut ring.get_conformer_mut().unwrap(), 0, 1, 1.5),
        Err(Error::Value { .. })
    ));
}

#[test]
fn transforms() {
    let mut mol = ethanol();
    let mut trans = vec![vec![0.0; 4]; 4];
    for (i, row) in trans.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    trans[0][3] = 1.0;
    trans[1][3] = 2.0;
    transform_conformer(&mut mol.get_conformer_mut().unwrap(), &trans).unwrap();
    assert_point(
        mol.get_conformer().unwrap().get_atom_position(0).unwrap(),
        [0.03600000000000003, 2.0151, 0.041],
    );
    // Python copies rows * rows values into a 4 x 4 transform, and crashes
    // for a larger array.
    let big = vec![vec![1.0; 200]; 200];
    assert!(matches!(
        transform_conformer(&mut mol.get_conformer_mut().unwrap(), &big),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        transform_conformer(&mut mol.get_conformer_mut().unwrap(), &trans[..3]),
        Err(Error::Value { .. })
    ));

    let mut mol = ethanol();
    canonicalize_conformer(&mut mol.get_conformer_mut().unwrap()).unwrap();
    let pos = mol.get_conformer().unwrap().get_atom_position(0).unwrap();
    assert_close(pos.x, -1.2017662434995573);
    assert_close(pos.y, 0.2604150690862285);
    assert!(pos.z.abs() < 1e-12);

    let mut mol = ethanol_two_confs();
    canonicalize_mol_with_params(
        &mut mol,
        &CanonicalizeMolParams {
            normalize_covar: true,
            ignore_hs: false,
        },
    )
    .unwrap();
    assert_point(
        mol.get_conformer_with_id(1)
            .unwrap()
            .get_atom_position(0)
            .unwrap(),
        [
            -0.9403362542875655,
            -0.2320832172234178,
            0.021914869411848587,
        ],
    );
    canonicalize_mol(&mut mol).unwrap();
    canonicalize_conformer_with_params(
        &mut mol.get_conformer_mut().unwrap(),
        &CanonicalizeConformerParams {
            center: Some(Point3D::new(0.0, 0.0, 0.0)),
            ..Default::default()
        },
    )
    .unwrap();
}

#[test]
fn conformers_outside_molecules() {
    let mol = ethanol();
    let mut free = Conformer::from(mol.get_conformer().unwrap());
    drop(mol);
    // Positions alone are enough here.
    assert_close(get_bond_length(&free, 0, 1).unwrap(), 1.5231098876968794);
    let mut trans = vec![vec![0.0; 4]; 4];
    for (i, row) in trans.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    transform_conformer(&mut free, &trans).unwrap();
    // The rest read the molecule. Python's copy still points at the freed
    // molecule, and crashes here.
    assert!(matches!(
        compute_centroid(&free),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        set_bond_length(&mut free, 0, 1, 1.6),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        canonicalize_conformer(&mut free),
        Err(Error::Runtime { .. })
    ));
}
