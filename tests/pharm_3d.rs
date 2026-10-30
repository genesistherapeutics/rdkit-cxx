// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.Pharm3D in RDKit 2026.03.6, and the tests
// of UnitTestEmbed.py (with its doctests), UnitTestExcludedVolume.py, and
// UnitTestPharmacophore.py. The CDK2 molecules are UnitTestEmbed.py's
// pickles as SMILES with their hydrogens; Python's values here are for
// those molecules and bounds computed from them, as its test3Embed and
// test4Search compute them. The Issue268 pharmacophore is its pickle's
// features and bounds.

use rdkit_cxx::Error;
use rdkit_cxx::chem::pharm_3d::embed_lib::*;
use rdkit_cxx::chem::pharm_3d::excluded_volume::{ExcludedVolume, ExcludedVolumeParams};
use rdkit_cxx::chem::pharm_3d::pharmacophore::{ExplicitPharmacophore, Pharmacophore};
use rdkit_cxx::chem::rd_chemical_features::FreeChemicalFeature;
use rdkit_cxx::chem::rd_dist_geom::get_molecule_bounds_matrix;
use rdkit_cxx::chem::rd_mol_chemical_features::{
    MolChemicalFeatureFactory, build_feature_factory, build_feature_factory_from_string,
};
use rdkit_cxx::chem::{
    Mol, SmilesParserParams, mol_from_mol_file, mol_from_smiles, mol_from_smiles_with_params,
};
use rdkit_cxx::distance_geometry::do_triangle_smoothing;
use rdkit_cxx::geometry::rd_geometry::Point3D;

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

const FDEF: &str = "
                   DefineFeature HAcceptor1 [N,O;H0]
                      Family HBondAcceptor
                      Weights 1.0
                   EndFeature
                   DefineFeature HDonor1 [N,O;!H0]
                      Family HBondDonor
                      Weights 1.0
                   EndFeature
                   DefineFeature Aromatic1 c1ccccc1
                      Family Aromatic
                      Weights 1.,1.,1.,1.,1.,1.
                   EndFeature
";

fn feat(family: &str, type_: &str, x: f64, y: f64) -> FreeChemicalFeature {
    FreeChemicalFeature::new(family, type_, Point3D::new(x, y, 0.0))
}

fn three_feats() -> Vec<FreeChemicalFeature> {
    vec![
        feat("HBondAcceptor", "HAcceptor1", 0.0, 0.0),
        feat("HBondDonor", "HDonor1", 2.65, 0.0),
        feat("Aromatic", "Aromatic1", 5.12, 0.908),
    ]
}

/// UnitTestEmbed.py's pharmacophore.
fn embed_pcophore() -> Pharmacophore {
    let mut p = Pharmacophore::new(&three_feats());
    for (i, j, low, up) in [(0, 1, 2.0, 3.3), (0, 2, 5.0, 5.4), (1, 2, 2.6, 3.0)] {
        p.set_lower_bound(i, j, low).unwrap();
        p.set_upper_bound(i, j, up).unwrap();
    }
    p
}

fn cdk2_mols() -> Vec<(String, Mol)> {
    let params = SmilesParserParams {
        remove_hs: false,
        ..Default::default()
    };
    std::fs::read_to_string(data("cdk2-syn-clip100.smi"))
        .unwrap()
        .lines()
        .map(|line| {
            let (smiles, name) = line.split_once(' ').unwrap();
            (
                name.to_owned(),
                mol_from_smiles_with_params(smiles, &params).unwrap(),
            )
        })
        .collect()
}

fn smoothed_bounds(mol: &Mol) -> Vec<Vec<f64>> {
    let mut bm = get_molecule_bounds_matrix(mol).unwrap();
    do_triangle_smoothing(&mut bm).unwrap();
    bm
}

fn fmt3(row: &[f64]) -> String {
    row.iter()
        .map(|x| format!("{x:.3}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn column(bm: &[Vec<f64>], j: usize) -> Vec<f64> {
    bm.iter().map(|row| row[j]).collect()
}

#[test]
fn pharmacophore_matches_python() {
    let mut pcophore = Pharmacophore::new(&three_feats());
    assert_eq!(pcophore.get_features().len(), 3);
    assert!(pcophore.get_feature(2).is_ok());
    assert!(matches!(pcophore.get_feature(3), Err(Error::Index { .. })));
    assert!(pcophore.to_string().contains("Aromatic"));

    pcophore.set_upper_bound(0, 1, 3.0).unwrap();
    assert_eq!(
        (
            pcophore.get_upper_bound(0, 1).unwrap(),
            pcophore.get_upper_bound(1, 0).unwrap()
        ),
        (3.0, 3.0)
    );
    pcophore.set_upper_bound(1, 0, 5.0).unwrap();
    assert_eq!(
        (
            pcophore.get_upper_bound(0, 1).unwrap(),
            pcophore.get_upper_bound(1, 0).unwrap()
        ),
        (5.0, 5.0)
    );
    pcophore.set_lower_bound(0, 1, 2.0).unwrap();
    assert_eq!(
        (
            pcophore.get_lower_bound(0, 1).unwrap(),
            pcophore.get_lower_bound(1, 0).unwrap()
        ),
        (2.0, 2.0)
    );
    pcophore.set_lower_bound(1, 0, 3.0).unwrap();
    assert_eq!(
        (
            pcophore.get_lower_bound(0, 1).unwrap(),
            pcophore.get_lower_bound(1, 0).unwrap()
        ),
        (3.0, 3.0)
    );
    for (i, j) in [(0, 3), (3, 0)] {
        assert!(matches!(
            pcophore.set_upper_bound(i, j, 2.0),
            Err(Error::Index { .. })
        ));
        assert!(matches!(
            pcophore.set_upper_bound_with_check_bounds(i, j, 2.0, true),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            pcophore.set_lower_bound(i, j, 2.0),
            Err(Error::Index { .. })
        ));
        assert!(matches!(
            pcophore.set_lower_bound_with_check_bounds(i, j, 2.0, true),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            pcophore.set_upper_bound_2d(i, j, 2),
            Err(Error::Index { .. })
        ));
        assert!(matches!(
            pcophore.set_upper_bound_2d_with_check_bounds(i, j, 2, true),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            pcophore.set_lower_bound_2d(i, j, 2),
            Err(Error::Index { .. })
        ));
        assert!(matches!(
            pcophore.set_lower_bound_2d_with_check_bounds(i, j, 2, true),
            Err(Error::Value { .. })
        ));
    }
    pcophore.set_upper_bound_2d(0, 1, 3).unwrap();
    assert_eq!(
        (
            pcophore.get_upper_bound_2d(0, 1).unwrap(),
            pcophore.get_upper_bound_2d(1, 0).unwrap()
        ),
        (3, 3)
    );
    pcophore.set_upper_bound_2d(1, 0, 5).unwrap();
    assert_eq!(pcophore.get_upper_bound_2d(0, 1).unwrap(), 5);
    pcophore.set_lower_bound_2d(0, 1, 3).unwrap();
    assert_eq!(
        (
            pcophore.get_lower_bound_2d(0, 1).unwrap(),
            pcophore.get_lower_bound_2d(1, 0).unwrap()
        ),
        (3, 3)
    );
    pcophore.set_lower_bound_2d(1, 0, 5).unwrap();
    assert_eq!(pcophore.get_lower_bound_2d(0, 1).unwrap(), 5);
    assert_eq!(
        Pharmacophore::new(&three_feats())
            .get_upper_bound_2d(0, 2)
            .unwrap(),
        1000
    );
    assert_eq!(
        Pharmacophore::new_with_init_mats(&three_feats(), false)
            .get_upper_bound(0, 1)
            .unwrap(),
        0.0
    );

    assert_eq!(
        embed_pcophore().to_string(),
        "              HBondAcceptor    HBondDonor      Aromatic \nHBondAcceptor         0.000         3.300         5.400 \n   HBondDonor         2.000         0.000         3.000 \n     Aromatic         5.000         2.600         0.000 \n"
    );
}

#[test]
fn explicit_pharmacophore_matches_python() {
    let mut p = ExplicitPharmacophore::default();
    // Python splits at the characters \\n, not at line ends.
    p.init_from_string("Acceptor 1.0 0.0 -2.5 1.25 # comment\\n\\n  Donor\t0 2 3 0.5")
        .unwrap();
    assert_eq!(p.get_features().len(), 2);
    assert_eq!(p.get_radii(), [1.25, 0.5]);
    assert_eq!(p.get_feature(1).unwrap().get_family(), "Donor");
    assert_eq!(p.get_feature(1).unwrap().get_type(), "Donor");
    assert_eq!(
        p.to_string(),
        "    Acceptor      1.0000   0.0000  -2.5000     1.25\n       Donor      0.0000   2.0000   3.0000     0.50\n"
    );
    p.set_radius(0, 2.0).unwrap();
    assert_eq!(p.get_radius(0).unwrap(), 2.0);
    assert!(matches!(p.get_radius(2), Err(Error::Index { .. })));
    p.init_from_file(&mut "Aromatic 1 1 1 1.5\n".as_bytes())
        .unwrap();
    assert_eq!(p.get_radii(), [1.5]);
    assert!(matches!(
        p.init_from_lines(&["Aromatic 1 1 1"]),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        p.init_from_lines(&["Aromatic 1 x 1 2"]),
        Err(Error::Parse { .. })
    ));
    assert_eq!(p.get_radii(), [1.5]);
    let p = ExplicitPharmacophore::new(&three_feats(), &[1.0, 2.0, 3.0]).unwrap();
    assert_eq!(p.get_radii(), [1.0, 2.0, 3.0]);
    assert!(matches!(
        ExplicitPharmacophore::new(&three_feats(), &[1.0]),
        Err(Error::Value { .. })
    ));
}

#[test]
#[allow(clippy::approx_constant)] // Python's test uses 3.14
fn excluded_volume_matches_python() {
    assert!(matches!(
        ExcludedVolume::new(vec![]),
        Err(Error::Value { .. })
    ));
    let feat_info = vec![(vec![0], 0.5, 1.0)];
    let ev = ExcludedVolume::new(feat_info.clone()).unwrap();
    assert_eq!(
        (ev.feat_info.clone(), ev.exclusion_dist, ev.index, ev.pos),
        (feat_info.clone(), 3.0, None, None)
    );
    let params = ExcludedVolumeParams {
        exclusion_dist: 3.14,
        ..Default::default()
    };
    assert_eq!(
        ExcludedVolume::new_with_params(feat_info, &params)
            .unwrap()
            .exclusion_dist,
        3.14
    );
}

#[test]
fn embed_lib_doctests() {
    let m = mol_from_smiles("CCO").unwrap();
    let idx = |atoms: Vec<rdkit_cxx::chem::AtomRef<'_>>| {
        atoms.iter().map(|a| a.get_idx()).collect::<Vec<_>>()
    };
    assert_eq!(
        idx(get_atom_heavy_neighbors(&m.get_atom_with_idx(0).unwrap())),
        [1]
    );
    assert_eq!(
        idx(get_atom_heavy_neighbors(&m.get_atom_with_idx(1).unwrap())),
        [0, 2]
    );

    let bounds = vec![
        vec![0.0, 2.0, 2.0],
        vec![1.0, 0.0, 2.0],
        vec![1.0, 1.0, 0.0],
    ];
    let params = ReplaceGroupParams {
        slop: 0.0,
        ..Default::default()
    };
    let (bm, idx) = replace_group_with_params(&[0, 1, 2], &bounds, &params).unwrap();
    assert_eq!((idx, bm.len(), bm[0].len()), (3, 4, 4));
    assert_eq!(fmt3(&bm[3]), "0.577, 0.577, 0.577, 0.000");
    assert_eq!(fmt3(&column(&bm, 3)), "1.155, 1.155, 1.155, 0.000");
    let (bm, _) = replace_group(&[0, 1, 2], &bounds).unwrap();
    assert_eq!(fmt3(&bm[3]), "0.572, 0.572, 0.572, 0.000");
    assert_eq!(fmt3(&column(&bm, 3)), "1.166, 1.166, 1.166, 0.000");
    let params = ReplaceGroupParams {
        use_dirs: true,
        ..Default::default()
    };
    let (bm, _) = replace_group_with_params(&[0, 1, 2], &bounds, &params).unwrap();
    assert_eq!(bm.len(), 5);
    assert_eq!(fmt3(&bm[4]), "2.070, 2.070, 2.070, 1.990, 0.000");
    assert_eq!(fmt3(&column(&bm, 4)), "2.324, 2.324, 2.324, 2.010, 0.000");

    let mut m = mol_from_smiles("c1ccccc1C").unwrap();
    let bounds7 = get_molecule_bounds_matrix(&m).unwrap();
    assert_eq!(m.get_num_conformers(), 0);
    embed_mol_with_params(
        &mut m,
        &bounds7,
        &EmbedMolParams {
            random_seed: 23,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(m.get_num_conformers(), 1);

    let mut evs = [ExcludedVolume::new_with_params(
        vec![(vec![0], 0.5, 1.0)],
        &ExcludedVolumeParams {
            exclusion_dist: 1.5,
            ..Default::default()
        },
    )
    .unwrap()];
    let bm = add_excluded_volumes(&bounds, &mut evs).unwrap();
    assert_eq!((bm.len(), evs[0].index), (4, Some(3)));
    assert_eq!(fmt3(&bm[3]), "0.500, 1.500, 1.500, 0.000");
    assert_eq!(fmt3(&column(&bm, 3)), "1.000, 3.000, 3.000, 0.000");

    let mut pcophore = Pharmacophore::new(&three_feats()[..2]);
    pcophore.set_lower_bound(0, 1, 1.0).unwrap();
    pcophore.set_upper_bound(0, 1, 2.0).unwrap();
    let bm = vec![
        vec![0.0, 3.0, 3.0],
        vec![2.0, 0.0, 3.0],
        vec![2.0, 2.0, 0.0],
    ];
    let bm = update_pharmacophore_bounds(bm, &[vec![0], vec![1]], &pcophore).unwrap();
    assert_eq!(
        [fmt3(&bm[0]), fmt3(&bm[1]), fmt3(&bm[2])],
        [
            "0.000, 2.000, 3.000",
            "1.000, 0.000, 3.000",
            "2.000, 2.000, 0.000"
        ]
    );

    let m = mol_from_smiles("OCCN").unwrap();
    let mut pcophore = Pharmacophore::new(&three_feats()[..2]);
    pcophore.set_lower_bound(0, 1, 2.5).unwrap();
    pcophore.set_upper_bound(0, 1, 3.5).unwrap();
    let atom_match = [vec![0], vec![3]];
    let params = EmbedPharmacophoreParams {
        random_seed: 23,
        silent: true,
        ..Default::default()
    };
    let (_, embeds, n_fail) =
        embed_pharmacophore_with_params(&m, &atom_match, &pcophore, &params).unwrap();
    assert_eq!((embeds.len(), n_fail), (10, 0));
    let mut impossible = Pharmacophore::new(&three_feats()[..2]);
    impossible.set_lower_bound(0, 1, 2.0).unwrap();
    impossible.set_upper_bound(0, 1, 2.1).unwrap();
    let (_, embeds, n_fail) =
        embed_pharmacophore_with_params(&m, &atom_match, &impossible, &params).unwrap();
    assert_eq!((embeds.len(), n_fail), (0, 10));

    assert!(!is_nan(0.0));

    pcophore.set_upper_bound(0, 1, 2.8).unwrap();
    let (bm, mut embeds, _) =
        embed_pharmacophore_with_params(&m, &atom_match, &pcophore, &params).unwrap();
    assert_eq!(embeds.len(), 10);
    let opt = OptimizeMolParams {
        atom_matches: Some(atom_match.to_vec()),
        ..Default::default()
    };
    let (e1, e2) = optimize_mol_with_params(&mut embeds[0], &bm, &opt).unwrap();
    assert!(e2 < e1);
    let conf = embeds[0].get_conformer().unwrap();
    let d03 = conf
        .get_atom_position(0)
        .unwrap()
        .distance(&conf.get_atom_position(3).unwrap());
    assert!((2.5 - 0.01..=2.8 + 0.01).contains(&d03), "{d03}");
    let (e1, e2) = optimize_mol(&mut embeds[1], &bm).unwrap();
    assert!(e2 < e1);
    let conf = embeds[1].get_conformer().unwrap();
    let d03 = conf
        .get_atom_position(0)
        .unwrap()
        .distance(&conf.get_atom_position(3).unwrap());
    assert!(d03 >= 2.5 - 0.01, "{d03}");
    assert!(d03 > 2.8 + 0.01, "{d03}");
    // A molecule without coordinates has no force field: Python's (-1, -1).
    assert_eq!(
        optimize_mol(&mut mol_from_smiles("OCCN").unwrap(), &bm).unwrap(),
        (-1.0, -1.0)
    );

    let factory = build_feature_factory(data("Pharm3D_BaseFeatures.fdef")).unwrap();
    let active = [
        FreeChemicalFeature::from_family("Acceptor", Point3D::new(0.0, 0.0, 0.0)),
        FreeChemicalFeature::from_family("Donor", Point3D::new(0.0, 0.0, 0.0)),
    ];
    let m = mol_from_smiles("FCCN").unwrap();
    let (matched, list) =
        match_pharmacophore_to_mol(&m, &factory, &Pharmacophore::new(&active)).unwrap();
    assert!(matched);
    let ids: Vec<Vec<Vec<u32>>> = list
        .unwrap()
        .iter()
        .map(|fs| fs.iter().map(|f| f.get_atom_ids()).collect())
        .collect();
    assert_eq!(ids, [vec![vec![0], vec![3]], vec![vec![3]]]);
    let (matched, list) = match_feats_to_mol(
        &m,
        &factory,
        &[FreeChemicalFeature::from_family(
            "Aromatic",
            Point3D::new(0.0, 0.0, 0.0),
        )],
    )
    .unwrap();
    assert!(!matched && list.is_none());

    assert_eq!(
        combi_enum(&[vec![1, 2], vec![10, 20]]),
        [vec![1, 10], vec![1, 20], vec![2, 10], vec![2, 20]]
    );
    assert_eq!(
        combi_enum(&[vec![1, 2], vec![10, 20], vec![100, 200]]).len(),
        8
    );
    assert_eq!(combi_enum::<u32>(&[]), [Vec::<u32>::new()]);

    let bounds = vec![
        vec![0.0, 3.0, 4.0],
        vec![2.0, 0.0, 3.0],
        vec![2.0, 2.0, 0.0],
    ];
    let bm = downsample_bounds_matrix_with_max_thresh(&bounds, &[0], 3.5).unwrap();
    assert_eq!(bm, [vec![0.0, 3.0], vec![2.0, 0.0]]);
    let bounds = vec![
        vec![0.0, 4.0, 3.0],
        vec![2.0, 0.0, 3.0],
        vec![2.0, 2.0, 0.0],
    ];
    assert_eq!(
        downsample_bounds_matrix_with_max_thresh(&bounds, &[0], 5.0)
            .unwrap()
            .len(),
        3
    );
    assert_eq!(
        downsample_bounds_matrix_with_max_thresh(&bounds, &[0, 1], 3.5)
            .unwrap()
            .len(),
        3
    );
    assert!(downsample_bounds_matrix(&bounds, &[]).unwrap().is_empty());

    let mut p = Pharmacophore::new(&three_feats());
    for (i, j, low, up) in [(0, 1, 1.1, 1.9), (0, 2, 2.1, 2.9), (1, 2, 2.1, 3.9)] {
        p.set_lower_bound(i, j, low).unwrap();
        p.set_upper_bound(i, j, up).unwrap();
    }
    let bounds = vec![
        vec![0.0, 2.0, 3.0],
        vec![1.0, 0.0, 4.0],
        vec![2.0, 3.0, 0.0],
    ];
    let screen = |m: &[u32]| {
        coarse_screen_pharmacophore(&m.iter().map(|&a| vec![a]).collect::<Vec<_>>(), &bounds, &p)
            .unwrap()
    };
    assert_eq!(
        [
            screen(&[0, 1]),
            screen(&[0, 2]),
            screen(&[1, 2]),
            screen(&[0, 1, 2]),
            screen(&[1, 0, 2]),
            screen(&[2, 1, 0])
        ],
        [true, false, false, true, false, false]
    );
    let mut four = three_feats();
    four.push(feat("HBondDonor", "HDonor1", 2.65, 0.0));
    let mut p = Pharmacophore::new(&four);
    for (i, j, low, up) in [
        (0, 1, 2.1, 2.9),
        (0, 2, 2.1, 2.9),
        (0, 3, 2.1, 2.9),
        (1, 2, 1.1, 1.9),
        (1, 3, 1.1, 1.9),
        (2, 3, 1.1, 1.9),
    ] {
        p.set_lower_bound(i, j, low).unwrap();
        p.set_upper_bound(i, j, up).unwrap();
    }
    let bounds = vec![
        vec![0.0, 3.0, 3.0, 3.0],
        vec![2.0, 0.0, 2.0, 2.0],
        vec![2.0, 1.0, 0.0, 2.0],
        vec![2.0, 1.0, 1.0, 0.0],
    ];
    let screen = |m: &[u32]| {
        coarse_screen_pharmacophore(&m.iter().map(|&a| vec![a]).collect::<Vec<_>>(), &bounds, &p)
            .unwrap()
    };
    assert_eq!(
        [
            screen(&[0, 1, 2, 3]),
            screen(&[0, 1, 3, 2]),
            screen(&[1, 0, 3, 2])
        ],
        [true, true, false]
    );

    let mut p = Pharmacophore::new(&active);
    p.set_upper_bound_2d(0, 1, 3).unwrap();
    let m = mol_from_smiles("FCC(N)CN").unwrap();
    assert!(check_2d_bounds(&[vec![0], vec![3]], &m, &p).unwrap());
    assert!(!check_2d_bounds(&[vec![0], vec![5]], &m, &p).unwrap());

    for (name, sign) in [
        ("mol-r.mol", -1.0),
        ("mol-s.mol", 1.0),
        ("mol-r-3.mol", -1.0),
        ("mol-s-3.mol", 1.0),
    ] {
        let mut mol = mol_from_mol_file(data(&format!("Pharm3D_{name}"))).unwrap();
        let vol = compute_chiral_volume(&mut mol, 1).unwrap();
        assert!(vol * sign > 3.0, "{name}: {vol}");
        assert_eq!(
            compute_chiral_volume_with_conf_id(&mut mol, 0, -1).unwrap(),
            0.0
        );
    }
}

/// A molecule's name, its bounds' sizes, and its features' atoms.
type Hit<S> = (S, (usize, usize), Vec<Vec<u32>>);

/// Each CDK2 molecule that matches the pharmacophore, with the sizes of its
/// bounds before and after downsampling and its features' atoms.
fn expected_hits(downsample: bool) -> Vec<Hit<&'static str>> {
    if downsample {
        vec![
            (
                "mol_102",
                (58, 29),
                vec![vec![5], vec![14], vec![40, 41, 43, 45, 47, 49]],
            ),
            (
                "mol_107",
                (49, 24),
                vec![vec![17], vec![34], vec![1, 2, 4, 6, 44, 46]],
            ),
            (
                "mol_126",
                (60, 32),
                vec![vec![17], vec![23], vec![1, 2, 4, 6, 55, 57]],
            ),
            (
                "mol_127",
                (47, 23),
                vec![vec![8], vec![10], vec![32, 33, 35, 37, 39, 41]],
            ),
            (
                "mol_129",
                (60, 30),
                vec![vec![5], vec![14], vec![42, 43, 45, 47, 49, 51]],
            ),
            (
                "mol_173",
                (39, 34),
                vec![vec![14], vec![1], vec![24, 25, 27, 29, 31, 33]],
            ),
            (
                "mol_174",
                (50, 36),
                vec![vec![8], vec![1], vec![29, 30, 32, 34, 38, 40]],
            ),
            (
                "mol_176",
                (40, 34),
                vec![vec![14], vec![1], vec![25, 26, 28, 30, 32, 34]],
            ),
            (
                "mol_178",
                (42, 33),
                vec![vec![14], vec![1], vec![27, 28, 30, 32, 34, 36]],
            ),
            (
                "mol_181",
                (49, 32),
                vec![vec![14], vec![1], vec![34, 35, 37, 39, 41, 43]],
            ),
            (
                "mol_197",
                (47, 32),
                vec![vec![9], vec![19], vec![1, 2, 4, 40, 42, 44]],
            ),
            (
                "mol_213",
                (44, 20),
                vec![vec![3], vec![1], vec![29, 30, 32, 34, 36, 38]],
            ),
            (
                "mol_214",
                (57, 24),
                vec![vec![3], vec![1], vec![7, 8, 10, 12, 52, 54]],
            ),
            (
                "mol_222",
                (56, 23),
                vec![vec![43], vec![1], vec![2, 3, 5, 7, 44, 46]],
            ),
            (
                "mol_223",
                (39, 22),
                vec![vec![3], vec![1], vec![13, 14, 16, 18, 20, 22]],
            ),
            (
                "mol_244",
                (51, 28),
                vec![vec![22], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_246",
                (61, 34),
                vec![vec![16], vec![30], vec![1, 2, 4, 6, 56, 58]],
            ),
            (
                "mol_247",
                (50, 27),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_248",
                (50, 27),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_249",
                (37, 26),
                vec![vec![19], vec![14], vec![3, 4, 6, 8, 10, 12]],
            ),
            (
                "mol_251",
                (39, 28),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_252",
                (40, 29),
                vec![vec![22], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_253",
                (39, 28),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_254",
                (50, 26),
                vec![vec![16], vec![41], vec![1, 2, 4, 6, 45, 47]],
            ),
            (
                "mol_255",
                (47, 30),
                vec![vec![22], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_257",
                (44, 27),
                vec![vec![19], vec![14], vec![3, 4, 6, 8, 10, 12]],
            ),
            (
                "mol_258",
                (46, 29),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_259",
                (46, 29),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_260",
                (57, 26),
                vec![vec![16], vec![48], vec![1, 2, 4, 6, 52, 54]],
            ),
            (
                "mol_261",
                (56, 24),
                vec![vec![17], vec![41], vec![1, 2, 4, 6, 51, 53]],
            ),
            (
                "mol_262",
                (56, 21),
                vec![vec![11], vec![13], vec![28, 29, 31, 33, 35, 37]],
            ),
            (
                "mol_263",
                (37, 17),
                vec![vec![35], vec![31], vec![20, 21, 23, 25, 27, 29]],
            ),
            (
                "mol_264",
                (41, 25),
                vec![vec![37], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_265",
                (40, 22),
                vec![vec![11], vec![34], vec![23, 24, 26, 28, 30, 32]],
            ),
            (
                "mol_266",
                (40, 22),
                vec![vec![11], vec![34], vec![23, 24, 26, 28, 30, 32]],
            ),
            (
                "mol_267",
                (50, 19),
                vec![vec![8], vec![10], vec![14, 15, 17, 19, 21, 23]],
            ),
            (
                "mol_268",
                (48, 27),
                vec![vec![25], vec![3], vec![6, 7, 9, 11, 13, 15]],
            ),
            (
                "mol_269",
                (37, 24),
                vec![vec![33], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_270",
                (38, 25),
                vec![vec![34], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_271",
                (37, 24),
                vec![vec![33], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_272",
                (35, 22),
                vec![vec![31], vec![14], vec![3, 4, 6, 8, 10, 12]],
            ),
            (
                "mol_274",
                (42, 29),
                vec![vec![16], vec![30], vec![1, 2, 4, 6, 37, 39]],
            ),
            (
                "mol_275",
                (31, 27),
                vec![vec![19], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_276",
                (32, 28),
                vec![vec![20], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_278",
                (32, 24),
                vec![vec![18], vec![5], vec![1, 2, 4, 8, 10, 12]],
            ),
            (
                "mol_279",
                (31, 27),
                vec![vec![19], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_282",
                (58, 28),
                vec![vec![3], vec![1], vec![17, 18, 20, 26, 28, 30]],
            ),
            (
                "mol_283",
                (54, 24),
                vec![vec![3], vec![1], vec![17, 18, 20, 22, 25, 27]],
            ),
            (
                "mol_285",
                (61, 34),
                vec![vec![19], vec![1], vec![29, 30, 32, 34, 37, 39]],
            ),
            (
                "mol_286",
                (62, 35),
                vec![vec![19], vec![1], vec![29, 30, 32, 34, 36, 40]],
            ),
        ]
    } else {
        vec![
            (
                "mol_102",
                (58, 58),
                vec![vec![5], vec![14], vec![40, 41, 43, 45, 47, 49]],
            ),
            (
                "mol_107",
                (49, 49),
                vec![vec![17], vec![34], vec![1, 2, 4, 6, 44, 46]],
            ),
            (
                "mol_126",
                (60, 60),
                vec![vec![17], vec![23], vec![1, 2, 4, 6, 55, 57]],
            ),
            (
                "mol_127",
                (47, 47),
                vec![vec![8], vec![10], vec![32, 33, 35, 37, 39, 41]],
            ),
            (
                "mol_129",
                (60, 60),
                vec![vec![5], vec![14], vec![42, 43, 45, 47, 49, 51]],
            ),
            (
                "mol_173",
                (39, 39),
                vec![vec![14], vec![1], vec![24, 25, 27, 29, 31, 33]],
            ),
            (
                "mol_174",
                (50, 50),
                vec![vec![8], vec![1], vec![29, 30, 32, 34, 38, 40]],
            ),
            (
                "mol_176",
                (40, 40),
                vec![vec![14], vec![1], vec![25, 26, 28, 30, 32, 34]],
            ),
            (
                "mol_178",
                (42, 42),
                vec![vec![14], vec![1], vec![27, 28, 30, 32, 34, 36]],
            ),
            (
                "mol_181",
                (49, 49),
                vec![vec![14], vec![1], vec![34, 35, 37, 39, 41, 43]],
            ),
            (
                "mol_197",
                (47, 47),
                vec![vec![9], vec![19], vec![1, 2, 4, 40, 42, 44]],
            ),
            (
                "mol_213",
                (44, 44),
                vec![vec![3], vec![1], vec![29, 30, 32, 34, 36, 38]],
            ),
            (
                "mol_214",
                (57, 57),
                vec![vec![3], vec![1], vec![7, 8, 10, 12, 52, 54]],
            ),
            (
                "mol_222",
                (56, 56),
                vec![vec![43], vec![1], vec![2, 3, 5, 7, 44, 46]],
            ),
            (
                "mol_223",
                (39, 39),
                vec![vec![3], vec![1], vec![13, 14, 16, 18, 20, 22]],
            ),
            (
                "mol_244",
                (51, 51),
                vec![vec![22], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_246",
                (61, 61),
                vec![vec![16], vec![30], vec![1, 2, 4, 6, 56, 58]],
            ),
            (
                "mol_247",
                (50, 50),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_248",
                (50, 50),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_249",
                (37, 37),
                vec![vec![19], vec![14], vec![3, 4, 6, 8, 10, 12]],
            ),
            (
                "mol_251",
                (39, 39),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_252",
                (40, 40),
                vec![vec![22], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_253",
                (39, 39),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_254",
                (50, 50),
                vec![vec![16], vec![41], vec![1, 2, 4, 6, 45, 47]],
            ),
            (
                "mol_255",
                (47, 47),
                vec![vec![22], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_257",
                (44, 44),
                vec![vec![19], vec![14], vec![3, 4, 6, 8, 10, 12]],
            ),
            (
                "mol_258",
                (46, 46),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_259",
                (46, 46),
                vec![vec![21], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_260",
                (57, 57),
                vec![vec![16], vec![48], vec![1, 2, 4, 6, 52, 54]],
            ),
            (
                "mol_261",
                (56, 56),
                vec![vec![17], vec![41], vec![1, 2, 4, 6, 51, 53]],
            ),
            (
                "mol_262",
                (56, 56),
                vec![vec![11], vec![13], vec![28, 29, 31, 33, 35, 37]],
            ),
            (
                "mol_263",
                (37, 37),
                vec![vec![35], vec![31], vec![20, 21, 23, 25, 27, 29]],
            ),
            (
                "mol_264",
                (41, 41),
                vec![vec![37], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_265",
                (40, 40),
                vec![vec![11], vec![34], vec![23, 24, 26, 28, 30, 32]],
            ),
            (
                "mol_266",
                (40, 40),
                vec![vec![11], vec![34], vec![23, 24, 26, 28, 30, 32]],
            ),
            (
                "mol_267",
                (50, 50),
                vec![vec![8], vec![10], vec![14, 15, 17, 19, 21, 23]],
            ),
            (
                "mol_268",
                (48, 48),
                vec![vec![25], vec![3], vec![6, 7, 9, 11, 13, 15]],
            ),
            (
                "mol_269",
                (37, 37),
                vec![vec![33], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_270",
                (38, 38),
                vec![vec![34], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_271",
                (37, 37),
                vec![vec![33], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_272",
                (35, 35),
                vec![vec![31], vec![14], vec![3, 4, 6, 8, 10, 12]],
            ),
            (
                "mol_274",
                (42, 42),
                vec![vec![16], vec![30], vec![1, 2, 4, 6, 37, 39]],
            ),
            (
                "mol_275",
                (31, 31),
                vec![vec![19], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_276",
                (32, 32),
                vec![vec![20], vec![17], vec![3, 4, 6, 8, 10, 15]],
            ),
            (
                "mol_278",
                (32, 32),
                vec![vec![18], vec![5], vec![1, 2, 4, 8, 10, 12]],
            ),
            (
                "mol_279",
                (31, 31),
                vec![vec![19], vec![16], vec![3, 4, 6, 8, 10, 14]],
            ),
            (
                "mol_282",
                (58, 58),
                vec![vec![3], vec![1], vec![17, 18, 20, 26, 28, 30]],
            ),
            (
                "mol_283",
                (54, 54),
                vec![vec![3], vec![1], vec![17, 18, 20, 22, 25, 27]],
            ),
            (
                "mol_285",
                (61, 61),
                vec![vec![19], vec![1], vec![29, 30, 32, 34, 37, 39]],
            ),
            (
                "mol_286",
                (62, 62),
                vec![vec![19], vec![1], vec![29, 30, 32, 34, 36, 40]],
            ),
        ]
    }
}

fn search(factory: &MolChemicalFeatureFactory, downsample: bool) -> Vec<Hit<String>> {
    let pcophore = embed_pcophore();
    let mut hits = Vec::new();
    for (name, mol) in cdk2_mols() {
        let bm = smoothed_bounds(&mol);
        let (matched, matches) = match_pharmacophore_to_mol(&mol, factory, &pcophore).unwrap();
        if !matched {
            continue;
        }
        let matches = matches.unwrap();
        let params = MatchPharmacophoreParams {
            use_downsampling: downsample,
            ..Default::default()
        };
        if let Some((_, m, sizes)) =
            match_pharmacophore_with_params(&matches, &bm, &pcophore, &params).unwrap()
        {
            hits.push((name, sizes, m.iter().map(|f| f.get_atom_ids()).collect()));
        }
    }
    hits
}

#[test]
fn search_matches_python() {
    // UnitTestEmbed.test1SearchFullMat and test2SearchDownsample.
    let factory = build_feature_factory_from_string(FDEF).unwrap();
    for downsample in [false, true] {
        let hits = search(&factory, downsample);
        let expected = expected_hits(downsample);
        assert_eq!(hits.len(), 50);
        for (hit, (name, sizes, atoms)) in hits.iter().zip(expected) {
            assert_eq!((hit.0.as_str(), hit.1, &hit.2), (name, sizes, &atoms));
        }
    }
}

#[test]
fn embed_one_matches_python() {
    // UnitTestEmbed.test3Embed, with its tolerance of 5: the energies come
    // from embedding and minimization, which round differently by platform.
    let factory = build_feature_factory_from_string(FDEF).unwrap();
    let pcophore = embed_pcophore();
    let expected = [
        (
            "mol_197",
            [
                177.59925973890168,
                19.410558035846492,
                93.20647801047038,
                3.931959508968619,
                93.1930950374198,
                3.9296131438347905,
                90.81816769073907,
                2.8295016789475893,
            ],
            0,
        ),
        (
            "mol_223",
            [
                246.38325509625798,
                9.653702341928259,
                134.7415636830309,
                0.24243871710543108,
                134.6678607320613,
                0.2429898998790186,
                85.15405602107276,
                0.2323699007677095,
            ],
            0,
        ),
        (
            "mol_269",
            [
                205.91953762892723,
                5.234913067327813,
                104.44755172369287,
                1.0113131534548931,
                103.2969649636922,
                1.1207796140214534,
                88.80465666793019,
                1.270281410963396,
            ],
            8,
        ),
    ];
    let mols = cdk2_mols();
    for (name, stats, n_failed) in expected {
        let mol = &mols.iter().find(|(n, _)| n == name).unwrap().1;
        let bm = smoothed_bounds(mol);
        let (_, matches) = match_pharmacophore_to_mol(mol, &factory, &pcophore).unwrap();
        let matches = matches.unwrap();
        let params = MatchPharmacophoreParams {
            use_downsampling: true,
            ..Default::default()
        };
        let (_, m, _) = match_pharmacophore_with_params(&matches, &bm, &pcophore, &params)
            .unwrap()
            .unwrap();
        // The features matched, owned again, as Python's match holds them.
        let owned: Vec<_> = m
            .iter()
            .map(|f| {
                let feats = (0..factory.get_num_mol_features(mol).unwrap())
                    .map(|i| factory.get_mol_feature(mol, i).unwrap())
                    .find(|g| {
                        g.get_family() == f.get_family() && g.get_atom_ids() == f.get_atom_ids()
                    });
                feats.unwrap()
            })
            .collect();
        let params = EmbedOneParams {
            count: 10,
            silent: true,
            embed: EmbedPharmacophoreParams {
                random_seed: 23,
                ..Default::default()
            },
        };
        let got = embed_one_with_params(mol, name, &owned, &pcophore, &params).unwrap();
        let got_values = [
            got.e1_mean,
            got.e2_mean,
            got.e2_dev,
            got.e3_mean,
            got.e3_dev,
            got.e4_mean,
            got.e4_dev,
        ];
        let expected_values = [
            stats[0], stats[2], stats[3], stats[4], stats[5], stats[6], stats[7],
        ];
        for (g, e) in got_values.iter().zip(expected_values) {
            assert!((g - e).abs() <= 5.0, "{name}: {got:?}");
        }
        assert_eq!(got.n_failed, n_failed, "{name}");
    }
}

#[test]
fn search_with_2d_limits_matches_python() {
    // UnitTestEmbed.test4Search.
    let factory = build_feature_factory(data("Pharm3D_BaseFeatures.fdef")).unwrap();
    let origin = Point3D::new(0.0, 0.0, 0.0);
    let active =
        ["Acceptor", "Donor", "Aromatic"].map(|f| FreeChemicalFeature::from_family(f, origin));
    let mut pcophore = Pharmacophore::new(&active);
    for (i, j, low, up, up_2d) in [
        (0, 1, 2.251, 2.451, 3),
        (0, 2, 4.970, 5.170, 6),
        (1, 2, 2.681, 2.881, 6),
    ] {
        pcophore.set_lower_bound(i, j, low).unwrap();
        pcophore.set_upper_bound(i, j, up).unwrap();
        pcophore.set_upper_bound_2d(i, j, up_2d).unwrap();
    }
    let (mut n_matches, mut n_hits) = (0, 0);
    for (_, mol) in cdk2_mols() {
        let bm = smoothed_bounds(&mol);
        let (can_match, matches) = match_pharmacophore_to_mol(&mol, &factory, &pcophore).unwrap();
        if can_match {
            n_matches += 1;
            let params = MatchPharmacophoreParams {
                use_downsampling: true,
                use_2d_limits: true,
                mol: Some(&mol),
                ..Default::default()
            };
            if match_pharmacophore_with_params(&matches.unwrap(), &bm, &pcophore, &params)
                .unwrap()
                .is_some()
            {
                n_hits += 1;
            }
        }
    }
    assert_eq!((n_matches, n_hits), (93, 67));
}

#[test]
fn issue268_matches_python() {
    let factory = build_feature_factory(data("Pharm3D_Issue268.fdef")).unwrap();
    let m1 = mol_from_mol_file(data("Pharm3D_Issue268_Mol1.mol")).unwrap();
    let m2 = mol_from_mol_file(data("Pharm3D_Issue268_Mol2.mol")).unwrap();
    // The pickle's features and bounds; its positions are not used.
    let origin = Point3D::new(0.0, 0.0, 0.0);
    let feats = [
        ("Acceptor", "SingleAtomAcceptor"),
        ("Donor", "SingleAtomDonor"),
        ("Aromatic", "Arom6"),
        ("PosIonizable", "BasicGroup"),
    ]
    .map(|(f, t)| FreeChemicalFeature::new(f, t, origin));
    let mut pcop = Pharmacophore::new_with_init_mats(&feats, false);
    let bm = [
        [0.0, 2.651, 5.37, 12.305],
        [2.051, 0.0, 3.081, 10.165],
        [4.77, 2.481, 0.0, 7.467],
        [11.705, 9.565, 6.867, 0.0],
    ];
    let bm_2d = [[0, 3, 6, 16], [0, 0, 6, 13], [0, 0, 0, 12], [0, 0, 0, 0]];
    for i in 0..4 {
        for j in i + 1..4 {
            pcop.set_upper_bound(i, j, bm[i][j]).unwrap();
            pcop.set_lower_bound(i, j, bm[j][i]).unwrap();
            pcop.set_upper_bound_2d(i, j, bm_2d[i][j]).unwrap();
            pcop.set_lower_bound_2d(i, j, bm_2d[j][i]).unwrap();
        }
    }
    let (_, list1) = match_feats_to_mol(&m1, &factory, pcop.get_features()).unwrap();
    let (_, list2) = match_feats_to_mol(&m2, &factory, pcop.get_features()).unwrap();
    let (list1, list2) = (list1.unwrap(), list2.unwrap());
    let mut b1 = get_molecule_bounds_matrix(&m1).unwrap();
    let mut b2 = get_molecule_bounds_matrix(&m2).unwrap();
    for smoothed in [false, true] {
        if smoothed {
            assert!(do_triangle_smoothing(&mut b1).unwrap());
            assert!(do_triangle_smoothing(&mut b2).unwrap());
        }
        for (list, b, m) in [(&list1, &b1, &m1), (&list2, &b2, &m2)] {
            assert_eq!(
                match_pharmacophore(list, b, &pcop)
                    .unwrap()
                    .unwrap()
                    .1
                    .len(),
                4
            );
            let params = MatchPharmacophoreParams {
                use_2d_limits: true,
                mol: Some(m),
                ..Default::default()
            };
            assert_eq!(
                match_pharmacophore_with_params(list, b, &pcop, &params)
                    .unwrap()
                    .unwrap()
                    .1
                    .len(),
                4
            );
            assert!(
                !get_all_pharmacophore_matches(list, b, &pcop)
                    .unwrap()
                    .is_empty()
            );
            assert!(!constrained_enum(list, &pcop, b).unwrap().is_empty());
        }
    }
}
