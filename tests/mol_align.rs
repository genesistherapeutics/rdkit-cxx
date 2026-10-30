// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{assert_close, assert_close_to};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_align::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::force_field::MMFFMolProperties;
use rdkit_cxx::geometry::Point3D;

// Expected values below come from Python's rdMolAlign on the same
// conformers. The RMSD alignments are direct linear algebra, compared to
// 1e-9; O3A's iterative search is compared loosely, and by what any
// alignment of a conformer onto itself gives.

/// Acetanilide with its hydrogens, conformer 0 of an ETKDG embedding.
const ACETANILIDE: &str = "
     RDKit          3D

 19 19  0  0  0  0  0  0  0  0999 V2000
    3.4099    0.0861    0.5661 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9310    0.3693    0.5700 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4474    1.3672    1.1595 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0885   -0.5435   -0.1256 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3159   -0.3197   -0.1512 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8006    0.4470   -1.1926 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1381    0.7281   -1.3188 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0207    0.2316   -0.3807 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5339   -0.5391    0.6664 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1844   -0.8221    0.7927 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9037    0.4564    1.4961 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8939    0.6237   -0.3014 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5004   -1.0172    0.4358 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5202   -1.3767   -0.6145 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1047    0.8463   -1.9440 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4388    1.3332   -2.1587 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0815    0.4600   -0.4889 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2444   -0.9038    1.3693 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8322   -1.4267    1.6206 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
 10  5  1  0
  1 11  1  0
  1 12  1  0
  1 13  1  0
  4 14  1  0
  6 15  1  0
  7 16  1  0
  8 17  1  0
  9 18  1  0
 10 19  1  0
M  END

";

/// Conformers 1 and 2 of the same embedding.
const ACETANILIDE_CONFS: [&[(f64, f64, f64)]; 2] = [
    &[
        (3.5254, 0.2624, 0.0131),
        (2.1456, -0.2946, 0.0528),
        (2.0344, -1.5387, 0.0722),
        (1.0231, 0.5589, 0.0673),
        (-0.3267, 0.1731, 0.1039),
        (-0.8110, -1.1063, 0.1324),
        (-2.1801, -1.3947, 0.1681),
        (-3.1184, -0.3835, 0.1763),
        (-2.6436, 0.9068, 0.1479),
        (-1.2738, 1.1902, 0.1122),
        (3.5035, 1.3508, -0.1345),
        (4.0182, -0.1962, -0.8662),
        (4.0608, -0.0605, 0.9152),
        (1.2272, 1.5942, 0.0485),
        (-0.1620, -1.9616, 0.1288),
        (-2.4886, -2.4376, 0.1894),
        (-4.1717, -0.6234, 0.2040),
        (-3.3411, 1.7262, 0.1530),
        (-1.0213, 2.2346, 0.0922),
    ],
    &[
        (3.4371, 0.2034, 0.1611),
        (2.0001, 0.0816, 0.5250),
        (1.6438, -0.1151, 1.7241),
        (1.0581, 0.1965, -0.5343),
        (-0.3591, 0.0934, -0.2896),
        (-1.1240, 1.1907, 0.0056),
        (-2.4751, 0.9731, 0.2220),
        (-3.0342, -0.2990, 0.1442),
        (-2.2413, -1.3761, -0.1533),
        (-0.8998, -1.1599, -0.3681),
        (3.7516, -0.7349, -0.3364),
        (3.9924, 0.3749, 1.1052),
        (3.6602, 1.0487, -0.5131),
        (1.3408, 0.3596, -1.5332),
        (-0.7214, 2.1921, 0.0736),
        (-3.0863, 1.8267, 0.4546),
        (-4.0829, -0.4773, 0.3109),
        (-2.5910, -2.3939, -0.2333),
        (-0.2687, -1.9845, -0.6002),
    ],
];

fn atoms_with(positions: &[(f64, f64, f64)]) -> Mol {
    let params = MolFromMolBlockParams {
        remove_hs: false,
        ..Default::default()
    };
    let mut mol = mol_from_mol_block_with_params(ACETANILIDE, &params).unwrap();
    let points: Vec<Point3D> = positions
        .iter()
        .map(|&(x, y, z)| Point3D::new(x, y, z))
        .collect();
    if !points.is_empty() {
        mol.get_conformer_mut()
            .unwrap()
            .set_positions(&points)
            .unwrap();
    }
    mol
}

/// Conformer 0 alone.
fn reference() -> Mol {
    atoms_with(&[])
}

/// Conformer `k` alone.
fn probe(k: usize) -> Mol {
    atoms_with(ACETANILIDE_CONFS[k - 1])
}

/// All three conformers.
fn acetanilide() -> Mol {
    let mut mol = reference();
    for positions in ACETANILIDE_CONFS {
        let mut conf = Conformer::new(19);
        let points: Vec<Point3D> = positions
            .iter()
            .map(|&(x, y, z)| Point3D::new(x, y, z))
            .collect();
        conf.set_positions(&points).unwrap();
        conf.set_3d(true);
        mol.add_conformer_with_assign_id(&conf, true).unwrap();
    }
    mol
}

fn heavy_map() -> Vec<(u32, u32)> {
    (0..10).map(|i| (i, i)).collect()
}

fn first_position(mol: &Mol, conf_id: i32) -> Point3D {
    mol.get_conformer_with_id(conf_id)
        .unwrap()
        .get_atom_position(0)
        .unwrap()
}

fn assert_point(actual: Point3D, expected: (f64, f64, f64)) {
    assert_close(actual.x, expected.0);
    assert_close(actual.y, expected.1);
    assert_close(actual.z, expected.2);
}

fn assert_row(actual: &[f64], expected: [f64; 4]) {
    for (a, e) in actual.iter().zip(expected) {
        assert_close(*a, e);
    }
}

fn value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn params_defaults() {
    assert_eq!(
        GetBestAlignmentTransformParams::default(),
        GetBestAlignmentTransformParams {
            prb_cid: -1,
            ref_cid: -1,
            map: Vec::new(),
            max_matches: 1_000_000,
            symmetrize_conjugated_terminal_groups: true,
            weights: Vec::new(),
            reflect: false,
            max_iters: 50,
            num_threads: 1,
        }
    );
    assert_eq!(
        GetBestRMSParams::default(),
        GetBestRMSParams {
            prb_id: -1,
            ref_id: -1,
            map: Vec::new(),
            max_matches: 1_000_000,
            symmetrize_conjugated_terminal_groups: true,
            weights: Vec::new(),
            num_threads: 1,
        }
    );
    assert_eq!(
        GetAllConformerBestRMSParams::default(),
        GetAllConformerBestRMSParams {
            num_threads: 1,
            map: Vec::new(),
            max_matches: 1_000_000,
            symmetrize_conjugated_terminal_groups: true,
            weights: Vec::new(),
        }
    );
    assert_eq!(
        GetAlignmentTransformParams::default(),
        GetAlignmentTransformParams {
            prb_cid: -1,
            ref_cid: -1,
            atom_map: Vec::new(),
            weights: Vec::new(),
            reflect: false,
            max_iters: 50,
        }
    );
    assert_eq!(
        CalcRMSParams::default(),
        CalcRMSParams {
            prb_id: -1,
            ref_id: -1,
            map: Vec::new(),
            max_matches: 1_000_000,
            symmetrize_conjugated_terminal_groups: true,
            weights: Vec::new(),
        }
    );
    assert_eq!(
        AlignMolConformersParams::default(),
        AlignMolConformersParams {
            atom_ids: Vec::new(),
            conf_ids: Vec::new(),
            weights: Vec::new(),
            reflect: false,
            max_iters: 50,
        }
    );
    assert_eq!(
        RandomTransformParams::default(),
        RandomTransformParams { cid: -1, seed: -1 }
    );
    let crippen = GetCrippenO3AParams::default();
    assert_eq!(
        crippen,
        GetCrippenO3AParams {
            prb_crippen_contribs: None,
            ref_crippen_contribs: None,
            prb_cid: -1,
            ref_cid: -1,
            reflect: false,
            max_iters: 50,
            options: 0,
            constraint_map: Vec::new(),
            constraint_weights: Vec::new(),
        }
    );
    assert_eq!(
        GetCrippenO3AForProbeConfsParams::default(),
        GetCrippenO3AForProbeConfsParams {
            num_threads: 1,
            prb_crippen_contribs: None,
            ref_crippen_contribs: None,
            ref_cid: -1,
            reflect: false,
            max_iters: 50,
            options: 0,
            constraint_map: Vec::new(),
            constraint_weights: Vec::new(),
        }
    );
    let mmff = GetO3AParams::default();
    assert!(mmff.prb_py_mmff_mol_properties.is_none() && mmff.ref_py_mmff_mol_properties.is_none());
    assert_eq!(
        (
            mmff.prb_cid,
            mmff.ref_cid,
            mmff.reflect,
            mmff.max_iters,
            mmff.options
        ),
        (-1, -1, false, 50, 0)
    );
    assert!(mmff.constraint_map.is_empty() && mmff.constraint_weights.is_empty());
    let confs = GetO3AForProbeConfsParams::default();
    assert_eq!(
        (
            confs.num_threads,
            confs.ref_cid,
            confs.reflect,
            confs.max_iters,
            confs.options
        ),
        (1, -1, false, 50, 0)
    );
    assert!(confs.prb_py_mmff_mol_properties.is_none() && confs.constraint_map.is_empty());
}

#[test]
fn alignment_transforms() {
    let (rmsd, trans) = get_alignment_transform(&probe(1), &reference()).unwrap();
    assert_close(rmsd, 0.8928731745517915);
    assert_row(
        &trans[0],
        [
            0.9986680328265547,
            -0.01706977523750763,
            0.0486906868269532,
            -0.0045580158124188305,
        ],
    );
    assert_row(
        &trans[2],
        [
            0.04466418756147247,
            0.758462821884721,
            -0.650184018696353,
            0.060797134092237304,
        ],
    );
    assert_eq!(trans[3], vec![0.0, 0.0, 0.0, 1.0]);
    let mut weights = vec![1.0; 9];
    weights.push(2.0);
    let params = GetAlignmentTransformParams {
        atom_map: heavy_map(),
        weights,
        reflect: true,
        max_iters: 30,
        ..Default::default()
    };
    let (rmsd, trans) =
        get_alignment_transform_with_params(&probe(1), &reference(), &params).unwrap();
    assert_close(rmsd, 0.7051496490696121);
    assert_row(
        &trans[0],
        [
            0.9773228804206964,
            -0.03763339338981556,
            -0.20838357686766187,
            -0.02209662297287096,
        ],
    );

    let mut prb = probe(1);
    assert_close(
        align_mol(&mut prb, &reference()).unwrap(),
        0.8928731745517915,
    );
    assert_point(
        first_position(&prb, 0),
        (
            3.5123050060894276,
            -0.018911733348943277,
            0.40875949473908096,
        ),
    );
    let mapped = GetAlignmentTransformParams {
        atom_map: heavy_map(),
        ..Default::default()
    };
    let mut prb = probe(1);
    assert_close(
        align_mol_with_params(&mut prb, &reference(), &mapped).unwrap(),
        0.6937072271379565,
    );

    // Python raises for these too.
    let wrong = GetAlignmentTransformParams {
        weights: vec![1.0; 3],
        ..Default::default()
    };
    assert!(value(get_alignment_transform_with_params(
        &probe(1),
        &reference(),
        &wrong
    )));
    let far = GetAlignmentTransformParams {
        atom_map: vec![(0, 19)],
        ..Default::default()
    };
    assert_eq!(
        align_mol_with_params(&mut probe(1), &reference(), &far),
        Err(Error::Index { index: 19, len: 19 })
    );
    let benzene = add_hs(&mol_from_smiles("c1ccccc1").unwrap()).unwrap();
    assert!(matches!(
        get_alignment_transform(&benzene, &reference()),
        Err(Error::Runtime { .. }) | Err(Error::Value { .. })
    ));
}

#[test]
fn best_alignments() {
    let (rmsd, trans, matched) = get_best_alignment_transform(&probe(1), &reference()).unwrap();
    assert_close(rmsd, 0.8029821548022692);
    assert_eq!(matched.len(), 19);
    assert_eq!(&matched[..4], &[(0, 0), (1, 1), (2, 2), (3, 3)]);
    assert_row(
        &trans[0],
        [
            0.9636550259883816,
            0.26145889313280746,
            -0.05484740731427951,
            0.005121690707780779,
        ],
    );
    // Any number of threads gives the same alignment.
    for num_threads in [1, 4, 0, i32::MIN] {
        let params = GetBestAlignmentTransformParams {
            num_threads,
            ..Default::default()
        };
        let (rmsd, _, matched) =
            get_best_alignment_transform_with_params(&probe(1), &reference(), &params).unwrap();
        assert_close(rmsd, 0.8029821548022692);
        assert_eq!(matched.len(), 19);
    }

    let mut prb = probe(1);
    assert_close(
        get_best_rms(&mut prb, &reference()).unwrap(),
        0.8029821548022692,
    );
    assert_point(
        first_position(&prb, 0),
        (3.4702794318494528, -0.12928838926107378, 0.666465338082832),
    );
    let mut prb = probe(1);
    assert_close(
        get_best_rms_with_params(&mut prb, &reference(), &Default::default()).unwrap(),
        0.8029821548022692,
    );
    let mapped = GetBestRMSParams {
        map: vec![heavy_map()],
        weights: vec![1.0; 10],
        ..Default::default()
    };
    assert_close(
        get_best_rms_with_params(&mut probe(1), &reference(), &mapped).unwrap(),
        0.6937072271379565,
    );
    let plain = GetBestRMSParams {
        symmetrize_conjugated_terminal_groups: false,
        max_matches: 10,
        ..Default::default()
    };
    assert_close(
        get_best_rms_with_params(&mut probe(1), &reference(), &plain).unwrap(),
        0.8029821548022692,
    );

    assert_close(
        calc_rms(&probe(1), &reference()).unwrap(),
        1.7092608408622902,
    );
    // One molecule can be both, read-only.
    let mol = acetanilide();
    let ids = CalcRMSParams {
        prb_id: 1,
        ref_id: 0,
        ..Default::default()
    };
    assert_close(
        calc_rms_with_params(&mol, &mol, &ids).unwrap(),
        1.7092608408622902,
    );
}

#[test]
fn all_conformers() {
    let mol = acetanilide();
    let expected = [0.8029821548022692, 0.35494074979983326, 0.7825745565387795];
    for num_threads in [1, 2] {
        let params = GetAllConformerBestRMSParams {
            num_threads,
            ..Default::default()
        };
        let rms = get_all_conformer_best_rms_with_params(&mol, &params).unwrap();
        assert_eq!(rms.len(), 3);
        for (a, e) in rms.iter().zip(expected) {
            assert_close(*a, e);
        }
    }
    assert_eq!(get_all_conformer_best_rms(&mol).unwrap().len(), 3);

    let mut mol = acetanilide();
    let rms = align_mol_conformers(&mut mol).unwrap();
    assert_eq!(rms.len(), 2);
    assert_close(rms[0], 0.8928731745517915);
    assert_close(rms[1], 0.49846097534395084);
    assert_point(
        first_position(&mol, 2),
        (3.390257128818704, 0.026768917387915044, 0.619296911216853),
    );
    let mut mol = acetanilide();
    let params = AlignMolConformersParams {
        atom_ids: (0..10).collect(),
        conf_ids: vec![2, 0],
        reflect: true,
        ..Default::default()
    };
    let rms = align_mol_conformers_with_params(&mut mol, &params).unwrap();
    assert_eq!(rms.len(), 1);
    assert_close(rms[0], 1.1342845356735207);
    let bad = AlignMolConformersParams {
        atom_ids: vec![19],
        ..Default::default()
    };
    assert_eq!(
        align_mol_conformers_with_params(&mut mol, &bad),
        Err(Error::Index { index: 19, len: 19 })
    );
    let missing = AlignMolConformersParams {
        conf_ids: vec![0, 7],
        ..Default::default()
    };
    assert!(value(align_mol_conformers_with_params(&mut mol, &missing)));
    assert_eq!(
        align_mol_conformers(&mut Mol::default()).unwrap(),
        Vec::<f64>::new()
    );
}

#[test]
fn random_transforms() {
    let mut mol = probe(1);
    let seeded = RandomTransformParams {
        seed: 42,
        ..Default::default()
    };
    random_transform_with_params(&mut mol, &seeded).unwrap();
    let pos = first_position(&mol, 0);
    for (a, e) in [
        (pos.x, -23.876328555449028),
        (pos.y, -5.256430404640382),
        (pos.z, 0.36720733970340547),
    ] {
        assert_close_to(a, e, 1e-6);
    }
    // The shape is unchanged: it aligns back.
    let mut again = probe(1);
    random_transform(&mut again).unwrap();
    assert!(calc_rms(&again, &probe(1)).unwrap() > 1e-3);
    assert!(align_mol(&mut again, &probe(1)).unwrap() < 1e-6);
    // RDKit reads the bounds of a conformer without atoms uninitialized.
    let mut empty = Mol::default();
    empty.add_conformer(&Conformer::new(0)).unwrap();
    random_transform(&mut empty).unwrap();
    assert!(value(random_transform(&mut reference_without_conformer())));
}

fn reference_without_conformer() -> Mol {
    mol_from_smiles("CC(=O)Nc1ccccc1").unwrap()
}

#[test]
fn best_alignment_inputs_rdkit_trusts() {
    // On RDKit's threads, Python calls std::terminate for these.
    for num_threads in [1, 2] {
        let far = GetBestRMSParams {
            map: vec![vec![(0, 100)]],
            num_threads,
            ..Default::default()
        };
        assert_eq!(
            get_best_rms_with_params(&mut probe(1), &reference(), &far),
            Err(Error::Index {
                index: 100,
                len: 19
            })
        );
        let short = GetAllConformerBestRMSParams {
            weights: vec![1.0],
            num_threads,
            ..Default::default()
        };
        assert!(value(get_all_conformer_best_rms_with_params(
            &acetanilide(),
            &short
        )));
        let short = GetBestAlignmentTransformParams {
            weights: vec![1.0],
            num_threads,
            ..Default::default()
        };
        assert!(value(get_best_alignment_transform_with_params(
            &probe(1),
            &reference(),
            &short
        )));
    }
    let uneven = CalcRMSParams {
        map: vec![heavy_map(), vec![(0, 0)]],
        weights: vec![1.0; 10],
        ..Default::default()
    };
    assert!(value(calc_rms_with_params(
        &probe(1),
        &reference(),
        &uneven
    )));
}

#[test]
fn o3a() {
    let reference = reference();
    let mut prb = probe(1);
    let mut o3a = get_crippen_o3a(&mut prb, &reference).unwrap();
    assert_close_to(o3a.score(), 53.131553536125054, 0.05);
    assert_eq!(o3a.matches().len(), 9);
    assert_eq!(o3a.weights().len(), 9);
    let (rmsd, trans) = o3a.trans().unwrap();
    assert_close_to(rmsd, 0.44116123180897954, 0.05);
    assert_eq!(trans.len(), 4);
    assert_close(o3a.align().unwrap(), rmsd);
    drop(o3a);
    // Aligning moved the probe.
    assert!(first_position(&prb, 0).distance(&first_position(&probe(1), 0)) > 1e-3);

    let mut prb = probe(1);
    let mut o3a = get_o3a(&mut prb, &reference).unwrap();
    assert_close_to(o3a.score(), 67.42485523915751, 0.05);
    assert_eq!(o3a.matches().len(), 9);
    assert_close_to(o3a.align().unwrap(), 0.41648855966236875, 0.05);

    // Given properties, as Python computes them.
    let mut prb = probe(1);
    let prb_props = MMFFMolProperties::new(&mut prb.clone()).unwrap();
    let ref_props = MMFFMolProperties::new(&mut reference.clone()).unwrap();
    let params = GetO3AParams {
        prb_py_mmff_mol_properties: Some(&prb_props),
        ref_py_mmff_mol_properties: Some(&ref_props),
        ..Default::default()
    };
    let o3a = get_o3a_with_params(&mut prb, &reference, &params).unwrap();
    assert_close_to(o3a.score(), 67.42485523915751, 0.05);
    drop(o3a);

    let mut prb = probe(1);
    let constrained = GetO3AParams {
        options: 4,
        constraint_map: vec![(0, 0), (1, 1)],
        constraint_weights: vec![50.0, 20.0],
        ..Default::default()
    };
    let o3a = get_o3a_with_params(&mut prb, &reference, &constrained).unwrap();
    assert_close_to(o3a.score(), 19.71178117885854, 0.05);
    assert_eq!(&o3a.matches()[..2], &[(0, 0), (1, 1)]);
    drop(o3a);

    let logp = [
        0.1441, -0.2783, -0.1526, -0.5188, 0.4619, 0.1581, 0.1581, 0.1581, 0.1581, 0.1581, 0.123,
        0.123, 0.123, 0.2142, 0.123, 0.123, 0.123, 0.123, 0.123,
    ];
    let contribs: Vec<(f64, f64)> = logp.iter().map(|&l| (l, 0.0)).collect();
    let mut prb = probe(1);
    let given = GetCrippenO3AParams {
        prb_crippen_contribs: Some(contribs.clone()),
        ref_crippen_contribs: Some(contribs),
        reflect: true,
        max_iters: 20,
        ..Default::default()
    };
    let o3a = get_crippen_o3a_with_params(&mut prb, &reference, &given).unwrap();
    assert_close_to(o3a.score(), 53.131553536125374, 0.05);
    drop(o3a);
    let mut prb = probe(1);
    let flat = GetCrippenO3AParams {
        prb_crippen_contribs: Some(vec![(0.5, 0.0); 19]),
        ref_crippen_contribs: Some(vec![(0.5, 0.0); 19]),
        ..Default::default()
    };
    let o3a = get_crippen_o3a_with_params(&mut prb, &reference, &flat).unwrap();
    assert_close_to(o3a.score(), 58.03295582788956, 0.05);
    drop(o3a);
    // Other than one per atom, as in Python, computes them.
    let mut prb = probe(1);
    let short = GetCrippenO3AParams {
        prb_crippen_contribs: Some(vec![(0.5, 0.0); 3]),
        ..Default::default()
    };
    let o3a = get_crippen_o3a_with_params(&mut prb, &reference, &short).unwrap();
    assert_close_to(o3a.score(), 53.131553536125054, 0.05);
}

#[test]
fn o3a_for_probe_conformers() {
    let reference = reference();
    let mut mol = acetanilide();
    let o3as = get_o3a_for_probe_confs(&mut mol, &reference).unwrap();
    let expected = [86.14, 67.42485523915751, 85.8189276652746];
    assert_eq!(o3as.len(), 3);
    for (o3a, e) in o3as.iter().zip(expected) {
        assert_close_to(o3a.score(), e, 0.05);
    }
    drop(o3as);

    let mut mol = acetanilide();
    let params = GetCrippenO3AForProbeConfsParams {
        num_threads: 2,
        ..Default::default()
    };
    let mut o3as =
        get_crippen_o3a_for_probe_confs_with_params(&mut mol, &reference, &params).unwrap();
    let scores = [64.6924, 53.131553536125054, 64.41494314002524];
    let rmsds = [0.0, 0.44116123180897954, 0.09336403231105261];
    for ((o3a, s), r) in o3as.iter_mut().zip(scores).zip(rmsds) {
        assert_close_to(o3a.score(), s, 0.05);
        let rmsd = o3a.align().unwrap();
        assert!((rmsd - r).abs() < 0.05, "{rmsd} != {r}");
    }
    // The first conformer is the reference's: every heavy atom matches
    // itself.
    let mut matched = o3as[0].matches();
    matched.sort();
    assert_eq!(matched, heavy_map());
    drop(o3as);
    assert_eq!(
        get_crippen_o3a_for_probe_confs(&mut mol, &reference)
            .unwrap()
            .len(),
        3
    );
}

#[test]
fn o3a_inputs_rdkit_trusts() {
    let reference = reference();
    // O3A counts dummy atoms as heavy in some places and not in others, and
    // writes past its tables; Python's heap breaks.
    let mut dummies = probe(1);
    for i in 10..19 {
        dummies.get_atom_with_idx_mut(i).unwrap().set_atomic_num(0);
    }
    assert!(value(get_crippen_o3a(&mut dummies, &reference)));
    // MolHistogram casts distances to unsigned int.
    for x in [f64::NAN, 1e12] {
        let mut far = probe(1);
        far.get_conformer_mut()
            .unwrap()
            .set_atom_position(3, Point3D::new(x, 0.0, 0.0))
            .unwrap();
        assert!(value(get_crippen_o3a(&mut far, &reference)));
        assert!(value(get_o3a_for_probe_confs(&mut far, &reference)));
    }
    // The cost function reads MMFF types by atom, and indexes a matrix by
    // type minus 1.
    let mut ethanol = add_hs(&mol_from_smiles("CCO").unwrap()).unwrap();
    let other = MMFFMolProperties::new(&mut ethanol).unwrap();
    let params = GetO3AParams {
        prb_py_mmff_mol_properties: Some(&other),
        ..Default::default()
    };
    assert!(value(get_o3a_with_params(
        &mut probe(1),
        &reference,
        &params
    )));
    // Constraints.
    let on_h = GetO3AParams {
        constraint_map: vec![(10, 10)],
        ..Default::default()
    };
    assert!(value(get_o3a_with_params(&mut probe(1), &reference, &on_h)));
    let far = GetCrippenO3AParams {
        constraint_map: vec![(0, 40)],
        ..Default::default()
    };
    assert_eq!(
        get_crippen_o3a_with_params(&mut probe(1), &reference, &far).err(),
        Some(Error::Index { index: 40, len: 19 })
    );
    for weights in [vec![1.0], vec![f64::NAN, 1.0]] {
        let params = GetCrippenO3AParams {
            constraint_map: vec![(0, 0), (1, 1)],
            constraint_weights: weights,
            ..Default::default()
        };
        assert!(value(get_crippen_o3a_with_params(
            &mut probe(1),
            &reference,
            &params
        )));
    }
    // The assignment solver sums costs made from the contributions in int.
    let huge = GetCrippenO3AParams {
        prb_crippen_contribs: Some(vec![(1e300, 0.0); 19]),
        ..Default::default()
    };
    assert!(value(get_crippen_o3a_with_params(
        &mut probe(1),
        &reference,
        &huge
    )));
    // A missing conformer.
    let missing = GetCrippenO3AParams {
        prb_cid: 3,
        ..Default::default()
    };
    assert!(value(get_crippen_o3a_with_params(
        &mut probe(1),
        &reference,
        &missing
    )));
}

#[test]
fn concurrent_alignments() {
    // Finding the matches symmetrizes terminal groups through a static RDKit
    // variable, which the embedder's pruning writes too.
    let handles: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                for _ in 0..3 {
                    let mut prb = probe(1);
                    assert_close(
                        get_best_rms_with_params(&mut prb, &reference(), &Default::default())
                            .unwrap(),
                        0.8029821548022692,
                    );
                    let mut mol = probe(1);
                    random_transform(&mut mol).unwrap();
                }
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
}
