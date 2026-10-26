// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_dist_geom::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values below come from Python's rdDistGeom. Embedded coordinates
// differ between platforms even with a seed, so the tests check what any
// embedding satisfies (bond lengths, bounds) and that a seed repeats on the
// platform at hand, not Python's coordinates.

fn with_hs(smiles: &str) -> Mol {
    add_hs(&mol_from_smiles(smiles).unwrap()).unwrap()
}

fn position(mol: &Mol, conf_id: i32, atom: u32) -> Point3D {
    mol.get_conformer_with_id(conf_id)
        .unwrap()
        .get_atom_position(atom)
        .unwrap()
}

fn positions(mol: &Mol, conf_id: i32) -> Vec<Point3D> {
    mol.get_conformer_with_id(conf_id).unwrap().get_positions()
}

/// Asserts that each bond of the conformer is as long as a bond between C,
/// O, and H atoms can be.
fn assert_sensible_bonds(mol: &Mol, conf_id: i32) {
    for bond in mol.get_bonds() {
        let length = position(mol, conf_id, bond.get_begin_atom_idx()).distance(&position(
            mol,
            conf_id,
            bond.get_end_atom_idx(),
        ));
        assert!((0.9..1.7).contains(&length), "bond length {length}");
    }
}

fn assert_same_positions(a: &[Point3D], b: &[Point3D]) {
    assert_eq!(a.len(), b.len());
    for (p, q) in a.iter().zip(b) {
        assert!(p.distance(q) < 1e-9, "{p:?} != {q:?}");
    }
}

/// Python's keyword defaults for `EmbedMolecule`.
fn keyword_params(random_seed: i32) -> EmbedParameters {
    EmbedParameters {
        use_exp_torsion_angle_prefs: true,
        use_basic_knowledge: true,
        use_macrocycle_torsions: true,
        use_macrocycle14config: true,
        random_seed,
        ..Default::default()
    }
}

#[test]
fn default_parameters() {
    let expected = EmbedParameters {
        max_iterations: 0,
        num_threads: 1,
        timeout: 0,
        random_seed: -1,
        clear_confs: true,
        use_random_coords: false,
        box_size_mult: 2.0,
        rand_neg_eig: true,
        num_zero_fail: 1,
        optimizer_force_tol: 1e-3,
        basin_thresh: 5.0,
        ignore_smoothing_failures: false,
        enforce_chirality: true,
        use_exp_torsion_angle_prefs: false,
        use_basic_knowledge: false,
        et_version: 2,
        verbose: false,
        prune_rms_thresh: -1.0,
        only_heavy_atoms_for_rms: true,
        embed_fragments_separately: true,
        use_small_ring_torsions: false,
        use_macrocycle_torsions: false,
        use_macrocycle14config: false,
        only_initial_embedding: false,
        use_legacy_implementation: true,
        initial_embedding_mode: InitialEmbeddingMode::DgEmbedding,
        bounds_mat_force_scaling: 1.0,
        use_symmetry_for_pruning: true,
        force_trans_amides: true,
        track_failures: false,
        enable_sequential_random_seeds: false,
        symmetrize_conjugated_terminal_groups_for_pruning: true,
        embed_force_field: EmbedFF::Uff,
        coord_map: BTreeMap::new(),
        cpci: BTreeMap::new(),
        bounds_mat: None,
        failures: Vec::new(),
    };
    assert_eq!(EmbedParameters::default(), expected);
    assert_eq!(
        GetMoleculeBoundsMatrixParams::default(),
        GetMoleculeBoundsMatrixParams {
            set15bounds: true,
            scale_vdw: false,
            do_triangle_smoothing: true,
            use_macrocycle14config: false,
            force_trans_amides: true,
            set14bounds: true,
            set13bounds: true,
        }
    );
    assert_eq!(
        GetMoleculeBoundsMatrixEmbedParametersParams::default(),
        GetMoleculeBoundsMatrixEmbedParametersParams {
            do_triangle_smoothing: true,
            scale_vdw: false,
            set15bounds: true,
            set14bounds: true,
            set13bounds: true,
        }
    );
}

#[test]
fn presets() {
    // (use_exp_torsion_angle_prefs, use_basic_knowledge, et_version,
    // use_small_ring_torsions, use_macrocycle_torsions, use_legacy_implementation)
    let summary = |p: EmbedParameters| {
        (
            p.use_exp_torsion_angle_prefs,
            p.use_basic_knowledge,
            p.et_version,
            p.use_small_ring_torsions,
            p.use_macrocycle_torsions,
            p.use_macrocycle14config,
            p.use_legacy_implementation,
        )
    };
    assert_eq!(summary(etkdg()), (true, true, 1, false, false, false, true));
    assert_eq!(
        summary(etkdg_v2()),
        (true, true, 2, false, false, false, true)
    );
    assert_eq!(
        summary(etkdg_v3()),
        (true, true, 2, false, true, true, true)
    );
    assert_eq!(
        summary(sr_etkdg_v3()),
        (true, true, 2, true, false, false, true)
    );
    assert_eq!(
        summary(etkdg_v4()),
        (true, true, 4, true, true, true, false)
    );
    assert_eq!(summary(etdg()), (true, false, 1, false, false, false, true));
    assert_eq!(
        summary(etdg_v2()),
        (true, false, 2, false, false, false, true)
    );
    assert_eq!(summary(kdg()), (false, true, 1, false, false, false, true));
    assert_eq!(summary(dg()), (false, false, 2, false, false, false, true));
    assert_eq!(
        embed_parameters_to_json(&etkdg_v3()).unwrap(),
        concat!(
            r#"{"basinThresh":"5","boundsMatForceScaling":"1","boxSizeMult":"2","#,
            r#""clearConfs":"true","embedForceField":"UFF","embedFragmentsSeparately":"true","#,
            r#""enableSequentialRandomSeeds":"false","enforceChirality":"true","ETversion":"2","#,
            r#""forceTransAmides":"true","ignoreSmoothingFailures":"false","maxIterations":"0","#,
            r#""numThreads":"1","numZeroFail":"1","onlyHeavyAtomsForRMS":"true","#,
            r#""optimizerForceTol":"0.001","pruneRmsThresh":"-1","randNegEig":"true","#,
            r#""randomSeed":"-1","symmetrizeConjugatedTerminalGroupsForPruning":"true","#,
            r#""timeout":"0","trackFailures":"false","useBasicKnowledge":"true","#,
            r#""useExpTorsionAnglePrefs":"true","useLegacyImplementation":"true","#,
            r#""useMacrocycle14config":"true","useMacrocycleTorsions":"true","#,
            r#""useRandomCoords":"false","useSmallRingTorsions":"false","#,
            r#""useSymmetryForPruning":"true","verbose":"false","onlyInitialEmbedding":"false","#,
            r#""initialEmbeddingMode":"DG_EMBEDDING"}"#
        )
    );
}

#[test]
fn embed_molecules() {
    let seeded = |params: &mut EmbedParameters| {
        let mut mol = with_hs("CCO");
        assert_eq!(embed_molecule_with_params(&mut mol, params).unwrap(), 0);
        assert_sensible_bonds(&mol, 0);
        positions(&mol, 0)
    };
    let first = seeded(&mut keyword_params(42));
    assert_same_positions(&first, &seeded(&mut keyword_params(42)));
    assert!(first[0].distance(&seeded(&mut keyword_params(7))[0]) > 1e-3);

    let mut mol = with_hs("CCO");
    let mut params = EmbedParameters {
        only_heavy_atoms_for_rms: false,
        ..keyword_params(42)
    };
    assert_eq!(
        embed_multiple_confs_with_params(&mut mol, 3, &mut params).unwrap(),
        vec![0, 1, 2]
    );
    for conf_id in 0..3 {
        assert_sensible_bonds(&mol, conf_id);
    }
    // The first conformer is the one EmbedMolecule makes.
    assert_same_positions(&positions(&mol, 0), &first);

    // Python's four threads give the conformers one thread gives; here the
    // calling thread makes them, and num_threads is left as it was.
    let mut by_threads = Vec::new();
    for num_threads in [1, 4, 0] {
        let mut mol = with_hs("CCO");
        let mut params = EmbedParameters {
            random_seed: 0xf00d,
            num_threads,
            ..etkdg_v3()
        };
        assert_eq!(
            embed_multiple_confs_with_params(&mut mol, 4, &mut params).unwrap(),
            vec![0, 1, 2, 3]
        );
        assert_eq!(params.num_threads, num_threads);
        by_threads.push(positions(&mol, 3));
    }
    assert_same_positions(&by_threads[0], &by_threads[1]);
    assert_same_positions(&by_threads[0], &by_threads[2]);

    // Every chair of cyclohexane is within 0.5 of the others.
    let mut mol = with_hs("C1CCCCC1");
    let mut params = EmbedParameters {
        random_seed: 42,
        prune_rms_thresh: 0.5,
        ..etkdg_v3()
    };
    let ids = embed_multiple_confs_with_params(&mut mol, 10, &mut params).unwrap();
    assert!(!ids.is_empty() && ids.len() < 10, "{ids:?}");

    // The keyword forms, unseeded.
    let mut mol = with_hs("CCO");
    assert_eq!(embed_molecule(&mut mol).unwrap(), 0);
    assert_eq!(mol.get_num_conformers(), 1);
    assert_eq!(embed_multiple_confs(&mut mol).unwrap().len(), 10);
    assert_eq!(mol.get_num_conformers(), 10);
}

#[test]
fn embedding_updates_the_parameters() {
    let mut mol = with_hs("CCO");
    let mut params = EmbedParameters {
        random_seed: 42,
        track_failures: true,
        use_random_coords: true,
        ..etkdg_v3()
    };
    assert_eq!(params.get_failure_counts(), Vec::<u32>::new());
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    assert_eq!(params.max_iterations, 90);
    assert_eq!(params.basin_thresh, 1e8);
    assert_eq!(
        params.initial_embedding_mode,
        InitialEmbeddingMode::RandomCoordinateEmbedding
    );
    assert_eq!(params.get_failure_counts().len(), 15);

    // maxIterations comes from the first fragment.
    let mut mol = with_hs("CC.O");
    let mut params = EmbedParameters {
        random_seed: 7,
        track_failures: true,
        ..Default::default()
    };
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    assert_eq!(params.max_iterations, 80);
}

#[test]
fn coordinate_maps_bounds_and_cpci() {
    let mut mol = with_hs("CCO");
    let mut params = EmbedParameters {
        random_seed: 42,
        ..etkdg_v3()
    };
    params.set_coord_map(&BTreeMap::from([
        (0, Point3D::new(0.0, 0.0, 0.0)),
        (2, Point3D::new(2.5, 0.0, 0.0)),
    ]));
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    // Python's is 2.4916 on x86-64.
    let distance = position(&mol, 0, 0).distance(&position(&mol, 0, 2));
    assert!((distance - 2.5).abs() < 0.05, "{distance}");
    // Python raises RDKit's ValueError for this.
    params.set_coord_map(&BTreeMap::from([(20, Point3D::new(0.0, 0.0, 0.0))]));
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Index { index: 20, len: 9 })
    );
    params.set_coord_map(&BTreeMap::from([(0, Point3D::new(f64::NAN, 0.0, 0.0))]));
    assert!(matches!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Value { .. })
    ));

    let mut mol = mol_from_smiles("CCO").unwrap();
    let bounds = get_molecule_bounds_matrix(&mol).unwrap();
    let mut params = EmbedParameters {
        random_seed: 42,
        ..etkdg_v3()
    };
    params.set_bounds_mat(&bounds).unwrap();
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    for i in 0..3 {
        for j in i + 1..3 {
            let d = position(&mol, 0, i).distance(&position(&mol, 0, j));
            let (i, j) = (i as usize, j as usize);
            assert!(bounds[j][i] - 0.1 < d && d < bounds[i][j] + 0.1, "{d}");
        }
    }
    assert!(matches!(
        params.set_bounds_mat(&[vec![0.0, 1.0, 2.0], vec![1.0, 0.0, 1.0]]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        params.set_bounds_mat(&[]),
        Err(Error::Value { .. })
    ));
    // RDKit checks the size against the molecule's.
    params.bounds_mat = Some(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
    assert!(matches!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Value { .. })
    ));
    params.bounds_mat = Some(vec![vec![0.0, 1.0], vec![1.0]]);
    assert!(matches!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Value { .. })
    ));

    let mut mol = with_hs("CCO");
    let mut params = EmbedParameters {
        random_seed: 42,
        ..Default::default()
    };
    params.set_cpci(&BTreeMap::from([((0, 2), 0.5)]));
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    assert_sensible_bonds(&mol, 0);
    params.set_cpci(&BTreeMap::from([((0, 20), 0.5)]));
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Index { index: 20, len: 9 })
    );
}

#[test]
fn embedding_inputs_rdkit_trusts() {
    assert!(matches!(
        embed_molecule(&mut Mol::default()),
        Err(Error::Value { .. })
    ));

    // computeInitialCoords casts the summed squares of the distances it
    // picks to int. With bounds this long, Python's seeded embedding seeds
    // its eigensolver from clock() instead, and differs from run to run.
    let n = 12;
    let bounds: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| match i.cmp(&j) {
                    std::cmp::Ordering::Less => 1e6,
                    std::cmp::Ordering::Greater => 1.0,
                    std::cmp::Ordering::Equal => 0.0,
                })
                .collect()
        })
        .collect();
    let mut mol = mol_from_smiles(&"C".repeat(n)).unwrap();
    let mut params = EmbedParameters {
        random_seed: 42,
        only_initial_embedding: true,
        ..Default::default()
    };
    params.set_bounds_mat(&bounds).unwrap();
    assert!(matches!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Value { .. })
    ));
    // Random starting coordinates make no such sum.
    params.use_random_coords = true;
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    // Nor do short bounds.
    let short: Vec<Vec<f64>> = bounds
        .iter()
        .map(|row| row.iter().map(|&b| b.min(10.0)).collect())
        .collect();
    params.use_random_coords = false;
    params.set_bounds_mat(&short).unwrap();
    assert_eq!(
        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
        0
    );
    let mut infinite = short.clone();
    infinite[0][1] = f64::INFINITY;
    params.bounds_mat = Some(infinite);
    assert!(matches!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Value { .. })
    ));
    // A long coordinate map makes long bounds too.
    let mut mol = with_hs("CCO");
    let mut params = EmbedParameters {
        random_seed: 42,
        ..Default::default()
    };
    params.set_coord_map(&BTreeMap::from([
        (0, Point3D::new(0.0, 0.0, 0.0)),
        (2, Point3D::new(1e9, 0.0, 0.0)),
    ]));
    assert!(matches!(
        embed_molecule_with_params(&mut mol, &mut params),
        Err(Error::Value { .. })
    ));

    // The minimizers work in four dimensions, past 65535 coordinates here.
    let mut big = mol_from_smiles(&"C".repeat(16384)).unwrap();
    assert!(matches!(embed_molecule(&mut big), Err(Error::Value { .. })));
}

#[test]
fn concurrent_embeddings() {
    // Every embedding swaps the process's SIGINT handler through RDKit's
    // static variables, and unseeded ones draw from RDKit's global random
    // generator; the wrapper runs one at a time.
    let handles: Vec<_> = (0..4)
        .map(|i| {
            std::thread::spawn(move || {
                for _ in 0..3 {
                    let mut mol = with_hs("CCCCO");
                    let mut params = EmbedParameters {
                        random_seed: if i % 2 == 0 { -1 } else { i },
                        ..etkdg_v3()
                    };
                    assert_eq!(
                        embed_molecule_with_params(&mut mol, &mut params).unwrap(),
                        0
                    );
                }
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
}

#[test]
fn bounds_matrices() {
    let bounds = get_molecule_bounds_matrix(&mol_from_smiles("CCO").unwrap()).unwrap();
    let expected = [
        [0.0, 1.524, 2.4211],
        [1.504, 0.0, 1.4038],
        [2.3299, 1.3838, 0.0],
    ];
    assert_eq!(bounds.len(), 3);
    for (row, expected) in bounds.iter().zip(expected) {
        for (a, e) in row.iter().zip(expected) {
            assert!((a - e).abs() < 1e-4, "{bounds:?}");
        }
    }
    let butane = mol_from_smiles("CCCC").unwrap();
    let unsmoothed = get_molecule_bounds_matrix_with_params(
        &butane,
        &GetMoleculeBoundsMatrixParams {
            do_triangle_smoothing: false,
            set15bounds: false,
            ..Default::default()
        },
    )
    .unwrap();
    let with_params = get_molecule_bounds_matrix_embed_parameters(&butane, &etkdg_v3()).unwrap();
    for row in [&unsmoothed[0], &with_params[0]] {
        for (a, e) in row.iter().zip([0.0, 1.524, 2.5195, 3.8754]) {
            assert!((a - e).abs() < 1e-4, "{row:?}");
        }
    }
    assert_eq!(
        get_molecule_bounds_matrix_embed_parameters_with_params(
            &butane,
            &etkdg_v3(),
            &GetMoleculeBoundsMatrixEmbedParametersParams {
                do_triangle_smoothing: false,
                ..Default::default()
            },
        )
        .unwrap()
        .len(),
        4
    );
    assert!(matches!(
        get_molecule_bounds_matrix(&Mol::default()),
        Err(Error::Value { .. })
    ));
}

#[test]
#[allow(clippy::approx_constant)] // RDKit's torsion table rounds pi to 3.142
fn experimental_torsions() {
    let mol = with_hs("CCOC");
    let torsions = get_experimental_torsions(&mol).unwrap();
    assert_eq!(torsions.len(), 2);
    let torsion = &torsions[0];
    assert_eq!(torsion.bond_index, 1);
    assert_eq!(torsion.torsion_index, 9);
    assert_eq!(torsion.smarts, "[C:1][CH2:2]!@;-[O:3][CX4:4]");
    assert_eq!(torsion.atom_indices, [0, 1, 2, 3]);
    assert_eq!(
        torsion.terms,
        TorsionTerms::Cosine {
            v: vec![0.0, 0.0, 2.5, 0.0, 0.0, 0.0],
            signs: vec![1; 6],
        }
    );
    let torsions = get_experimental_torsions_with_params(&mol, &etkdg_v4()).unwrap();
    assert_eq!(torsions.len(), 2);
    assert_eq!(torsions[0].torsion_index, 175);
    assert_eq!(
        torsions[0].terms,
        TorsionTerms::Gaussian {
            positions: vec![1.396, 3.142],
            widths: vec![0.349, 0.278],
            heights: vec![0.084, 0.87],
        }
    );
    let v3 = EmbedParameters {
        et_version: 3,
        ..Default::default()
    };
    assert!(matches!(
        get_experimental_torsions_with_params(&mol, &v3),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        get_experimental_torsions(&Mol::default()),
        Err(Error::Value { .. })
    ));
}

#[test]
fn ordered_failure_causes() {
    use EmbedFailureCauses as E;
    assert_eq!(
        ordered_embed_failure_causes(),
        vec![
            (E::InitialCoords, 0),
            (E::FirstMinimization, 1),
            (E::CheckTetrahedralCenters, 2),
            (E::CheckChiralCenters, 3),
            (E::MinimizeFourthDimension, 4),
            (E::EtkMinimization, 5),
            (E::LinearDoubleBond, 8),
            (E::CheckChiralCenters2, 10),
            (E::FinalChiralBounds, 6),
            (E::FinalCenterInVolume, 7),
            (E::BadDoubleBondStereo, 9),
            (E::ExceededTimeout, 11),
        ]
    );
    assert_eq!(
        ordered_embed_failure_causes_with_legacy_implementation(false),
        vec![
            (E::InitialCoords, 0),
            (E::Minimization, 12),
            (E::KtermViolation, 13),
            (E::CheckTetrahedralCenters, 2),
            (E::CheckChiralCenters, 3),
            (E::FinalChiralBounds, 6),
            (E::FinalCenterInVolume, 7),
            (E::BadDoubleBondStereo, 9),
            (E::Clash, 14),
            (E::LinearDoubleBond, 9),
            (E::ExceededTimeout, 11),
        ]
    );
}
