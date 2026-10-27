// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_fingerprint_generator::*;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn morgan() -> FingerprintGenerator64 {
    get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 1,
        fp_size: 64,
        ..Default::default()
    })
    .unwrap()
}

// Expected values come from Python's rdFingerprintGenerator.

#[test]
fn morgan_fingerprints() {
    let generator = morgan();
    let mol = smiles("CCO");
    assert_eq!(
        generator.get_fingerprint(&mol).unwrap().get_on_bits(),
        vec![2, 16, 30, 33, 38, 39]
    );
    assert_eq!(
        generator
            .get_count_fingerprint(&mol)
            .unwrap()
            .get_nonzero_elements(),
        BTreeMap::from([(2, 1), (16, 1), (30, 1), (33, 1), (38, 1), (39, 1)])
    );
    let unfolded: [u32; 6] = [
        864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
    ];
    assert_eq!(
        generator
            .get_sparse_count_fingerprint(&mol)
            .unwrap()
            .get_nonzero_elements()
            .into_keys()
            .collect::<Vec<_>>(),
        unfolded.map(u64::from)
    );
    // Python shows the bits past i32::MAX as negative ints.
    let sparse = generator.get_sparse_fingerprint(&mol).unwrap();
    assert_eq!(sparse.get_num_bits(), u32::MAX);
    assert_eq!(sparse.get_on_bits(), unfolded);
    assert!(sparse.get_bit(4018048386).unwrap());
    let from_atoms = GetFingerprintParams {
        from_atoms: Some(vec![0]),
        ..Default::default()
    };
    assert_eq!(
        generator
            .get_fingerprint_with_params(&mol, &from_atoms)
            .unwrap()
            .get_on_bits(),
        vec![33, 38]
    );
    // An empty list is no list, as in Python.
    let empty = GetFingerprintParams {
        from_atoms: Some(vec![]),
        ..Default::default()
    };
    assert_eq!(
        generator
            .get_fingerprint_with_params(&mol, &empty)
            .unwrap()
            .get_num_on_bits(),
        6
    );
    let bits = generator.get_fingerprint_as_num_py(&mol).unwrap();
    assert_eq!(bits.len(), 64);
    assert_eq!(bits.iter().filter(|&&b| b == 1).count(), 6);
    let counts = generator.get_count_fingerprint_as_num_py(&mol).unwrap();
    assert_eq!(counts.iter().sum::<u32>(), 6);
}

#[test]
fn oversized_sparse_fingerprints() {
    // RDKit's pickle of a SparseBitVect subtracts its int-stored bits in int.
    let sparse = morgan().get_sparse_fingerprint(&smiles("CCO")).unwrap();
    assert!(matches!(sparse.to_binary(), Err(Error::Value { .. })));
    assert!(matches!(sparse.to_base64(), Err(Error::Value { .. })));
    assert!(matches!(
        rdkit_cxx::data_structs::convert_to_explicit(&sparse),
        Err(Error::Value { .. })
    ));
}

#[test]
fn additional_output() {
    let generator = morgan();
    let mut output = AdditionalOutput::default();
    output.allocate_bit_info_map();
    output.collect_atom_counts();
    output.allocate_atom_to_bits();
    let params = GetFingerprintParams {
        additional_output: Some(output),
        ..Default::default()
    };
    let (fp, output) = generator
        .get_fingerprint_with_outputs(&smiles("CCO"), &params)
        .unwrap();
    assert_eq!(fp.get_num_on_bits(), 6);
    assert_eq!(
        output.get_bit_info_map().unwrap(),
        &BTreeMap::from([
            (2, vec![(1, 1)]),
            (16, vec![(1, 0)]),
            (30, vec![(2, 1)]),
            (33, vec![(0, 0)]),
            (38, vec![(0, 1)]),
            (39, vec![(2, 0)]),
        ])
    );
    assert_eq!(output.get_atom_counts().unwrap(), &[2, 2, 2]);
    assert_eq!(
        output.get_atom_to_bits().unwrap(),
        &[vec![33, 38], vec![16, 2], vec![39, 30]]
    );
    assert!(output.get_bit_paths().is_none());

    let rdkit = get_rdkit_fp_generator_with_params(&GetRDKitFPGeneratorParams {
        fp_size: 64,
        max_path: 2,
        ..Default::default()
    })
    .unwrap();
    let mut output = AdditionalOutput::default();
    output.collect_bit_paths();
    let params = GetFingerprintParams {
        additional_output: Some(output),
        ..Default::default()
    };
    let (fp, output) = rdkit
        .get_fingerprint_with_outputs(&smiles("CCO"), &params)
        .unwrap();
    assert_eq!(fp.get_on_bits(), vec![0, 28, 31, 44, 50, 59]);
    assert_eq!(output.get_bit_paths().unwrap()[&0], vec![vec![0, 1]]);
}

#[test]
fn other_generators() {
    let mol = smiles("CCO");
    let atom_pairs = get_atom_pair_generator_with_params(&GetAtomPairGeneratorParams {
        fp_size: 64,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(
        atom_pairs
            .get_count_fingerprint(&mol)
            .unwrap()
            .get_nonzero_elements(),
        BTreeMap::from([(15, 1), (31, 2)])
    );
    assert_eq!(
        atom_pairs.get_fingerprint(&mol).unwrap().get_on_bits(),
        vec![60, 61]
    );
    let torsions =
        get_topological_torsion_generator_with_params(&GetTopologicalTorsionGeneratorParams {
            fp_size: 64,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        torsions
            .get_count_fingerprint(&smiles("CCCCO"))
            .unwrap()
            .get_nonzero_elements(),
        BTreeMap::from([(0, 2)])
    );
    let features = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 1,
        fp_size: 64,
        atom_invariants_generator: Some(get_morgan_feature_atom_inv_gen()),
        ..Default::default()
    })
    .unwrap();
    assert_eq!(
        features.get_fingerprint(&mol).unwrap().get_on_bits(),
        vec![0, 3, 37, 57, 61]
    );
    let simulated = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 1,
        fp_size: 64,
        count_simulation: true,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(
        simulated.get_fingerprint(&mol).unwrap().get_on_bits(),
        vec![0, 4, 8, 24, 28, 56]
    );
    // The default generators and invariant generators build. A torsion
    // needs four atoms.
    for generator in [
        get_morgan_generator(),
        get_atom_pair_generator(),
        get_rdkit_fp_generator(),
        get_topological_torsion_generator(),
    ] {
        assert!(
            generator
                .get_fingerprint(&smiles("CCCCO"))
                .unwrap()
                .get_num_on_bits()
                > 0
        );
    }
    let _ = (
        get_morgan_atom_inv_gen(true),
        get_morgan_feature_atom_inv_gen_with_patterns(&[mol_from_smarts("[#8]").unwrap()]),
        get_morgan_bond_inv_gen(),
        get_morgan_bond_inv_gen_with_params(&GetMorganBondInvGenParams::default()),
        get_atom_pair_atom_inv_gen(),
        get_atom_pair_atom_inv_gen_with_include_chirality(true),
        get_rdkit_atom_inv_gen(),
    );
}

#[test]
fn bulk_fingerprints() {
    let mols = vec![smiles("CCO"), smiles("c1ccccc1"), smiles("CN")];
    let generator = morgan();
    // RDKit's thread count negates a negative count in int, undefined for
    // i32::MIN; more threads than molecules have nothing to do.
    for threads in [1, 2, 0, i32::MIN, i32::MAX] {
        let fps = generator
            .get_fingerprints_with_num_threads(&mols, threads)
            .unwrap();
        assert_eq!(
            fps.iter().map(|f| f.get_num_on_bits()).collect::<Vec<_>>(),
            vec![6, 2, 3]
        );
    }
    assert_eq!(generator.get_count_fingerprints(&mols).unwrap().len(), 3);
    assert_eq!(generator.get_sparse_fingerprints(&mols).unwrap().len(), 3);
    assert_eq!(
        generator
            .get_sparse_count_fingerprints_with_num_threads(&mols, 3)
            .unwrap()
            .len(),
        3
    );
    assert_eq!(
        get_fps(&mols)
            .unwrap()
            .iter()
            .map(|f| f.get_num_on_bits())
            .collect::<Vec<_>>(),
        vec![6, 3, 3]
    );
    assert_eq!(
        get_sparse_count_fps_with_fp_type(&mols, FPType::AtomPairFp)
            .unwrap()
            .iter()
            .map(|f| f.get_nonzero_elements().len())
            .collect::<Vec<_>>(),
        vec![3, 3, 1]
    );
    assert!(get_count_fps(&[]).unwrap().is_empty());
    assert_eq!(
        get_sparse_fps_with_fp_type(&mols, FPType::RdkitFp)
            .unwrap()
            .len(),
        3
    );
}

#[test]
fn bulk_errors_do_not_abort() {
    // RDKit's threads do not catch: an exception in one calls std::terminate.
    let generator = get_atom_pair_generator_with_params(&GetAtomPairGeneratorParams {
        use_2d: false,
        ..Default::default()
    })
    .unwrap();
    let mols = vec![smiles("CCO"), smiles("CN")];
    assert!(matches!(
        generator.get_fingerprints_with_num_threads(&mols, 2),
        Err(Error::Value { .. })
    ));
}

#[test]
fn options() {
    let mut generator = morgan();
    let options = generator.get_options();
    assert_eq!(options.fp_size, 64);
    assert!(!options.count_simulation);
    assert_eq!(options.count_bounds, vec![1, 2, 4, 8]);
    assert_eq!(
        options.kind,
        FingerprintKindOptions::Morgan(MorganFingerprintOptions {
            radius: 1,
            only_nonzero_invariants: false,
            include_redundant_environments: false,
        })
    );
    assert!(generator.get_info_string().unwrap().contains("radius=1"));
    // The default options are those of GetMorganGenerator().
    let default = FingerprintOptions::default();
    assert_eq!(default, get_morgan_generator().get_options());
    assert_eq!(default.fp_size, 2048);
    assert!(matches!(
        default.kind,
        FingerprintKindOptions::Morgan(MorganFingerprintOptions { radius: 3, .. })
    ));

    let mut bigger = options.clone();
    bigger.fp_size = 128;
    bigger.set_count_bounds(&[1, 3]);
    generator.set_options(&bigger).unwrap();
    assert_eq!(generator.get_options().fp_size, 128);
    assert_eq!(generator.get_options().count_bounds, vec![1, 3]);

    let mut wrong_kind = bigger.clone();
    wrong_kind.kind = FingerprintKindOptions::AtomPair(AtomPairFingerprintOptions::default());
    assert!(matches!(
        generator.set_options(&wrong_kind),
        Err(Error::Value { .. })
    ));
    let mut no_bounds = bigger.clone();
    no_bounds.count_simulation = true;
    no_bounds.count_bounds.clear();
    assert!(matches!(
        generator.set_options(&no_bounds),
        Err(Error::Value { .. })
    ));
    // A refused change leaves the options alone.
    assert_eq!(generator.get_options(), bigger);
}

#[test]
fn generator_parameters_that_rdkit_trusts() {
    // getSparseFingerprint divides by the number of count bounds.
    let no_bounds = GetMorganGeneratorParams {
        count_simulation: true,
        count_bounds: Some(vec![]),
        ..Default::default()
    };
    assert!(matches!(
        get_morgan_generator_with_params(&no_bounds),
        Err(Error::Value { .. })
    ));
    // Bit vectors hold at most i32::MAX bits (U3).
    let huge = GetRDKitFPGeneratorParams {
        fp_size: u32::MAX,
        ..Default::default()
    };
    assert!(matches!(
        get_rdkit_fp_generator_with_params(&huge),
        Err(Error::Value { .. })
    ));
    let no_bits = GetRDKitFPGeneratorParams {
        num_bits_per_feature: 0,
        ..Default::default()
    };
    assert!(matches!(
        get_rdkit_fp_generator_with_params(&no_bits),
        Err(Error::Value { .. })
    ));
    // Torsion codes are shifted into 64 bits, 9 bits an atom (11 with
    // chirality).
    for (count, chirality, ok) in [
        (7, false, true),
        (8, false, false),
        (5, true, true),
        (6, true, false),
    ] {
        let params = GetTopologicalTorsionGeneratorParams {
            torsion_atom_count: count,
            include_chirality: chirality,
            ..Default::default()
        };
        assert_eq!(
            get_topological_torsion_generator_with_params(&params).is_ok(),
            ok,
            "{count} {chirality}"
        );
    }
    let mut tt = get_topological_torsion_generator();
    let mut options = tt.get_options();
    options.kind =
        FingerprintKindOptions::TopologicalTorsion(TopologicalTorsionFingerprintOptions {
            torsion_atom_count: 20,
            only_shortest_paths: false,
        });
    assert!(matches!(tt.set_options(&options), Err(Error::Value { .. })));
}

#[test]
fn call_arguments_that_rdkit_trusts() {
    let mol = smiles("CCO");
    for generator in [
        morgan(),
        get_rdkit_fp_generator(),
        get_topological_torsion_generator(),
        get_atom_pair_generator(),
    ] {
        // Morgan, RDKit, and torsion set these in unchecked bitsets.
        for params in [
            GetFingerprintParams {
                from_atoms: Some(vec![3]),
                ..Default::default()
            },
            GetFingerprintParams {
                ignore_atoms: Some(vec![99]),
                ..Default::default()
            },
        ] {
            assert!(matches!(
                generator.get_fingerprint_with_params(&mol, &params),
                Err(Error::Index { .. })
            ));
        }
        // The custom invariants are indexed by atom and bond.
        for params in [
            GetFingerprintParams {
                custom_atom_invariants: Some(vec![1, 2]),
                ..Default::default()
            },
            GetFingerprintParams {
                custom_bond_invariants: Some(vec![1]),
                ..Default::default()
            },
        ] {
            assert!(matches!(
                generator.get_count_fingerprint_with_params(&mol, &params),
                Err(Error::Value { .. })
            ));
        }
        let custom = GetFingerprintParams {
            custom_atom_invariants: Some(vec![1, 2, 3]),
            custom_bond_invariants: Some(vec![1, 2]),
            ..Default::default()
        };
        assert!(generator.get_fingerprint_with_params(&mol, &custom).is_ok());
    }
}

#[test]
fn atom_pairs_from_3d_distances() {
    let generator = get_atom_pair_generator_with_params(&GetAtomPairGeneratorParams {
        use_2d: false,
        ..Default::default()
    })
    .unwrap();
    // No conformer.
    assert!(matches!(
        generator.get_fingerprint(&smiles("CCO")),
        Err(Error::Value { .. })
    ));
    let block = |x: f64| {
        format!(
            "
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
{x:>10.4}    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"
        )
    };
    let near = mol_from_mol_block(&block(1.5)).unwrap();
    assert!(generator.get_fingerprint(&near).unwrap().get_num_on_bits() > 0);
    // RDKit casts each distance to unsigned int.
    let mut far = near.clone();
    far.get_conformer_mut()
        .unwrap()
        .set_atom_position(1, rdkit_cxx::geometry::Point3D::new(f64::NAN, 0.0, 0.0))
        .unwrap();
    assert!(matches!(
        generator.get_fingerprint(&far),
        Err(Error::Value { .. })
    ));
    far.get_conformer_mut()
        .unwrap()
        .set_atom_position(1, rdkit_cxx::geometry::Point3D::new(1e300, 0.0, 0.0))
        .unwrap();
    assert!(matches!(
        generator.get_fingerprint(&far),
        Err(Error::Value { .. })
    ));
}

#[test]
fn params_defaults() {
    let morgan = GetMorganGeneratorParams::default();
    assert_eq!(
        (
            morgan.radius,
            morgan.fp_size,
            morgan.count_simulation,
            morgan.use_bond_types
        ),
        (3, 2048, false, true)
    );
    let pairs = GetAtomPairGeneratorParams::default();
    assert_eq!(
        (
            pairs.min_distance,
            pairs.max_distance,
            pairs.use_2d,
            pairs.count_simulation
        ),
        (1, 30, true, true)
    );
    let rdkit = GetRDKitFPGeneratorParams::default();
    assert_eq!(
        (
            rdkit.min_path,
            rdkit.max_path,
            rdkit.num_bits_per_feature,
            rdkit.use_hs
        ),
        (1, 7, 2, true)
    );
    let tt = GetTopologicalTorsionGeneratorParams::default();
    assert_eq!((tt.torsion_atom_count, tt.count_simulation), (4, true));
    assert_eq!(GetFingerprintParams::default().conf_id, -1);
    // The generators' own options agree with Python's.
    let options = get_atom_pair_generator().get_options();
    assert_eq!(
        options.kind,
        FingerprintKindOptions::AtomPair(AtomPairFingerprintOptions {
            min_distance: 1,
            max_distance: 30,
            use_2d: true,
        })
    );
    assert_eq!(
        get_rdkit_fp_generator().get_options().num_bits_per_feature,
        2
    );
}
