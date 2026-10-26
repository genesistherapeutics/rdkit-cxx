// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.Pharm2D in RDKit 2026.09.1, and the tests
// of UnitTestGenerate.py, UnitTestGobbi.py, UnitTestMatcher.py,
// UnitTestSignature.py, and UnitTestUtils.py with the doctests.

use rdkit_cxx::Error;
use rdkit_cxx::chem::pharm_2d::generate::{
    Gen2DFingerprintParams, gen_2d_fingerprint, gen_2d_fingerprint_with_outputs,
    gen_2d_fingerprint_with_params,
};
use rdkit_cxx::chem::pharm_2d::matcher::{
    GetAtomsMatchingBitParams, get_atoms_matching_bit, get_atoms_matching_bit_with_params,
};
use rdkit_cxx::chem::pharm_2d::sig_factory::{SigFactory, SigFactoryParams, Signature};
use rdkit_cxx::chem::pharm_2d::{gobbi_pharm_2d, utils};
use rdkit_cxx::chem::rd_mol_chemical_features::build_feature_factory;
use rdkit_cxx::chem::{Mol, SDMolSupplier, mol_from_smiles, mol_to_smiles};

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn on_bits(sig: &Signature) -> Vec<u32> {
    match sig {
        Signature::Bits(bits) => bits.get_on_bits(),
        other => panic!("not bits: {other:?}"),
    }
}

fn sig_len(sig: &Signature) -> u64 {
    match sig {
        Signature::Bits(bits) => bits.get_num_bits().into(),
        Signature::Counts(counts) => counts.get_length() as u64,
        Signature::LongCounts(counts) => counts.get_length() as u64,
        other => panic!("unknown signature: {other:?}"),
    }
}

fn smiles(s: &str) -> Mol {
    mol_from_smiles(s).unwrap()
}

/// UnitTestSignature.py's factory.
fn base_factory(triangle_prune_bins: bool) -> SigFactory {
    let feat_factory = build_feature_factory(data("Pharm2D_BaseFeatures.fdef")).unwrap();
    let params = SigFactoryParams {
        min_point_count: 2,
        max_point_count: 3,
        triangle_prune_bins,
        ..Default::default()
    };
    let mut factory = SigFactory::new_with_params(feat_factory, &params);
    factory.set_bins(&[(0, 2), (2, 5), (5, 8)]).unwrap();
    factory.init().unwrap();
    factory
}

#[test]
fn utils_match_python() {
    for (combo, tgt) in [
        ([0, 1, 2], 4),
        ([0, 0, 0], 0),
        ([2, 2, 2], 9),
        ([1, 1, 2], 7),
        ([1, 2, 2], 8),
    ] {
        assert_eq!(utils::count_up_to(3, 3, &combo).unwrap(), tgt, "{combo:?}");
    }
    assert_eq!(utils::get_triangles(2), []);
    assert_eq!(utils::get_triangles(3), [(0, 1, 2)]);
    assert_eq!(utils::get_triangles(4), [(0, 1, 3), (1, 2, 4)]);
    assert_eq!(utils::get_triangles(5), [(0, 1, 4), (1, 2, 5), (2, 3, 6)]);

    let bins = [(1, 2), (2, 3), (5, 6)];
    for (ds, tgt) in [
        ([0, 0, 0], true),
        ([0, 0, 1], true),
        ([0, 0, 2], false),
        ([0, 1, 2], true),
        ([1, 1, 2], true),
    ] {
        assert_eq!(
            utils::bins_triangle_inequality(bins[ds[0]], bins[ds[1]], bins[ds[2]]),
            tgt
        );
        assert_eq!(utils::scaffold_passes(&ds, &bins).unwrap(), tgt);
    }
    assert!(matches!(
        utils::scaffold_passes(&[0, 0], &bins),
        Err(Error::Key { .. })
    ));
    assert_eq!(utils::get_possible_scaffolds(2, &bins).unwrap().len(), 3);
    assert_eq!(utils::get_possible_scaffolds(3, &bins).unwrap().len(), 24);
    assert_eq!(utils::get_possible_scaffolds(1, &bins).unwrap().len(), 0);
    assert_eq!(
        utils::get_possible_scaffolds_with_use_triangle_inequality(3, &bins, false)
            .unwrap()
            .len(),
        27
    );

    assert!(matches!(
        utils::order_triangle(&[0, 2], &[1, 2, 3]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        utils::order_triangle(&[0, 2, 4], &[1, 2]),
        Err(Error::Value { .. })
    ));
    let ot = |f: [usize; 3], d: [u32; 3]| utils::order_triangle(&f, &d).unwrap();
    assert_eq!(ot([1, 3, 1], [2, 3, 4]), (vec![1, 3, 1], vec![4, 3, 2]));
    assert_eq!(ot([1, 3, 1], [4, 3, 2]), (vec![1, 3, 1], vec![4, 3, 2]));
    for d in [
        [1, 2, 3],
        [1, 3, 2],
        [2, 1, 3],
        [2, 3, 1],
        [3, 1, 2],
        [3, 2, 1],
    ] {
        assert_eq!(ot([1, 1, 1], d), (vec![1, 1, 1], vec![3, 2, 1]));
        assert_eq!(ot([0, 0, 0], d), (vec![0, 0, 0], vec![3, 2, 1]));
    }
    assert_eq!(ot([0, 2, 4], [1, 2, 3]), (vec![0, 2, 4], vec![1, 2, 3]));
    assert_eq!(ot([0, 0, 1], [3, 2, 1]), (vec![0, 0, 1], vec![3, 2, 1]));
    assert_eq!(ot([0, 0, 1], [1, 3, 2]), (vec![0, 0, 1], vec![1, 3, 2]));
    assert_eq!(ot([0, 0, 1], [1, 2, 3]), (vec![0, 0, 1], vec![1, 3, 2]));

    assert_eq!(
        utils::uniquify_combinations(&[vec![1, 2, 3], vec![3, 2, 1]]),
        [vec![3, 2, 1]]
    );
    assert_eq!(
        utils::uniquify_combinations(&[vec![1], vec![1], vec![2]]),
        [vec![1], vec![2]]
    );

    assert_eq!(
        utils::get_all_combinations(&[vec![0], vec![1], vec![2]]),
        [vec![0, 1, 2]]
    );
    assert_eq!(
        utils::get_all_combinations(&[vec![0], vec![1, 3], vec![2]]),
        [vec![0, 1, 2], vec![0, 3, 2]]
    );
    assert_eq!(
        utils::get_all_combinations(&[vec![0, 1], vec![1, 3], vec![2]]),
        [vec![0, 1, 2], vec![0, 3, 2], vec![1, 3, 2]]
    );
    assert_eq!(
        utils::get_all_combinations_with_no_dups(&[vec![0, 1], vec![1, 3], vec![2]], false),
        [vec![0, 1, 2], vec![0, 3, 2], vec![1, 1, 2], vec![1, 3, 2]]
    );
    assert_eq!(
        utils::get_index_combinations(3, 2),
        [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]
    );
    assert_eq!(utils::num_combinations(5, 2), 15);
    assert_eq!(utils::num_combinations(8, 3), 120);

    let n_point = utils::n_point_dist_dict();
    assert_eq!(n_point[&3], [(0, 1), (0, 2), (1, 2)]);
    assert_eq!(n_point[&10].len(), 17);
    assert_eq!(utils::n_dist_point_dict()[&17], 10);
    assert_eq!(utils::n_dist_point_dict().len(), 9);
}

#[test]
fn unique_combinations_match_python() {
    let c = |xs: &[&[u32]]| -> Vec<Vec<u32>> { xs.iter().map(|x| x.to_vec()).collect() };
    let choices = vec![c(&[&[11], &[12]]), c(&[&[11], &[12]]), c(&[&[31], &[32]])];
    assert_eq!(
        utils::get_unique_combinations(&choices, &[1, 1, 2]).unwrap(),
        [
            vec![(1, vec![11]), (1, vec![12]), (2, vec![31])],
            vec![(1, vec![11]), (1, vec![12]), (2, vec![32])]
        ]
    );
    let choices = vec![
        c(&[&[11, 111], &[12], &[14]]),
        c(&[&[11, 111], &[12], &[13]]),
        c(&[&[31], &[32]]),
    ];
    let combos = utils::get_unique_combinations(&choices, &[1, 1, 2]).unwrap();
    assert_eq!(combos.len(), 12);
    assert_eq!(
        combos[0],
        [(1, vec![11, 111]), (1, vec![12]), (2, vec![31])]
    );
    assert_eq!(combos[11], [(1, vec![13]), (1, vec![14]), (2, vec![32])]);
    let choices = vec![c(&[&[11], &[12]]), c(&[&[11], &[13]]), c(&[&[11], &[14]])];
    assert_eq!(
        utils::get_unique_combinations(&choices, &[1, 1, 1])
            .unwrap()
            .len(),
        4
    );
    let choices = vec![c(&[&[0], &[4]]), c(&[&[0]])];
    assert_eq!(
        utils::get_unique_combinations(&choices, &[0, 1]).unwrap(),
        [vec![(0, vec![4]), (1, vec![0])]]
    );
    let empty: Vec<Vec<Vec<u32>>> = vec![vec![], vec![], vec![]];
    assert!(
        utils::get_unique_combinations(&empty, &[0, 1, 1])
            .unwrap()
            .is_empty()
    );
    let choices = vec![c(&[&[0], &[4]]), c(&[&[0], &[4]]), c(&[&[0], &[4]])];
    assert!(
        utils::get_unique_combinations(&choices, &[0, 0, 0])
            .unwrap()
            .is_empty()
    );
    let choices = vec![
        c(&[&[4], &[6], &[7], &[10]]),
        c(&[&[2], &[3], &[5], &[13], &[14]]),
        c(&[&[2], &[3], &[5], &[13], &[14]]),
    ];
    assert_eq!(
        utils::get_unique_combinations(&choices, &[0, 1, 1])
            .unwrap()
            .len(),
        40
    );
    assert!(matches!(
        utils::get_unique_combinations(&choices, &[0, 1]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn signature_sizes_and_bits_match_python() {
    // UnitTestSignature.test1Sizes.
    let mut factory = base_factory(false);
    for (max, size) in [(2, 45), (3, 990), (4, 18000)] {
        factory.max_point_count = max;
        factory.init().unwrap();
        assert_eq!(sig_len(&factory.get_signature().unwrap()), size);
    }
    // test2BitIdx.
    let factory = base_factory(false);
    for (patts, dists, bit) in [
        (&[0, 0][..], &[0][..], 0),
        (&[0, 0], &[2], 1),
        (&[0, 0], &[5], 2),
        (&[0, 1], &[5], 5),
        (&[1, 1], &[4], 16),
        (&[1, 1], &[7], 17),
        (&[0, 0, 0], &[1, 1, 1], 45),
        (&[0, 0, 1], &[1, 1, 1], 72),
        (&[0, 0, 1], &[1, 1, 3], 75),
        (&[0, 0, 1], &[3, 1, 1], 81),
        (&[0, 0, 1], &[3, 3, 1], 84),
    ] {
        assert_eq!(factory.get_bit_idx(patts, dists).unwrap(), bit);
        let (cnt, feats, _) = factory.get_bit_info(bit).unwrap();
        assert_eq!((cnt, feats.as_slice()), (patts.len(), patts));
    }
    // test3BitIdx.
    let mut factory = base_factory(false);
    factory.set_bins(&[(0, 2), (2, 4), (4, 8)]).unwrap();
    assert_eq!(factory.get_sig_size(), 990);
    for (patts, dists, ans) in [
        ([0, 0, 0], [1, 3, 1], 54),
        ([0, 0, 0], [3, 1, 1], 54),
        ([0, 0, 0], [1, 1, 3], 54),
        ([0, 0, 0], [1, 3, 3], 57),
        ([0, 0, 1], [1, 3, 1], 75),
    ] {
        assert_eq!(factory.get_bit_idx(&patts, &dists).unwrap(), ans);
        assert_eq!(factory.get_bit_info(ans).unwrap().1, patts);
    }
    // test4BitIdx.
    let mut factory = base_factory(false);
    factory.triangle_prune_bins = true;
    factory.init().unwrap();
    assert_eq!(sig_len(&factory.get_signature().unwrap()), 885);
    for (patts, dists, ans) in [
        ([0, 0, 0], [1, 3, 1], 52),
        ([0, 0, 0], [1, 1, 3], 52),
        ([0, 0, 0], [3, 1, 1], 52),
        ([0, 0, 0], [1, 3, 3], 55),
        ([0, 0, 1], [1, 3, 1], 71),
    ] {
        assert_eq!(factory.get_bit_idx(&patts, &dists).unwrap(), ans);
        assert_eq!(factory.get_bit_info(ans).unwrap().1, patts);
    }
    assert!(matches!(
        factory.get_bit_idx(&[0, 0, 0, 0], &[1, 3, 1, 2]),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        factory.get_bit_idx(&[0], &[1]),
        Err(Error::Index { .. })
    ));
    factory.max_point_count = 2;
    assert!(matches!(
        factory.get_bit_idx(&[0, 0, 0], &[1, 3, 1]),
        Err(Error::Index { .. })
    ));
    factory.max_point_count = 3;
    assert!(matches!(
        factory.get_bit_idx(&[0, 5, 0], &[1, 3, 1]),
        Err(Error::Index { .. })
    ));
    assert!(matches!(
        factory.get_bit_idx(&[0, 0], &[9]),
        Err(Error::Index { .. })
    ));
    assert!(matches!(
        factory.get_bit_info(885),
        Err(Error::Index { .. })
    ));
}

#[test]
fn simple_signatures_match_python() {
    // UnitTestSignature.test5SimpleSig to test9BondOrderSigs.
    let mut factory = base_factory(false);
    factory.set_bins(&[(1, 3), (3, 7), (7, 10)]).unwrap();
    for (smi, bits) in [
        ("O=CCC=O", vec![1]),
        ("O=CC(CC=O)CCC=O", vec![1, 2, 67]),
        ("O=Cc1ccccc1", vec![3]),
        ("O=CCCCCCCCCc1ccccc1", vec![]),
    ] {
        let sig = gen_2d_fingerprint(&smiles(smi), &factory).unwrap();
        assert_eq!(sig_len(&sig), 990);
        assert_eq!(on_bits(&sig), bits, "{smi}");
    }

    factory.use_counts = true;
    factory.init().unwrap();
    for (smi, counts) in [
        ("O=CCC=O", vec![(1, 1)]),
        ("O=CC(CC=O)CCC=O", vec![(1, 2), (2, 1), (67, 1)]),
    ] {
        let Signature::Counts(sig) = gen_2d_fingerprint(&smiles(smi), &factory).unwrap() else {
            panic!("not counts");
        };
        assert_eq!(sig.get_length(), 990);
        assert_eq!(
            sig.get_nonzero_elements().into_iter().collect::<Vec<_>>(),
            counts
        );
    }

    factory.use_counts = false;
    factory.skip_feats = vec!["Acceptor".to_owned()];
    factory.init().unwrap();
    let sig = gen_2d_fingerprint(&smiles("O=CCC=O"), &factory).unwrap();
    assert_eq!(sig_len(&sig), 570);
    assert!(on_bits(&sig).is_empty());

    let mut factory = base_factory(false);
    factory.set_bins(&[(1, 4), (4, 7), (7, 10)]).unwrap();
    let mol = smiles("[O-]CCC(=O)");
    assert_eq!(on_bits(&gen_2d_fingerprint(&mol, &factory).unwrap()), [1]);
    factory.include_bond_order = true;
    assert_eq!(on_bits(&gen_2d_fingerprint(&mol, &factory).unwrap()), [0]);

    // Python's perms and dMat arguments.
    let mut factory = base_factory(false);
    factory.set_bins(&[(1, 3), (3, 7), (7, 10)]).unwrap();
    let mol = smiles("O=CC(CC=O)CCC=O");
    let params = Gen2DFingerprintParams {
        perms: Some(vec![vec![0, 0, 0]]),
        d_mat: None,
    };
    assert_eq!(
        on_bits(&gen_2d_fingerprint_with_params(&mol, &factory, &params).unwrap()),
        [67]
    );
    let n = mol.get_num_atoms() as usize;
    let params = Gen2DFingerprintParams {
        perms: None,
        d_mat: Some(vec![vec![2.0; n]; n]),
    };
    assert_eq!(
        on_bits(&gen_2d_fingerprint_with_params(&mol, &factory, &params).unwrap()),
        [0, 45]
    );
}

#[test]
fn default_factory_matches_python() {
    let factory = rdkit_cxx::chem::pharm_2d::default_sig_factory().unwrap();
    assert_eq!(factory.get_num_bins(), 7);
    let sig = gen_2d_fingerprint(&smiles("OCCC(=O)"), &factory).unwrap();
    assert_eq!(sig_len(&sig), 19355);
    assert_eq!(on_bits(&sig), [2, 16, 21, 84, 1274, 4361]);
    assert_eq!(factory.get_bit_info(21).unwrap(), (2, vec![0, 3], vec![0]));
    assert_eq!(
        factory.get_bit_description(21).unwrap(),
        "Acceptor Hydrophobe |0 0|0 0|"
    );
    assert_eq!(
        factory.get_bit_description(2).unwrap(),
        "Acceptor Acceptor |0 2|2 0|"
    );
    assert_eq!(
        factory.get_bit_description(4361).unwrap(),
        "Acceptor Donor Hydrophobe |0 2 0|2 0 0|0 0 0|"
    );
    assert_eq!(
        factory.get_feat_families(),
        [
            "Acceptor",
            "Aromatic",
            "Donor",
            "Hydrophobe",
            "NegIonizable",
            "PosIonizable"
        ]
    );
    let with_file = rdkit_cxx::chem::pharm_2d::default_sig_factory_with_params(
        &rdkit_cxx::chem::pharm_2d::DefaultSigFactoryParams {
            fdef_file: Some(data("Pharm2D_BaseFeatures.fdef").into()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(with_file.get_feat_families().len(), 5);
}

#[test]
fn gobbi_matches_python() {
    let factory = gobbi_pharm_2d::factory();
    assert_eq!(factory.get_sig_size(), 39972);
    assert_eq!(
        gobbi_pharm_2d::default_bins(),
        [(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 100)]
    );
    assert!(gobbi_pharm_2d::fdef().contains("DefineFeature AcidicGroup"));
    assert_eq!(
        factory.get_feat_families(),
        ["AG", "AR", "BG", "HA", "HD", "LH", "RR", "X"]
    );

    // UnitTestGobbi.test1Sigs: each family's features' first atoms.
    let first_atoms = |smi: &str| -> Vec<Vec<u32>> {
        factory
            .get_mol_feats(&smiles(smi))
            .unwrap()
            .into_iter()
            .map(|feats| {
                let mut atoms: Vec<u32> = feats.iter().map(|f| f[0]).collect();
                atoms.sort();
                atoms
            })
            .collect()
    };
    //                                           AG     AR     BG     HA      HD      LH      RR     X
    assert_eq!(
        first_atoms("OCCC=O"),
        [
            vec![],
            vec![],
            vec![],
            vec![0, 4],
            vec![0],
            vec![],
            vec![],
            vec![]
        ]
    );
    assert_eq!(
        first_atoms("OCCC(=O)O"),
        [
            vec![3],
            vec![],
            vec![],
            vec![0, 4],
            vec![0, 5],
            vec![],
            vec![],
            vec![]
        ]
    );
    assert_eq!(
        first_atoms("CCCN"),
        [
            vec![],
            vec![],
            vec![3],
            vec![3],
            vec![3],
            vec![],
            vec![],
            vec![]
        ]
    );
    assert_eq!(
        first_atoms("CCCCC"),
        [
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![1, 3],
            vec![],
            vec![]
        ]
    );
    assert_eq!(
        first_atoms("CC1CCC1"),
        [
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![1, 3],
            vec![1],
            vec![]
        ]
    );
    assert_eq!(
        first_atoms("[SiH3]C1CCC1"),
        [
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![1],
            vec![1],
            vec![0]
        ]
    );
    assert_eq!(
        first_atoms("[SiH3]c1ccccc1"),
        [
            vec![],
            vec![1],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![0]
        ]
    );

    // test2Sigs.
    for (smi, tgt) in [
        ("O=CCC=O", vec![149]),
        ("OCCC=O", vec![149, 156]),
        ("OCCC(=O)O", vec![22, 29, 149, 154, 156, 184, 28822, 30134]),
    ] {
        let sig = gen_2d_fingerprint(&smiles(smi), &factory).unwrap();
        assert_eq!(sig_len(&sig), 39972);
        assert_eq!(on_bits(&sig), tgt);
    }

    // testOrderBug and testOrderBug2.
    let mut suppl = SDMolSupplier::new(data("Pharm2D_orderBug.sdf")).unwrap();
    let m1 = suppl.get(0).unwrap();
    let m2 = suppl.get(1).unwrap();
    assert_eq!(
        on_bits(&gen_2d_fingerprint(&m1, &factory).unwrap()),
        on_bits(&gen_2d_fingerprint(&m2, &factory).unwrap())
    );
    let m1 = smiles("Oc1nc(Oc2ncccc2)ccc1");
    let m2 = smiles(&mol_to_smiles(&m1).unwrap());
    assert_eq!(
        on_bits(&gen_2d_fingerprint(&m1, &factory).unwrap()),
        on_bits(&gen_2d_fingerprint(&m2, &factory).unwrap())
    );

    // testBitInfo.
    let m = smiles("OCC=CC(=O)O");
    let (sig, bi) = gen_2d_fingerprint_with_outputs(&m, &factory, &Default::default()).unwrap();
    assert_eq!(
        bi.keys().map(|&k| k as u32).collect::<Vec<_>>(),
        on_bits(&sig)
    );
    assert_eq!(
        bi.keys().copied().collect::<Vec<_>>(),
        [23, 30, 150, 154, 157, 185, 28878, 30184]
    );
    assert_eq!(bi[&28878], [vec![vec![0], vec![5], vec![6]]]);
    let mut bit157 = bi[&157].clone();
    bit157.sort();
    assert_eq!(bit157, [vec![vec![0], vec![6]], vec![vec![5], vec![0]]]);
}

#[test]
fn matcher_matches_python() {
    // UnitTestMatcher.test1_simple.
    let mut factory = base_factory(true);
    let mol = smiles("OCC(=O)CCCN");
    factory.skip_feats = vec!["Donor".to_owned()];
    factory.init().unwrap();
    assert_eq!(factory.get_sig_size(), 510);
    let sig = gen_2d_fingerprint(&mol, &factory).unwrap();
    assert_eq!(on_bits(&sig), [1, 2, 11, 52, 117]);
    let matches: Vec<_> = on_bits(&sig)
        .into_iter()
        .map(|bit| get_atoms_matching_bit(&factory, bit.into(), &mol).unwrap())
        .collect();
    let pair = |a: u32, b: u32| vec![vec![a], vec![b]];
    let ac_n = vec![pair(0, 7), pair(3, 7)];
    let triple = vec![vec![vec![0], vec![3], vec![7]]];
    assert_eq!(
        matches,
        [vec![pair(0, 3)], ac_n.clone(), ac_n, triple.clone(), triple]
    );

    // test2Bug28: every bit matches.
    let mut factory = gobbi_pharm_2d::factory();
    factory
        .set_bins(&[(2, 3), (3, 4), (4, 5), (5, 8), (8, 100)])
        .unwrap();
    let mol = smiles(r"Cc([s]1)nnc1SCC(\CS2)=C(/C([O-])=O)N3C(=O)[C@H]([C@@H]23)NC(=O)C[n]4cnnn4");
    let sig = gen_2d_fingerprint(&mol, &factory).unwrap();
    let params = GetAtomsMatchingBitParams {
        just_one: true,
        ..Default::default()
    };
    for bit in on_bits(&sig) {
        let atoms =
            get_atoms_matching_bit_with_params(&factory, bit.into(), &mol, &params).unwrap();
        assert_eq!(atoms.len(), 1, "{bit}");
    }

    // Python's _exampleCode.
    let feat_factory = build_feature_factory(data("Pharm2D_BaseFeatures.fdef")).unwrap();
    let mut factory = SigFactory::new(feat_factory);
    factory.set_bins(&[(1, 2), (2, 5), (5, 8)]).unwrap();
    let mol = smiles("OCC(=O)CCCN");
    let sig = gen_2d_fingerprint(&mol, &factory).unwrap();
    assert!(!on_bits(&sig).is_empty());
    let mut not_shortest = SigFactory::new_with_params(
        build_feature_factory(data("Pharm2D_BaseFeatures.fdef")).unwrap(),
        &SigFactoryParams {
            shortest_paths_only: false,
            ..Default::default()
        },
    );
    not_shortest.set_bins(&[(1, 2), (2, 5), (5, 8)]).unwrap();
    assert!(matches!(
        get_atoms_matching_bit(&not_shortest, 0, &mol),
        Err(Error::Value { .. })
    ));
}

/// Each molecule's Gobbi bits, default-factory bits, Gobbi features, and
/// for each Gobbi bit its description, its bit info (sorted), and its
/// matching atoms, as text; Python's lists of integers print as Rust's.
fn pharm2d_text(mol: &Mol) -> (usize, usize, String) {
    let gobbi = gobbi_pharm_2d::factory();
    let default = rdkit_cxx::chem::pharm_2d::default_sig_factory().unwrap();
    let (sig, info) = gen_2d_fingerprint_with_outputs(mol, &gobbi, &Default::default()).unwrap();
    let bits = on_bits(&sig);
    let default_bits = on_bits(&gen_2d_fingerprint(mol, &default).unwrap());
    let mut lines = vec![
        format!("{bits:?}"),
        format!("{default_bits:?}"),
        format!("{:?}", gobbi.get_mol_feats(mol).unwrap()),
    ];
    for &b in &bits {
        let mut bit_info = info[&u64::from(b)].clone();
        bit_info.sort();
        let matches = get_atoms_matching_bit(&gobbi, b.into(), mol).unwrap();
        lines.push(format!(
            "{b} {} {bit_info:?} {matches:?}",
            gobbi.get_bit_description(b.into()).unwrap()
        ));
    }
    (bits.len(), default_bits.len(), lines.join("\n"))
}

#[test]
fn fingerprints_match_python() {
    for (smi, n_gobbi, n_default, digest) in [
        (
            "OCC=CC(=O)O",
            8,
            26,
            "227f0cfe3b34ba880462a0f1955654c0ed87025e",
        ),
        (
            "CC(=O)Nc1ccc(O)cc1",
            21,
            8,
            "3b1b04d56a69bae993e12947f95976171a9f6aa6",
        ),
        (
            "CN1CCC[C@H]1c2cccnc2",
            28,
            9,
            "6be9b2bf89313c033ad287bf9d4a67bee3c3213d",
        ),
        (
            "O=C([O-])CC[NH3+]",
            3,
            9,
            "c1cdd76247302994a2208e16b04fb8d0722fdddc",
        ),
        (
            "c1ccc2[nH]ccc2c1",
            2,
            1,
            "0d1b9d3ba50bf317b0c1776af21ca00d48b94d71",
        ),
        (
            "OCCC(=O)O",
            8,
            18,
            "25c869751934a8bb635742b10a9ffc3316d32c0e",
        ),
        (
            "Cc([s]1)nnc1SCC(\\CS2)=C(/C([O-])=O)N3C(=O)[C@H]([C@@H]23)NC(=O)C[n]4cnnn4",
            762,
            257,
            "c9991f8f463fda12b9fb082af46b851ccf8f6b8a",
        ),
        (
            "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
            26,
            60,
            "36ce1042f58eadd8cc987dcfa408ff167a5c001d",
        ),
        (
            "Oc1nc(Oc2ncccc2)ccc1",
            50,
            17,
            "efbdc5453a16f5d01d5db6f4c18b164903e789cf",
        ),
        (
            "NCCCCC(N)C(=O)O",
            149,
            104,
            "3e18f4de7bec4ae9b41c7ceb31199f8141897e63",
        ),
    ] {
        let (gobbi, default, text) = pharm2d_text(&smiles(smi));
        assert_eq!((gobbi, default), (n_gobbi, n_default), "{smi}");
        assert_eq!(
            sha1_smol::Sha1::from(text).digest().to_string(),
            digest,
            "{smi}"
        );
    }
}
