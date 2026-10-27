// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/DataStructs/UnitTestBitEnsemble.py (but its database tests),
// UnitTestDocTests.py (with the doctests of BitUtils, LazySignature, and
// VectCollection), UnitTestFlatMatrix.py, and UnitTestTopNContainer.py,
// with values from Python's rdkit.DataStructs in RDKit 2025.03.6.

use std::cell::Cell;

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rdmolops::{
    RDKFingerprintParams, rdk_fingerprint, rdk_fingerprint_with_params,
};
use rdkit_cxx::data_structs::bit_ensemble::BitEnsemble;
use rdkit_cxx::data_structs::bit_utils::construct_ensemble_bv;
use rdkit_cxx::data_structs::lazy_signature::LazySig;
use rdkit_cxx::data_structs::top_n_container::TopNContainer;
use rdkit_cxx::data_structs::vect_collection::VectCollection;
use rdkit_cxx::data_structs::*;
use rdkit_cxx::rd_random::Random;

fn bits(size: u32, on: &[u32]) -> ExplicitBitVect {
    let mut bv = ExplicitBitVect::new(size);
    bv.set_bits_from_list(on).unwrap();
    bv
}

#[test]
fn bit_ensemble() {
    let mut bv = SparseBitVect::new(100);
    for b in [1, 11, 13] {
        bv.set_bit(b).unwrap();
    }
    let mut set = BitEnsemble::default();
    set.set_bits(vec![1, 11, 21, 31]);
    assert_eq!(set.get_num_bits(), 4);
    let given = BitEnsemble::new_with_bits(vec![1, 11, 21, 31]);
    let mut added = BitEnsemble::default();
    for bit in [1, 11, 21, 31] {
        added.add_bit(bit);
    }
    for ensemble in [set, given, added] {
        assert_eq!(ensemble.get_bits(), [1, 11, 21, 31]);
        assert_eq!(ensemble.score_with_on_bits(&bv), 2);
        assert_eq!(ensemble.score_with_index(&bv).unwrap(), 2);
    }
    let past = BitEnsemble::new_with_bits(vec![100]);
    assert!(matches!(
        past.score_with_index(&bv),
        Err(Error::Index { .. })
    ));
}

#[test]
fn construct_ensemble() {
    // BitUtils' doctest.
    let bv = bits(128, &[1, 5, 47, 99, 120]);
    let r = construct_ensemble_bv(&bv, &[0, 1, 2, 3, 45, 46, 47, 48, 49]).unwrap();
    assert_eq!(r.get_num_bits(), 9);
    assert_eq!(r.get_on_bits(), [1, 6]);
    let mut sparse = SparseBitVect::new(128);
    sparse.set_bits_from_list(&[1, 5, 47]).unwrap();
    let r: SparseBitVect = construct_ensemble_bv(&sparse, &[47, 2, 1]).unwrap();
    assert_eq!(r.get_on_bits(), [0, 2]);
    assert!(matches!(
        construct_ensemble_bv(&bv, &[128]),
        Err(Error::Index { .. })
    ));
}

#[test]
fn lazy_signature() {
    // LazySignature's doctests.
    let calls = Cell::new(0);
    let mut obj = LazySig::new(
        |x: usize| {
            calls.set(calls.get() + 1);
            x
        },
        10,
    )
    .unwrap();
    assert_eq!(obj.len(), 10);
    assert_eq!(obj.get(1).unwrap(), 1);
    assert_eq!(obj.get(-1).unwrap(), 9);
    assert!(matches!(obj.get(10), Err(Error::Index { .. })));
    assert!(matches!(obj.get(-10), Err(Error::Index { .. })));
    // Python refuses index 0 too.
    assert!(matches!(obj.get(0), Err(Error::Index { .. })));
    assert_eq!(obj.get(8).unwrap(), 8);
    assert_eq!(obj.get(8).unwrap(), 8);
    assert_eq!(calls.get(), 3);
    assert!(matches!(
        LazySig::new(|_: usize| 1, 0),
        Err(Error::Value { .. })
    ));
}

#[test]
fn top_n_container() {
    let mut cont = TopNContainer::new(4);
    for foo in 0..10 {
        cont.insert_with_extra(foo as f64, foo.to_string()).unwrap();
    }
    assert_eq!(cont.get_pts(), [6.0, 7.0, 8.0, 9.0]);
    let extras: Vec<String> = cont
        .get_extras()
        .iter()
        .map(|e| e.clone().unwrap())
        .collect();
    assert_eq!(extras, ["6", "7", "8", "9"]);

    let mut rng = Random::new(23);
    let mut cont: TopNContainer = TopNContainer::new(50);
    for _ in 0..1000 {
        cont.insert(rng.random()).unwrap();
    }
    assert!(cont.get_pts().windows(2).all(|w| w[0] <= w[1]));

    let mut cont = TopNContainer::new(10);
    for _ in 0..1000 {
        let v = rng.random();
        cont.insert_with_extra(v, v + 1.0).unwrap();
    }
    for i in 1..10 {
        let (v, e) = cont.get(i).unwrap();
        let (last_v, last_e) = cont.get(i - 1).unwrap();
        assert!(v >= last_v && e >= last_e);
    }
    cont.reverse();
    for i in 1..10 {
        let (v, e) = cont.get(i).unwrap();
        let (last_v, last_e) = cont.get(i - 1).unwrap();
        assert!(v <= last_v && e <= last_e);
    }

    let mut cont = TopNContainer::new(-1);
    for i in 0..10 {
        cont.insert_with_extra((9 - i) as f64, (9 - i).to_string())
            .unwrap();
        assert_eq!(cont.len(), i + 1);
    }
    let expected: Vec<f64> = (0..10).map(f64::from).collect();
    assert_eq!(cont.get_pts(), expected);

    // Values from Python.
    let mut cont = TopNContainer::new(4);
    for v in [5, 1, 7, 3, 9, 2, 8] {
        cont.insert_with_extra(v as f64, v.to_string()).unwrap();
    }
    assert_eq!(cont.get_pts(), [5.0, 7.0, 8.0, 9.0]);
    assert_eq!(cont.get(-1).unwrap(), (9.0, Some(&"9".to_owned())));
    assert!(matches!(cont.get(4), Err(Error::Index { .. })));
    let mut cont = TopNContainer::new_with_most_neg(3, 0.0);
    cont.insert(-1.0).unwrap();
    cont.insert(2.0).unwrap();
    cont.insert_with_extra(2.0, 'b').unwrap();
    assert_eq!(cont.get_pts(), [0.0, 2.0, 2.0]);
    assert_eq!(cont.get_extras(), [None, None, Some('b')]);
    let mut empty: TopNContainer = TopNContainer::new(0);
    assert!(matches!(empty.insert(1.0), Err(Error::Index { .. })));
    // Once a NaN is kept, the list is unsorted, and only Python's bisect
    // probes give Python's order.
    let mut cont = TopNContainer::new(-1);
    for (i, v) in [0.0, 1.0, 2.0, 1.0, f64::NAN, 1.0].into_iter().enumerate() {
        cont.insert_with_extra(v, i).unwrap();
    }
    let extras: Vec<usize> = cont.get_extras().iter().map(|e| e.unwrap()).collect();
    assert_eq!(extras, [0, 1, 3, 2, 4, 5]);
}

#[test]
fn vect_collection() {
    // VectCollection's doctests.
    let mut vc = VectCollection::default();
    assert_eq!(vc.get_or_vect().unwrap(), None);
    assert_eq!(vc.get_num_bits().unwrap(), -1);
    vc.add_vect(1, bits(10, &[1, 3, 5]));
    vc.add_vect(2, bits(10, &[6, 8]));
    assert_eq!(vc.get_num_bits().unwrap(), 10);
    assert!(!vc.get_bit(0).unwrap());
    assert!(vc.get_bit(1).unwrap());
    assert!(!vc.get_bit(9).unwrap());
    assert!(vc.get_bit(6).unwrap());
    assert_eq!(vc.get_on_bits().unwrap(), [1, 3, 5, 6, 8]);
    // Keys are unique: adding one again replaces its vector.
    vc.add_vect(1, bits(10, &[7, 9]));
    assert!(!vc.get_bit(1).unwrap());
    assert!(vc.get_bit(9).unwrap());
    assert_eq!(vc.num_children(), 2);
    let cs = vc.get_children();
    assert_eq!(cs[0].0, 1);
    assert_eq!(cs[0].1.get_on_bits(), [7, 9]);
    assert_eq!(cs[1].0, 2);
    assert_eq!(cs[1].1.get_on_bits(), [6, 8]);
    vc.add_vect(3, bits(10, &[5, 6]));
    assert_eq!(vc.num_children(), 3);
    assert_eq!(vc.get_on_bits().unwrap(), [5, 6, 7, 8, 9]);
    vc.detach_vects_not_matching_bit(6).unwrap();
    assert_eq!(vc.num_children(), 2);
    assert_eq!(vc.get_on_bits().unwrap(), [5, 6, 8]);
    vc.add_vect(1, bits(10, &[7, 9]));
    assert_eq!(vc.num_children(), 3);
    assert_eq!(vc.get_on_bits().unwrap(), [5, 6, 7, 8, 9]);
    vc.detach_vects_matching_bit(6).unwrap();
    assert_eq!(vc.num_children(), 1);
    assert_eq!(vc.get_on_bits().unwrap(), [7, 9]);
    // Copies are independent.
    vc.add_vect(3, bits(10, &[5, 6]));
    assert_eq!(vc.get_on_bits().unwrap(), [5, 6, 7, 9]);
    let mut vc2 = vc.clone();
    vc.detach_vects_not_matching_bit(6).unwrap();
    assert_eq!(vc.get_on_bits().unwrap(), [5, 6]);
    assert_eq!(vc2.get_on_bits().unwrap(), [5, 6, 7, 9]);
    assert_eq!(vc2.or_vect().unwrap().unwrap().get_on_bits(), [5, 6, 7, 9]);
    // Uniquify removes duplicates.
    let mut vc = VectCollection::default();
    vc.add_vect(1, bits(10, &[7, 9]));
    vc.add_vect(2, bits(10, &[7, 9]));
    vc.add_vect(3, bits(10, &[2, 3, 5]));
    assert_eq!(vc.num_children(), 3);
    vc.uniquify();
    assert_eq!(vc.num_children(), 2);
    assert_eq!(vc.get_children()[1].0, 3);
    vc.reset().unwrap();
    assert!(matches!(vc.get_bit(10), Err(Error::Index { .. })));
    vc.add_vect(4, bits(11, &[]));
    assert!(matches!(vc.get_on_bits(), Err(Error::Value { .. })));
    // With every vector gone, the size is the last one seen.
    let mut vc = VectCollection::default();
    vc.add_vect(1, bits(10, &[1]));
    vc.get_on_bits().unwrap();
    vc.detach_vects_matching_bit(1).unwrap();
    assert_eq!(vc.get_num_bits().unwrap(), 10);
    assert!(matches!(vc.get_bit(1), Err(Error::Value { .. })));
}

fn fp(smiles: &str) -> ExplicitBitVect {
    rdk_fingerprint(&mol_from_smiles(smiles).unwrap()).unwrap()
}

fn fp_sized(smiles: &str, fp_size: u32, tgt_density: f64) -> ExplicitBitVect {
    rdk_fingerprint_with_params(
        &mol_from_smiles(smiles).unwrap(),
        &RDKFingerprintParams {
            fp_size,
            tgt_density,
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn fingerprint_similarity_matches_python() {
    let fps = [fp("CCOC"), fp("CCO"), fp("COC")];
    assert_eq!(fingerprint_similarity(&fps[0], &fps[1]).unwrap(), 0.6);
    let folded: Vec<ExplicitBitVect> = ["CCOC", "CCO", "COC"]
        .iter()
        .enumerate()
        .map(|(i, s)| fp_sized(s, 16 * (i as u32 + 1), 0.0))
        .collect();
    assert_eq!(
        fingerprint_similarity(&folded[0], &folded[1]).unwrap(),
        0.5555555555555556
    );
    assert_eq!(
        fingerprint_similarity(&folded[1], &folded[0]).unwrap(),
        0.5555555555555556
    );
    assert_eq!(
        fingerprint_similarity_with_metric(&folded[2], &folded[0], dice_similarity).unwrap(),
        0.6153846153846154
    );
    assert!(matches!(
        fingerprint_similarity(&fps[0], &ExplicitBitVect::new(0)),
        Err(Error::Value { .. })
    ));

    let dense = fp_sized("CCOC", 1024, 0.8);
    assert_eq!(dense.get_num_bits(), 128);
    let same = fold_to_target_density_with_params(
        &dense,
        &FoldToTargetDensityParams {
            density: 0.1,
            min_length: 2,
        },
    )
    .unwrap();
    assert_eq!(
        same.get_on_bits(),
        [21, 28, 31, 40, 50, 52, 59, 64, 108, 121]
    );
    let big = fp("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
    let folded = fold_to_target_density(&big).unwrap();
    assert_eq!(
        (folded.get_num_bits(), folded.get_num_on_bits()),
        (2048, 228)
    );
    let sparse_fold = fold_to_target_density_with_params(
        &bits(256, &(0..200).collect::<Vec<u32>>()),
        &FoldToTargetDensityParams::default(),
    )
    .unwrap();
    assert_eq!(sparse_fold.get_num_bits(), 128);

    let names: Vec<&str> = similarity_functions::<ExplicitBitVect>()
        .iter()
        .map(|f| f.0)
        .collect();
    assert_eq!(
        names,
        [
            "Tanimoto",
            "Dice",
            "Cosine",
            "Sokal",
            "Russel",
            "RogotGoldberg",
            "AllBit",
            "Kulczynski",
            "McConnaughey",
            "Asymmetric",
            "BraunBlanquet"
        ]
    );
    let vals: Vec<f64> = similarity_functions()
        .iter()
        .map(|f| (f.1)(&fps[0], &fps[1]).unwrap())
        .collect();
    assert_eq!(
        vals,
        [
            0.6,
            0.75,
            0.7745966692414834,
            0.42857142857142855,
            0.0029296875,
            0.8745098039215686,
            0.998046875,
            0.8,
            0.6,
            1.0,
            0.6
        ]
    );
}
