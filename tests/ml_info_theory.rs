// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports the BitClusterer test of rdkit/ML/InfoTheory/UnitTestCorrMatGen.py,
// with values from Python's rdkit.ML.InfoTheory modules in RDKit 2026.09.1.
// Gains come from RDKit's C++ `InfoGain`, whose rounding differs across
// platforms (arm64's build differs in the last bits), so they are compared
// to within rounding.
// Python's BitRank functions fail with NumPy 2, so their values come from
// Python's code with `numpy.integer` read as `numpy.int64`.

mod common;

use common::assert_close;
use rdkit_cxx::Error;
use rdkit_cxx::data_structs::{ExplicitBitVect, SparseBitVect};
use rdkit_cxx::ml::data::data_utils::init_random_numbers;
use rdkit_cxx::ml::info_theory::bit_clusterer::BitClusterer;
use rdkit_cxx::ml::info_theory::bit_rank::*;
use rdkit_cxx::ml::info_theory::entropy::*;
use rdkit_cxx::ml::info_theory::{BitCorrMatGenerator, info_gain};
use rdkit_cxx::rd_random;

fn assert_all_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &e) in actual.iter().zip(expected) {
        assert_close(a, e);
    }
}

#[test]
fn entropy_matches_python() {
    // NumPy's log can round differently from the C library's.
    assert_close(py_info_entropy(&[5.0, 6.0, 1.0]), 1.3250112108241772);
    assert_eq!(py_info_entropy(&[0.0, 0.0]), 0.0);
    assert_close(py_info_entropy(&[3.0, 0.0, 7.0, 2.5]), 1.4269608142719166);
    let table = vec![vec![5.0, 1.0], vec![2.0, 7.0], vec![0.0, 3.0]];
    assert_close(py_info_gain(&table).unwrap(), 0.3653023710044675);
    assert_close(info_gain(&table).unwrap(), 0.3653023710044674);
    assert_eq!(
        py_info_gain(&[vec![0.0, 0.0], vec![0.0, 0.0]]).unwrap(),
        0.0
    );
    assert!(matches!(
        py_info_gain(&[vec![1.0], vec![1.0, 2.0]]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn bit_rank_matches_python() {
    let vects: Vec<Vec<u32>> = (0..12)
        .map(|i| (0..6).map(|j| (i * j * 7 + j) % 3 % 2).collect())
        .collect();
    let acts: Vec<u32> = (0..12).map(|i| i % 3 % 2).collect();
    assert_eq!(form_counts(&vects, &acts, 2, 2).unwrap(), [[8, 0], [0, 4]]);
    let gains = [
        0.0,
        0.2516291673878228,
        0.9182958340544894,
        0.0,
        0.2516291673878228,
        0.9182958340544894,
    ];
    assert_all_close(&calc_info_gains(&vects, &acts, 2).unwrap(), &gains);
    let (order, metrics) = rank_bits(&vects, &acts).unwrap();
    assert_eq!(order, [5, 2, 4, 1, 3, 0]);
    assert_all_close(&metrics, &gains);
    assert!(matches!(
        form_counts(&vects, &acts[1..], 2, 2),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        form_counts(&vects, &acts, 2, 1),
        Err(Error::Index { .. })
    ));

    let sv: Vec<SparseBitVect> = (0..10)
        .map(|i| {
            let mut s = SparseBitVect::new(20);
            for b in 0..20 {
                if (i * 3 + b * 5) % 7 < 2 {
                    s.set_bit(b).unwrap();
                }
            }
            s
        })
        .collect();
    let sacts: Vec<u32> = (0..10).map(|i| i % 2).collect();
    let (res, gains) = analyze_sparse_vects(&sv, &sacts).unwrap();
    let (a, b, c) = (0.01437846047807878, 0.10803154614559995, 0.0);
    let expected_gains = [a, a, b, c, b, a, a, a, a, b, c, b, a, a, a, a, b, c, b, a];
    assert_all_close(&gains, &expected_gains);
    for (bit, gain, n_act, n_inact) in [(2usize, b, 0, 2), (4, b, 2, 0), (5, a, 1, 2)] {
        assert_eq!(
            (res[bit].0, res[bit].2, res[bit].3),
            (bit as u32, n_act, n_inact)
        );
        assert_close(res[bit].1, gain);
    }
    let (order, info) = sparse_rank_bits(&sv, &sacts).unwrap();
    assert_eq!(info, res);
    // Python's ranking, whose ties NumPy orders as it likes:
    // [18, 16, 11, 9, 2, 4, 5, 19, 12, 13, 15, 14, 0, 8, 6, 7, 1, 17, 10, 3].
    // The six tied gains can round apart in RDKit's C++ on some platforms,
    // so only their set is pinned.
    let mut top: Vec<usize> = order[..6].to_vec();
    top.sort_unstable();
    assert_eq!(top, [2, 4, 9, 11, 16, 18]);
    assert!(order.windows(2).all(|w| gains[w[0]] >= gains[w[1]]));
}

#[test]
fn bit_clusterer() {
    // Python's setUp: the second half of each fingerprint copies the first.
    let (nbits, d, nfp) = (200u32, 40, 1000);
    init_random_numbers((100, 23));
    let fps: Vec<ExplicitBitVect> = (0..nfp)
        .map(|_| {
            let mut fp = ExplicitBitVect::new(nbits);
            let mut obits: Vec<u32> = (0..nbits / 2).collect();
            rd_random::shuffle_with_random(&mut obits, rd_random::random).unwrap();
            for &bit in &obits[..d] {
                fp.set_bit(bit).unwrap();
                fp.set_bit(bit + nbits / 2).unwrap();
            }
            fp
        })
        .collect();
    let blist: Vec<u32> = (0..nbits).collect();
    let mut cmg = BitCorrMatGenerator::default();
    cmg.set_bit_list(&blist).unwrap();
    for fp in &fps {
        cmg.collect_votes(fp).unwrap();
    }
    let corr_mat = cmg.get_corr_matrix();

    // Python's test0CorrMat, which needs Python's random numbers.
    let val = |i: usize, j: usize| corr_mat[j * (j - 1) / 2 + i];
    let half = nbits as usize / 2;
    let avr: f64 = (0..half).map(|i| val(i, i + half)).sum();
    let navr: f64 = (0..half).map(|i| val(i, i + 1)).sum();
    assert_eq!(2.0 * avr / nbits as f64, 400.0);
    assert_eq!(2.0 * navr / nbits as f64, 158.3);

    let mut bcl = BitClusterer::new(blist, half);
    bcl.cluster_bits(&corr_mat).unwrap();
    let cls = bcl.get_clusters().to_vec();
    for cl in &cls {
        assert_eq!(cl.len(), 2);
        assert_eq!(cl[0] + nbits / 2, cl[1]);
    }
    bcl.set_clusters(cls.clone()).unwrap();
    assert!(matches!(
        bcl.set_clusters(cls[..cls.len() - 1].to_vec()),
        Err(Error::Value { .. })
    ));

    let mut tfp = ExplicitBitVect::new(nbits);
    let obits: Vec<u32> = (0..nbits / 4).chain(nbits / 2..3 * nbits / 4).collect();
    tfp.set_bits_from_list(&obits).unwrap();
    let rvc = bcl.map_to_cluster_scores(&tfp).unwrap();
    assert_eq!(rvc.len(), half);
    for (i, &score) in rvc.iter().enumerate() {
        assert_eq!(score, if i < nbits as usize / 4 { 2 } else { 0 });
    }
    let nfp = bcl.map_to_cluster_fp(&tfp).unwrap();
    assert_eq!(nfp.get_num_bits() as usize, half);
    for i in 0..half {
        assert_eq!(nfp.get_bit(i as u32).unwrap(), i < nbits as usize / 4);
    }
}
