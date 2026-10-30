// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::data_manip::metric::get_euclidean_dist_mat;
use rdkit_cxx::data_structs::ExplicitBitVect;
use rdkit_cxx::ml::cluster::clustering::murtagh_dist_cluster;
use rdkit_cxx::sim_div_filters::*;

// Expected values come from Python's rdSimDivPickers.

fn points() -> Vec<Vec<f64>> {
    [
        [0.0, 0.0],
        [1.0, 0.2],
        [5.0, 5.0],
        [5.5, 4.0],
        [9.0, 1.0],
        [0.5, 8.0],
        [2.0, 2.5],
        [7.0, 7.5],
        [3.0, 9.0],
        [8.5, 3.0],
        [4.0, 0.5],
        [6.0, 2.0],
    ]
    .iter()
    .map(|p| p.to_vec())
    .collect()
}

fn matrix() -> Vec<f64> {
    get_euclidean_dist_mat(&points()).unwrap()
}

fn lookup(d: &[f64], i: u32, j: u32) -> f64 {
    let (i, j) = (i.max(j) as usize, i.min(j) as usize);
    if i == j { 0.0 } else { d[i * (i - 1) / 2 + j] }
}

fn ebv(bits: &[u32]) -> ExplicitBitVect {
    let mut bv = ExplicitBitVect::new(64);
    bv.set_bits_from_list(bits).unwrap();
    bv
}

fn fingerprints() -> Vec<ExplicitBitVect> {
    let bits: [&[u32]; 20] = [
        &[3, 4, 6, 19, 23, 25, 34, 41, 52],
        &[4, 11, 27, 55],
        &[3, 5, 7, 8, 14, 15, 27, 35, 36, 52],
        &[5, 6, 28, 50],
        &[7, 9, 11, 17, 18, 19, 26, 34, 35, 36, 43, 52],
        &[7, 8, 12, 24, 47],
        &[20, 27, 29, 34, 43, 49, 63],
        &[5, 11, 15, 19, 36, 44, 46, 49, 50, 61, 62],
        &[4, 7, 10, 18, 21, 26, 28, 32, 38, 46, 56, 63],
        &[2, 4, 19, 26, 31, 35, 42, 48, 59],
        &[4, 22, 29, 31, 37, 38, 43, 44, 51],
        &[7, 8, 34, 39, 60],
        &[1, 10, 22, 24, 29, 36, 39, 42, 45, 56, 58],
        &[7, 16, 27, 36, 63],
        &[5, 10, 25, 31, 50, 55, 58],
        &[8, 17, 27, 35, 45, 51, 52, 55, 56, 61, 62],
        &[5, 9, 11, 14, 24, 43, 45, 56, 58, 59],
        &[1, 11, 16, 18, 31, 37, 53],
        &[18, 40, 47, 53],
        &[6, 29, 49, 55, 57, 60],
    ];
    bits.iter().map(|b| ebv(b)).collect()
}

type Clusters<'a> = &'a [&'a [u32]];

fn seeded(seed: i32) -> MaxMinPickParams {
    MaxMinPickParams {
        seed,
        ..Default::default()
    }
}

#[test]
fn defaults() {
    assert_eq!(
        MaxMinPickParams::default(),
        MaxMinPickParams {
            first_picks: vec![],
            seed: -1,
        }
    );
    assert_eq!(
        LeaderPickParams::default(),
        LeaderPickParams {
            pick_size: 0,
            first_picks: vec![],
        }
    );
}

#[test]
fn max_min() {
    let picker = MaxMinPicker::default();
    let d = matrix();
    assert_eq!(
        picker.pick_with_params(&d, 12, 5, &seeded(42)).unwrap(),
        [4, 5, 0, 7, 10]
    );
    let first = MaxMinPickParams {
        first_picks: vec![3, 7],
        ..Default::default()
    };
    assert_eq!(
        picker.pick_with_params(&d, 12, 5, &first).unwrap(),
        [3, 7, 0, 5, 4]
    );
    assert_eq!(picker.pick(&d, 12, 5).unwrap().len(), 5);
    let dist = |i, j| lookup(&d, i, j);
    assert_eq!(
        picker
            .lazy_pick_with_params(dist, 12, 5, &seeded(42))
            .unwrap(),
        [4, 5, 0, 7, 10]
    );
    // pick_size 0 picks one item, as in Python.
    assert_eq!(
        picker
            .lazy_pick_with_params(dist, 12, 0, &seeded(42))
            .unwrap(),
        [4]
    );
    let (picks, threshold) = picker
        .lazy_pick_with_threshold_with_params(dist, 12, 10, 3.0, &seeded(42))
        .unwrap();
    assert_eq!(picks, [4, 5, 0, 7, 10, 3]);
    assert!((threshold - 3.8078865529319543).abs() < 1e-12);
    let first = MaxMinPickParams {
        first_picks: vec![1, 2],
        ..Default::default()
    };
    assert_eq!(
        picker
            .lazy_pick_with_threshold_with_params(dist, 12, 2, 3.0, &first)
            .unwrap(),
        (vec![1, 2], -1.0)
    );
    assert_eq!(picker.lazy_pick(dist, 12, 4).unwrap().len(), 4);
    assert_eq!(
        picker
            .lazy_pick_with_threshold(dist, 12, 4, 0.5)
            .unwrap()
            .0
            .len(),
        4
    );
}

#[test]
fn max_min_bit_vectors() {
    let picker = MaxMinPicker::default();
    let fps = fingerprints();
    assert_eq!(
        picker
            .lazy_bit_vector_pick_with_params(&fps, 20, 6, &seeded(3))
            .unwrap(),
        [11, 1, 3, 18, 8, 16]
    );
    let first = MaxMinPickParams {
        first_picks: vec![0],
        seed: 3,
    };
    assert_eq!(
        picker
            .lazy_bit_vector_pick_with_params(&fps, 15, 6, &first)
            .unwrap(),
        [0, 5, 7, 10, 1, 8]
    );
    let (picks, threshold) = picker
        .lazy_bit_vector_pick_with_threshold_with_params(&fps, 20, 20, 0.7, &seeded(3))
        .unwrap();
    assert_eq!(
        picks,
        [
            11, 1, 3, 18, 8, 16, 0, 6, 9, 10, 2, 7, 19, 4, 14, 17, 12, 15, 5, 13
        ]
    );
    assert_eq!(threshold, 0.75);
    assert_eq!(picker.lazy_bit_vector_pick(&fps, 20, 3).unwrap().len(), 3);
    assert_eq!(
        picker
            .lazy_bit_vector_pick_with_threshold(&fps, 20, 3, 0.1)
            .unwrap()
            .0
            .len(),
        3
    );
}

#[test]
fn leader() {
    let picker = LeaderPicker::default();
    let fps = fingerprints();
    assert_eq!(
        picker.lazy_bit_vector_pick(&fps, 20, 0.75).unwrap(),
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 19]
    );
    let params = LeaderPickParams {
        pick_size: 3,
        first_picks: vec![5],
    };
    assert_eq!(
        picker
            .lazy_bit_vector_pick_with_params(&fps, 20, 0.75, &params)
            .unwrap(),
        [5, 0, 1]
    );
    let d = matrix();
    let dist = |i, j| lookup(&d, i, j);
    assert_eq!(
        picker.lazy_pick(dist, 12, 3.0).unwrap(),
        [0, 2, 4, 5, 6, 7, 11]
    );
    let params = LeaderPickParams {
        pick_size: 2,
        first_picks: vec![4],
    };
    assert_eq!(
        picker
            .lazy_pick_with_params(dist, 12, 3.0, &params)
            .unwrap(),
        [4, 0]
    );
    // With a negative threshold no item removes itself. RDKit's state for
    // threaded Unix builds (Python's on Linux) then picks the first one again
    // and again; the other, every item in turn.
    let picks = picker.lazy_bit_vector_pick(&fps[..10], 10, -1.0).unwrap();
    assert!(
        picks == [0; 10] || picks == (0..10).collect::<Vec<_>>(),
        "{picks:?}"
    );
}

#[test]
fn hierarchical() {
    let d = matrix();
    let expected: [(ClusterMethod, &[u32], Clusters); 7] = [
        (
            ClusterMethod::Ward,
            &[1, 3, 11, 5],
            &[&[0, 1, 6], &[2, 3, 7], &[4, 9, 10, 11], &[5, 8]],
        ),
        (
            ClusterMethod::Slink,
            &[0, 2, 5, 7],
            &[&[0, 1, 6], &[2, 3, 11, 10, 4, 9], &[5, 8], &[7]],
        ),
        (
            ClusterMethod::Clink,
            &[1, 2, 11, 5],
            &[&[0, 1, 6], &[2, 3, 7], &[4, 9, 10, 11], &[5, 8]],
        ),
        (
            ClusterMethod::Upgma,
            &[1, 11, 5, 7],
            &[&[0, 1, 6], &[2, 3, 10, 11, 4, 9], &[5, 8], &[7]],
        ),
        (
            ClusterMethod::Mcquitty,
            &[1, 11, 5, 7],
            &[&[0, 1, 6], &[2, 3, 10, 11, 4, 9], &[5, 8], &[7]],
        ),
        (
            ClusterMethod::Gower,
            &[0, 11, 5, 7],
            &[&[0, 1, 6, 10], &[2, 3, 11, 4, 9], &[5, 8], &[7]],
        ),
        (
            ClusterMethod::Centroid,
            &[0, 11, 5, 7],
            &[&[0, 1, 6, 10], &[2, 3, 11, 4, 9], &[5, 8], &[7]],
        ),
    ];
    for (method, picks, clusters) in expected {
        let picker = HierarchicalClusterPicker::new(method);
        assert_eq!(picker.pick(&d, 12, 4).unwrap(), picks, "{method:?}");
        assert_eq!(picker.cluster(&d, 12, 4).unwrap(), clusters, "{method:?}");
    }
    let all = HierarchicalClusterPicker::new(ClusterMethod::Ward)
        .cluster(&d, 12, 12)
        .unwrap();
    assert_eq!(all, (0..12).map(|i| vec![i]).collect::<Vec<_>>());
    // The matrix is not changed, though hc_ changes its copy.
    assert_eq!(d, matrix());
}

#[test]
fn rdkit_errors() {
    let d = matrix();
    let dist = |i, j| lookup(&d, i, j);
    let mm = MaxMinPicker::default();
    let lp = LeaderPicker::default();
    let value = |r: Result<Vec<u32>, Error>| matches!(r, Err(Error::Value { .. }));
    assert!(value(mm.pick(&d, 12, 12)));
    assert!(value(mm.lazy_pick(dist, 0, 0)));
    assert!(value(mm.lazy_pick(dist, 5, 6)));
    let past = MaxMinPickParams {
        first_picks: vec![5],
        ..Default::default()
    };
    assert!(value(mm.lazy_pick_with_params(dist, 5, 2, &past)));
    let big = LeaderPickParams {
        pick_size: 6,
        ..Default::default()
    };
    assert!(value(lp.lazy_pick_with_params(dist, 5, 0.1, &big)));
    let past = LeaderPickParams {
        first_picks: vec![9],
        ..Default::default()
    };
    assert!(value(lp.lazy_pick_with_params(dist, 5, 0.1, &past)));
    assert!(value(lp.lazy_pick(dist, 0, 0.1)));
    let ward = HierarchicalClusterPicker::new(ClusterMethod::Ward);
    assert!(value(ward.pick(&d, 12, 12)));
    assert!(matches!(ward.cluster(&d, 3, 4), Err(Error::Runtime { .. })));
    // Python reads past the list (IndexError).
    let fps = fingerprints();
    assert!(matches!(
        mm.lazy_bit_vector_pick(&fps, 21, 3),
        Err(Error::Index { index: 20, len: 20 })
    ));
    assert!(matches!(
        lp.lazy_bit_vector_pick(&fps, 21, 0.5),
        Err(Error::Index { .. })
    ));
    // Bit vectors of different sizes, which on RDKit's threads abort in
    // Python.
    let mut mixed = fingerprints();
    mixed.push(ExplicitBitVect::new(128));
    assert!(value(lp.lazy_bit_vector_pick(&mixed, 21, 0.1)));
    assert!(value(mm.lazy_bit_vector_pick_with_params(
        &mixed,
        21,
        21,
        &seeded(1)
    )));
}

#[test]
fn max_min_guards() {
    let mm = MaxMinPicker::default();
    // Distances at or below -1, or NaN, leave RDKit with no pick (a null
    // dereference, a segfault in Python).
    for bad in [-2.0, -1.0, f64::NAN] {
        assert!(matches!(
            mm.lazy_pick(|_, _| bad, 10, 3),
            Err(Error::Value { .. })
        ));
        let d = vec![bad; 45];
        assert!(matches!(mm.pick(&d, 10, 3), Err(Error::Value { .. })));
    }
    // Distances above -1 are picked by, as in Python.
    assert_eq!(
        mm.lazy_pick_with_params(|i, j| -0.5 + 0.01 * (i + j) as f64, 10, 3, &seeded(1))
            .unwrap()
            .len(),
        3
    );
    // A matrix shorter than the pool's (read past its end in Python).
    assert!(matches!(
        mm.pick(&[1.0], 30000, 5),
        Err(Error::Value { .. })
    ));
    // getDistFromLTM's unsigned index wraps past 65537 items.
    assert!(matches!(mm.pick(&[], 65538, 1), Err(Error::Value { .. })));
    // The work array's size wraps in unsigned int.
    assert!(matches!(
        mm.lazy_pick(|_, _| 1.0, 1 << 28, 2),
        Err(Error::Value { .. })
    ));
}

#[test]
fn leader_guards() {
    // RDKit numbers the items in int.
    assert!(matches!(
        LeaderPicker::default().lazy_pick(|_, _| 1.0, u32::MAX, 0.5),
        Err(Error::Value { .. })
    ));
    // Any distance compares, as in Python.
    assert_eq!(
        LeaderPicker::default()
            .lazy_pick(|_, _| f64::NAN, 5, 0.5)
            .unwrap(),
        [0]
    );
}

#[test]
fn panics() {
    let calls = std::cell::Cell::new(0);
    let panicking = |_, _| -> f64 {
        calls.set(calls.get() + 1);
        if calls.get() == 3 {
            panic!("no distance");
        }
        1.0
    };
    let err = MaxMinPicker::default()
        .lazy_pick_with_params(panicking, 10, 5, &seeded(1))
        .unwrap_err();
    assert!(
        matches!(&err, Error::Runtime { message } if message.contains("no distance")),
        "{err:?}"
    );
    assert_eq!(calls.get(), 3);
    calls.set(0);
    let err = LeaderPicker::default()
        .lazy_pick(panicking, 10, 0.5)
        .unwrap_err();
    assert!(matches!(&err, Error::Runtime { message } if message.contains("no distance")));
    assert_eq!(calls.get(), 3);
}

#[test]
fn hierarchical_guards() {
    let d = matrix();
    let ward = HierarchicalClusterPicker::new(ClusterMethod::Ward);
    // No picks: RDKit reads cluster -1 (a segfault in Python).
    assert!(matches!(ward.cluster(&d, 4, 0), Err(Error::Value { .. })));
    assert!(matches!(ward.pick(&d, 1, 0), Err(Error::Value { .. })));
    // Fewer than 2 items, which Murtagh's code cannot cluster.
    assert!(matches!(ward.cluster(&[], 1, 1), Err(Error::Value { .. })));
    // A short matrix.
    assert!(matches!(ward.cluster(&d, 13, 2), Err(Error::Value { .. })));
    assert!(matches!(ward.pick(&d, 13, 2), Err(Error::Value { .. })));
    // Distances that are not finite.
    for bad in [f64::INFINITY, f64::NAN] {
        let mut d = matrix();
        d[5] = bad;
        assert!(matches!(ward.pick(&d, 12, 4), Err(Error::Value { .. })));
        assert!(matches!(ward.cluster(&d, 12, 4), Err(Error::Value { .. })));
    }
    // Distances past Murtagh's 1e20 sentinel cluster as scaled ones.
    let scale = 2f64.powi(100);
    let large: Vec<f64> = d.iter().map(|x| x * scale).collect();
    for method in [
        ClusterMethod::Ward,
        ClusterMethod::Centroid,
        ClusterMethod::Slink,
    ] {
        let picker = HierarchicalClusterPicker::new(method);
        assert_eq!(
            picker.cluster(&large, 12, 4).unwrap(),
            picker.cluster(&d, 12, 4).unwrap()
        );
        assert_eq!(
            picker.pick(&large, 12, 4).unwrap(),
            picker.pick(&d, 12, 4).unwrap()
        );
    }
}

#[test]
fn threads() {
    // The hierarchical picker shares Murtagh's lock with Clustering.
    let d = matrix();
    let picks = HierarchicalClusterPicker::new(ClusterMethod::Ward)
        .pick(&d, 12, 4)
        .unwrap();
    let merges = murtagh_dist_cluster(&d, 12, 2).unwrap();
    std::thread::scope(|scope| {
        for t in 0..8 {
            let (d, picks, merges) = (&d, &picks, &merges);
            scope.spawn(move || {
                for _ in 0..100 {
                    if t % 2 == 0 {
                        let picker = HierarchicalClusterPicker::new(ClusterMethod::Ward);
                        assert_eq!(&picker.pick(d, 12, 4).unwrap(), picks);
                    } else {
                        assert_eq!(&murtagh_dist_cluster(d, 12, 2).unwrap(), merges);
                    }
                }
            });
        }
    });
}
