// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/ML/Cluster/UnitTestButina.py, UnitTestCluster.py, and
// murtagh_test.py, with values from Python's rdkit.ML.Cluster modules in
// RDKit 2026.09.1.

// Python's test data, 3.14, is not pi.
#![allow(clippy::approx_constant)]

use rdkit_cxx::Error;
use rdkit_cxx::ml::cluster::butina::{self, ClusterDataParams};
use rdkit_cxx::ml::cluster::cluster_utils::*;
use rdkit_cxx::ml::cluster::clusters::{Cluster, ClusterParams, PrintParams, cmp, cmptol};
use rdkit_cxx::ml::cluster::{murtagh, resemblance, standardize};

fn flat(data: &[f64], n_pts: usize, thresh: f64, reordering: bool) -> Vec<Vec<usize>> {
    butina::cluster_data_flat_with_reordering(data, n_pts, thresh, reordering).unwrap()
}

#[test]
fn butina_unit_tests() {
    let dists = [1., 2., 1., 4., 3., 2., 6., 5., 4., 2., 7., 6., 5., 3., 1.];
    assert_eq!(
        butina::cluster_data_flat(&dists, 6, 1.1).unwrap(),
        [vec![1, 0, 2], vec![5, 4], vec![3]]
    );
    let dists = [
        0.5, 1., 0.5, 2., 1.5, 1., 3., 2.5, 2., 1., 5., 4.5, 4., 3., 2., 8., 7.5, 7., 6., 5., 3.,
        9., 8.5, 8., 7., 6., 4., 1.,
    ];
    assert_eq!(
        flat(&dists, 8, 2.1, false),
        [vec![3, 0, 1, 2, 4], vec![7, 6], vec![5]]
    );
    // Everything a singleton.
    assert_eq!(
        flat(&[1., 2., 1.], 3, 0.9, false),
        [vec![2], vec![1], vec![0]]
    );
    // Everything in one cluster.
    assert_eq!(
        flat(&[1., 2., 1., 3., 2., 1.], 4, 2.0, false),
        [vec![2, 0, 1, 3]]
    );
    // One in the middle leaves the edges lonely.
    let dists = [1.5, 2.5, 1., 3.5, 2., 1., 5., 3.5, 2.5, 1.5];
    assert_eq!(
        flat(&dists, 5, 1.1, false),
        [vec![2, 1, 3], vec![4], vec![0]]
    );
    // Zero distances.
    let dists = [1., 2., 0., 2., 0., 0., 4., 2., 2., 2.];
    assert_eq!(
        flat(&dists, 5, 0.9, false),
        [vec![3, 1, 2], vec![4], vec![0]]
    );
    // Reordering, without changes.
    let dists = [1., 2., 1., 4., 3., 2., 6., 5., 4., 2., 7., 6., 5., 3., 1.];
    assert_eq!(
        flat(&dists, 6, 1.1, true),
        [vec![1, 0, 2], vec![5, 4], vec![3]]
    );
    // Reordering, with changes.
    let dists = [
        2., 3.5, 1.5, 5., 3., 1.5, 7., 5., 3.5, 2., 8., 6., 4.5, 3., 1., 9., 7., 5.5, 4., 2., 1.,
    ];
    assert_eq!(
        flat(&dists, 7, 2.1, false),
        [vec![4, 3, 5, 6], vec![2, 1], vec![0]]
    );
    assert_eq!(
        flat(&dists, 7, 2.1, true),
        [vec![4, 3, 5, 6], vec![1, 0, 2]]
    );
    // Empty input.
    assert!(flat(&[], 0, 2.0, true).is_empty());
    // Mismatched input.
    let err = butina::cluster_data_flat(&[1., 2., 3.], 1, 2.0).unwrap_err();
    assert_eq!(
        err,
        Error::Value {
            message: "Mismatched input data dimension and nPts".to_owned()
        }
    );
}

fn points() -> Vec<Vec<f64>> {
    (0..30)
        .map(|i| {
            vec![
                (i * 37 % 17) as f64 / 3.0,
                (i * 11 % 13) as f64 / 2.0,
                (i * 5 % 7) as f64 * 1.25,
            ]
        })
        .collect()
}

#[test]
fn butina_points_match_python() {
    let pts = points();
    let expected: Vec<Vec<usize>> = vec![
        vec![24, 3, 6, 13, 17],
        vec![21, 10, 14, 28],
        vec![9, 2, 16, 20, 27],
        vec![4, 11, 15, 22, 25],
        vec![29, 1, 8, 18],
        vec![26, 5, 19],
        vec![23, 12],
        vec![7],
        vec![0],
    ];
    assert_eq!(butina::cluster_data(&pts, 30, 2.5).unwrap(), expected);
    let reordered = butina::cluster_data_with_params(
        &pts,
        30,
        2.5,
        &ClusterDataParams {
            reordering: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        reordered,
        [
            vec![24, 3, 6, 13, 17],
            vec![9, 2, 16, 20, 27],
            vec![4, 11, 15, 22, 25],
            vec![29, 1, 8, 18],
            vec![21, 10, 14, 28],
            vec![26, 5, 19],
            vec![23, 12],
            vec![7],
            vec![0],
        ]
    );
    let m = butina::compute_distance_matrix(&pts, 30, |a, b| butina::euclidean_dist(a, b)).unwrap();
    assert_eq!(m[7][3], 3.2542706982944396);
    assert_eq!(m[29][0], 7.194229941032219);
    assert_eq!(
        butina::euclidean_dist(&[1.5, 2.0], &[0.25, -3.0]).unwrap(),
        5.153882032022076
    );
    let from_matrix = butina::cluster_data_with_params(
        &m,
        30,
        2.5,
        &ClusterDataParams {
            is_dist_data: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(from_matrix, expected);
    let with_func = butina::cluster_data_with_dist_func(
        &pts,
        30,
        2.5,
        |a: &Vec<f64>, b: &Vec<f64>| butina::euclidean_dist(a, b),
        false,
    )
    .unwrap();
    assert_eq!(with_func, expected);
    assert!(matches!(
        butina::cluster_data(&pts, 31, 2.5),
        Err(Error::Index { .. })
    ));
    assert!(matches!(
        butina::euclidean_dist(&[1.0, 2.0], &[1.0, 2.0, 3.0]),
        Err(Error::Value { .. })
    ));
}

fn leaf(index: i64) -> Cluster {
    Cluster::new_with_params(&ClusterParams {
        index,
        ..Default::default()
    })
}

fn node(metric: f64, children: Vec<Cluster>, index: i64) -> Cluster {
    Cluster::new_with_params(&ClusterParams {
        metric,
        children,
        index,
        ..Default::default()
    })
}

fn indices(clusters: &[&Cluster]) -> Vec<i64> {
    clusters.iter().map(|c| c.get_index()).collect()
}

#[test]
fn divide() {
    let tree = || {
        let c1 = node(10.0, vec![leaf(1), leaf(2)], 7);
        let c2 = node(15.0, vec![leaf(3), leaf(4)], 8);
        let c3 = node(20.0, vec![leaf(5), leaf(6)], 9);
        let c4 = node(25.0, vec![c2, c3], 10);
        node(30.0, vec![c4, c1], 11)
    };
    let mut c5 = tree();
    let cs = split_into_n_clusters(&mut c5, 4).unwrap();
    assert_eq!(indices(&cs), [9, 8, 1, 2]);
    let cs2 = split_into_n_clusters_with_breadth_first(&mut c5, 4, false).unwrap();
    assert_eq!(indices(&cs2), [8, 7, 5, 6]);

    for breadth_first in [true, false] {
        let mut c5 = tree();
        let n = c5.len();
        assert!(matches!(
            split_into_n_clusters_with_breadth_first(&mut c5, n + 1, breadth_first),
            Err(Error::Value { .. })
        ));
        let points = indices(&c5.get_points());
        let all = split_into_n_clusters_with_breadth_first(&mut c5, n, breadth_first).unwrap();
        assert_eq!(indices(&all), points);
        let none = split_into_n_clusters_with_breadth_first(&mut c5, 0, breadth_first).unwrap();
        assert_eq!(indices(&none), [11]);
        for n in 0..c5.len() {
            let res = split_into_n_clusters_with_breadth_first(&mut c5, n, breadth_first);
            // Python's code fails (an assertion) for 7 and above.
            assert_eq!(res.is_err(), n >= 7, "{n}");
        }
    }
}

#[test]
fn cluster_class() {
    let mut root = Cluster::new_with_params(&ClusterParams {
        index: 1,
        position: vec![1.0],
        ..Default::default()
    });
    let mut c1 = Cluster::new_with_params(&ClusterParams {
        index: 10,
        position: vec![10.0],
        ..Default::default()
    });
    for i in [30, 31] {
        c1.add_child(Cluster::new_with_params(&ClusterParams {
            index: i,
            position: vec![i as f64],
            ..Default::default()
        }));
    }
    c1.add_child(Cluster::new_with_params(&ClusterParams {
        index: 32,
        position: vec![32.0],
        ..Default::default()
    }));
    assert!(c1.find_subtree(32).unwrap().is_terminal());
    let mut c2 = leaf(11);
    c2.add_children(vec![leaf(40), leaf(41)]);
    root.add_child(c1);
    root.add_child(c2);
    assert_eq!(
        indices(&get_node_list(&mut root)),
        [30, 31, 32, 10, 40, 41, 11, 1]
    );
    assert_eq!(
        indices(&get_nodes_down_to_centroids(&mut root)),
        [30, 31, 32, 10, 40, 41, 11, 1]
    );
    let mut subtree = root.find_subtree(11).unwrap().clone();
    assert_eq!(indices(&get_node_list(&mut subtree)), [40, 41, 11]);
    assert!(!root.is_terminal());
    assert!(root.find_subtree(99).is_none());

    assert_eq!(root.get_data(), None);
    root.set_data(Some(3.14));
    assert_eq!(root.get_data(), Some(3.14));
    assert_eq!(root.get_metric(), 0.0);
    root.set_metric(0.1);
    assert_eq!(root.get_metric(), 0.1);
    assert_eq!(root.get_index(), 1);
    root.set_index(100);
    assert_eq!(root.get_index(), 100);
    // Python gives [30, 31, 32, []]: cluster 11, given its children in one
    // AddChildren call while a leaf, caches itself as its only point. The
    // port gathers the points afresh, so 40 and 41 come back.
    let empty: &[f64] = &[];
    assert_eq!(
        root.get_points_positions(),
        [&[30.0][..], &[31.0], &[32.0], empty, empty]
    );

    let removed = root.remove_child(0).unwrap();
    assert_eq!(removed.get_index(), 10);
    assert_eq!(indices(&get_node_list(&mut root)), [40, 41, 11, 100]);
    assert!(matches!(root.remove_child(1), Err(Error::Index { .. })));
    assert_eq!(root.get_name(), "Cluster(100)");
    root.set_name(Some("abc".to_owned()));
    assert_eq!(root.get_name(), "abc");
    let text = root.print_with_params(&PrintParams {
        show_data: true,
        ..Default::default()
    });
    assert!(text.contains("abc"));
    assert!(text.contains("Cluster(41)"));
    assert!(text.contains("Metric"));
}

#[test]
fn print_matches_python() {
    let mut root = Cluster::new_with_params(&ClusterParams {
        index: 1,
        ..Default::default()
    });
    let mut c1 = leaf(10);
    c1.add_child(Cluster::new_with_params(&ClusterParams {
        index: 30,
        metric: 2.0,
        ..Default::default()
    }));
    c1.add_child(Cluster::new_with_params(&ClusterParams {
        index: 31,
        metric: 1.0,
        ..Default::default()
    }));
    c1.add_child(leaf(32));
    let mut c2 = node(0.5, Vec::new(), 11);
    c2.add_children(vec![leaf(40), leaf(41)]);
    root.add_child(c1);
    root.add_child(c2);
    root.set_data(Some(3.14));
    let text = root.print_with_params(&PrintParams {
        show_data: true,
        ..Default::default()
    }) + &root.print_with_params(&PrintParams {
        level: 1,
        offset: "--".to_owned(),
        ..Default::default()
    });
    assert_eq!(
        text,
        "Cluster(1)\t Data: 3.140000\t Metric: 0.000000\n  Cluster(10)\t Metric: 0.000000\n    \
         Cluster(32)\t Metric: 0.000000\n    Cluster(31)\t Metric: 1.000000\n    Cluster(30)\t \
         Metric: 2.000000\n  Cluster(11)\t Metric: 0.500000\n    Cluster(40)\t Metric: 0.000000\n    \
         Cluster(41)\t Metric: 0.000000\n  Cluster(1)-- Metric: 0.000000\n    Cluster(10)-- Metric: \
         0.000000\n      Cluster(32)-- Metric: 0.000000\n      Cluster(31)-- Metric: 1.000000\n      \
         Cluster(30)-- Metric: 2.000000\n    Cluster(11)-- Metric: 0.500000\n      Cluster(40)-- \
         Metric: 0.000000\n      Cluster(41)-- Metric: 0.000000\n"
    );
    let children: Vec<i64> = root.get_children().iter().map(Cluster::get_index).collect();
    assert_eq!(children, [10, 11]);
}

/// Romesburg's data from "Cluster Analysis for Researchers".
fn romesburg() -> (Vec<Vec<f64>>, Vec<f64>) {
    let d = vec![
        vec![10., 5.],
        vec![20., 20.],
        vec![30., 10.],
        vec![30., 15.],
        vec![5., 10.],
    ];
    let mut dists = Vec::new();
    for i in 0..d.len() {
        for j in 0..i {
            dists.push(
                d[i].iter()
                    .zip(&d[j])
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>(),
            );
        }
    }
    (d, dists)
}

/// Python's `Print(showData=True)` of each method's tree, which is the same
/// from the points and from their squared distances; the trees differ only
/// in their metrics.
#[test]
fn murtagh_matches_python() {
    let (d, dists) = romesburg();
    let template = |top: &str, left: &str, right: &str| {
        format!(
            "Cluster(9)\t Metric: {top}\n  Cluster(7)\t Metric: 50.000000\n    Cluster(1)\t Data: \
             0.000000\t Metric: 0.000000\n    Cluster(5)\t Data: 4.000000\t Metric: 0.000000\n  \
             Cluster(8)\t Metric: {left}\n    Cluster(2)\t Data: 1.000000\t Metric: 0.000000\n    \
             Cluster(6)\t Metric: {right}\n      Cluster(3)\t Data: 2.000000\t Metric: 0.000000\n      \
             Cluster(4)\t Data: 3.000000\t Metric: 0.000000\n"
        )
    };
    let cases = [
        (
            murtagh::slink(),
            template("325.000000", "125.000000", "25.000000"),
        ),
        (
            murtagh::clink(),
            template("650.000000", "200.000000", "25.000000"),
        ),
        (
            murtagh::upgma(),
            template("475.000000", "162.500000", "25.000000"),
        ),
        (
            murtagh::centroid(),
            template("423.611111", "156.250000", "25.000000"),
        ),
        (
            murtagh::mcquitty(),
            template("437.500000", "162.500000", "25.000000"),
        ),
        (
            murtagh::gower(),
            template("382.812500", "156.250000", "25.000000"),
        ),
        (
            murtagh::wards(),
            template("508.333333", "104.166667", "12.500000").replace("50.000000", "25.000000"),
        ),
    ];
    let show_data = PrintParams {
        show_data: true,
        ..Default::default()
    };
    for (method, expected) in cases {
        let mut from_dists = murtagh::cluster_data_dists(&dists, 5, method).unwrap();
        assert_eq!(from_dists.len(), 1);
        assert_eq!(
            from_dists[0].print_with_params(&show_data),
            expected,
            "{method}"
        );
        let mut from_points = murtagh::cluster_data(&d, 5, method).unwrap();
        assert_eq!(
            from_points[0].print_with_params(&show_data),
            expected,
            "{method}"
        );
        assert_eq!(
            from_dists[0].compare_with_ignore_extras(&mut from_points[0], false),
            0
        );
    }
    assert_eq!(
        murtagh::methods().iter().map(|m| m.1).collect::<Vec<_>>(),
        [1, 4, 2, 3, 7]
    );
}

#[test]
fn murtagh_upgma_from_points_and_dists() {
    // Python's testMurtaghUPGMA, with fixed points instead of random ones.
    let pts: Vec<Vec<f64>> = (0..5)
        .map(|i| {
            (0..5)
                .map(|j| ((i * 7 + j * 3) % 11) as f64 / 11.0)
                .collect()
        })
        .collect();
    let mut ds = Vec::new();
    for i in 0..5 {
        for j in 0..i {
            let d: Vec<f64> = pts[i].iter().zip(&pts[j]).map(|(a, b)| a - b).collect();
            ds.push(d.iter().map(|x| x * x).fold(0.0, |a, b| a + b));
        }
    }
    let mut c1 = murtagh::cluster_data(&pts, 5, murtagh::upgma())
        .unwrap()
        .remove(0);
    let mut c2 = murtagh::cluster_data_dists(&ds, 5, murtagh::upgma())
        .unwrap()
        .remove(0);
    assert_eq!(c1.len(), c2.len());
    assert_eq!(c1.compare_with_ignore_extras(&mut c2, false), 0);
}

#[test]
fn cluster_utils_match_python() {
    let (_, dists) = romesburg();
    let mut c = murtagh::cluster_data_dists(&dists, 5, murtagh::upgma())
        .unwrap()
        .remove(0);
    assert_eq!(find_cluster_centroid_from_dists(&mut c, &dists).unwrap(), 1);
    assert_eq!(indices(&get_node_list(&mut c)), [3, 4, 6, 2, 8, 1, 5, 7, 9]);
    assert_eq!(indices(&c.get_points()), [3, 4, 2, 1, 5]);
    // Python's calls, in its order, on one tree that each call reorders.
    let expected: [&[i64]; 6] = [
        &[9],
        &[9],
        &[8, 7],
        &[7, 6, 2],
        &[6, 2, 1, 5],
        &[2, 1, 5, 3, 4],
    ];
    for (n, want) in expected.iter().enumerate() {
        for breadth_first in [true, false] {
            let got = split_into_n_clusters_with_breadth_first(&mut c, n, breadth_first).unwrap();
            assert_eq!(&indices(&got)[..], *want, "{n} {breadth_first}");
        }
    }
    for n in 6..9 {
        assert!(split_into_n_clusters(&mut c, n).is_err());
    }
    let mut no_data = node(1.0, vec![leaf(1), leaf(2)], 3);
    assert!(matches!(
        find_cluster_centroid_from_dists(&mut no_data, &dists),
        Err(Error::Value { .. })
    ));
}

#[test]
fn compare_and_cmp() {
    assert_eq!(cmp(&1, &2), -1);
    assert_eq!(cmp(&2.0, &1.0), 1);
    assert_eq!(cmp("a", "a"), 0);
    assert_eq!(cmp(&f64::NAN, &1.0), 0);
    assert_eq!(cmptol(), 1e-6);
    let mut a = node(1.0, vec![leaf(1), leaf(2)], 3);
    let mut b = node(1.0 + 1e-7, vec![leaf(1), leaf(2)], 3);
    assert_eq!(a.compare_with_ignore_extras(&mut b, false), 0);
    let mut c = node(2.0, vec![leaf(1), leaf(2)], 3);
    assert_eq!(a.compare_with_ignore_extras(&mut c, false), -1);
    assert_eq!(a.compare(&mut c), 0);
    let mut d = node(2.0, vec![leaf(1)], 3);
    assert_eq!(a.compare(&mut d), 1);
}

#[test]
fn resemblance_matches_python() {
    let m = [0.4, 0.2, 0.3, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    assert_eq!(
        resemblance::find_min_val_in_list(&m, 5).unwrap(),
        (0, 3, 0.1)
    );
    assert_eq!(
        resemblance::find_min_val_in_list_with_min_idx(&m, 5, 7).unwrap(),
        (1, 4, 0.8)
    );
    assert!(matches!(
        resemblance::find_min_val_in_list(&m, 4),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        resemblance::show_metric_mat(&[0.4, 0.2, 0.3, 1.0 / 3.0, 0.5, 2.0], 4).unwrap(),
        "   ---      0.400000  0.200000  0.333333\n   ---       ---      0.300000  0.500000\n   \
         ---       ---       ---      2.000000\n   ---       ---       ---       ---    \n"
    );
    let pts = vec![vec![0.0, 0.0], vec![3.0, 4.0], vec![6.0, 8.0]];
    assert_eq!(
        resemblance::euclidean_distance(&pts).unwrap(),
        [5.0, 10.0, 5.0]
    );
    let clusters: Vec<Cluster> = pts
        .iter()
        .map(|p| {
            Cluster::new_with_params(&ClusterParams {
                position: p.clone(),
                ..Default::default()
            })
        })
        .collect();
    assert_eq!(
        resemblance::calc_metric_matrix(&clusters, resemblance::euclidean_distance).unwrap(),
        [5.0, 10.0, 5.0]
    );
    // The metric may be a closure that captures.
    let scale = 2.0;
    let scaled = |points: &[Vec<f64>]| {
        let distances = resemblance::euclidean_distance(points)?;
        Ok(distances.into_iter().map(|d| d * scale).collect())
    };
    assert_eq!(
        resemblance::calc_metric_matrix(&clusters, scaled).unwrap(),
        [10.0, 20.0, 10.0]
    );
    assert_eq!(resemblance::methods()[0].0, "Euclidean");
}

#[test]
fn standardize_methods() {
    let mat = vec![vec![1.0], vec![3.0], vec![4.0]];
    let methods = standardize::methods();
    assert_eq!(methods[0].0, "None");
    assert_eq!((methods[0].1)(&mat).unwrap(), mat);
    // Python's StandardizeMatrix of the vector [1., 3., 4.].
    let expected = [-1.0910894511799618, 0.2182178902359925, 0.8728715609439697];
    let got: Vec<f64> = standardize::std_dev(&mat)
        .unwrap()
        .into_iter()
        .map(|r| r[0])
        .collect();
    assert_eq!(got, expected);
}

/// Python's ClusterData clusters the first nPts rows and ignores the rest.
#[test]
fn murtagh_clusters_the_first_n_pts_rows() {
    let d: Vec<Vec<f64>> = (0..5)
        .map(|i| vec![f64::from(i), f64::from(i * i)])
        .collect();
    let mut first = murtagh::cluster_data(&d, 3, murtagh::wards())
        .unwrap()
        .remove(0);
    let mut only = murtagh::cluster_data(&d[..3], 3, murtagh::wards())
        .unwrap()
        .remove(0);
    let mut points: Vec<i64> = first.get_points().iter().map(|c| c.get_index()).collect();
    points.sort();
    assert_eq!(points, [1, 2, 3]);
    assert_eq!(first.compare_with_ignore_extras(&mut only, false), 0);
}
