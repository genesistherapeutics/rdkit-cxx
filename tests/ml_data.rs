// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/ML/Data/UnitTestQuantize.py, UnitTestStats.py,
// UnitTestMLData.py, UnitTestFilter.py, and UnitTestDoctests.py (with the
// doctests of DataUtils and SplitData), with values from Python's
// rdkit.ML.Data modules in RDKit 2025.09.6.

mod common;

use std::sync::{Mutex, MutexGuard};

use common::{assert_close, temp_path};
use rdkit_cxx::Error;
use rdkit_cxx::ml::data::data_utils::*;
use rdkit_cxx::ml::data::ml_data::{MLDataSet, MLDataSetParams, MLQuantDataSet, Value};
use rdkit_cxx::ml::data::quantize::*;
use rdkit_cxx::ml::data::split_data::*;
use rdkit_cxx::ml::data::stats::*;
use rdkit_cxx::ml::data::transforms::get_avail_transforms;

/// The tests that seed the shared generator take turns.
fn seeded() -> MutexGuard<'static, ()> {
    static LOCK: Mutex<()> = Mutex::new(());
    let guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    init_random_numbers((23, 42));
    guard
}

fn quant(pairs: &[(f64, usize)], n_bounds: usize, n_possible_res: usize) -> (Vec<f64>, f64) {
    let (vals, res): (Vec<f64>, Vec<usize>) = pairs.iter().copied().unzip();
    find_var_mult_quant_bounds(&vals, n_bounds, &res, n_possible_res).unwrap()
}

#[test]
fn quantize_matches_python() {
    let one1 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 0),
        (1.4, 0),
        (1.6, 0),
        (2., 1),
        (2.1, 1),
        (2.2, 1),
        (2.3, 1),
    ];
    assert_eq!(quant(&one1, 1, 2), (vec![1.8], 0.9709505944546688));
    let (vals, res): (Vec<f64>, Vec<usize>) = one1.iter().copied().unzip();
    assert_eq!(
        find_var_quant_bound(&vals, &res, 2).unwrap(),
        (1.8, 0.9709505944546688)
    );
    let one2 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 0),
        (1.4, 1),
        (1.6, 0),
        (2., 1),
        (2.1, 1),
        (2.2, 1),
        (2.3, 1),
    ];
    assert_eq!(quant(&one2, 1, 2), (vec![1.8], 0.6099865470109875));
    let one3 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 2),
        (1.4, 2),
        (1.6, 2),
        (2., 2),
        (2.1, 1),
        (2.2, 1),
        (2.3, 1),
    ];
    assert_eq!(
        quant(&one3, 1, 3),
        (vec![1.2999999999999998], 0.8812908992306928)
    );
    let one4 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.2, 1),
        (1.4, 0),
        (1.4, 0),
        (1.6, 0),
        (2., 1),
        (2.1, 1),
        (2.1, 1),
        (2.1, 1),
        (2.1, 1),
        (2.2, 1),
        (2.3, 1),
    ];
    assert_eq!(quant(&one4, 1, 2), (vec![1.8], 0.6893917467430879));
    let one5 = [
        (1., 0),
        (1.1, 0),
        (2.2, 1),
        (1.2, 0),
        (1.6, 0),
        (1.4, 0),
        (2., 1),
        (2.1, 1),
        (1.4, 0),
        (2.3, 1),
    ];
    assert_eq!(quant(&one5, 1, 2), (vec![1.8], 0.9709505944546688));
    let mult1 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 2),
        (1.4, 2),
        (1.6, 2),
        (2., 2),
        (2.1, 1),
        (2.1, 1),
        (2.1, 1),
        (2.2, 1),
        (2.3, 1),
    ];
    assert_eq!(
        quant(&mult1, 2, 3),
        (vec![1.2999999999999998, 2.05], 1.5545851693377992)
    );
    let mult2 = [
        (1., 0),
        (2.1, 1),
        (1.1, 0),
        (1.2, 0),
        (1.4, 2),
        (1.6, 2),
        (2., 2),
        (1.4, 2),
        (2.1, 1),
        (2.2, 1),
        (2.1, 1),
        (2.3, 1),
    ];
    assert_eq!(
        quant(&mult2, 2, 3),
        (vec![1.2999999999999998, 2.05], 1.5545851693377992)
    );
    let mult3 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 2),
        (1.4, 2),
        (1.6, 2),
        (2., 2),
        (2.1, 1),
        (2.1, 1),
        (2.1, 1),
        (2.2, 1),
        (2.3, 1),
        (3.0, 3),
        (3.1, 3),
        (3.2, 3),
        (3.3, 3),
    ];
    assert_eq!(
        quant(&mult3, 3, 4),
        (vec![1.2999999999999998, 2.05, 2.65], 1.9772170014624826)
    );
    let mult4 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 1),
        (1.4, 1),
        (1.6, 1),
        (2., 1),
        (2.1, 0),
        (2.1, 0),
        (2.1, 0),
        (2.2, 0),
        (2.3, 0),
    ];
    assert_eq!(
        quant(&mult4, 2, 2),
        (vec![1.2999999999999998, 2.05], 0.9182958340544894)
    );
    let mult5 = [
        (1., 0),
        (1.1, 0),
        (1.2, 0),
        (1.4, 1),
        (1.4, 0),
        (1.6, 1),
        (2., 1),
        (2.1, 0),
        (2.1, 0),
        (2.1, 0),
        (2.2, 1),
        (2.3, 0),
    ];
    assert_eq!(
        quant(&mult5, 2, 2),
        (vec![1.2999999999999998, 2.05], 0.3470664196983776)
    );
    let many: Vec<(f64, usize)> = (0..40)
        .map(|i| ((i * 13 % 29) as f64 / 7.0, i * 5 % 3))
        .collect();
    assert_eq!(
        quant(&many, 4, 3),
        (
            vec![
                0.5,
                0.9285714285714286,
                1.3571428571428572,
                3.3571428571428568
            ],
            0.36432252022542366
        )
    );
    // Github issue 18.
    let issue18 = [(0., 0), (1., 0), (2., 1), (3., 1), (4., 1)];
    assert_eq!(quant(&issue18, 1, 2), (vec![1.5], 0.9709505944546688));
    assert_eq!(quant(&[(1.0, 1)], 2, 2), (vec![0.0], 0.0));
    assert_eq!(
        quant(&[(1.0, 1), (1.0, 0), (1.0, 1)], 2, 2),
        (vec![0.0], 0.0)
    );
    assert_eq!(quant(&[], 2, 2), (vec![], -1e8));
    assert!(matches!(
        find_var_mult_quant_bounds(&[1.0, 2.0], 1, &[0], 2),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        find_var_mult_quant_bounds(&[1.0, 2.0], 1, &[0, 2], 2),
        Err(Error::Index { .. })
    ));
    // Python gives ([1.0, nan], 0.5409...) here, from where its sort happens
    // to leave the NaNs.
    let nan = f64::NAN;
    assert!(matches!(
        find_var_mult_quant_bounds(&[1.0, nan, 0.0, 2.0, nan, 3.0], 2, &[0, 1, 0, 1, 1, 0], 2),
        Err(Error::Value { .. })
    ));
    assert!(feq(1.0, 1.0 + 1e-9));
    assert!(!feq(1.0, 1.0 + 1e-7));
    assert!(feq_with_tol(1.0, 1.05, 0.1));
}

fn qdat(name: &str) -> std::path::PathBuf {
    let path = temp_path(&format!("{name}.qdat"));
    std::fs::write(
        &path,
        "# Testing data\n# ----------\n# Variable Table\n# foo1 []\n# foo2 []\n# foo3 []\n# foo4 \
         []\n# res []\n# ---------\np1 0 0 0 0 0\np2 0 0 0 1 0\np3 1 0 0 0 1\np4 2 1 0 0 1\np5 2 2 \
         1 0 1\n",
    )
    .unwrap();
    path
}

/// RDKit's `test.dat`, whose lines end in spaces, in a file of the test's
/// own.
fn dat(name: &str) -> std::path::PathBuf {
    let path = temp_path(&format!("{name}.dat"));
    std::fs::write(
        &path,
        "# Testing data\n# ----------\n# Variable Table\n# foo1 []\n# foo2 []\n# foo3 []\n# foo4 \
         []\n# res [2,4]\n# ---------\np1 foo 1 1.0 1 1.1   \np2 foo 2 1.0 1 2.1   \np3 foo 3 1.2 \
         1.1 3.1 \np4 foo 4 1.0 1 4.1   \np5 foo 5 1.1 1 5.1   \n\n\n",
    )
    .unwrap();
    path
}

fn names(names: &[&str]) -> Vec<String> {
    names.iter().map(|s| s.to_string()).collect()
}

fn row(values: &[Value]) -> Vec<Value> {
    values.to_vec()
}

use Value::{Float as F, Int as I, Str as S};

#[test]
fn quant_data_set() {
    let d = build_quant_data_set(qdat("quant_data_set")).unwrap();
    assert_eq!(d.get_n_pts(), 5);
    assert_eq!(d.get_n_vars(), 4);
    assert_eq!(d.get_n_results(), 1);
    assert_eq!(
        d.get_var_names(),
        names(&["foo1", "foo2", "foo3", "foo4", "res"])
    );
    assert_eq!(d.get_pt_names(), names(&["p1", "p2", "p3", "p4", "p5"]));
    assert_eq!(d.get_n_possible_vals(), [3, 3, 2, 2, 2]);
    assert_eq!(d.get_quant_bounds(), vec![Vec::<f64>::new(); 5]);
    assert_eq!(
        d.get_results(),
        [vec![0], vec![0], vec![1], vec![1], vec![1]]
    );
    assert_eq!(d.get_all_data()[1], [0, 0, 0, 1, 0]);
    assert_eq!(d.get_input_data()[3], [2, 1, 0, 0]);
    assert_eq!(
        d.get_named_data().unwrap()[2],
        ("p3".to_owned(), vec![1, 0, 0, 0, 1])
    );
    let mut out = Vec::new();
    let all: Vec<Vec<Value>> = d
        .get_all_data()
        .into_iter()
        .map(|r| r.into_iter().map(Value::Int).collect())
        .collect();
    write_data(&mut out, d.get_var_names(), d.get_quant_bounds(), &all).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("DataUtils"));
    assert!(text.contains("foo1"));
    assert!(text.contains("2 2 1 0 1"));
    assert_eq!(
        calc_n_possible_using_map(&all, &[0, 1, 2, 3, 4], d.get_quant_bounds()).unwrap(),
        [3, 3, 2, 2, 2]
    );
    let again = MLQuantDataSet::new_with_params(
        d.get_all_data(),
        &MLDataSetParams {
            var_names: Some(d.get_var_names().to_vec()),
            q_bounds: Some(d.get_quant_bounds().to_vec()),
            pt_names: Some(d.get_pt_names().to_vec()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(again, d);
}

#[test]
fn general_data_set() {
    let mut d = build_data_set(dat("general_data_set")).unwrap();
    assert_eq!(d.get_n_pts(), 5);
    assert_eq!(d.get_n_vars(), 4);
    assert_eq!(d.get_n_results(), 1);
    assert_eq!(
        d.get_var_names(),
        names(&["foo1", "foo2", "foo3", "foo4", "res"])
    );
    assert_eq!(d.get_pt_names(), names(&["p1", "p2", "p3", "p4", "p5"]));
    assert_eq!(d.get_n_possible_vals(), [0, 6, 0, 0, 3]);
    assert_eq!(
        d.get_quant_bounds(),
        [vec![], vec![], vec![], vec![], vec![2.0, 4.0]]
    );
    let results: Vec<Value> = d.get_results().into_iter().map(|r| r[0].clone()).collect();
    assert_eq!(results, [F(1.1), F(2.1), F(3.1), F(4.1), F(5.1)]);
    let foo = || S("foo".to_owned());
    assert_eq!(
        d.get_all_data()[1],
        row(&[foo(), I(2), F(1.0), I(1), F(2.1)])
    );
    assert_eq!(d.get_input_data()[4], row(&[foo(), I(5), F(1.1), I(1)]));
    assert_eq!(
        d.get_named_data().unwrap()[3],
        row(&[S("p4".to_owned()), foo(), I(4), F(1.0), I(1), F(4.1)])
    );
    assert_eq!(
        d.get(1).unwrap(),
        row(&[S("p2".to_owned()), foo(), I(2), F(1.0), I(1), F(2.1)])
    );
    let p2 = d.get(1).unwrap();
    d.set(3, p2.clone()).unwrap();
    assert_eq!(d.get(3).unwrap(), p2);
    assert!(matches!(
        d.set(3, p2[1..].to_vec()),
        Err(Error::Value { .. })
    ));
    d.add_point(row(&[
        S("p6".to_owned()),
        foo(),
        I(6),
        F(1.0),
        I(1),
        F(6.1),
    ]))
    .unwrap();
    assert_eq!(d.get_n_pts(), 6);
    assert!(matches!(
        d.add_points(vec![vec![]], vec![]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn write_data_matches_python() {
    let d = build_data_set(dat("write_data_matches_python")).unwrap();
    let mut out = Vec::new();
    write_data(
        &mut out,
        d.get_var_names(),
        &[
            vec![],
            vec![1.5, 2.0],
            vec![],
            vec![1e-5, 1e16],
            vec![2.0, 4.0],
        ],
        &d.get_all_data(),
    )
    .unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "# Quantized data from DataUtils\n# ----------\n# Variable Table\n# foo1 []\n# foo2 [1.5, \
         2.0]\n# foo3 []\n# foo4 [1e-05, 1e+16]\n# res [2.0, 4.0]\n# ----------\nfoo 1 1.0 1 \
         1.1\nfoo 2 1.0 1 2.1\nfoo 3 1.2 1.1 3.1\nfoo 4 1.0 1 4.1\nfoo 5 1.1 1 5.1\n"
    );
}

#[test]
fn randomize_activities_matches_python() {
    let _guard = seeded();
    let mut d = build_data_set(dat("randomize_activities_matches_python")).unwrap();
    randomize_activities_with_shuffle(&mut d, true).unwrap();
    let last: Vec<Value> = d.get_all_data().into_iter().map(|r| r[4].clone()).collect();
    assert_eq!(last, [F(2.1), F(1.1), F(3.1), F(4.1), F(5.1)]);
    // RDKit 2025.09.6's Python raises NameError without shuffle.
    let mut d = build_data_set(dat("randomize_activities_matches_python")).unwrap();
    assert!(matches!(
        randomize_activities_with_shuffle(&mut d, false),
        Err(Error::Value { .. })
    ));
}

#[test]
fn random_choices_match_python() {
    let _guard = seeded();
    assert_eq!(permutation(10), [1, 5, 6, 4, 2, 3, 0, 7, 8, 9]);
    assert_eq!(permutation(7), [4, 5, 6, 0, 2, 1, 3]);

    // SplitData's doctests.
    init_random_numbers((23, 42));
    assert_eq!(
        split_indices(10, 0.5).unwrap(),
        (vec![1, 5, 6, 4, 2], vec![3, 0, 7, 8, 9])
    );
    assert_eq!(
        split_indices(10, 0.5).unwrap(),
        (vec![5, 2, 9, 8, 7], vec![6, 0, 3, 1, 4])
    );
    let legacy = SplitIndicesParams {
        legacy: true,
        ..Default::default()
    };
    init_random_numbers((23, 42));
    assert_eq!(
        split_indices_with_params(10, 0.5, &legacy).unwrap(),
        (vec![3, 5, 7, 8, 9], vec![0, 1, 2, 4, 6])
    );
    assert_eq!(
        split_indices_with_params(10, 0.5, &legacy).unwrap(),
        (vec![0, 1, 2, 3, 5, 8, 9], vec![4, 6, 7])
    );
    let replacement = SplitIndicesParams {
        replacement: true,
        ..Default::default()
    };
    init_random_numbers((23, 42));
    assert_eq!(
        split_indices_with_params(10, 0.5, &replacement).unwrap(),
        (vec![9, 9, 8, 0, 5], vec![1, 2, 3, 4, 6, 7])
    );
    assert_eq!(
        split_indices_with_params(10, 0.5, &replacement).unwrap(),
        (vec![4, 5, 1, 1, 4], vec![0, 2, 3, 6, 7, 8, 9])
    );
    assert!(matches!(split_indices(10, -0.1), Err(Error::Value { .. })));
    assert!(matches!(split_indices(10, 1.1), Err(Error::Value { .. })));

    init_random_numbers((23, 42));
    let data: Vec<char> = "abcdefghij".chars().collect();
    assert_eq!(
        split_data_set(&data, 0.3).unwrap(),
        (vec!['b', 'f', 'g'], vec!['e', 'c', 'd', 'a', 'h', 'i', 'j'])
    );
    assert!(matches!(
        split_data_set(&data, 1.1),
        Err(Error::Value { .. })
    ));
}

fn d1() -> Vec<Vec<f64>> {
    [
        [0., 0., 0.],
        [1., 0., 1.],
        [2., 0., 0.],
        [3., 0., 1.],
        [4., 0., 0.],
        [5., 0., 0.],
    ]
    .iter()
    .map(|r| r.to_vec())
    .collect()
}

#[test]
fn filter_data_matches_python() {
    let _guard = seeded();
    let d = d1();
    let rows = |idx: &[usize]| -> Vec<Vec<f64>> { idx.iter().map(|&i| d[i].clone()).collect() };
    // The fraction, then the kept and rejected indices, without and with a
    // subset of the rows.
    type Case<'a> = (f64, &'a [usize], &'a [usize], &'a [usize], &'a [usize]);
    let cases: [Case; 5] = [
        (0.5, &[2, 3, 0, 1], &[4, 5], &[3, 4, 5, 1], &[]),
        (0.7, &[3, 0, 1], &[2, 4, 5], &[3, 5, 1], &[4]),
        (0.75, &[3, 0, 1], &[2, 4, 5], &[3, 5, 1], &[4]),
        (0.333, &[2, 3, 0, 1, 4, 5], &[], &[3, 4, 5], &[1]),
        (0.25, &[2, 3, 0, 4], &[1, 5], &[3, 4, 5], &[1]),
    ];
    let subset = FilterDataParams {
        indices_to_use: Some(vec![5, 4, 3, 1]),
        ..Default::default()
    };
    for (frac, keep, rej, keep_subset, rej_subset) in cases {
        init_random_numbers((23, 42));
        assert_eq!(filter_data(&d, 1.0, frac).unwrap(), (rows(keep), rows(rej)));
        init_random_numbers((23, 42));
        assert_eq!(
            filter_data_indices_only(&d, 1.0, frac).unwrap(),
            (keep.to_vec(), rej.to_vec())
        );
        init_random_numbers((23, 42));
        assert_eq!(
            filter_data_indices_only_with_params(&d, 1.0, frac, &subset).unwrap(),
            (keep_subset.to_vec(), rej_subset.to_vec())
        );
        // Python's test4 and test5: all the indices give the same rows.
        let all = FilterDataParams {
            indices_to_use: Some((0..6).collect()),
            ..Default::default()
        };
        init_random_numbers((23, 42));
        assert_eq!(
            filter_data_with_params(&d, 1.0, frac, &all).unwrap(),
            (rows(keep), rows(rej))
        );
    }
    assert!(matches!(
        filter_data(&d, 1.0, -1.0),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        filter_data(&d, 1.0, 2.0),
        Err(Error::Value { .. })
    ));
    let col3 = FilterDataParams {
        col: 3,
        ..Default::default()
    };
    assert!(matches!(
        filter_data_with_params(&d, 1.0, 0.5, &col3),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        filter_data(&d, 7.0, 0.5),
        Err(Error::Value { .. })
    ));
}

#[test]
fn count_results_matches_python() {
    let d = d1();
    assert_eq!(count_results(&d).unwrap(), [(0.0, 4), (1.0, 2)]);
    assert_eq!(
        count_results_with_params(
            &d,
            &CountResultsParams {
                col: 0,
                bounds: Some(vec![1.5, 3.5])
            }
        )
        .unwrap(),
        [(0.0, 2), (1.0, 2), (2.0, 2)]
    );
}

#[test]
fn text_file_to_data_matches_python() {
    let path = temp_path("t.csv");
    std::fs::write(
        &path,
        "id,a,\"b c\",res\r\nm1,1,2.5,x\n\"m,2\",3,\"4\"\"\",0\n\nm3,-5,1e3,  7\n",
    )
    .unwrap();
    let t = text_file_to_data(&path).unwrap();
    assert_eq!(t.get_var_names(), names(&["id", "a", "b c", "res"]));
    assert_eq!(t.get_pt_names(), names(&["m1", "m,2", "m3"]));
    assert_eq!(
        t.get_all_data(),
        [
            row(&[I(1), F(2.5), S("x".to_owned())]),
            row(&[I(3), S("4\"".to_owned()), I(0)]),
            row(&[I(-5), F(1000.0), I(7)]),
        ]
    );
    assert_eq!(t.get_n_possible_vals(), [4, 0, 0]);
    let only =
        text_file_to_data_with_only_cols(&path, Some(names(&["id", "res", "zzz", "a"]))).unwrap();
    assert_eq!(only.get_var_names(), names(&["id", "res", "res", "a"]));
    assert_eq!(
        only.get_all_data(),
        [
            row(&[S("x".to_owned()), S("x".to_owned()), I(1)]),
            row(&[I(0), I(0), I(3)]),
            row(&[I(7), I(7), I(-5)]),
        ]
    );
    let rows = vec![names(&["n", "x"]), names(&["a", "1"]), names(&["b"])];
    assert!(matches!(text_to_data(rows), Err(Error::Value { .. })));
}

#[test]
fn calc_n_possible_and_take_ensemble() {
    let data = vec![
        row(&[I(1), F(2.0), S("a".to_owned()), F(3.5)]),
        row(&[I(4), I(0), S("b".to_owned()), I(1)]),
    ];
    let none = vec![Vec::new(); 4];
    assert_eq!(
        calc_n_possible_using_map(&data, &[0, 1, 2, 3], &none).unwrap(),
        [5, 3, 0, 0]
    );
    assert_eq!(
        calc_n_possible_using_map_with_n_q_bounds(&data, &[1, 0], &none[..2], Some(&[0, 1]))
            .unwrap(),
        [3, 0]
    );
    // DataUtils' doctests.
    assert_eq!(
        take_ensemble(&[10, 20, 30, 40, 50], &[1, 2, 3]).unwrap(),
        [20, 30, 40]
    );
    let v = ["foo", "10", "20", "30", "40", "50", "1"];
    assert_eq!(
        take_ensemble_with_is_data_vect(&v, &[1, 2, 3], true).unwrap(),
        ["foo", "20", "30", "40", "1"]
    );
    assert!(matches!(
        take_ensemble(&[1], &[1]),
        Err(Error::Index { .. })
    ));
}

#[test]
fn ml_data_set_defaults() {
    let d = MLDataSet::new(vec![row(&[I(1), F(2.5), I(0)]), row(&[I(3), F(1.0), I(1)])]).unwrap();
    assert_eq!(d.get_n_vars(), 2);
    assert_eq!(d.get_n_possible_vals(), [4, 0, 2]);
    assert_eq!(d.get_var_names(), names(&["", ""]));
    assert!(matches!(MLDataSet::new(vec![]), Err(Error::Index { .. })));
    let q = MLQuantDataSet::new(vec![vec![1, 2, 0], vec![3, 1, 1]]).unwrap();
    assert_eq!(q.get_n_possible_vals(), [4, 3, 2]);
    assert_eq!(q.get_input_data(), [vec![1, 2], vec![3, 1]]);
    assert!(matches!(
        MLQuantDataSet::new(vec![vec![1, 2], vec![3]]),
        Err(Error::Value { .. })
    ));
    let names: Vec<&str> = get_avail_transforms().iter().map(|t| t.0).collect();
    assert_eq!(names, ["Center", "Normalize", "Standardize"]);
}

fn stats_data() -> Vec<Vec<f64>> {
    [
        [7., 4., 3.],
        [4., 1., 8.],
        [6., 3., 5.],
        [8., 6., 1.],
        [8., 5., 7.],
        [7., 2., 9.],
        [5., 3., 3.],
        [9., 5., 8.],
        [7., 4., 5.],
        [8., 2., 2.],
    ]
    .iter()
    .map(|r| r.to_vec())
    .collect()
}

#[test]
fn stats_match_python() {
    let d = stats_data();
    assert_eq!(
        form_correlation_matrix(&d).unwrap(),
        [
            [1.0, 0.6686573150705388, -0.10131373970455891],
            [0.6686573150705388, 1.0, -0.28792771037232723],
            [-0.10131373970455891, -0.28792771037232723, 1.0]
        ]
    );
    let cov = vec![
        vec![4., 2., 0.6],
        vec![4.2, 2.1, 0.59],
        vec![3.9, 2.0, 0.58],
        vec![4.3, 2.1, 0.62],
        vec![4.1, 2.2, 0.63],
    ];
    // Python's matrix product is BLAS's, which may round differently.
    let expected = [
        [0.025, 0.007499999999999998, 0.0017499999999999985],
        [
            0.007499999999999998,
            0.007000000000000012,
            0.0013500000000000025,
        ],
        [
            0.0017499999999999985,
            0.0013500000000000025,
            0.00043000000000000075,
        ],
    ];
    for (got, want) in form_covariance_matrix(&cov).unwrap().iter().zip(expected) {
        for (g, w) in got.iter().zip(want) {
            assert_close(*g, w);
        }
    }

    let (vals, vects) = principal_components(&d).unwrap();
    for (g, w) in vals
        .iter()
        .zip([1.7687741361034233, 0.927075916898894, 0.30414994699768233])
    {
        assert_close(*g, w);
    }
    // Python's second value is rows of LAPACK's eigenvector matrix, not
    // eigenvectors; the port's are eigenvectors of the correlation matrix.
    let corr = form_correlation_matrix(&d).unwrap();
    for (val, vect) in vals.iter().zip(&vects) {
        let norm: f64 = vect.iter().map(|x| x * x).sum();
        assert_close(norm, 1.0);
        for i in 0..3 {
            let cv: f64 = (0..3).map(|j| corr[i][j] * vect[j]).sum();
            assert_close(cv, val * vect[i]);
        }
    }
    let py_vects = [
        [-0.6420045763498946, -0.6632174243435959, 0.3846722916884436],
        [0.3416691692479823, 0.20166618906061493, 0.9179286066874492],
        [-0.6863616413605773, 0.7207450285897329, 0.09713033013430539],
    ];
    let (vals, _) = principal_components_with_reverse_order(&d, false).unwrap();
    assert!(vals[0] < vals[1] && vals[1] < vals[2]);

    let py_vects: Vec<Vec<f64>> = py_vects.iter().map(|r| r.to_vec()).collect();
    let pts = transform_points(&py_vects, &d).unwrap();
    for (got, want) in pts.iter().zip([
        [
            -1.2036209823525186,
            -1.7926500625885375,
            0.08776265687676767,
        ],
        [4.635406478170171, 1.166986895922917, 0.47026414586082793],
    ]) {
        for (g, w) in got.iter().zip(want) {
            assert_close(*g, w);
        }
    }

    let v = [1.5, 2.25, 3.0, 7.125];
    assert_eq!(mean_and_dev(&v), (3.46875, 2.5132461578603875));
    assert_eq!(
        mean_and_dev_with_sample_sd(&v, false),
        (3.46875, 2.176535018670731)
    );
    assert_eq!(mean_and_dev(&[]), (0.0, 0.0));
    assert_eq!(mean_and_dev(&[4.0]), (4.0, 0.0));
    assert_eq!(r2(&v, 0.75), 0.9604205318491033);
    assert_eq!(
        get_confidence_interval(1.5, 10).unwrap(),
        1.0729608100951311
    );
    assert_eq!(
        get_confidence_interval_with_level(1.5, 31, 99.0).unwrap(),
        0.7408718708604465
    );
    assert_eq!(
        get_confidence_interval_with_level(2.0, 1000, 80.0).unwrap(),
        0.08133378141953071
    );
    assert!(matches!(
        get_confidence_interval_with_level(1.0, 3, 50.0),
        Err(Error::Key { .. })
    ));
    assert_eq!(t_table().len(), 48);
    assert_eq!(t_confs()[2], (95.0, 3));
    let column: Vec<Vec<f64>> = v.iter().map(|&x| vec![x]).collect();
    let std: Vec<f64> = standardize_matrix(&column)
        .unwrap()
        .into_iter()
        .map(|r| r[0])
        .collect();
    assert_eq!(
        std,
        [
            -0.7833494518006403,
            -0.48493061301944396,
            -0.1865117742382477,
            1.4547918390583319
        ]
    );
}

#[test]
fn transform_preserves_dot_products() {
    // Python's testTransform2: rotating points into the PCA frame keeps
    // their angles.
    let d: Vec<Vec<f64>> = [
        [2.068703704, 2.040555556, 2.068703704],
        [1.48537037, -0.186756425, 1.48537037],
        [1.917469136, 0.785465797, 1.917469136],
        [2.068703704, 1.125743575, 2.068703704],
        [2.138703704, 1.283243575, 2.138703704],
        [1.730740741, 1.457222222, -0.179901738],
        [1.973796296, 1.889320988, 0.792320484],
        [1.390462963, -0.37502803, 0.171950113],
        [1.822561728, 0.597194192, 0.604048879],
        [2.189074074, 1.359539872, 1.366394558],
    ]
    .iter()
    .map(|r| r.to_vec())
    .collect();
    let (_, vects) = principal_components(&d).unwrap();
    let pts = transform_points(&vects, &d).unwrap();
    let n = d.len() as f64;
    let avg: Vec<f64> = (0..3)
        .map(|j| d.iter().map(|r| r[j]).sum::<f64>() / n)
        .collect();
    let unit = |v: Vec<f64>| {
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        v.into_iter().map(|x| x / norm).collect::<Vec<f64>>()
    };
    let centered: Vec<Vec<f64>> = d
        .iter()
        .map(|r| unit(r.iter().zip(&avg).map(|(x, a)| x - a).collect()))
        .collect();
    let pts: Vec<Vec<f64>> = pts.into_iter().map(unit).collect();
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    for i in 0..d.len() {
        for j in 0..d.len() {
            let (a, b) = (dot(&centered[i], &centered[j]), dot(&pts[i], &pts[j]));
            assert!((a - b).abs() < 1e-4, "{a} {b}");
        }
    }
}
