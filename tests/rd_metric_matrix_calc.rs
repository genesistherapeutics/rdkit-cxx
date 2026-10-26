// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::data_manip::metric::*;
use rdkit_cxx::data_structs::{ExplicitBitVect, SparseBitVect};

// Expected values come from Python's rdMetricMatrixCalc.

fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} vs {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1e-12, "{actual:?} vs {expected:?}");
    }
}

fn ebv(size: u32, bits: &[u32]) -> ExplicitBitVect {
    let mut bv = ExplicitBitVect::new(size);
    bv.set_bits_from_list(bits).unwrap();
    bv
}

fn sbv(size: u32, bits: &[u32]) -> SparseBitVect {
    let mut bv = SparseBitVect::new(size);
    bv.set_bits_from_list(bits).unwrap();
    bv
}

#[test]
fn euclidean() {
    let rows = vec![
        vec![0.5, 1.25, -2.0],
        vec![3.0, 0.0, 1.5],
        vec![-1.0, 2.0, 2.0],
        vec![0.0, 0.0, 0.0],
    ];
    assert_close(
        &get_euclidean_dist_mat(&rows).unwrap(),
        &[
            4.479118216792229,
            4.337337893224368,
            4.5,
            2.4109126902482387,
            3.3541019662496847,
            3.0,
        ],
    );
    assert!(
        get_euclidean_dist_mat(&[vec![1.0, 2.0]])
            .unwrap()
            .is_empty()
    );
    assert_eq!(get_euclidean_dist_mat(&[vec![], vec![]]).unwrap(), [0.0]);
    assert!(matches!(
        get_euclidean_dist_mat(&[]),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        get_euclidean_dist_mat(&[vec![1.0], vec![1.0, 2.0]]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn tanimoto() {
    let e = [
        ebv(16, &[0, 1, 2, 5]),
        ebv(16, &[1, 2, 3, 4, 9]),
        ebv(16, &[0, 5, 15]),
        ebv(16, &[]),
    ];
    assert_close(
        &get_tanimoto_dist_mat(&e).unwrap(),
        &[0.7142857142857143, 0.6, 1.0, 1.0, 1.0, 1.0],
    );
    assert_close(
        &get_tanimoto_sim_mat(&e).unwrap(),
        &[0.2857142857142857, 0.4, 0.0, 0.0, 0.0, 0.0],
    );
    let s = [
        sbv(1000, &[0, 10, 999]),
        sbv(1000, &[10, 500]),
        sbv(1000, &[0, 10, 999, 7]),
    ];
    assert_close(&get_tanimoto_dist_mat(&s).unwrap(), &[0.75, 0.25, 0.8]);
    assert_close(&get_tanimoto_sim_mat(&s).unwrap(), &[0.25, 0.75, 0.2]);
    // RDKit folds the larger vector to the smaller's size.
    let folded = [
        ebv(8, &[0, 1, 2]),
        ebv(8, &[1, 2, 3, 4]),
        ebv(16, &[1, 9, 3]),
    ];
    assert_close(&get_tanimoto_sim_mat(&folded).unwrap(), &[0.4, 0.25, 0.5]);
}

#[test]
fn tanimoto_errors() {
    assert!(matches!(
        get_tanimoto_sim_mat(&[ebv(8, &[1])]),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        get_tanimoto_dist_mat::<ExplicitBitVect>([]),
        Err(Error::Runtime { .. })
    ));
    // Sizes that do not fold, and a fold factor RDKit refuses.
    for other in [ebv(12, &[1]), ebv(1, &[0])] {
        assert!(matches!(
            get_tanimoto_sim_mat([&ebv(8, &[1]), &other]),
            Err(Error::Value { .. })
        ));
    }
}

#[test]
fn tanimoto_guards() {
    // RDKit divides by the smaller size (SIGFPE in Python).
    for list in [[ebv(8, &[1]), ebv(0, &[])], [ebv(0, &[]), ebv(8, &[1])]] {
        assert!(matches!(
            get_tanimoto_dist_mat(&list),
            Err(Error::Value { .. })
        ));
    }
    // Vectors of no bits compare as RDKit compares them when no folding is
    // needed.
    assert_eq!(
        get_tanimoto_sim_mat(&[ebv(0, &[]), ebv(0, &[])]).unwrap(),
        [0.0]
    );
    // Folding an unfolded fingerprint (2^32 - 1 bits) computes its size in
    // int; comparing two of them folds nothing.
    let generator = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 1,
        ..Default::default()
    })
    .unwrap();
    let unfolded = generator
        .get_sparse_fingerprint(&mol_from_smiles("CCO").unwrap())
        .unwrap();
    let other = generator
        .get_sparse_fingerprint(&mol_from_smiles("CCN").unwrap())
        .unwrap();
    assert!(matches!(
        get_tanimoto_sim_mat([&unfolded, &sbv(1024, &[1])]),
        Err(Error::Value { .. })
    ));
    let sim = get_tanimoto_sim_mat([&unfolded, &other]).unwrap();
    assert_eq!(sim.len(), 1);
    assert!(sim[0] > 0.0 && sim[0] < 1.0);
}

#[test]
fn item_limit() {
    // RDKit's unsigned index into the triangle wraps past 65537 items.
    let rows = vec![Vec::new(); 65538];
    assert!(matches!(
        get_euclidean_dist_mat(&rows),
        Err(Error::Value { .. })
    ));
    let bvs: Vec<ExplicitBitVect> = (0..65538).map(|_| ebv(1, &[])).collect();
    assert!(matches!(
        get_tanimoto_sim_mat(&bvs),
        Err(Error::Value { .. })
    ));
}
