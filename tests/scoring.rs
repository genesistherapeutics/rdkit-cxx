// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/ML/Scoring/UnitTestScoring.py, with values from Python's
// rdkit.ML.Scoring.Scoring in RDKit 2025.03.6, which match to the last bit.

use rdkit_cxx::Error;
use rdkit_cxx::ml::scoring::scoring::*;

const ACC: f64 = 1e-4;

type Scores = Vec<Vec<f64>>;

/// 400 molecules, 20 actives, best first: Python's setUp.
fn cases() -> (Scores, Scores, Scores, Scores) {
    let act = vec![vec![1.0]; 20];
    let dcy = vec![vec![0.0]; 380];
    let best = [act.clone(), dcy.clone()].concat();
    let worst = [dcy, act].concat();
    (best, worst, vec![vec![1.0]; 400], vec![vec![0.0]; 400])
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < ACC
}

#[test]
fn enrichment_factor() {
    let (best, worst, all_actives, all_decoys) = cases();
    let fractions = [20.0 / 400.0];
    assert!(close(
        calc_enrichment(&best, 0, &fractions).unwrap()[0],
        20.0
    ));
    assert!(close(
        calc_enrichment(&worst, 0, &fractions).unwrap()[0],
        0.0
    ));
    let empty: Vec<Vec<f64>> = Vec::new();
    assert!(matches!(
        calc_enrichment(&empty, 0, &fractions),
        Err(Error::Value { .. })
    ));
    assert!(close(
        calc_enrichment(&all_actives, 0, &fractions).unwrap()[0],
        1.0
    ));
    assert_eq!(calc_enrichment(&all_decoys, 0, &fractions).unwrap()[0], 0.0);
    let small = [fractions[0] / 100.0];
    assert!(close(calc_enrichment(&best, 0, &small).unwrap()[0], 20.0));
    assert!(matches!(
        calc_enrichment(&best, 0, &[]),
        Err(Error::Value { .. })
    ));
    assert!(close(calc_enrichment(&best, 0, &[0.0]).unwrap()[0], 20.0));
    assert!(matches!(
        calc_enrichment(&best, 0, &[-0.05]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        calc_enrichment(&best, 0, &[1.5]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn rie() {
    let (best, worst, _, all_decoys) = cases();
    let (alpha, ratio) = (100.0f64, 20.0 / 400.0);
    let rie_max = ((1.0 - (-alpha * ratio).exp()) / (1.0 - (-alpha).exp())) / ratio;
    assert!(close(calc_rie(&best, 0, alpha).unwrap(), rie_max));
    let rie_min = ((1.0 - (alpha * ratio).exp()) / (1.0 - alpha.exp())) / ratio;
    assert!(close(calc_rie(&worst, 0, alpha).unwrap(), rie_min));
    let empty: Vec<Vec<f64>> = Vec::new();
    assert!(matches!(
        calc_rie(&empty, 0, alpha),
        Err(Error::Value { .. })
    ));
    assert!(matches!(calc_rie(&best, 0, 0.0), Err(Error::Value { .. })));
    assert_eq!(calc_rie(&all_decoys, 0, alpha).unwrap(), 0.0);
}

#[test]
fn auc() {
    let (best, worst, all_actives, all_decoys) = cases();
    assert!(close(calc_auc(&best, 0).unwrap(), 1.0));
    assert!(close(calc_auc(&worst, 0).unwrap(), 0.0));
    let empty: Vec<Vec<f64>> = Vec::new();
    assert!(matches!(calc_auc(&empty, 0), Err(Error::Value { .. })));
    assert!(close(calc_auc(&all_actives, 0).unwrap(), 0.0));
    assert!(close(calc_auc(&all_decoys, 0).unwrap(), 0.0));
}

#[test]
fn bedroc() {
    let (best, worst, all_actives, all_decoys) = cases();
    assert!(close(calc_bedroc(&best, 0, 100.0).unwrap(), 1.0));
    assert!(close(calc_bedroc(&worst, 0, 100.0).unwrap(), 0.0));
    let empty: Vec<Vec<f64>> = Vec::new();
    assert!(matches!(
        calc_bedroc(&empty, 0, 100.0),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        calc_bedroc(&best, 0, 0.0),
        Err(Error::Value { .. })
    ));
    assert_eq!(calc_bedroc(&all_actives, 0, 100.0).unwrap(), 1.0);
    assert_eq!(calc_bedroc(&all_decoys, 0, 100.0).unwrap(), 0.0);
}

/// Scores with actives spread through the list; every value is Python's,
/// exactly.
fn spread() -> Vec<Vec<f64>> {
    (0..400)
        .map(|i| {
            let active = if (i * 7) % 11 < 3 { 1.0 } else { 0.0 };
            vec![i as f64 * 0.5, active]
        })
        .collect()
}

#[test]
fn matches_python_bit_for_bit() {
    let scores = spread();
    assert_eq!(calc_auc(&scores, 1).unwrap(), 0.5028847063274379);
    for (alpha, bedroc, rie) in [
        (20.0, 0.28356038148442175, 1.0361191284577835),
        (100.0, 0.3359889684581409, 1.2329870402115561),
        (1.5, 0.4472212953437237, 1.0029951048513608),
    ] {
        assert_eq!(calc_bedroc(&scores, 1, alpha).unwrap(), bedroc);
        assert_eq!(calc_rie(&scores, 1, alpha).unwrap(), rie);
    }
    assert_eq!(
        calc_enrichment(&scores, 1, &[0.01, 0.05, 0.1, 0.5, 0.0, 1.0]).unwrap(),
        [
            0.9174311926605505,
            1.1009174311926606,
            1.0091743119266054,
            1.0091743119266054,
            1.0041535442055776
        ]
    );
    let roc = calc_roc(&scores, 1).unwrap();
    let at = |v: &[f64]| [0, 1, 2, 57, 399].map(|i| v[i]);
    assert_eq!(
        at(&roc.fpr),
        [
            0.0,
            0.003436426116838488,
            0.006872852233676976,
            0.14432989690721648,
            1.0
        ]
    );
    assert_eq!(
        at(&roc.tpr),
        [
            0.009174311926605505,
            0.009174311926605505,
            0.009174311926605505,
            0.14678899082568808,
            1.0
        ]
    );

    let small: Vec<Vec<f64>> = [1., 0., 1., 1., 0., 0., 1., 0., 0., 0.]
        .iter()
        .map(|&a| vec![a])
        .collect();
    assert_eq!(calc_auc(&small, 0).unwrap(), 0.7916666666666666);
    assert_eq!(calc_bedroc(&small, 0, 20.0).unwrap(), 0.8829456820213197);
    assert_eq!(calc_rie(&small, 0, 20.0).unwrap(), 2.206625518817827);
    assert_eq!(
        calc_enrichment(&small, 0, &[0.1, 0.3, 0.35]).unwrap(),
        [2.5, 1.6666666666666667, 1.875]
    );
}

#[test]
fn short_rows_are_index_errors() {
    let scores = vec![vec![1.0], vec![]];
    assert!(matches!(calc_auc(&scores, 0), Err(Error::Index { .. })));
}

/// Python's `math.exp` raises OverflowError, and its division
/// ZeroDivisionError, where these would return inf or NaN.
#[test]
fn rie_and_bedroc_raise_where_python_does() {
    let one = vec![vec![2.0, 1.0]];
    assert!(matches!(calc_rie(&one, 1, 710.0), Err(Error::Value { .. })));
    assert!(matches!(
        calc_bedroc(&one, 1, 710.0),
        Err(Error::Value { .. })
    ));
    let many: Scores = (0..500).map(|i| vec![f64::from(i % 7 == 0)]).collect();
    assert!(calc_rie(&many, 0, 1000.0).is_ok());
    assert!(matches!(
        calc_bedroc(&many, 0, 1000.0),
        Err(Error::Value { .. })
    ));
    let (best, ..) = cases();
    for result in [calc_rie(&best, 0, 1e-300), calc_bedroc(&best, 0, 1e-300)] {
        assert!(
            matches!(&result, Err(Error::Value { message }) if message == "float division by zero"),
            "{result:?}"
        );
    }
    assert!(calc_bedroc(&best, 0, 709.0).is_ok());
}
