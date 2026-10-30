// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/ML/SLT/UnitTestRisk.py, with values from Python's
// rdkit.ML.SLT.Risk in RDKit 2025.09.6.

use rdkit_cxx::ml::slt::risk::*;

const D_LIST: [(f64, f64); 7] = [
    (1., 40.),
    (2., 40.),
    (3., 21.),
    (7., 16.),
    (8., 12.),
    (9., 11.),
    (10., 11.),
];
const N_PTS: f64 = 95.0;

fn check(res: impl Iterator<Item = f64>, target: [f64; 7]) {
    for (r, t) in res.zip(target) {
        assert!((r - t).abs() < 1e-4, "{r} {t}");
    }
}

#[test]
fn burges() {
    check(
        D_LIST
            .iter()
            .map(|&(h, w)| burges_risk_bound(h, N_PTS, w, 0.1)),
        [0.7445, 0.8157, 0.6698, 0.7649, 0.7506, 0.7658, 0.7896],
    );
}

#[test]
fn cherkassky() {
    let eps = (4.0 / N_PTS.sqrt()).min(1.0);
    check(
        D_LIST
            .iter()
            .map(|&(h, w)| cherkassky_risk_bound(h, N_PTS, w, eps)),
        [0.6654, 0.7450, 0.5378, 0.6329, 0.6010, 0.6175, 0.6501],
    );
}

#[test]
fn cristiani() {
    let eps = (4.0 / N_PTS.sqrt()).min(1.0);
    check(
        D_LIST
            .iter()
            .map(|&(h, w)| cristiani_risk_bound(h, N_PTS, w, eps)),
        [1.5617, 1.7438, 1.4797, 1.7394, 1.7235, 1.7653, 1.8235],
    );
}

#[test]
fn matches_python_bit_for_bit() {
    assert_eq!(burges_risk_bound(4.0, 100.0, 5.0, 0.05), 0.540205249425039);
    assert_eq!(
        cristiani_risk_bound(4.0, 100.0, 5.0, 0.05),
        1.2775930199783305
    );
    assert_eq!(
        cherkassky_risk_bound(4.0, 100.0, 5.0, 0.05),
        0.33278891716100323
    );
    let params = CherkasskyRiskBoundParams { a1: 2.0, a2: 1.5 };
    assert_eq!(
        cherkassky_risk_bound_with_params(4.0, 100.0, 5.0, 0.05, &params),
        0.5530675639418245
    );
    assert_eq!(
        CherkasskyRiskBoundParams::default(),
        CherkasskyRiskBoundParams { a1: 1.0, a2: 2.0 }
    );
}

#[test]
fn log2_values() {
    assert_eq!(log2(10.0), 3.3219280948873626);
    for (x, l) in [
        (1.0, 0.0),
        (2.0, 1.0),
        (4.0, 2.0),
        (0.5, -1.0),
        (0.25, -2.0),
        (0.125, -3.0),
    ] {
        assert_eq!(log2(x), l);
    }
}
