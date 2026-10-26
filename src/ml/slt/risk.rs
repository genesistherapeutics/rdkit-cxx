// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.SLT.Risk`, a Python module: bounds on the risk of
//! a classifier from statistical learning theory.
//!
//! Each bound takes the VC dimension, the number of data points, the
//! number misclassified, and the confidence, and computes as Python does,
//! so the results match Python's to the last bit.

/// The base-2 logarithm, as `ln(x) / ln(2)`. Mirrors `Risk.log2`.
pub fn log2(x: f64) -> f64 {
    x.ln() / 2f64.ln()
}

/// Burges's risk bound (Eqn. 3 of "A Tutorial on Support Vector Machines
/// for Pattern Recognition", Data Mining and Knowledge Discovery 2, 1998).
/// Mirrors `Risk.BurgesRiskBound`.
pub fn burges_risk_bound(vc_dim: f64, n_data: f64, n_wrong: f64, conf: f64) -> f64 {
    let (h, l, eta) = (vc_dim, n_data, conf);
    let numerator = h * ((2.0 * l / h).ln() + 1.0) - (eta / 4.0).ln();
    let struct_risk = (numerator / l).sqrt();
    let r_emp = n_wrong / l;
    r_emp + struct_risk
}

/// Cristianini and Shawe-Taylor's risk bound (Theorem 4.6 of "An
/// Introduction to Support Vector Machines", 2000). Mirrors
/// `Risk.CristianiRiskBound`.
pub fn cristiani_risk_bound(vc_dim: f64, n_data: f64, n_wrong: f64, conf: f64) -> f64 {
    let (d, delta, l, k) = (vc_dim, conf, n_data, n_wrong);
    let struct_risk = ((4.0 / n_data)
        * (d * log2((2.0 * std::f64::consts::E * l) / d) + log2(4.0 / delta)))
    .sqrt();
    let r_emp = 2.0 * k / l;
    r_emp + struct_risk
}

/// Parameters for [`cherkassky_risk_bound_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CherkasskyRiskBoundParams {
    /// A constant of the bound, from 0 to 4.
    pub a1: f64,
    /// A constant of the bound, from 0 to 2.
    pub a2: f64,
}

impl Default for CherkasskyRiskBoundParams {
    fn default() -> CherkasskyRiskBoundParams {
        CherkasskyRiskBoundParams { a1: 1.0, a2: 2.0 }
    }
}

/// Cherkassky and Mulier's risk bound (Eqns. 4.22 and 4.23 of "Learning
/// From Data", 1998). Mirrors `Risk.CherkasskyRiskBound(VCDim, nData,
/// nWrong, conf)`.
pub fn cherkassky_risk_bound(vc_dim: f64, n_data: f64, n_wrong: f64, conf: f64) -> f64 {
    cherkassky_risk_bound_with_params(
        vc_dim,
        n_data,
        n_wrong,
        conf,
        &CherkasskyRiskBoundParams::default(),
    )
}

/// Mirrors `Risk.CherkasskyRiskBound(VCDim, nData, nWrong, conf, a1, a2)`.
pub fn cherkassky_risk_bound_with_params(
    vc_dim: f64,
    n_data: f64,
    n_wrong: f64,
    conf: f64,
    params: &CherkasskyRiskBoundParams,
) -> f64 {
    let (h, n, eta) = (vc_dim, n_data, conf);
    let r_emp = n_wrong / n_data;
    let numerator = h * ((params.a2 * n / h).ln() + 1.0) - (eta / 4.0).ln();
    let eps = params.a1 * numerator / n;
    let struct_risk = eps / 2.0 * (1.0 + (1.0 + (4.0 * r_emp / eps)).sqrt());
    r_emp + struct_risk
}
