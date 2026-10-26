// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Scoring.Scoring`, a pure-Python module: rank statistics
//! for virtual screening (ROC curves and their area, BEDROC, RIE, and
//! enrichment factors).
//!
//! `scores` lists the samples best first; column `col` of each row is
//! nonzero for an active, as Python tests the column's truth. The
//! arithmetic follows Python's, operation by operation, so the results
//! match Python's to the last bit.
//!
//! ```no_run
//! use rdkit_cxx::ml::scoring::scoring::calc_bedroc;
//!
//! let scores = vec![vec![1.0], vec![0.0], vec![1.0], vec![0.0]];
//! let bedroc = calc_bedroc(&scores, 0, 20.0).unwrap();
//! ```

use crate::{Error, Result};

/// A ROC curve: the false- and true-positive rates after each sample.
/// Mirrors the `RocCurve` named tuple `Scoring.CalcROC` returns.
#[derive(Clone, Debug, PartialEq)]
pub struct RocCurve {
    /// False-positive rates, `FP / (TN + FP)`.
    pub fpr: Vec<f64>,
    /// True-positive rates, `TP / (TP + FN)`.
    pub tpr: Vec<f64>,
}

fn value_error(message: &str) -> Error {
    Error::Value {
        message: message.to_owned(),
    }
}

/// Whether row `i` is an active: its column `col` is nonzero (Python's truth
/// test, so NaN counts as active). `Error::Index` for a short row.
fn is_active<S: AsRef<[f64]>>(scores: &[S], i: usize, col: usize) -> Result<bool> {
    let row = scores[i].as_ref();
    match row.get(col) {
        Some(&value) => Ok(value != 0.0),
        None => Err(Error::Index {
            index: u32::try_from(col).unwrap_or(u32::MAX),
            len: u32::try_from(row.len()).unwrap_or(u32::MAX),
        }),
    }
}

/// The ROC curve of the ranked samples. Mirrors `Scoring.CalcROC`.
///
/// `Error::Value` for no samples, as Python raises; `Error::Index` for a row
/// without column `col`.
pub fn calc_roc<S: AsRef<[f64]>>(scores: &[S], col: usize) -> Result<RocCurve> {
    let num_mol = scores.len();
    if num_mol == 0 {
        return Err(value_error("score list is empty"));
    }
    let mut tp = Vec::with_capacity(num_mol);
    let mut fp = Vec::with_capacity(num_mol);
    let (mut num_actives, mut num_inactives) = (0u64, 0u64);
    for i in 0..num_mol {
        if is_active(scores, i, col)? {
            num_actives += 1;
        } else {
            num_inactives += 1;
        }
        tp.push(num_actives);
        fp.push(num_inactives);
    }
    // Python divides only when the denominator is positive; the counts are
    // then all zero, which as floats are the rates.
    let rates = |counts: Vec<u64>, total: u64| -> Vec<f64> {
        if total > 0 {
            counts
                .into_iter()
                .map(|n| 1.0 * n as f64 / total as f64)
                .collect()
        } else {
            counts.into_iter().map(|n| n as f64).collect()
        }
    };
    Ok(RocCurve {
        fpr: rates(fp, num_inactives),
        tpr: rates(tp, num_actives),
    })
}

/// The area under the ROC curve, by the trapezoidal rule. Mirrors
/// `Scoring.CalcAUC`.
///
/// Errors as [`calc_roc`].
pub fn calc_auc<S: AsRef<[f64]>>(scores: &[S], col: usize) -> Result<f64> {
    let roc = calc_roc(scores, col)?;
    let mut auc = 0.0;
    for i in 0..scores.len() - 1 {
        auc += (roc.fpr[i + 1] - roc.fpr[i]) * (roc.tpr[i + 1] + roc.tpr[i]);
    }
    Ok(0.5 * auc)
}

/// RIE and the number of actives. Python's `_RIEHelper`.
fn rie_helper<S: AsRef<[f64]>>(scores: &[S], col: usize, alpha: f64) -> Result<(f64, u64)> {
    let num_mol = scores.len();
    if num_mol == 0 {
        return Err(value_error("score list is empty"));
    }
    if alpha <= 0.0 {
        return Err(value_error("alpha must be greater than zero"));
    }
    let n = num_mol as f64;
    let denom = 1.0 / n * ((1.0 - (-alpha).exp()) / ((alpha / n).exp() - 1.0));
    let mut num_actives = 0u64;
    let mut sum_exp = 0.0;
    for i in 0..num_mol {
        if is_active(scores, i, col)? {
            num_actives += 1;
            sum_exp += (-(alpha * (i + 1) as f64) / n).exp();
        }
    }
    let rie = if num_actives > 0 {
        sum_exp / (num_actives as f64 * denom)
    } else {
        0.0
    };
    Ok((rie, num_actives))
}

/// The robust initial enhancement of the ranked samples, with exponential
/// weight `alpha` (Sheridan et al., J. Chem. Inf. Comput. Sci. 41,
/// 1395-1406, 2001). Mirrors `Scoring.CalcRIE`.
///
/// `Error::Value` for no samples or an `alpha` that is not positive, as
/// Python raises; `Error::Index` for a row without column `col`.
pub fn calc_rie<S: AsRef<[f64]>>(scores: &[S], col: usize, alpha: f64) -> Result<f64> {
    Ok(rie_helper(scores, col, alpha)?.0)
}

/// The Boltzmann-enhanced discrimination of the ROC curve, with exponential
/// weight `alpha` (Truchon and Bayly, J. Chem. Inf. Model. 47, 488-508,
/// 2007). Mirrors `Scoring.CalcBEDROC`.
///
/// Errors as [`calc_rie`].
pub fn calc_bedroc<S: AsRef<[f64]>>(scores: &[S], col: usize, alpha: f64) -> Result<f64> {
    let (rie, num_actives) = rie_helper(scores, col, alpha)?;
    if num_actives == 0 {
        return Ok(0.0);
    }
    let ratio = 1.0 * num_actives as f64 / scores.len() as f64;
    let rie_max = (1.0 - (-alpha * ratio).exp()) / (ratio * (1.0 - (-alpha).exp()));
    let rie_min = (1.0 - (alpha * ratio).exp()) / (ratio * (1.0 - alpha.exp()));
    Ok(if rie_max != rie_min {
        (rie - rie_min) / (rie_max - rie_min)
    } else {
        // Every sample is active.
        1.0
    })
}

/// The enrichment factor at each fraction of the ranked samples. Mirrors
/// `Scoring.CalcEnrichment`.
///
/// As in Python, a fraction whose cutoff (`ceil(n * fraction)` samples)
/// reaches the last sample adds no value, and every value is zero when no
/// sample is active. `Error::Value` for no samples, no fractions, or a
/// fraction outside 0 to 1, as Python raises; `Error::Index` for a row
/// without column `col`.
pub fn calc_enrichment<S: AsRef<[f64]>>(
    scores: &[S],
    col: usize,
    fractions: &[f64],
) -> Result<Vec<f64>> {
    let num_mol = scores.len();
    if num_mol == 0 {
        return Err(value_error("score list is empty"));
    }
    if fractions.is_empty() {
        return Err(value_error("fraction list is empty"));
    }
    // Python's check passes NaN, for which its `math.ceil` then raises.
    if fractions.iter().any(|f| !(0.0..=1.0).contains(f)) {
        return Err(value_error("fractions must be between [0,1]"));
    }
    let mut num_per_frac: Vec<f64> = fractions
        .iter()
        .map(|&f| (num_mol as f64 * f).ceil())
        .collect();
    num_per_frac.push(num_mol as f64);
    // Python pops the front; an index into the cutoffs does the same.
    let mut next = 0;
    let mut num_actives = 0u64;
    let mut enrich = Vec::new();
    for i in 0..num_mol {
        if i as f64 > num_per_frac[next] - 1.0 && i > 0 {
            enrich.push(1.0 * num_actives as f64 * num_mol as f64 / i as f64);
            next += 1;
        }
        if is_active(scores, i, col)? {
            num_actives += 1;
        }
    }
    Ok(if num_actives > 0 {
        enrich.into_iter().map(|e| e / num_actives as f64).collect()
    } else {
        vec![0.0; fractions.len()]
    })
}
