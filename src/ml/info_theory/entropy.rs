// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.InfoTheory.entropy`, a pure-Python module: informational
//! entropy and gain, as Mitchell's "Machine Learning" defines them.
//!
//! Python's module also exports `InfoEntropy` and `InfoGain`, which are
//! `rdInfoTheory`'s ([`info_entropy`](super::info_entropy),
//! [`info_gain`](super::info_gain)); these are its pure-Python versions.

use crate::ml::info_theory::rd_info_theory;
use crate::{Error, Result};

/// Python's `sum(values)`.
fn seq_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0, |acc, x| acc + x)
}

/// The informational entropy, in bits, of the counts of each result; 0 for
/// no instances. Mirrors `entropy.PyInfoEntropy`.
///
/// NumPy's `log` may round differently from the C library's on some
/// processors, so results can differ from Python's in the last bit.
pub fn py_info_entropy(results: &[f64]) -> f64 {
    let n_instances = seq_sum(results.iter().copied());
    if n_instances == 0.0 {
        return 0.0;
    }
    let log2 = 2f64.ln();
    seq_sum(results.iter().map(|&r| {
        let p = r / n_instances;
        let t = if p > 0.0 { p } else { 1.0 };
        -p * t.ln() / log2
    }))
}

/// The expected information gain of a variable, from the counts of each
/// result (columns) for each of its values (rows); 0 for no instances.
/// Mirrors `entropy.PyInfoGain`, which takes the entropies with
/// `rdInfoTheory.InfoEntropy`.
///
/// `Error::Value` for rows of different lengths; errors from
/// [`rd_info_theory::info_entropy`] pass through.
pub fn py_info_gain(var_mat: &[Vec<f64>]) -> Result<f64> {
    let ncols = var_mat.first().map_or(0, Vec::len);
    if var_mat.iter().any(|row| row.len() != ncols) {
        return Err(Error::Value {
            message: "the rows have different lengths".to_owned(),
        });
    }
    let mut overall_res = vec![0.0; ncols];
    for row in var_mat {
        for (o, &x) in overall_res.iter_mut().zip(row) {
            *o += x;
        }
    }
    let mut term2 = 0.0;
    for row in var_mat {
        let variable_res = crate::chem::py_math::np_sum(row);
        term2 += variable_res * rd_info_theory::info_entropy(row)?;
    }
    let t_sum = seq_sum(overall_res.iter().copied());
    if t_sum != 0.0 {
        term2 *= 1.0 / t_sum;
        Ok(rd_info_theory::info_entropy(&overall_res)? - term2)
    } else {
        Ok(0.0)
    }
}
