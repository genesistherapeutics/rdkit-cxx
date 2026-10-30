// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.Resemblance`, a Python module: resemblance
//! (metric) matrices.
//!
//! A metric matrix of `n` objects is stored as its upper triangle, column
//! by column: the value for `row < col` is at `col * (col - 1) / 2 + row`,
//! so the order is (0, 1), (0, 2), (1, 2), (0, 3), and so on.

use crate::ml::cluster::clusters::Cluster;
use crate::py::format_f;
use crate::{Error, Result};

/// A metric: it takes points and returns their metric matrix. A Rust
/// addition, for the entries of [`methods`]; `calc_metric_matrix` takes any
/// function or closure.
pub type Metric = fn(&[Vec<f64>]) -> Result<Vec<f64>>;

/// The Euclidean distances between the points, as a metric matrix. Mirrors
/// `Resemblance.EuclideanDistance`.
///
/// Python's raises `TypeError` in RDKit 2024.09.6 (it sizes its array with
/// a float); this computes what it means to, summing the squared
/// differences of each pair as Python's `sum` does.
///
/// `Error::Value` for points of different lengths.
pub fn euclidean_distance(in_data: &[Vec<f64>]) -> Result<Vec<f64>> {
    let n_objs = in_data.len();
    let mut res = Vec::with_capacity(n_objs * n_objs.saturating_sub(1) / 2);
    for col in 1..n_objs {
        for row in 0..col {
            let (a, b) = (&in_data[row], &in_data[col]);
            if a.len() != b.len() {
                return Err(Error::Value {
                    message: "the points have different lengths".to_owned(),
                });
            }
            let ss = a
                .iter()
                .zip(b)
                .map(|(x, y)| (x - y) * (x - y))
                .fold(0.0, |acc, t| acc + t);
            res.push(ss.sqrt());
        }
    }
    Ok(res)
}

/// The metric matrix of the clusters' positions. Mirrors
/// `Resemblance.CalcMetricMatrix`.
///
/// Python's raises `TypeError` in RDKit 2024.09.6 (it hands the metric a
/// `map` object, which has no length); this hands it the positions.
pub fn calc_metric_matrix(
    in_data: &[Cluster],
    metric_func: impl Fn(&[Vec<f64>]) -> Result<Vec<f64>>,
) -> Result<Vec<f64>> {
    let positions: Vec<Vec<f64>> = in_data.iter().map(|c| c.get_position().to_vec()).collect();
    metric_func(&positions)
}

fn check_length(mat: &[f64], n_objs: usize) -> Result<()> {
    if Some(mat.len()) != n_objs.checked_mul(n_objs.saturating_sub(1)).map(|n| n / 2) {
        return Err(Error::Value {
            message: "bad matrix length in FindMinValInList".to_owned(),
        });
    }
    Ok(())
}

/// The smallest value in a metric matrix of `n_objs` objects, as `(row,
/// col, value)`: the first smallest, with NaN smallest of all, as NumPy's
/// `argmin` finds it. Mirrors `Resemblance.FindMinValInList(mat, nObjs)`.
///
/// `Error::Value` when `mat` is not `n_objs * (n_objs - 1) / 2` long (Python's
/// assertion), or empty.
pub fn find_min_val_in_list(mat: &[f64], n_objs: usize) -> Result<(usize, usize, f64)> {
    check_length(mat, n_objs)?;
    let mut min_idx = None;
    for (i, &x) in mat.iter().enumerate() {
        match min_idx {
            None => min_idx = Some(i),
            Some(m) => {
                let best: f64 = mat[m];
                if !best.is_nan() && (x.is_nan() || x < best) {
                    min_idx = Some(i);
                }
            }
        }
    }
    let min_idx = min_idx.ok_or_else(|| Error::Value {
        message: "attempt to get argmin of an empty sequence".to_owned(),
    })?;
    find_min_val_in_list_with_min_idx(mat, n_objs, min_idx)
}

/// The row, column, and value of the item at `min_idx` of a metric matrix
/// of `n_objs` objects. Mirrors `Resemblance.FindMinValInList(mat, nObjs,
/// minIdx)`.
///
/// `Error::Value` as [`find_min_val_in_list`]; `Error::Index` for an index
/// past the matrix.
pub fn find_min_val_in_list_with_min_idx(
    mat: &[f64],
    n_objs: usize,
    min_idx: usize,
) -> Result<(usize, usize, f64)> {
    check_length(mat, n_objs)?;
    let value = *mat.get(min_idx).ok_or(Error::index(min_idx, mat.len()))?;
    let mut n_so_far = 0;
    let mut col = 0;
    while n_so_far <= min_idx {
        col += 1;
        n_so_far += col;
    }
    let row = min_idx + col - n_so_far;
    Ok((row, col, value))
}

/// The metric matrix as a table, each value `%10.6f` above the diagonal and
/// `   ---    ` on and below it, one line per row. Mirrors
/// `Resemblance.ShowMetricMat`, which prints the text this returns.
///
/// Python's raises `TypeError` in RDKit 2024.09.6 for more than one object
/// (it indexes the matrix with a float); this computes what it means to.
///
/// `Error::Value` as [`find_min_val_in_list`].
pub fn show_metric_mat(metric_mat: &[f64], n_objs: usize) -> Result<String> {
    check_length(metric_mat, n_objs)?;
    let mut out = String::new();
    for row in 0..n_objs {
        for col in 0..n_objs {
            if col <= row {
                out.push_str("   ---    ");
            } else {
                out.push_str(&format!(
                    "{:>10}",
                    format_f(metric_mat[col * (col - 1) / 2 + row])
                ));
            }
        }
        out.push('\n');
    }
    Ok(out)
}

/// The metrics offered to users, as `(name, function, description)`.
/// Mirrors `Resemblance.methods`.
pub fn methods() -> Vec<(&'static str, Metric, &'static str)> {
    vec![("Euclidean", euclidean_distance, "Euclidean Distance")]
}
