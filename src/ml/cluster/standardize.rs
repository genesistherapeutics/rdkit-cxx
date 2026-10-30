// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Cluster.Standardize`, a Python module: the ways of
//! standardizing a data matrix before clustering it.

use crate::Result;
use crate::ml::data::stats;

/// A standardization: it takes a data matrix (rows of points) and returns
/// the standardized one. A Rust addition, naming the functions of `methods`.
pub type Standardization = fn(&[Vec<f64>]) -> Result<Vec<Vec<f64>>>;

/// Standardizes each column by its mean and sample standard deviation, with
/// [`stats::standardize_matrix`]. Mirrors `Standardize.StdDev`.
pub fn std_dev(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    stats::standardize_matrix(mat)
}

fn no_standardization(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    Ok(mat.to_vec())
}

/// The standardizations offered to users, as `(name, function,
/// description)`. Mirrors `Standardize.methods`.
pub fn methods() -> Vec<(&'static str, Standardization, &'static str)> {
    vec![
        ("None", no_standardization, "No Standardization"),
        ("Standard Deviation", std_dev, "Use the standard deviation"),
    ]
}
