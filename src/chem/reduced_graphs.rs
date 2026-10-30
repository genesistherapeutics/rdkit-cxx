// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.ReducedGraphs`, a Python module that re-exports
//! `rdkit.Chem.rdReducedGraphs`, as this does, and adds a similarity for
//! the ErG fingerprints it makes.

pub use crate::chem::rd_reduced_graphs::*;

use crate::{Error, Result};

/// The Tanimoto similarity of two real-valued vectors, such as ErG
/// fingerprints: their dot product over the sum of their squared norms
/// less it, or 0 when either is 0. `Error::Value` for vectors of different
/// lengths, as NumPy raises. Mirrors `ReducedGraphs.TanimotoSimilarity`.
pub fn tanimoto_similarity(arr1: &[f64], arr2: &[f64]) -> Result<f64> {
    if arr1.len() != arr2.len() {
        return Err(Error::Value {
            message: format!(
                "the vectors' lengths differ: {} and {}",
                arr1.len(),
                arr2.len()
            ),
        });
    }
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let numer = dot(arr1, arr2);
    if numer == 0.0 {
        return Ok(0.0);
    }
    let denom = dot(arr1, arr1) + dot(arr2, arr2) - numer;
    if denom == 0.0 {
        return Ok(0.0);
    }
    Ok(numer / denom)
}
