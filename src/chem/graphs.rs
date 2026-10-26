// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Graphs`, a pure-Python module of functions on
//! molecular graphs.

use crate::chem::Mol;
use crate::{Error, Result};

/// The coefficients of the characteristic polynomial of the matrix `mat`,
/// an `n` by `n` matrix for a molecule of `n` atoms, highest power first,
/// by the Le Verrier-Faddeev-Frame method (Trinajstić, _Chemical Graph
/// Theory_, 2nd edition, p. 76). Python's `mat` is optional, but without it
/// Python fails (`A` is never set), so it is required here. A matrix of
/// another shape is `Error::Value`. Mirrors `Graphs.CharacteristicPolynomial`.
pub fn characteristic_polynomial(mol: &Mol, mat: &[Vec<f64>]) -> Result<Vec<f64>> {
    let n_atoms = mol.get_num_atoms() as usize;
    if mat.len() != n_atoms || mat.iter().any(|row| row.len() != n_atoms) {
        return Err(Error::Value {
            message: format!("the matrix must be {n_atoms} by {n_atoms}, one row per atom"),
        });
    }
    let a = mat;
    let mut an: Vec<Vec<f64>> = a.to_vec();
    let mut res = vec![0.0; n_atoms + 1];
    res[0] = 1.0;
    for n in 1..=n_atoms {
        let trace: f64 = (0..n_atoms).map(|i| an[i][i]).sum();
        res[n] = 1.0 / n as f64 * trace;
        let mut bn = an.clone();
        for (i, row) in bn.iter_mut().enumerate() {
            row[i] -= res[n];
        }
        an = (0..n_atoms)
            .map(|i| {
                (0..n_atoms)
                    .map(|j| (0..n_atoms).map(|k| a[i][k] * bn[k][j]).sum())
                    .collect()
            })
            .collect();
    }
    for x in &mut res[1..] {
        *x *= -1.0;
    }
    Ok(res)
}
