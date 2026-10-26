// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.EState.EState`, a pure-Python module: Kier and Hall's
//! electrotopological state (EState) indices (Hall, Mohney, and Kier, JCICS
//! 31, 76-81 (1991)).
//!
//! Python keeps the indices it computes on the molecule and reads them back
//! unless forced; these never do, so their `force` arguments are not
//! ported and each call computes the indices afresh.

use crate::Result;
use crate::chem::descriptors::atomic_number;
use crate::chem::py_util::{py_max, py_min};
use crate::chem::{
    GetDistanceMatrixParams, Mol, get_distance_matrix_with_params, get_periodic_table,
};

/// The principal quantum number of an element's outer shell. Mirrors
/// `EState.GetPrincipleQuantumNumber`.
pub fn get_principle_quantum_number(at_num: i32) -> u32 {
    [2, 10, 18, 36, 54, 86]
        .iter()
        .position(|&last| at_num <= last)
        .map_or(7, |row| row as u32 + 1)
}

/// Each atom's EState index: its intrinsic state, from its valence
/// electrons, principal quantum number, and degree, plus the perturbations
/// of the atoms in its fragment. Mirrors `EState.EStateIndices`.
pub fn e_state_indices(mol: &Mol) -> Result<Vec<f64>> {
    let tbl = get_periodic_table();
    let n_atoms = mol.get_num_atoms() as usize;
    let mut is = vec![0.0; n_atoms];
    for (i, at) in mol.get_atoms().enumerate() {
        let d = at.get_degree();
        if d > 0 {
            let at_num = at.get_atomic_num();
            let dv = i64::from(tbl.get_n_outer_elecs(atomic_number(at_num))?)
                - i64::from(at.get_total_num_hs()?);
            let n = f64::from(get_principle_quantum_number(at_num));
            is[i] = (4.0 / (n * n) * dv as f64 + 1.0) / f64::from(d);
        }
    }
    let dists = get_distance_matrix_with_params(mol, &GetDistanceMatrixParams::default())?;
    let mut accum = vec![0.0; n_atoms];
    for i in 0..n_atoms {
        for j in i + 1..n_atoms {
            let p = dists[i][j] + 1.0;
            if p < 1e6 {
                let tmp = (is[i] - is[j]) / (p * p);
                accum[i] += tmp;
                accum[j] -= tmp;
            }
        }
    }
    Ok(accum.iter().zip(&is).map(|(a, i)| a + i).collect())
}

/// The largest EState index; `Error::Value` for a molecule without atoms,
/// as Python raises. Mirrors `EState.MaxEStateIndex`.
pub fn max_e_state_index(mol: &Mol) -> Result<f64> {
    py_max(e_state_indices(mol)?)
}

/// The smallest EState index. Mirrors `EState.MinEStateIndex`.
pub fn min_e_state_index(mol: &Mol) -> Result<f64> {
    py_min(e_state_indices(mol)?)
}

/// The largest absolute EState index. Mirrors `EState.MaxAbsEStateIndex`.
pub fn max_abs_e_state_index(mol: &Mol) -> Result<f64> {
    py_max(e_state_indices(mol)?.into_iter().map(f64::abs))
}

/// The smallest absolute EState index. Mirrors `EState.MinAbsEStateIndex`.
pub fn min_abs_e_state_index(mol: &Mol) -> Result<f64> {
    py_min(e_state_indices(mol)?.into_iter().map(f64::abs))
}
