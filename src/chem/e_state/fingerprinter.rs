// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.EState.Fingerprinter`, a pure-Python module: EState
//! fingerprints (Hall and Kier, JCICS 35, 1039-1045 (1995)).

use crate::Result;
use crate::chem::Mol;
use crate::chem::e_state::atom_types::es_patterns;
use crate::chem::e_state::e_state::e_state_indices;

/// The molecule's EState fingerprints: for each atom type of
/// [`es_patterns`], in order, the number of its unique matches, and the sum
/// of the EState indices of the atoms that start them. Mirrors
/// `Fingerprinter.FingerprintMol`.
pub fn fingerprint_mol(mol: &Mol) -> Result<(Vec<i64>, Vec<f64>)> {
    let es_indices = e_state_indices(mol)?;
    let patterns = es_patterns();
    let mut counts = vec![0; patterns.len()];
    let mut sums = vec![0.0; patterns.len()];
    for (i, (_, pattern)) in patterns.iter().enumerate() {
        let matches = mol.get_substruct_matches(pattern)?;
        counts[i] = matches.len() as i64;
        for m in &matches {
            sums[i] += es_indices[m[0] as usize];
        }
    }
    Ok((counts, sums))
}
