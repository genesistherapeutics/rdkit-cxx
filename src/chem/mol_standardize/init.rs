// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/Chem/MolStandardize/__init__.py`,
//! re-exported from `rdkit_cxx::chem::mol_standardize` as Python's is from
//! `rdkit.Chem.MolStandardize`.

use crate::Result;
use crate::chem::mol_standardize::rd_mol_standardize::TautomerEnumerator;
use crate::chem::{Mol, mol_to_smiles};

/// The molecule's tautomers, the canonical one, as a default
/// `TautomerEnumerator` scores them, first and the others after it in the
/// order of their SMILES. Mirrors `MolStandardize.ReorderTautomers`.
pub fn reorder_tautomers(molecule: &Mol) -> Result<Vec<Mol>> {
    let enumerator = TautomerEnumerator::default();
    let canon = enumerator.canonicalize(molecule)?;
    let csmi = mol_to_smiles(&canon)?;
    let mut stpl = Vec::new();
    for taut in enumerator.enumerate(molecule)?.tautomers() {
        let smi = mol_to_smiles(&taut)?;
        if smi != csmi {
            stpl.push((smi, taut));
        }
    }
    // Python sorts (SMILES, molecule) pairs, whose SMILES differ.
    stpl.sort_by(|a, b| a.0.cmp(&b.0));
    let mut res = vec![canon];
    res.extend(stpl.into_iter().map(|(_, taut)| taut));
    Ok(res)
}
