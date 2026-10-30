// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Wiener Index (RDKitCB_12): the sum of the topological distances between
//! all pairs of atoms.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#wiener-index),
//! by Greg Landrum.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{Mol, get_distance_matrix, mol_from_smiles};

fn wiener_index(m: &Mol) -> Result<f64> {
    let amat = get_distance_matrix(m)?;
    // The distances above the diagonal: each pair once.
    let mut res = 0.0;
    for (i, row) in amat.iter().enumerate() {
        res += row[i + 1..].iter().sum::<f64>();
    }
    Ok(res)
}

fn main() -> Result<()> {
    let butane = wiener_index(&mol_from_smiles("CCCC")?)?;
    println!("{butane:?}");
    assert_eq!(butane, 10.0);

    let isobutane = wiener_index(&mol_from_smiles("CC(C)C")?)?;
    println!("{isobutane:?}");
    assert_eq!(isobutane, 9.0);
    Ok(())
}
