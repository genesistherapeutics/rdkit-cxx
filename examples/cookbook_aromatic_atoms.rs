// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Identify Aromatic Atoms (RDKitCB_9): tell aromatic carbons from olefinic
//! ones, with SMARTS or with `rdqueries`.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#identify-aromatic-atoms),
//! by Paolo Tosco.

use rdkit_cxx::Result;
use rdkit_cxx::chem::rdqueries::{
    hybridization_equals_query_atom, is_aliphatic_query_atom, is_aromatic_query_atom,
};
use rdkit_cxx::chem::{HybridizationType, mol_from_smarts, mol_from_smiles};

fn main() -> Result<()> {
    let mol = mol_from_smiles("c1ccccc1C=CCC")?;
    let aromatic_carbon = mol_from_smarts("c")?;
    let matches = mol.get_substruct_matches(&aromatic_carbon)?;
    println!("{matches:?}");
    assert_eq!(matches, [[0], [1], [2], [3], [4], [5]]);

    // The RDKit includes a SMARTS extension that allows hybridization
    // queries; here we query for SP2 aliphatic carbons:
    let olefinic_carbon = mol_from_smarts("[C^2]")?;
    let matches = mol.get_substruct_matches(&olefinic_carbon)?;
    println!("{matches:?}");
    assert_eq!(matches, [[6], [7]]);

    // There is also an alternative, more efficient approach, using the
    // rdqueries module:
    let q = is_aromatic_query_atom();
    let idxs: Vec<u32> = mol
        .get_atoms_matching_query(&q)?
        .iter()
        .map(|x| x.get_idx())
        .collect();
    println!("{idxs:?}");
    assert_eq!(idxs, [0, 1, 2, 3, 4, 5]);

    // HybridizationEqualsQueryAtom takes the HybridizationType value.
    let q = hybridization_equals_query_atom(HybridizationType::Sp2 as i32)?;
    let idxs: Vec<u32> = mol
        .get_atoms_matching_query(&q)?
        .iter()
        .map(|x| x.get_idx())
        .collect();
    println!("{idxs:?}");
    assert_eq!(idxs, [0, 1, 2, 3, 4, 5, 6, 7]);

    let mut qcombined = is_aliphatic_query_atom();
    qcombined.expand_query(&q)?;
    let idxs: Vec<u32> = mol
        .get_atoms_matching_query(&qcombined)?
        .iter()
        .map(|x| x.get_idx())
        .collect();
    println!("{idxs:?}");
    assert_eq!(idxs, [6, 7]);
    Ok(())
}
