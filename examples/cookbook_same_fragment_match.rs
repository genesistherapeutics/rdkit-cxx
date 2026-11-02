// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Within the Same Fragment (RDKitCB_20): match a pattern only within one
//! fragment.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#within-the-same-fragment),
//! by Greg Landrum.

use std::collections::BTreeSet;

use rdkit_cxx::Result;
use rdkit_cxx::chem::{Mol, get_mol_frags, mol_from_smarts, mol_from_smiles};

/// The first match that lies within a single fragment, if any.
fn fragsearch(m: &Mol, p: &Mol) -> Result<Option<BTreeSet<u32>>> {
    let matches: Vec<BTreeSet<u32>> = m
        .get_substruct_matches(p)?
        .into_iter()
        .map(|x| x.into_iter().collect())
        .collect();
    let frags: Vec<BTreeSet<u32>> = get_mol_frags(m)
        .into_iter()
        .map(|y| y.into_iter().collect())
        .collect();
    for frag in &frags {
        for m in &matches {
            if m.is_subset(frag) {
                return Ok(Some(m.clone()));
            }
        }
    }
    Ok(None)
}

fn main() -> Result<()> {
    let p = mol_from_smarts("O.N")?;
    let m1 = mol_from_smiles("OCCCN.CCC")?;
    let m2 = mol_from_smiles("OCCC.CCCN")?;

    println!("{}", m1.has_substruct_match(&p)?);
    println!("{}", m2.has_substruct_match(&p)?);
    assert!(m1.has_substruct_match(&p)?);
    assert!(m2.has_substruct_match(&p)?);

    let found1 = fragsearch(&m1, &p)?;
    let found2 = fragsearch(&m2, &p)?;
    println!("{found1:?}");
    println!("{found2:?}");
    assert_eq!(found1, Some(BTreeSet::from([0, 4])));
    assert_eq!(found2, None);
    Ok(())
}
