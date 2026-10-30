// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Identify Aromatic Rings (RDKitCB_8): flag a ring as aromatic when all of
//! its bonds are.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#identify-aromatic-rings),
//! by Benjamin Datko and Greg Landrum.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

fn is_ring_aromatic(mol: &Mol, bond_ring: &[u32]) -> Result<bool> {
    for &id in bond_ring {
        if !mol.get_bond_with_idx(id)?.get_is_aromatic() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn main() -> Result<()> {
    let m = mol_from_smiles("c1cccc2c1CCCC2")?;
    let ri = m.get_ring_info();

    // You can interrogate the RingInfo object to tell you the atoms that make
    // up each ring:
    let atom_rings = ri.atom_rings();
    println!("{atom_rings:?}");
    assert_eq!(atom_rings, [[0, 1, 2, 3, 4, 5], [4, 5, 6, 7, 8, 9]]);

    // or the bonds that make up each ring:
    let bond_rings = ri.bond_rings();
    println!("{bond_rings:?}");
    assert_eq!(bond_rings, [[0, 1, 2, 3, 4, 9], [4, 5, 6, 7, 8, 10]]);

    // To detect aromatic rings, loop over the bonds in each ring and flag the
    // ring as aromatic if all bonds are aromatic:
    let first = is_ring_aromatic(&m, &bond_rings[0])?;
    let second = is_ring_aromatic(&m, &bond_rings[1])?;
    println!("{first}");
    println!("{second}");
    assert!(first);
    assert!(!second);
    Ok(())
}
