// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Hybridization Type and Count (RDKitCB_26): get each atom's hybridization
//! and count the SP3 atoms.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#hybridization-type-and-count),
//! by Jean-Marc Nuzillard and Andrew Dalke.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{HybridizationType, mol_from_smiles};

fn main() -> Result<()> {
    let m = mol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")?;
    let hybridizations: Vec<HybridizationType> =
        m.get_atoms().map(|x| x.get_hybridization()).collect();
    for (idx, hybridization) in hybridizations.iter().enumerate() {
        println!("{idx} {hybridization:?}");
    }
    use HybridizationType::{Sp2, Sp3};
    assert_eq!(
        hybridizations,
        [
            Sp3, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp3, Sp3
        ]
    );

    // If you want to count hybridization type (e.g., SP3):
    let num_sp3 = m
        .get_atoms()
        .filter(|x| x.get_hybridization() == Sp3)
        .count();
    println!("{num_sp3}");
    assert_eq!(num_sp3, 3);
    Ok(())
}
