// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Bonds and Bonding" chapter.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{HybridizationType, mol_from_smiles};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[("RDKitCB_26", "Hybridization Type and Count", hybridization)];

/// Lists each atom's hybridization and counts the sp3 atoms. By Jean-Marc
/// Nuzillard and Andrew Dalke.
fn hybridization() -> Result<()> {
    let m = mol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")?;
    let hybridizations: Vec<HybridizationType> =
        m.get_atoms().map(|atom| atom.get_hybridization()).collect();
    for (idx, hybridization) in hybridizations.iter().enumerate() {
        println!("{idx} {hybridization:?}");
    }
    use HybridizationType::{Sp2, Sp3};
    assert_eq!(
        hybridizations,
        [Sp3, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp2, Sp3, Sp3]
    );

    let num_sp3 = m
        .get_atoms()
        .filter(|atom| atom.get_hybridization() == Sp3)
        .count();
    println!("{num_sp3}");
    assert_eq!(num_sp3, 3);
    Ok(())
}
