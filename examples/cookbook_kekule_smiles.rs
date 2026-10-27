// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Kekule SMILES (RDKitCB_4): kekulize a molecule and write Kekule SMILES.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#kekule-smiles),
//! by Paolo Tosco.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    SmilesWriteParams, kekulize, mol_from_smiles, mol_to_smiles, mol_to_smiles_with_params,
};

fn main() -> Result<()> {
    let mut mol = mol_from_smiles("CN1C(NC2=NC=CC=C2)=CC=C1")?;
    println!("{}", mol_to_smiles(&mol)?);
    assert_eq!(mol_to_smiles(&mol)?, "Cn1cccc1Nc1ccccn1");

    kekulize(&mut mol)?;
    let params = SmilesWriteParams {
        do_kekule: true,
        ..Default::default()
    };
    let kekule = mol_to_smiles_with_params(&mol, &params)?;
    println!("{kekule}");
    assert_eq!(kekule, "CN1C=CC=C1NC1=NC=CC=C1");
    Ok(())
}
