// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Sidechain-Core Enumeration (RDKitCB_29): replace sidechains on a core and
//! enumerate the combinations.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#sidechain-core-enumeration),
//! by Chris Earnshaw, Stephen Roughley, Greg Landrum, and Vincent Scalfani.

use rdkit_cxx::Result;
use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles, mol_to_smiles, replace_substructs};

fn main() -> Result<()> {
    // core is '*c1c(C)cccc1(O)'
    // chain is 'CN*'
    let mut rxn = reaction_from_smarts("[c:1][#0].[#0][*:2]>>[c:1]-[*:2]")?;
    let reacts = [mol_from_smiles("*c1c(C)cccc1(O)")?, mol_from_smiles("CN*")?];
    let products = rxn.run_reactants(&reacts)?;
    println!("{}", products.len());
    println!("{}", products[0].len());
    println!("{}", mol_to_smiles(&products[0][0])?);
    assert_eq!(products.len(), 1);
    assert_eq!(products[0].len(), 1);
    assert_eq!(mol_to_smiles(&products[0][0])?, "CNc1c(C)cccc1O");

    // The reaction-based approach is flexible. If you can generate your
    // sidechains so that the atom to attach to the core is atom zero, there
    // is an easier way to do this kind of simple replacement:
    let core = mol_from_smiles("*c1c(C)cccc1(O)")?;
    let chain = mol_from_smiles("NC")?;
    let dummy = mol_from_smarts("[#0]")?;
    let products = replace_substructs(&core, &dummy, &chain)?;
    println!("{}", mol_to_smiles(&products[0])?);
    assert_eq!(mol_to_smiles(&products[0])?, "CNc1c(C)cccc1O");

    // Here is an example in a loop for an imidazolium core with alkyl chains
    let core = mol_from_smiles("*[n+]1cc[nH]c1")?;
    let chains = ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"];
    let mut product_smi = Vec::new();
    for chain in chains {
        let chain_mol = mol_from_smiles(chain)?;
        let product_mol = replace_substructs(&core, &dummy, &chain_mol)?;
        product_smi.push(mol_to_smiles(&product_mol[0])?);
    }
    println!("{product_smi:?}");
    assert_eq!(
        product_smi,
        [
            "C[n+]1cc[nH]c1",
            "CC[n+]1cc[nH]c1",
            "CCC[n+]1cc[nH]c1",
            "CCCC[n+]1cc[nH]c1",
            "CCCCC[n+]1cc[nH]c1",
            "CCCCCC[n+]1cc[nH]c1",
        ]
    );
    Ok(())
}
