// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reversing Reactions (RDKitCB_6): decompose a reaction product into its
//! reactants.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#reversing-reactions),
//! by Greg Landrum, with the Pictet-Spengler reaction from Hartenfeller et
//! al., J. Chem. Inf. Model. 2011, 51, 3093. The recipe draws the molecules;
//! this prints their SMILES.
//!
//! The approach isn't perfect and won't work for every reaction: reactions
//! with extensive query information in their reactants are likely to be
//! problematic.

use rdkit_cxx::Result;
use rdkit_cxx::chem::rd_chem_reactions::{ChemicalReaction, reaction_from_smarts};
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles, sanitize_mol};

fn smiles(mols: &[Mol]) -> Result<Vec<String>> {
    mols.iter().map(mol_to_smiles).collect()
}

fn main() -> Result<()> {
    // Pictet-Spengler rxn
    let mut rxn = reaction_from_smarts(
        "[cH1:1]1:[c:2](-[CH2:7]-[CH2:8]-[NH2:9]):[c:3]:[c:4]:[c:5]:[c:6]:1.[#6:11]-[CH1;R0:10]=[OD1]>>[c:1]12:[c:2](-[CH2:7]-[CH2:8]-[NH1:9]-[C:10]-2(-[#6:11])):[c:3]:[c:4]:[c:5]:[c:6]:1",
    )?;

    // The reverse reaction swaps the reactant and product templates.
    let mut rxn2 = ChemicalReaction::default();
    for i in 0..rxn.get_num_reactant_templates() {
        rxn2.add_product_template(&*rxn.get_reactant_template(i)?);
    }
    for i in 0..rxn.get_num_product_templates() {
        rxn2.add_reactant_template(&*rxn.get_product_template(i)?);
    }
    rxn2.initialize()?;

    let reacts = [
        mol_from_smiles("NCCc1ccccc1")?,
        mol_from_smiles("C1CC1C(=O)")?,
    ];
    let mut ps = rxn.run_reactants(&reacts)?;
    let mut ps0 = ps.swap_remove(0);
    for p in &mut ps0 {
        sanitize_mol(p)?;
    }
    println!("{:?}", smiles(&ps0)?);
    assert_eq!(smiles(&ps0)?, ["c1ccc2c(c1)CCNC2C1CC1"]);

    let reacts = ps0;
    let mut rps = rxn2.run_reactants(&reacts)?;
    let mut rps0 = rps.swap_remove(0);
    for rp in &mut rps0 {
        sanitize_mol(rp)?;
    }
    println!("{:?}", smiles(&rps0)?);
    assert_eq!(smiles(&rps0)?, ["NCCc1ccccc1", "O=CC1CC1"]);
    Ok(())
}
