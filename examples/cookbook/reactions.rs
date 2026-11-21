// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Reactions" chapter.

use rdkit_cxx::Result;
use rdkit_cxx::chem::rd_chem_reactions::{
    ChemicalReaction, create_difference_fingerprint_for_reaction, reaction_from_smarts,
};
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles, sanitize_mol};
use rdkit_cxx::data_structs::c_data_structs::tanimoto_similarity_sparse_int_vect;

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[
    ("RDKitCB_6", "Reversing Reactions", reversing_reactions),
    ("RDKitCB_27", "Reaction Fingerprints and Similarity", reaction_fingerprints),
];

/// Sanitizes each product of the first outcome and writes its SMILES.
fn first_outcome(outcomes: Vec<Vec<Mol>>) -> Result<Vec<String>> {
    let mut smiles = Vec::new();
    for mut product in outcomes.into_iter().next().unwrap_or_default() {
        sanitize_mol(&mut product)?;
        smiles.push(mol_to_smiles(&product)?);
    }
    Ok(smiles)
}

/// Runs a reaction backward by swapping its reactant and product templates.
/// By Greg Landrum, with a reaction from Hartenfeller et al. This works poorly
/// for reactions whose reactants carry many queries.
fn reversing_reactions() -> Result<()> {
    // The Pictet-Spengler reaction.
    let mut rxn = reaction_from_smarts(
        "[cH1:1]1:[c:2](-[CH2:7]-[CH2:8]-[NH2:9]):[c:3]:[c:4]:[c:5]:[c:6]:1.\
         [#6:11]-[CH1;R0:10]=[OD1]>>\
         [c:1]12:[c:2](-[CH2:7]-[CH2:8]-[NH1:9]-[C:10]-2(-[#6:11])):[c:3]:[c:4]:[c:5]:[c:6]:1",
    )?;
    let mut reverse = ChemicalReaction::default();
    for idx in 0..rxn.get_num_reactant_templates() {
        reverse.add_product_template(&rxn.get_reactant_template(idx)?);
    }
    for idx in 0..rxn.get_num_product_templates() {
        reverse.add_reactant_template(&rxn.get_product_template(idx)?);
    }
    reverse.initialize()?;

    let reactants = [mol_from_smiles("NCCc1ccccc1")?, mol_from_smiles("C1CC1C(=O)")?];
    let products = first_outcome(rxn.run_reactants(&reactants)?)?;
    println!("{products:?}");
    assert_eq!(products, ["c1ccc2c(c1)CCNC2C1CC1"]);

    let products: Vec<Mol> = products
        .iter()
        .map(|smi| mol_from_smiles(smi))
        .collect::<Result<_>>()?;
    let recovered = first_outcome(reverse.run_reactants(&products)?)?;
    println!("{recovered:?}");
    assert_eq!(recovered, ["NCCc1ccccc1", "O=CC1CC1"]);
    Ok(())
}

/// Compares reactions by their difference fingerprints, after Schneider et
/// al. By Greg Landrum.
fn reaction_fingerprints() -> Result<()> {
    let rxn1 = reaction_from_smarts("CCCO>>CCC=O")?;
    let rxn2 = reaction_from_smarts("CC(O)C>>CC(=O)C")?;
    let rxn3 = reaction_from_smarts("NCCO>>NCC=O")?;
    // Each fingerprint is the products' less the reactants'.
    let fp1 = create_difference_fingerprint_for_reaction(&rxn1)?;
    let fp2 = create_difference_fingerprint_for_reaction(&rxn2)?;
    let fp3 = create_difference_fingerprint_for_reaction(&rxn3)?;

    // What changes in the first two reactions has nothing in common.
    let similarity = tanimoto_similarity_sparse_int_vect(&fp1, &fp2)?;
    println!("{similarity}");
    assert_eq!(similarity, 0.0);
    // The first and third share some of it.
    let similarity = tanimoto_similarity_sparse_int_vect(&fp1, &fp3)?;
    println!("{similarity}");
    assert_eq!(similarity, 0.42857142857142855);
    Ok(())
}
