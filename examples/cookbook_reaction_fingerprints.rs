// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reaction Fingerprints and Similarity (RDKitCB_27): build difference
//! fingerprints of reactions and compare them.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#reaction-fingerprints-and-similarity),
//! by Greg Landrum, with the fingerprint of Schneider et al., J. Chem. Inf.
//! Model. 2015, 55, 39.

use rdkit_cxx::Result;
use rdkit_cxx::chem::rd_chem_reactions::{
    create_difference_fingerprint_for_reaction, reaction_from_smarts,
};
use rdkit_cxx::data_structs::c_data_structs::tanimoto_similarity_sparse_int_vect;

fn main() -> Result<()> {
    let rxn1 = reaction_from_smarts("CCCO>>CCC=O")?;
    let rxn2 = reaction_from_smarts("CC(O)C>>CC(=O)C")?;
    let rxn3 = reaction_from_smarts("NCCO>>NCC=O")?;

    // A difference fingerprint subtracts the reactants' fingerprint from the
    // products'.
    let fp1 = create_difference_fingerprint_for_reaction(&rxn1)?;
    let fp2 = create_difference_fingerprint_for_reaction(&rxn2)?;
    let fp3 = create_difference_fingerprint_for_reaction(&rxn3)?;

    // The parts that change in rxn1 and rxn2 have nothing in common, as far
    // as the fingerprint is concerned:
    let sim12 = tanimoto_similarity_sparse_int_vect(&fp1, &fp2)?;
    println!("{sim12}");
    assert_eq!(sim12, 0.0);

    // whereas rxn1 and rxn3 share some:
    let sim13 = tanimoto_similarity_sparse_int_vect(&fp1, &fp3)?;
    println!("{sim13}");
    assert_eq!(sim13, 0.42857142857142855);
    Ok(())
}
