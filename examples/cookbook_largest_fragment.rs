// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Largest Fragment (RDKitCB_31): select the largest fragment of a molecule.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#largest-fragment),
//! by Andrew Dalke and Susan Leung.

use rdkit_cxx::Result;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::LargestFragmentChooser;
use rdkit_cxx::chem::{get_mol_frags_as_mols, mol_from_smiles, mol_to_smiles};

const SMILES: &str = "CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl.CO.C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N";
const LARGEST: &str = "Nc1nc2ncc(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)nc2c(=O)[nH]1";

fn main() -> Result<()> {
    let mol = mol_from_smiles(SMILES)?;
    let mol_frags = get_mol_frags_as_mols(&mol)?;
    // Python's max returns the first of equal fragments; Rust's max_by_key
    // returns the last, so search in reverse to match.
    let largest_mol = mol_frags
        .iter()
        .rev()
        .max_by_key(|m| m.get_num_atoms())
        .unwrap_or(&mol);
    println!("{}", mol_to_smiles(largest_mol)?);
    assert_eq!(mol_to_smiles(largest_mol)?, LARGEST);

    // The same result can also be achieved with MolStandardize:
    let largest_fragment = LargestFragmentChooser::default();
    let largest_mol = largest_fragment.choose(&mol)?;
    println!("{}", mol_to_smiles(&largest_mol)?);
    assert_eq!(mol_to_smiles(&largest_mol)?, LARGEST);
    Ok(())
}
