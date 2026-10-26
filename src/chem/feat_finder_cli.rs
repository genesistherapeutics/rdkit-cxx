// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FeatFinderCLI`, a Python command-line program
//! that lists the pharmacophore features of the molecules in a SMILES
//! file. Its reusable function is ported; the command line (`main`,
//! `initParser`, `processArgs`, `existingFile`, `splitExpr`) is not.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rd_mol_chemical_features::MolChemicalFeatureFactory;

/// The `family-type` names of the features each atom is in, in the order
/// of the factory's features; `None` for an atom in none. Mirrors
/// `FeatFinderCLI.GetAtomFeatInfo`.
pub fn get_atom_feat_info(
    factory: &MolChemicalFeatureFactory,
    mol: &Mol,
) -> Result<Vec<Option<Vec<String>>>> {
    let mut res: Vec<Option<Vec<String>>> = vec![None; mol.get_num_atoms() as usize];
    for i in 0..factory.get_num_mol_features(mol)? {
        let feat = factory.get_mol_feature(mol, i)?;
        let feature = format!("{}-{}", feat.get_family(), feat.get_type());
        for id in feat.get_atom_ids() {
            res[id as usize]
                .get_or_insert_with(Vec::new)
                .push(feature.clone());
        }
    }
    Ok(res)
}
