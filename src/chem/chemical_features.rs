// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.ChemicalFeatures`, a Python module that
//! re-exports `rdChemicalFeatures` and `rdMolChemicalFeatures`, as this one
//! does, and adds `MolChemicalFeatureFactory.GetFeaturesForMol`.

pub use crate::chem::rd_chemical_features::*;
pub use crate::chem::rd_mol_chemical_features::*;

use crate::Result;
use crate::chem::Mol;

/// Every feature the factory finds in the molecule. Python installs it as
/// `MolChemicalFeatureFactory.GetFeaturesForMol`. Mirrors
/// `ChemicalFeatures.MCFF_GetFeaturesForMol(self, mol)`.
pub fn mcff_get_features_for_mol<'a>(
    factory: &'a MolChemicalFeatureFactory,
    mol: &'a Mol,
) -> Result<Vec<MolChemicalFeature<'a>>> {
    mcff_get_features_for_mol_with_params(factory, mol, &GetMolFeatureParams::default())
}

/// As [`mcff_get_features_for_mol`], only for the family
/// `params.include_only` when it is not empty, placed by the conformer
/// `params.conf_id`. Mirrors
/// `ChemicalFeatures.MCFF_GetFeaturesForMol(self, mol, includeOnly, confId)`.
pub fn mcff_get_features_for_mol_with_params<'a>(
    factory: &'a MolChemicalFeatureFactory,
    mol: &'a Mol,
    params: &GetMolFeatureParams,
) -> Result<Vec<MolChemicalFeature<'a>>> {
    let count = factory.get_num_mol_features_with_include_only(mol, &params.include_only)?;
    (0..count)
        .map(|i| factory.get_mol_feature_with_params(mol, i, params))
        .collect()
}
