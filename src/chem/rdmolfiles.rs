//! Mirrors `rdkit.Chem.rdmolfiles`.

use crate::Result;
use crate::chem::Mol;
use crate::error::to_string;
use crate::ffi::molfiles::ffi;

/// Parses a SMILES string, sanitizing the result. Mirrors
/// `rdmolfiles.MolFromSmiles`, which returns `None` where this returns
/// `Error::Parse` or a sanitization error.
pub fn mol_from_smiles(smiles: &str) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_smiles(smiles)?))
}

/// The canonical SMILES for a molecule. Mirrors `rdmolfiles.MolToSmiles`.
pub fn mol_to_smiles(mol: &Mol) -> Result<String> {
    to_string(ffi::mol_to_smiles(mol.as_ffi())?)
}
