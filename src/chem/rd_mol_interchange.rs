// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolInterchange`, which `rdkit.Chem` re-exports:
//! molecules to and from JSON, in CommonChem's format or RDKit's extension
//! of it.
//!
//! RDKit's JSON parser trusts the indices, counts, and query trees it reads,
//! so the wrapper reads the JSON first and refuses what RDKit would trust
//! (CONVENTIONS.md U2): fewer than two stereo atoms for a bond, an empty atom
//! ring, a substance group bracket of more than three points, cstates and
//! attachment points for atoms or bonds the molecule lacks, a query on a
//! bond a stereo group holds (RDKit leaves the group pointing at the bond it
//! deletes), and queries RDKit does not make: a description with a class of
//! query RDKit never gives it, an AND, OR, or XOR query with fewer than two
//! children, a recursive query without atoms, and query values beyond 2^30
//! (U7).

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdmolops::{mol_list_into_vec, with_borrowed_mols_shared};
use crate::error::{expect_ok, to_string};
use crate::ffi::mol_interchange::ffi;

/// Parameters for the JSON parser. Mirrors
/// `rdMolInterchange.JSONParseParameters`; `Default` is a default-constructed
/// RDKit `JSONParseParameters` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct JSONParseParameters {
    /// Gives bonds the JSON marks aromatic the aromatic bond type.
    pub set_aromatic_bonds: bool,
    /// Checks valences strictly (`Error::AtomValence`).
    pub strict_valence_check: bool,
    /// Reads the molecules' properties.
    pub parse_properties: bool,
    /// Reads the molecules' conformers.
    pub parse_conformers: bool,
    /// Uses the atoms' hydrogen counts from the JSON; turn it off when
    /// reading queries.
    pub use_h_counts: bool,
}

impl Default for JSONParseParameters {
    fn default() -> JSONParseParameters {
        let data = expect_ok(ffi::json_parse_parameters_default());
        JSONParseParameters {
            set_aromatic_bonds: data.set_aromatic_bonds,
            strict_valence_check: data.strict_valence_check,
            parse_properties: data.parse_properties,
            parse_conformers: data.parse_conformers,
            use_h_counts: data.use_h_counts,
        }
    }
}

impl JSONParseParameters {
    fn to_ffi(&self) -> ffi::JSONParseParametersData {
        ffi::JSONParseParametersData {
            set_aromatic_bonds: self.set_aromatic_bonds,
            strict_valence_check: self.strict_valence_check,
            parse_properties: self.parse_properties,
            parse_conformers: self.parse_conformers,
            use_h_counts: self.use_h_counts,
        }
    }
}

/// Parameters for the JSON writer. Mirrors
/// `rdMolInterchange.JSONWriteParameters`; `Default` is a default-constructed
/// RDKit `JSONWriteParameters` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct JSONWriteParameters {
    /// Writes RDKit's extensions to CommonChem (`rdkitjson`): aromaticity,
    /// rings, CIP labels, partial charges, queries, and stereo and substance
    /// groups.
    pub use_rdkit_extensions: bool,
}

impl Default for JSONWriteParameters {
    fn default() -> JSONWriteParameters {
        let data = expect_ok(ffi::json_write_parameters_default());
        JSONWriteParameters {
            use_rdkit_extensions: data.use_rdkit_extensions,
        }
    }
}

impl JSONWriteParameters {
    fn to_ffi(&self) -> ffi::JSONWriteParametersData {
        ffi::JSONWriteParametersData {
            use_rdkit_extensions: self.use_rdkit_extensions,
        }
    }
}

/// The molecules in a JSON document. Mirrors `rdMolInterchange.JSONToMols`.
///
/// # Examples
///
/// ```no_run
/// use rdkit_cxx::chem::{json_to_mols, mol_from_smiles, mol_to_json, mol_to_smiles};
///
/// let json = mol_to_json(&mol_from_smiles("c1ccccc1O")?)?;
/// let mols = json_to_mols(&json)?;
/// assert_eq!(mol_to_smiles(&mols[0])?, "Oc1ccccc1");
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub fn json_to_mols(json_block: &str) -> Result<Vec<Mol>> {
    json_to_mols_with_params(json_block, &JSONParseParameters::default())
}

/// The molecules in a JSON document, read as `params` says. Text that is not
/// JSON, or JSON RDKit cannot read, is `Error::Parse`, and RDKit's
/// invariant violations (an aromatic atom or ring bond the molecule lacks)
/// `Error::Runtime`, as Python raises them; what RDKit would read without a
/// check (see the module documentation) is `Error::Value`. Mirrors
/// `rdMolInterchange.JSONToMols(jsonBlock, params)`.
pub fn json_to_mols_with_params(
    json_block: &str,
    params: &JSONParseParameters,
) -> Result<Vec<Mol>> {
    Ok(mol_list_into_vec(ffi::json_to_mols(
        json_block,
        &params.to_ffi(),
    )?))
}

/// A molecule as JSON. `Error::Runtime` for a molecule without atoms, as in
/// Python. Mirrors `rdMolInterchange.MolToJSON`.
pub fn mol_to_json(mol: &Mol) -> Result<String> {
    mol_to_json_with_params(mol, &JSONWriteParameters::default())
}

/// Mirrors `rdMolInterchange.MolToJSON(mol, params)`.
pub fn mol_to_json_with_params(mol: &Mol, params: &JSONWriteParameters) -> Result<String> {
    to_string(ffi::mol_to_json(mol.as_ffi(), &params.to_ffi())?)
}

/// Molecules as one JSON document. `Error::Runtime` if one has no atoms.
/// Mirrors `rdMolInterchange.MolsToJSON`.
pub fn mols_to_json(mols: &[Mol]) -> Result<String> {
    mols_to_json_with_params(mols, &JSONWriteParameters::default())
}

/// Mirrors `rdMolInterchange.MolsToJSON(mols, params)`.
pub fn mols_to_json_with_params(mols: &[Mol], params: &JSONWriteParameters) -> Result<String> {
    let json = with_borrowed_mols_shared(mols.iter().map(|mol| ("", mol)), |list| {
        ffi::mols_to_json(list, &params.to_ffi())
    })?;
    to_string(json)
}
