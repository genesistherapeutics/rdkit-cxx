// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolEnumerator`: the molecules a molecule's position
//! variation bonds, link nodes, and repeat units (SRU substance groups) stand
//! for, as mol files and CXSMILES write them.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smiles;
//! use rdkit_cxx::chem::rd_mol_enumerator::enumerate;
//!
//! let mol = mol_from_smiles("CCO |LN:1:1.3|")?;
//! let bundle = enumerate(&mol)?;
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```
//!
//! RDKit reads these features from properties and substance groups that any
//! caller can set, and trusts some of what it reads, so the functions refuse
//! with `Error::Value`:
//! - an operation with more than 4000 variation points, each of which RDKit
//!   enumerates one recursion deeper;
//! - link nodes RDKit's reactions might not find, after which it reads a
//!   product that is not there: on atoms past index 64 (RDKit labels them
//!   with 16-bit isotopes), across bonds that are not single, with both bonds
//!   to one atom, or two on one atom;
//! - SRU groups with atoms out of range, and what RDKit would misread when
//!   it removes atoms: a bond's `_MolFileBondEndPts`, or substance groups
//!   whose `PARENT`s form a cycle.

use crate::Result;
use crate::chem::MolBundle;
use crate::chem::rdchem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::rdmolops::mol_list_into_vec;
use crate::error::expect_ok;
use crate::ffi::mol_enumerator::ffi;

rdkit_enum! {
    /// An enumeration operation. Mirrors `rdMolEnumerator.EnumeratorType`.
    EnumeratorType: ffi::EnumeratorType {
        LinkNode = LinkNode,
        PositionVariation = PositionVariation,
        RepeatUnit = RepeatUnit,
    }
}

/// Parameters for `enumerate_with_params`, with RDKit's defaults (A3).
/// Mirrors `rdMolEnumerator.MolEnumeratorParams`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MolEnumeratorParams {
    /// Unused by RDKit's enumeration. Mirrors `MolEnumeratorParams.sanitize`.
    pub sanitize: bool,
    /// The most molecules to make. Mirrors
    /// `MolEnumeratorParams.maxToEnumerate`.
    pub max_to_enumerate: usize,
    /// Random enumeration, which RDKit has not implemented: set, an
    /// enumeration with variation points is `Error::Runtime`. Mirrors
    /// `MolEnumeratorParams.doRandom`.
    pub do_random: bool,
    /// Unused by RDKit's enumeration. Mirrors
    /// `MolEnumeratorParams.randomSeed`.
    pub random_seed: i32,
    /// The operation; `None`, as in `MolEnumeratorParams()`, makes
    /// `enumerate_with_params` fail with `Error::Value`. Python sets it
    /// with the constructor or `SetEnumerationOperator`, and cannot read it.
    pub operation: Option<EnumeratorType>,
}

/// Parameters without an operation. Mirrors `MolEnumeratorParams()`.
impl Default for MolEnumeratorParams {
    fn default() -> MolEnumeratorParams {
        let data = expect_ok(ffi::mol_enumerator_params_default());
        MolEnumeratorParams {
            sanitize: data.sanitize,
            max_to_enumerate: data.max_to_enumerate,
            do_random: data.do_random,
            random_seed: data.random_seed,
            operation: None,
        }
    }
}

impl MolEnumeratorParams {
    /// Parameters with the operation `typ`. Mirrors
    /// `MolEnumeratorParams(typ)`.
    pub fn new(typ: EnumeratorType) -> MolEnumeratorParams {
        MolEnumeratorParams {
            operation: Some(typ),
            ..Default::default()
        }
    }

    /// Mirrors `MolEnumeratorParams.SetEnumerationOperator`.
    pub fn set_enumeration_operator(&mut self, typ: EnumeratorType) {
        self.operation = Some(typ);
    }

    fn to_ffi(&self) -> ffi::MolEnumeratorParamsData {
        ffi::MolEnumeratorParamsData {
            sanitize: self.sanitize,
            max_to_enumerate: self.max_to_enumerate,
            do_random: self.do_random,
            random_seed: self.random_seed,
            has_operation: self.operation.is_some(),
            operation: self.operation.unwrap_or(EnumeratorType::LinkNode).to_ffi(),
        }
    }
}

/// The molecules the molecule stands for: its position variation bonds are
/// enumerated first, then its repeat units, then its link nodes (which RDKit
/// supports only in molecules of one fragment, and not with repeat units),
/// up to 1000 molecules for each operation. An empty bundle if the molecule
/// has none of these. Mirrors `rdMolEnumerator.Enumerate(mol)`.
pub fn enumerate(mol: &Mol) -> Result<MolBundle> {
    enumerate_with_max_per_operation(mol, 0)
}

/// As `enumerate`, with at most `max_per_operation` molecules from each
/// operation, and repeat units repeated at most that many times; 0 for
/// RDKit's defaults. Mirrors `rdMolEnumerator.Enumerate(mol,
/// maxPerOperation)`.
pub fn enumerate_with_max_per_operation(mol: &Mol, max_per_operation: u32) -> Result<MolBundle> {
    let list = ffi::mol_enumerator_enumerate(mol.as_ffi(), max_per_operation)?;
    Ok(MolBundle::from_mols(mol_list_into_vec(list)))
}

/// The molecules one operation makes; see `MolEnumeratorParams`. Mirrors
/// `rdMolEnumerator.Enumerate(mol, enumParams)`.
pub fn enumerate_with_params(mol: &Mol, enum_params: &MolEnumeratorParams) -> Result<MolBundle> {
    let list = ffi::mol_enumerator_enumerate_with_params(mol.as_ffi(), &enum_params.to_ffi())?;
    Ok(MolBundle::from_mols(mol_list_into_vec(list)))
}
