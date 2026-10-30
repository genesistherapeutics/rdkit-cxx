// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolHash`: NextMove's molecular hashes.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::to_string;
use crate::ffi::mol_hash::ffi;

rdkit_enum! {
    /// A hash function. Mirrors `rdMolHash.HashFunction`.
    HashFunction: ffi::HashFunction {
        AnonymousGraph = AnonymousGraph,
        ElementGraph = ElementGraph,
        CanonicalSmiles = CanonicalSmiles,
        MurckoScaffold = MurckoScaffold,
        ExtendedMurcko = ExtendedMurcko,
        MolFormula = MolFormula,
        AtomBondCounts = AtomBondCounts,
        DegreeVector = DegreeVector,
        Mesomer = Mesomer,
        HetAtomTautomer = HetAtomTautomer,
        HetAtomProtomer = HetAtomProtomer,
        RedoxPair = RedoxPair,
        Regioisomer = Regioisomer,
        NetCharge = NetCharge,
        SmallWorldIndexBr = SmallWorldIndexBR,
        SmallWorldIndexBrl = SmallWorldIndexBRL,
        ArthorSubstructureOrder = ArthorSubstructureOrder,
        HetAtomTautomerv2 = HetAtomTautomerv2,
        HetAtomProtomerv2 = HetAtomProtomerv2,
    }
}

/// Parameters for `mol_hash_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MolHashParams {
    /// Writes CXSMILES extensions.
    pub use_cx_smiles: bool,
    /// The CXSMILES fields to leave out, as `rdmolfiles.CXSmilesFields`
    /// flags.
    pub cx_flags_to_skip: u32,
}

/// The molecule's hash. Mirrors `rdMolHash.MolHash`, which hashes a copy.
/// `Error::Runtime` when the molecule's property cache is not up to date.
/// The formula, mesomer, redox-pair, net-charge, tautomer, and protomer
/// hashes sum formal charges in `int`, and the tautomer and protomer hashes
/// hydrogen counts too; a molecule whose sums could overflow it (only
/// possible with millions of atoms) is `Error::Value`.
pub fn mol_hash(mol: &Mol, func: HashFunction) -> Result<String> {
    mol_hash_with_params(mol, func, &MolHashParams::default())
}

/// Mirrors `rdMolHash.MolHash(mol, func, useCxSmiles, cxFlagsToSkip)`.
pub fn mol_hash_with_params(
    mol: &Mol,
    func: HashFunction,
    params: &MolHashParams,
) -> Result<String> {
    to_string(ffi::mol_hash(
        mol.as_ffi(),
        func.to_ffi(),
        params.use_cx_smiles,
        params.cx_flags_to_skip,
    )?)
}
