// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.SpacialScore`, a pure-Python module: the spacial score
//! (SPS), a measure of a molecule's spatial complexity from its atoms'
//! hybridization, stereo, ring membership, and heavy neighbors
//! (Krzyzanowski, Pahl, Grigalunas, and Waldmann, J. Med. Chem. (2023),
//! <https://doi.org/10.1021/acs.jmedchem.3c00689>).
//!
//! Each atom scores h·s·r·n², for hybridization h (1 for sp, 2 for sp2, 3
//! for sp3, 4 otherwise), stereo s (2 for a stereocenter, specified or not,
//! or an atom of a double bond with stereo, else 1), ring r (2 for an atom
//! in a ring that is not aromatic, else 1), and degree n. SPS is the sum,
//! and nSPS, the default, divides it by the number of heavy atoms.

use crate::chem::{
    BondStereo, BondType, FindMolChiralCentersParams, HybridizationType, Mol,
    find_mol_chiral_centers_with_params, find_potential_stereo_bonds,
};
use crate::{Error, Result};

/// The normalized spacial score (nSPS). Mirrors `SpacialScore.SPS`.
pub fn sps(mol: &Mol) -> Result<f64> {
    sps_with_normalize(mol, true)
}

/// The spacial score, normalized by the number of heavy atoms when
/// `normalize`; `Error::Value` then for a molecule without heavy atoms,
/// where Python divides by zero. Python returns the unnormalized score as an
/// `int`. Mirrors `SpacialScore.SPS(mol, normalize)`.
pub fn sps_with_normalize(mol: &Mol, normalize: bool) -> Result<f64> {
    let mut mol = mol.clone();
    find_potential_stereo_bonds(&mut mol)?;
    let chiral_idxs: Vec<u32> = find_mol_chiral_centers_with_params(
        &mut mol,
        &FindMolChiralCentersParams {
            include_unassigned: true,
            include_cip: false,
            use_legacy_implementation: Some(false),
            ..FindMolChiralCentersParams::default()
        },
    )?
    .into_iter()
    .map(|(idx, _)| idx)
    .collect();
    let stereo_bond_atoms: Vec<u32> = mol
        .get_bonds()
        .filter(|bond| {
            bond.get_bond_type() == BondType::Double && bond.get_stereo() != BondStereo::Stereonone
        })
        .flat_map(|bond| [bond.get_begin_atom_idx(), bond.get_end_atom_idx()])
        .collect();
    let mut score: u64 = 0;
    for atom in mol.get_atoms() {
        let idx = atom.get_idx();
        let hybridization = match atom.get_hybridization() {
            HybridizationType::Sp => 1,
            HybridizationType::Sp2 => 2,
            HybridizationType::Sp3 => 3,
            _ => 4,
        };
        let stereo = if chiral_idxs.contains(&idx) || stereo_bond_atoms.contains(&idx) {
            2
        } else {
            1
        };
        let ring = if !atom.get_is_aromatic() && atom.is_in_ring() {
            2
        } else {
            1
        };
        let degree = u64::from(atom.get_degree());
        score += hybridization * stereo * ring * degree * degree;
    }
    if !normalize {
        return Ok(score as f64);
    }
    match mol.get_num_heavy_atoms() {
        0 => Err(Error::Value {
            message: "division by zero: the molecule has no heavy atoms".to_owned(),
        }),
        heavy => Ok(score as f64 / f64::from(heavy)),
    }
}
