// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Randomize`, a Python module: molecules with their
//! atoms in a random order.
//!
//! Python draws from the `random` module's global generator; the ports take
//! a [`Random`], which draws what Python's would from the same seed.

use crate::Result;
use crate::chem::{
    Mol, MolFromMolBlockParams, mol_from_mol_block_with_params, mol_to_mol_block, renumber_atoms,
};
use crate::rd_random::Random;

/// A copy of the molecule with its atoms in a random order, one that
/// differs from theirs when it has two atoms or more. The empty molecule
/// comes back as it is, where Python's binding refuses its empty atom
/// order. Mirrors `Randomize.RandomizeMol`.
pub fn randomize_mol(mol: &Mol, random: &mut Random) -> Result<Mol> {
    let atom_indices: Vec<u32> = (0..mol.get_num_atoms()).collect();
    let mut randomized = random.sample(&atom_indices, atom_indices.len());
    if atom_indices.len() > 1 {
        while randomized == atom_indices {
            randomized = random.sample(&atom_indices, atom_indices.len());
        }
    }
    renumber_atoms(mol, &randomized)
}

/// The mol block, read without sanitizing or removing hydrogens, with its
/// atoms in a random order. Mirrors `Randomize.RandomizeMolBlock`.
pub fn randomize_mol_block(molblock: &str, random: &mut Random) -> Result<String> {
    let mol = mol_from_mol_block_with_params(
        molblock,
        &MolFromMolBlockParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )?;
    mol_to_mol_block(&randomize_mol(&mol, random)?)
}
