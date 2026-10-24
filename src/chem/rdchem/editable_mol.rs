// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Result;
use crate::chem::rdchem::{
    Atom, BondRef, BondType, Mol, RWMol, ReplaceAtomParams, ReplaceBondParams,
};

/// A copy of a molecule to add and remove atoms and bonds in. Mirrors
/// `rdchem.EditableMol`; `RWMol` does all this and more.
#[derive(Clone, Debug)]
pub struct EditableMol {
    mol: RWMol,
}

impl EditableMol {
    /// An editable copy of `m`. Mirrors `EditableMol(m)`.
    pub fn new(m: &Mol) -> EditableMol {
        EditableMol {
            mol: RWMol::from(m.clone()),
        }
    }

    /// Mirrors `EditableMol.RemoveAtom`.
    pub fn remove_atom(&mut self, idx: u32) -> Result<()> {
        self.mol.remove_atom(idx)
    }

    /// Removes the bond between two atoms, if there is one. Mirrors
    /// `EditableMol.RemoveBond`.
    pub fn remove_bond(&mut self, idx1: u32, idx2: u32) -> Result<()> {
        self.mol.remove_bond(idx1, idx2)
    }

    /// Adds a bond of unspecified type, returning the new number of bonds.
    /// Mirrors `EditableMol.AddBond`.
    pub fn add_bond(&mut self, begin_atom_idx: u32, end_atom_idx: u32) -> Result<u32> {
        self.mol.add_bond(begin_atom_idx, end_atom_idx)
    }

    /// Mirrors `EditableMol.AddBond(beginAtomIdx, endAtomIdx, order)`.
    pub fn add_bond_with_order(
        &mut self,
        begin_atom_idx: u32,
        end_atom_idx: u32,
        order: BondType,
    ) -> Result<u32> {
        self.mol
            .add_bond_with_order(begin_atom_idx, end_atom_idx, order)
    }

    /// Adds a copy of `atom`, returning its index. Mirrors
    /// `EditableMol.AddAtom`.
    pub fn add_atom(&mut self, atom: &Atom) -> u32 {
        self.mol.add_atom(atom)
    }

    /// Replaces an atom with a copy of `new_atom`. Mirrors
    /// `EditableMol.ReplaceAtom`.
    pub fn replace_atom(&mut self, index: u32, new_atom: &Atom) -> Result<()> {
        self.mol.replace_atom(index, new_atom)
    }

    /// Mirrors `EditableMol.ReplaceAtom(index, newAtom, updateLabel,
    /// preserveProps)`.
    pub fn replace_atom_with_params(
        &mut self,
        index: u32,
        new_atom: &Atom,
        params: &ReplaceAtomParams,
    ) -> Result<()> {
        self.mol.replace_atom_with_params(index, new_atom, params)
    }

    /// Replaces a bond with a copy of `new_bond`; the bond keeps its atoms.
    /// Mirrors `EditableMol.ReplaceBond`.
    pub fn replace_bond(&mut self, index: u32, new_bond: &BondRef<'_>) -> Result<()> {
        self.replace_bond_with_preserve_props(index, new_bond, false)
    }

    /// Keeps the replaced bond's properties unless the new bond sets them.
    /// Mirrors `EditableMol.ReplaceBond(index, newBond, preserveProps)`.
    pub fn replace_bond_with_preserve_props(
        &mut self,
        index: u32,
        new_bond: &BondRef<'_>,
        preserve_props: bool,
    ) -> Result<()> {
        let params = ReplaceBondParams {
            preserve_props,
            ..Default::default()
        };
        self.mol.replace_bond_with_params(index, new_bond, &params)
    }

    /// Defers removals until `commit_batch_edit`. Mirrors
    /// `EditableMol.BeginBatchEdit`.
    pub fn begin_batch_edit(&mut self) -> Result<()> {
        self.mol.begin_batch_edit()
    }

    /// Mirrors `EditableMol.CommitBatchEdit`.
    pub fn commit_batch_edit(&mut self) {
        self.mol.commit_batch_edit()
    }

    /// Mirrors `EditableMol.RollbackBatchEdit`.
    pub fn rollback_batch_edit(&mut self) {
        self.mol.rollback_batch_edit()
    }

    /// A copy of the molecule. Mirrors `EditableMol.GetMol`.
    pub fn get_mol(&self) -> Mol {
        self.mol.get_mol()
    }
}
