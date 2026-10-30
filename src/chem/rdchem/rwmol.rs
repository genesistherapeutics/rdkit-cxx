// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::ops::{Deref, DerefMut};

use crate::Result;
use crate::chem::rdchem::{AsBond, Atom, BondType, Mol};
use crate::error::expect_ok;
use crate::ffi::mol::ffi;

/// An editable molecule. Mirrors `rdchem.RWMol`.
///
/// It dereferences to `Mol`, so every `Mol` method works on it, and converts
/// to and from `Mol` without copying: both own an `RDKit::RWMol`
/// (CONVENTIONS.md W1).
#[derive(Clone, Debug, Default)]
pub struct RWMol {
    mol: Mol,
}

/// Parameters for `RWMol::replace_atom_with_params`, with Python's defaults.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReplaceAtomParams {
    /// Makes the new atom the active atom.
    pub update_label: bool,
    /// Keeps the replaced atom's properties unless the new atom sets them.
    pub preserve_props: bool,
}

/// Parameters for `RWMol::replace_bond_with_params`, with Python's defaults.
#[derive(Clone, Debug, PartialEq)]
pub struct ReplaceBondParams {
    /// Keeps the replaced bond's properties unless the new bond sets them.
    pub preserve_props: bool,
    /// Keeps substance groups that refer to the bond.
    pub keep_s_groups: bool,
}

impl Default for ReplaceBondParams {
    fn default() -> ReplaceBondParams {
        ReplaceBondParams {
            preserve_props: false,
            keep_s_groups: true,
        }
    }
}

impl RWMol {
    /// A copy as a `Mol`. Mirrors `RWMol.GetMol`.
    pub fn get_mol(&self) -> Mol {
        self.mol.clone()
    }

    /// Adds a copy of `atom`, returning its index. Mirrors `RWMol.AddAtom`.
    pub fn add_atom(&mut self, atom: &Atom) -> u32 {
        expect_ok(ffi::rwmol_add_atom(self.mol.pin_mut(), atom.raw()))
    }

    /// Adds a bond of unspecified type, returning the new number of bonds.
    /// Mirrors `RWMol.AddBond`.
    pub fn add_bond(&mut self, begin_atom_idx: u32, end_atom_idx: u32) -> Result<u32> {
        self.add_bond_with_order(begin_atom_idx, end_atom_idx, BondType::Unspecified)
    }

    /// Adds a bond, returning the new number of bonds. Mirrors
    /// `RWMol.AddBond(beginAtomIdx, endAtomIdx, order)`.
    pub fn add_bond_with_order(
        &mut self,
        begin_atom_idx: u32,
        end_atom_idx: u32,
        order: BondType,
    ) -> Result<u32> {
        Ok(ffi::rwmol_add_bond(
            self.mol.pin_mut(),
            begin_atom_idx,
            end_atom_idx,
            order.to_ffi(),
        )?)
    }

    /// Mirrors `RWMol.RemoveAtom`. `Error::Value`, with nothing removed, for
    /// a molecule whose `_MolFileBondEndPts` bond properties RDKit would
    /// misread (one that is empty, or in parentheses holds anything but
    /// numbers below 2^32 between spaces) or whose substance groups' `PARENT`
    /// chains form a cycle, which RDKit would follow forever; so do every
    /// removal and the functions that remove atoms or bonds.
    ///
    /// The atom and its bonds leave the stereo groups that held them (every
    /// copy), and a group left with neither atoms nor bonds goes, as in RDKit
    /// 2026.03; RDKit 2024.03 leaves a group pointing at a removed bond, or at
    /// a repeated atom's second copy. The functions whose RDKit call removes
    /// atoms or bonds itself cannot do that, so they return `Error::Value`
    /// for a molecule whose stereo groups hold a bond, repeat an atom, or
    /// share an atom between two ABS groups (which RDKit merges when it
    /// combines molecules).
    pub fn remove_atom(&mut self, idx: u32) -> Result<()> {
        Ok(ffi::rwmol_remove_atom(self.mol.pin_mut(), idx)?)
    }

    /// Removes the bond between two atoms, if there is one. Mirrors
    /// `RWMol.RemoveBond`, with the errors of `remove_atom`. The bond leaves
    /// the stereo groups that held it, as for `remove_atom`.
    pub fn remove_bond(&mut self, idx1: u32, idx2: u32) -> Result<()> {
        Ok(ffi::rwmol_remove_bond(self.mol.pin_mut(), idx1, idx2)?)
    }

    /// Replaces an atom with a copy of `new_atom`. Mirrors `RWMol.ReplaceAtom`.
    pub fn replace_atom(&mut self, index: u32, new_atom: &Atom) -> Result<()> {
        self.replace_atom_with_params(index, new_atom, &ReplaceAtomParams::default())
    }

    /// Mirrors `RWMol.ReplaceAtom(index, newAtom, updateLabel, preserveProps)`.
    pub fn replace_atom_with_params(
        &mut self,
        index: u32,
        new_atom: &Atom,
        params: &ReplaceAtomParams,
    ) -> Result<()> {
        Ok(ffi::rwmol_replace_atom(
            self.mol.pin_mut(),
            index,
            new_atom.raw(),
            params.update_label,
            params.preserve_props,
        )?)
    }
}

impl RWMol {
    /// Replaces a bond with a copy of `new_bond`, which may come from another
    /// molecule; the bond keeps its atoms. Mirrors `RWMol.ReplaceBond`.
    pub fn replace_bond(&mut self, index: u32, new_bond: &dyn AsBond) -> Result<()> {
        self.replace_bond_with_params(index, new_bond, &ReplaceBondParams::default())
    }

    /// Mirrors `RWMol.ReplaceBond(index, newBond, preserveProps, keepSGroups)`.
    /// Without `keep_s_groups`, it removes the substance groups that hold the
    /// bond as `remove_bond` does, with its errors.
    pub fn replace_bond_with_params(
        &mut self,
        index: u32,
        new_bond: &dyn AsBond,
        params: &ReplaceBondParams,
    ) -> Result<()> {
        Ok(ffi::rwmol_replace_bond(
            self.mol.pin_mut(),
            index,
            new_bond.bond(),
            params.preserve_props,
            params.keep_s_groups,
        )?)
    }

    /// Adds a copy of `mol`'s atoms and bonds. Mirrors `RWMol.InsertMol`.
    pub fn insert_mol(&mut self, mol: &Mol) {
        expect_ok(ffi::rwmol_insert_mol(self.mol.pin_mut(), mol.as_ffi()))
    }

    /// Defers atom and bond removals until `commit_batch_edit`, so indices
    /// stay valid meanwhile. Fails if a batch is already open. Mirrors
    /// `RWMol.BeginBatchEdit`; Python's `with` block runs the batch.
    pub fn begin_batch_edit(&mut self) -> Result<()> {
        Ok(ffi::rwmol_begin_batch_edit(self.mol.pin_mut())?)
    }

    /// Applies the batch's removals. Mirrors `RWMol.CommitBatchEdit`.
    /// `Error::Value`, with the batch still open and nothing removed, for a
    /// molecule whose `_MolFileBondEndPts` bond properties or substance
    /// groups' `PARENT` chains RDKit's removals would misread (a malformed
    /// value, or a cycle), as for `remove_atom` and `remove_bond`. The removed
    /// atoms and bonds leave the stereo groups, as for `remove_atom`.
    pub fn commit_batch_edit(&mut self) -> Result<()> {
        Ok(ffi::rwmol_commit_batch_edit(self.mol.pin_mut())?)
    }

    /// Discards the batch's removals. Mirrors `RWMol.RollbackBatchEdit`.
    pub fn rollback_batch_edit(&mut self) {
        expect_ok(ffi::rwmol_rollback_batch_edit(self.mol.pin_mut()))
    }
}

impl Deref for RWMol {
    type Target = Mol;

    fn deref(&self) -> &Mol {
        &self.mol
    }
}

impl DerefMut for RWMol {
    fn deref_mut(&mut self) -> &mut Mol {
        &mut self.mol
    }
}

/// Mirrors `RWMol(mol)`, without copying.
impl From<Mol> for RWMol {
    fn from(mol: Mol) -> RWMol {
        RWMol { mol }
    }
}

impl From<RWMol> for Mol {
    fn from(rw: RWMol) -> Mol {
        rw.mol
    }
}
