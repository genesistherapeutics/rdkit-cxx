use std::ops::{Deref, DerefMut};

use crate::Result;
use crate::chem::rdchem::{Atom, BondType, Mol};
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

    /// Mirrors `RWMol.RemoveAtom`.
    pub fn remove_atom(&mut self, idx: u32) -> Result<()> {
        Ok(ffi::rwmol_remove_atom(self.mol.pin_mut(), idx)?)
    }

    /// Removes the bond between two atoms, if there is one. Mirrors
    /// `RWMol.RemoveBond`.
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
