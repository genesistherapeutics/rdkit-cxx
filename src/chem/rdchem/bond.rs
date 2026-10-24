use std::fmt;
use std::pin::Pin;

use crate::Result;
use crate::chem::rdchem::atom::expect_ok_index;
use crate::chem::rdchem::{AtomRef, BondType, Mol};
use crate::error::expect_ok;
use crate::ffi::bond::ffi;

/// Stamps the getters that need no owning molecule onto a type with
/// `fn raw(&self) -> &ffi::Bond` (CONVENTIONS.md W3).
macro_rules! bond_getters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// The bond's index in its molecule. Mirrors `Bond.GetIdx`.
            pub fn get_idx(&self) -> u32 {
                expect_ok(ffi::bond_get_idx(self.raw()))
            }

            /// Mirrors `Bond.GetBondType`.
            pub fn get_bond_type(&self) -> BondType {
                BondType::from_ffi(expect_ok(ffi::bond_get_bond_type(self.raw())))
            }

            /// The bond order (1.5 for aromatic). RDKit has no order for some
            /// types, such as `Other`. Mirrors `Bond.GetBondTypeAsDouble`.
            pub fn get_bond_type_as_double(&self) -> Result<f64> {
                Ok(ffi::bond_get_bond_type_as_double(self.raw())?)
            }

            /// Mirrors `Bond.GetBeginAtomIdx`.
            pub fn get_begin_atom_idx(&self) -> u32 {
                expect_ok(ffi::bond_get_begin_atom_idx(self.raw()))
            }

            /// Mirrors `Bond.GetEndAtomIdx`.
            pub fn get_end_atom_idx(&self) -> u32 {
                expect_ok(ffi::bond_get_end_atom_idx(self.raw()))
            }

            /// The index of the bond's other atom, given one of its atoms.
            /// Mirrors `Bond.GetOtherAtomIdx`.
            pub fn get_other_atom_idx(&self, this_idx: u32) -> Result<u32> {
                Ok(ffi::bond_get_other_atom_idx(self.raw(), this_idx)?)
            }

            /// Mirrors `Bond.GetIsAromatic`.
            pub fn get_is_aromatic(&self) -> bool {
                expect_ok(ffi::bond_get_is_aromatic(self.raw()))
            }

            /// Mirrors `Bond.GetIsConjugated`.
            pub fn get_is_conjugated(&self) -> bool {
                expect_ok(ffi::bond_get_is_conjugated(self.raw()))
            }
        }
    };
}

/// Stamps the setters onto a type with
/// `fn raw_mut(&mut self) -> Pin<&mut ffi::Bond>` (W3).
macro_rules! bond_setters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `Bond.SetBondType`.
            pub fn set_bond_type(&mut self, b_t: BondType) {
                expect_ok(ffi::bond_set_bond_type(self.raw_mut(), b_t.to_ffi()))
            }

            /// Mirrors `Bond.SetIsAromatic`.
            pub fn set_is_aromatic(&mut self, what: bool) {
                expect_ok(ffi::bond_set_is_aromatic(self.raw_mut(), what))
            }

            /// Mirrors `Bond.SetIsConjugated`.
            pub fn set_is_conjugated(&mut self, what: bool) {
                expect_ok(ffi::bond_set_is_conjugated(self.raw_mut(), what))
            }
        }
    };
}

/// A bond in a molecule, borrowed. Mirrors `rdchem.Bond`.
#[derive(Clone, Copy)]
pub struct BondRef<'m> {
    mol: &'m Mol,
    bond: &'m ffi::Bond,
}

impl<'m> BondRef<'m> {
    pub(crate) fn new(mol: &'m Mol, bond: &'m ffi::Bond) -> BondRef<'m> {
        BondRef { mol, bond }
    }

    fn raw(&self) -> &'m ffi::Bond {
        self.bond
    }

    /// Mirrors `Bond.GetBeginAtom`.
    pub fn get_begin_atom(&self) -> AtomRef<'m> {
        expect_ok_index(self.mol.get_atom_with_idx(self.get_begin_atom_idx()))
    }

    /// Mirrors `Bond.GetEndAtom`.
    pub fn get_end_atom(&self) -> AtomRef<'m> {
        expect_ok_index(self.mol.get_atom_with_idx(self.get_end_atom_idx()))
    }

    /// Mirrors `Bond.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &'m Mol {
        self.mol
    }
}

bond_getters!(BondRef<'m>);

impl fmt::Debug for BondRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BondRef")
            .field("idx", &self.get_idx())
            .finish()
    }
}

/// A bond in a molecule, borrowed exclusively for editing. Mirrors
/// `rdchem.Bond`.
///
/// It holds an index and finds the bond on each call (CONVENTIONS.md W2).
pub struct BondMut<'m> {
    mol: &'m mut Mol,
    idx: u32,
}

impl<'m> BondMut<'m> {
    /// `idx` has been checked against the molecule.
    pub(crate) fn new(mol: &'m mut Mol, idx: u32) -> BondMut<'m> {
        BondMut { mol, idx }
    }

    fn raw(&self) -> &ffi::Bond {
        expect_ok_index(self.mol.get_bond_with_idx(self.idx)).raw()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::Bond> {
        expect_ok_index(self.mol.bond_ffi_mut(self.idx))
    }

    /// Mirrors `Bond.GetBeginAtom`.
    pub fn get_begin_atom(&self) -> AtomRef<'_> {
        expect_ok_index(self.mol.get_atom_with_idx(self.get_begin_atom_idx()))
    }

    /// Mirrors `Bond.GetEndAtom`.
    pub fn get_end_atom(&self) -> AtomRef<'_> {
        expect_ok_index(self.mol.get_atom_with_idx(self.get_end_atom_idx()))
    }

    /// Mirrors `Bond.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &Mol {
        self.mol
    }
}

bond_getters!(BondMut<'m>);
bond_setters!(BondMut<'m>);

impl fmt::Debug for BondMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BondMut").field("idx", &self.idx).finish()
    }
}
