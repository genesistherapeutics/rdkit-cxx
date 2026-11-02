// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::atom::expect_ok_index;
use crate::chem::rdchem::props::{prop_getters, prop_setters};
use crate::chem::rdchem::{AtomRef, BondDir, BondStereo, BondType, Mol};
use crate::error::{expect_ok, to_string};
use crate::ffi::bond::ffi;
use crate::ffi::props::ffi::{RDProps, bond_as_props};

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

            /// Mirrors `Bond.GetBondDir`.
            pub fn get_bond_dir(&self) -> BondDir {
                BondDir::from_ffi(expect_ok(ffi::bond_get_bond_dir(self.raw())))
            }

            /// Mirrors `Bond.GetStereo`.
            pub fn get_stereo(&self) -> BondStereo {
                BondStereo::from_ffi(expect_ok(ffi::bond_get_stereo(self.raw())))
            }

            /// The atoms that set the bond's stereochemistry. Mirrors
            /// `Bond.GetStereoAtoms`.
            pub fn get_stereo_atoms(&self) -> Vec<i32> {
                expect_ok(ffi::bond_get_stereo_atoms(self.raw()))
            }

            /// Mirrors `Bond.HasOwningMol`.
            pub fn has_owning_mol(&self) -> bool {
                expect_ok(ffi::bond_has_owning_mol(self.raw()))
            }

            /// Mirrors `Bond.HasQuery`.
            pub fn has_query(&self) -> bool {
                expect_ok(ffi::bond_has_query(self.raw()))
            }

            /// The bond as SMARTS, or as SMILES when it has no query. Mirrors
            /// `Bond.GetSmarts`.
            pub fn get_smarts(&self) -> Result<String> {
                self.get_smarts_with_all_bonds_explicit(false)
            }

            /// Writes single and aromatic bond symbols too when
            /// `all_bonds_explicit`. Mirrors `Bond.GetSmarts(allBondsExplicit)`.
            pub fn get_smarts_with_all_bonds_explicit(
                &self,
                all_bonds_explicit: bool,
            ) -> Result<String> {
                to_string(ffi::bond_get_smarts(self.raw(), all_bonds_explicit)?)
            }
        }
    };
}

/// Stamps the getters that need the owning molecule onto a view type with
/// `fn raw(&self) -> &ffi::Bond` (W4).
macro_rules! bond_view_getters {
    ($type:ident) => {
        impl $type<'_> {
            /// Whether the bond is in a ring, finding the molecule's rings if
            /// needed. Mirrors `Bond.IsInRing`.
            pub fn is_in_ring(&self) -> bool {
                expect_ok(ffi::bond_is_in_ring(self.raw()))
            }

            /// Whether the bond is in a ring of `size` bonds, finding the
            /// molecule's rings if needed. Mirrors `Bond.IsInRingSize`.
            pub fn is_in_ring_size(&self, size: i32) -> bool {
                expect_ok(ffi::bond_is_in_ring_size(self.raw(), size))
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

            /// Mirrors `Bond.SetBondDir`.
            pub fn set_bond_dir(&mut self, what: BondDir) {
                expect_ok(ffi::bond_set_bond_dir(self.raw_mut(), what.to_ffi()))
            }

            /// Mirrors `Bond.SetStereo`.
            pub fn set_stereo(&mut self, what: BondStereo) {
                expect_ok(ffi::bond_set_stereo(self.raw_mut(), what.to_ffi()))
            }

            /// Sets the atoms that set the bond's stereochemistry: a neighbor
            /// of each end. Mirrors `Bond.SetStereoAtoms`.
            pub fn set_stereo_atoms(&mut self, bgn_idx: u32, end_idx: u32) -> Result<()> {
                Ok(ffi::bond_set_stereo_atoms(self.raw_mut(), bgn_idx, end_idx)?)
            }

            /// Inverts atropisomer stereochemistry, returning whether there
            /// was any. Mirrors `Bond.InvertChirality`.
            pub fn invert_chirality(&mut self) -> bool {
                expect_ok(ffi::bond_invert_chirality(self.raw_mut()))
            }
        }
    };
}

/// A bond outside any molecule: a query bond from `rdqueries`. Mirrors
/// `rdchem.Bond`, which Python only constructs as a `QueryBond`;
/// `RWMol::replace_bond` copies it into a molecule.
pub struct Bond {
    /// Never null.
    ptr: UniquePtr<ffi::Bond>,
}

impl Bond {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::Bond>) -> Bond {
        assert!(!ptr.is_null(), "RDKit returned a null bond");
        Bond { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::Bond {
        self.ptr.as_ref().unwrap()
    }

    pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::Bond> {
        self.ptr.pin_mut()
    }

    fn props(&self) -> &RDProps {
        expect_ok(bond_as_props(self.raw()))
    }
}

/// A copy that keeps the bond's query. Mirrors `copy.copy` of a `Bond`.
impl Clone for Bond {
    fn clone(&self) -> Bond {
        Bond::from_ffi(expect_ok(ffi::bond_clone(self.raw())))
    }
}

impl fmt::Debug for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Bond")
            .field("bond_type", &self.get_bond_type())
            .finish()
    }
}

bond_getters!(Bond);
bond_setters!(Bond);
prop_getters!(Bond);
prop_setters!(Bond);

mod sealed {
    use crate::ffi::bond::ffi;

    pub trait BondLike {
        fn bond(&self) -> &ffi::Bond;
    }
}

/// A bond that methods taking another bond read: `Bond`, `BondRef`, or
/// `BondMut` (CONVENTIONS.md S3). A Rust addition; sealed.
pub trait AsBond: sealed::BondLike {}

impl sealed::BondLike for Bond {
    fn bond(&self) -> &ffi::Bond {
        self.raw()
    }
}

impl sealed::BondLike for BondRef<'_> {
    fn bond(&self) -> &ffi::Bond {
        self.raw()
    }
}

impl sealed::BondLike for BondMut<'_> {
    fn bond(&self) -> &ffi::Bond {
        self.raw()
    }
}

impl AsBond for Bond {}
impl AsBond for BondRef<'_> {}
impl AsBond for BondMut<'_> {}

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

    pub(crate) fn raw(&self) -> &'m ffi::Bond {
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

    /// The bond's other atom, given one of its atoms. Mirrors
    /// `Bond.GetOtherAtom`.
    pub fn get_other_atom(&self, what: &AtomRef<'_>) -> Result<AtomRef<'m>> {
        let idx = self.get_other_atom_idx(what.get_idx())?;
        Ok(expect_ok_index(self.mol.get_atom_with_idx(idx)))
    }

    /// The bond's contribution to an atom's valence; 0 for an atom not in the
    /// bond. Mirrors `Bond.GetValenceContrib`.
    pub fn get_valence_contrib(&self, at: &AtomRef<'_>) -> f64 {
        expect_ok(ffi::bond_get_valence_contrib(self.bond, at.raw()))
    }

    /// Mirrors `Bond.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &'m Mol {
        self.mol
    }
}

impl<'m> BondRef<'m> {
    fn props(&self) -> &'m RDProps {
        expect_ok(bond_as_props(self.bond))
    }
}

bond_getters!(BondRef<'m>);
bond_view_getters!(BondRef);
prop_getters!(BondRef<'m>);

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

    pub(crate) fn raw(&self) -> &ffi::Bond {
        expect_ok_index(self.mol.get_bond_with_idx(self.idx)).raw()
    }

    pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::Bond> {
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

impl BondMut<'_> {
    fn props(&self) -> &RDProps {
        expect_ok(bond_as_props(self.raw()))
    }
}

bond_getters!(BondMut<'m>);
bond_view_getters!(BondMut);
bond_setters!(BondMut<'m>);
prop_getters!(BondMut<'m>);
prop_setters!(BondMut<'m>);

impl fmt::Debug for BondMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BondMut").field("idx", &self.idx).finish()
    }
}
