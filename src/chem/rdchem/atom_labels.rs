// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! `rdchem`'s module-level atom functions, which take any atom: owned or in a
//! molecule.

use std::pin::Pin;

use crate::Result;
use crate::chem::rdchem::{Atom, AtomMut, AtomRef};
use crate::error::{expect_ok, to_string};
use crate::ffi::atom::ffi;

mod sealed {
    use std::pin::Pin;

    use crate::ffi::atom::ffi;

    pub trait AtomLike {
        fn atom(&self) -> &ffi::Atom;
    }

    pub trait AtomLikeMut: AtomLike {
        fn atom_mut(&mut self) -> Pin<&mut ffi::Atom>;
    }
}

/// An atom the module-level functions read: `Atom`, `AtomRef`, or `AtomMut`.
pub trait AsAtom: sealed::AtomLike {}

/// An atom the module-level functions write: `Atom` or `AtomMut`.
pub trait AsAtomMut: sealed::AtomLikeMut {}

impl sealed::AtomLike for Atom {
    fn atom(&self) -> &ffi::Atom {
        self.raw()
    }
}

impl sealed::AtomLikeMut for Atom {
    fn atom_mut(&mut self) -> Pin<&mut ffi::Atom> {
        self.raw_mut()
    }
}

impl sealed::AtomLike for AtomRef<'_> {
    fn atom(&self) -> &ffi::Atom {
        self.raw()
    }
}

impl sealed::AtomLike for AtomMut<'_> {
    fn atom(&self) -> &ffi::Atom {
        self.raw()
    }
}

impl sealed::AtomLikeMut for AtomMut<'_> {
    fn atom_mut(&mut self) -> Pin<&mut ffi::Atom> {
        self.raw_mut()
    }
}

impl AsAtom for Atom {}
impl AsAtom for AtomRef<'_> {}
impl AsAtom for AtomMut<'_> {}
impl AsAtomMut for Atom {}
impl AsAtomMut for AtomMut<'_> {}

/// The atom's MDL alias text, or empty. Mirrors `rdchem.GetAtomAlias`.
pub fn get_atom_alias(atom: &impl AsAtom) -> Result<String> {
    to_string(ffi::atom_get_alias(atom.atom())?)
}

/// Sets the MDL alias text; empty clears it. Mirrors `rdchem.SetAtomAlias`.
pub fn set_atom_alias(atom: &mut impl AsAtomMut, rlabel: &str) {
    expect_ok(ffi::atom_set_alias(atom.atom_mut(), rlabel))
}

/// The atom's MDL R label (0 to 99), or 0. Mirrors `rdchem.GetAtomRLabel`.
pub fn get_atom_r_label(atom: &impl AsAtom) -> i32 {
    expect_ok(ffi::atom_get_r_label(atom.atom()))
}

/// Sets the MDL R label; 0 clears it. Mirrors `rdchem.SetAtomRLabel`.
pub fn set_atom_r_label(atom: &mut impl AsAtomMut, rlabel: i32) -> Result<()> {
    Ok(ffi::atom_set_r_label(atom.atom_mut(), rlabel)?)
}

/// The atom's MDL value text, or empty. Mirrors `rdchem.GetAtomValue`.
pub fn get_atom_value(atom: &impl AsAtom) -> Result<String> {
    to_string(ffi::atom_get_value(atom.atom())?)
}

/// Sets the MDL value text; empty clears it. Mirrors `rdchem.SetAtomValue`.
pub fn set_atom_value(atom: &mut impl AsAtomMut, rlabel: &str) {
    expect_ok(ffi::atom_set_value(atom.atom_mut(), rlabel))
}

/// Text written after the atom in SMILES, or empty. Mirrors
/// `rdchem.GetSupplementalSmilesLabel`.
pub fn get_supplemental_smiles_label(atom: &impl AsAtom) -> Result<String> {
    to_string(ffi::atom_get_supplemental_smiles_label(atom.atom())?)
}

/// Mirrors `rdchem.SetSupplementalSmilesLabel`.
pub fn set_supplemental_smiles_label(atom: &mut impl AsAtomMut, label: &str) {
    expect_ok(ffi::atom_set_supplemental_smiles_label(
        atom.atom_mut(),
        label,
    ))
}

/// The electrons the atom uses for pi bonding; needs a molecule. Mirrors
/// `rdchem.GetNumPiElectrons`.
pub fn get_num_pi_electrons(atom: &impl AsAtom) -> Result<u32> {
    Ok(ffi::atom_get_num_pi_electrons(atom.atom())?)
}
