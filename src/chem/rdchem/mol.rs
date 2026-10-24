use std::fmt;
use std::ops::Range;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::atom::expect_ok_index;
use crate::chem::rdchem::props::{prop_getters, prop_setters};
use crate::chem::rdchem::{AtomMut, AtomRef, BondMut, BondRef};
use crate::error::expect_ok;
use crate::ffi::mol::ffi;
use crate::ffi::props::ffi::{RDProps, mol_as_props, mol_clear_computed_props};

/// A molecule. Mirrors `rdchem.Mol`.
///
/// `Mol` is `Send` but not `Sync`: RDKit writes computed properties through
/// read-only molecules, so parallel code gives each thread its own molecule
/// (CONVENTIONS.md S1).
///
/// ```compile_fail
/// fn assert_sync<T: Sync>() {}
/// assert_sync::<rdkit_cxx::chem::Mol>();
/// ```
pub struct Mol {
    /// Never null.
    ptr: UniquePtr<ffi::RWMol>,
}

impl Mol {
    /// Takes ownership of a molecule from the C++ wrappers.
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::RWMol>) -> Mol {
        assert!(!ptr.is_null(), "RDKit returned a null molecule");
        Mol { ptr }
    }

    pub(crate) fn as_ffi(&self) -> &ffi::RWMol {
        self.ptr.as_ref().unwrap()
    }

    pub(crate) fn pin_mut(&mut self) -> std::pin::Pin<&mut ffi::RWMol> {
        self.ptr.pin_mut()
    }

    /// The number of atoms in the graph. Mirrors `Mol.GetNumAtoms`.
    pub fn get_num_atoms(&self) -> u32 {
        expect_ok(ffi::mol_get_num_atoms(self.as_ffi()))
    }

    /// The number of atoms, counting implicit hydrogens unless `only_explicit`.
    /// Mirrors `Mol.GetNumAtoms(onlyExplicit)`. Counting implicit hydrogens
    /// needs the property cache.
    pub fn get_num_atoms_with_only_explicit(&self, only_explicit: bool) -> Result<u32> {
        Ok(ffi::mol_get_num_atoms_with_only_explicit(
            self.as_ffi(),
            only_explicit,
        )?)
    }

    /// The number of atoms with atomic number above 1. Mirrors
    /// `Mol.GetNumHeavyAtoms`.
    pub fn get_num_heavy_atoms(&self) -> u32 {
        expect_ok(ffi::mol_get_num_heavy_atoms(self.as_ffi()))
    }

    /// The number of bonds in the graph. Mirrors `Mol.GetNumBonds`.
    pub fn get_num_bonds(&self) -> u32 {
        expect_ok(ffi::mol_get_num_bonds(self.as_ffi(), true))
    }

    /// The number of bonds, counting bonds to implicit hydrogens unless
    /// `only_heavy`. Mirrors `Mol.GetNumBonds(onlyHeavy)`. Counting implicit
    /// hydrogens needs the property cache.
    pub fn get_num_bonds_with_only_heavy(&self, only_heavy: bool) -> Result<u32> {
        Ok(ffi::mol_get_num_bonds(self.as_ffi(), only_heavy)?)
    }

    /// Mirrors `Mol.GetAtomWithIdx`.
    pub fn get_atom_with_idx(&self, idx: u32) -> Result<AtomRef<'_>> {
        let atom = ffi::mol_get_atom_with_idx(self.as_ffi(), idx)?;
        Ok(AtomRef::new(self, atom))
    }

    /// Mirrors `Mol.GetBondWithIdx`.
    pub fn get_bond_with_idx(&self, idx: u32) -> Result<BondRef<'_>> {
        let bond = ffi::mol_get_bond_with_idx(self.as_ffi(), idx)?;
        Ok(BondRef::new(self, bond))
    }

    /// The bond between two atoms, if they are bonded. Mirrors
    /// `Mol.GetBondBetweenAtoms`.
    pub fn get_bond_between_atoms(&self, idx1: u32, idx2: u32) -> Result<Option<BondRef<'_>>> {
        let idx = ffi::mol_get_bond_between_atoms(self.as_ffi(), idx1, idx2)?;
        Ok(u32::try_from(idx)
            .ok()
            .map(|idx| expect_ok_index(self.get_bond_with_idx(idx))))
    }

    /// Mirrors `Mol.GetAtomWithIdx`, for editing.
    pub fn get_atom_with_idx_mut(&mut self, idx: u32) -> Result<AtomMut<'_>> {
        self.atom_ffi_mut(idx)?;
        Ok(AtomMut::new(self, idx))
    }

    /// Mirrors `Mol.GetBondWithIdx`, for editing.
    pub fn get_bond_with_idx_mut(&mut self, idx: u32) -> Result<BondMut<'_>> {
        self.bond_ffi_mut(idx)?;
        Ok(BondMut::new(self, idx))
    }

    pub(crate) fn atom_ffi_mut(
        &mut self,
        idx: u32,
    ) -> Result<std::pin::Pin<&mut crate::ffi::atom::ffi::Atom>> {
        Ok(ffi::mol_get_atom_with_idx_mut(self.pin_mut(), idx)?)
    }

    pub(crate) fn bond_ffi_mut(
        &mut self,
        idx: u32,
    ) -> Result<std::pin::Pin<&mut crate::ffi::bond::ffi::Bond>> {
        Ok(ffi::mol_get_bond_with_idx_mut(self.pin_mut(), idx)?)
    }

    /// Recomputes implicit valences and related properties, failing on a bad
    /// valence. Mirrors `Mol.UpdatePropertyCache`.
    pub fn update_property_cache(&mut self) -> Result<()> {
        self.update_property_cache_with_strict(true)
    }

    /// Recomputes implicit valences and related properties, failing on a bad
    /// valence only when `strict`. Mirrors `Mol.UpdatePropertyCache(strict)`.
    pub fn update_property_cache_with_strict(&mut self, strict: bool) -> Result<()> {
        Ok(ffi::mol_update_property_cache(self.pin_mut(), strict)?)
    }

    /// Mirrors `Mol.GetAtoms`.
    pub fn get_atoms(&self) -> Atoms<'_> {
        Atoms {
            mol: self,
            range: 0..self.get_num_atoms(),
        }
    }

    /// Mirrors `Mol.GetBonds`.
    pub fn get_bonds(&self) -> Bonds<'_> {
        Bonds {
            mol: self,
            range: 0..self.get_num_bonds(),
        }
    }
}

/// Defines an iterator over a molecule's items by index (CONVENTIONS.md W5).
macro_rules! index_iterator {
    ($(#[$meta:meta])* $name:ident, $item:ident, $get:ident) => {
        $(#[$meta])*
        #[derive(Clone, Debug)]
        pub struct $name<'m> {
            mol: &'m Mol,
            range: Range<u32>,
        }

        impl<'m> Iterator for $name<'m> {
            type Item = $item<'m>;

            fn next(&mut self) -> Option<$item<'m>> {
                let idx = self.range.next()?;
                Some(expect_ok_index(self.mol.$get(idx)))
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                self.range.size_hint()
            }
        }

        impl<'m> DoubleEndedIterator for $name<'m> {
            fn next_back(&mut self) -> Option<$item<'m>> {
                let idx = self.range.next_back()?;
                Some(expect_ok_index(self.mol.$get(idx)))
            }
        }

        impl ExactSizeIterator for $name<'_> {}
    };
}

index_iterator!(
    /// The atoms of a molecule. Mirrors the sequence `Mol.GetAtoms` returns.
    Atoms,
    AtomRef,
    get_atom_with_idx
);
index_iterator!(
    /// The bonds of a molecule. Mirrors the sequence `Mol.GetBonds` returns.
    Bonds,
    BondRef,
    get_bond_with_idx
);

impl Mol {
    fn props(&self) -> &RDProps {
        expect_ok(mol_as_props(self.as_ffi()))
    }

    /// Removes the properties RDKit computed, ring information included.
    /// Mirrors `Mol.ClearComputedProps`.
    pub fn clear_computed_props(&mut self) {
        self.clear_computed_props_with_include_rings(true)
    }

    /// Mirrors `Mol.ClearComputedProps(includeRings)`.
    pub fn clear_computed_props_with_include_rings(&mut self, include_rings: bool) {
        expect_ok(mol_clear_computed_props(self.as_ffi(), include_rings))
    }
}

prop_getters!(Mol);
prop_setters!(Mol, mol);

/// An empty molecule. Mirrors `Mol()`.
impl Default for Mol {
    fn default() -> Mol {
        Mol::from_ffi(expect_ok(ffi::mol_new()))
    }
}

/// A deep copy. Mirrors `Mol.__copy__` and `Mol.__deepcopy__`.
impl Clone for Mol {
    fn clone(&self) -> Mol {
        Mol::from_ffi(expect_ok(ffi::mol_clone(self.as_ffi())))
    }
}

impl fmt::Debug for Mol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Mol")
            .field("num_atoms", &self.get_num_atoms())
            .field("num_bonds", &self.get_num_bonds())
            .finish()
    }
}
