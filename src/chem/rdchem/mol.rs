use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::error::expect_ok;
use crate::ffi::mol::ffi;

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

    #[allow(dead_code)] // used by later wrappers
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
}

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
