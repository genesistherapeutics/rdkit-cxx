use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::Mol;
use crate::chem::rdchem::atom::expect_ok_index;
use crate::chem::rdchem::props::{prop_getters, prop_setters};
use crate::error::expect_ok;
use crate::ffi::conformer::ffi;
use crate::ffi::props::ffi::{RDProps, conformer_as_props};
use crate::geometry::Point3D;

/// Stamps the conformer getters onto a type with
/// `fn raw(&self) -> &ffi::Conformer` (CONVENTIONS.md W3).
macro_rules! conformer_getters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `Conformer.GetId`.
            pub fn get_id(&self) -> u32 {
                expect_ok(ffi::conformer_get_id(self.raw()))
            }

            /// Mirrors `Conformer.GetNumAtoms`.
            pub fn get_num_atoms(&self) -> u32 {
                expect_ok(ffi::conformer_get_num_atoms(self.raw()))
            }

            /// Mirrors `Conformer.Is3D`.
            pub fn is_3d(&self) -> bool {
                expect_ok(ffi::conformer_is_3d(self.raw()))
            }

            /// Mirrors `Conformer.HasOwningMol`.
            pub fn has_owning_mol(&self) -> bool {
                expect_ok(ffi::conformer_has_owning_mol(self.raw()))
            }

            /// Mirrors `Conformer.GetAtomPosition`.
            pub fn get_atom_position(&self, aid: u32) -> Result<Point3D> {
                Ok(ffi::conformer_get_atom_position(self.raw(), aid)?)
            }

            /// Every atom's position. Mirrors `Conformer.GetPositions`, which
            /// returns an array of shape (atoms, 3).
            pub fn get_positions(&self) -> Vec<Point3D> {
                expect_ok(ffi::conformer_get_positions(self.raw()))
            }
        }
    };
}

/// Stamps the conformer setters onto a type with
/// `fn raw_mut(&mut self) -> Pin<&mut ffi::Conformer>` (W3).
macro_rules! conformer_setters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `Conformer.SetAtomPosition`.
            pub fn set_atom_position(&mut self, aid: u32, loc: Point3D) -> Result<()> {
                Ok(ffi::conformer_set_atom_position(self.raw_mut(), aid, &loc)?)
            }

            /// Sets every atom's position; `positions` has one point per atom.
            /// Mirrors `Conformer.SetPositions`.
            pub fn set_positions(&mut self, positions: &[Point3D]) -> Result<()> {
                Ok(ffi::conformer_set_positions(self.raw_mut(), positions)?)
            }

            /// Mirrors `Conformer.Set3D`.
            pub fn set_3d(&mut self, v: bool) {
                expect_ok(ffi::conformer_set_3d(self.raw_mut(), v))
            }

            /// Mirrors `Conformer.SetId`.
            pub fn set_id(&mut self, id: u32) {
                expect_ok(ffi::conformer_set_id(self.raw_mut(), id))
            }
        }
    };
}

/// A conformer outside any molecule. Mirrors `rdchem.Conformer()`;
/// `Mol::add_conformer` copies it into a molecule.
pub struct Conformer {
    /// Never null.
    ptr: UniquePtr<ffi::Conformer>,
}

impl Conformer {
    fn from_ffi(ptr: UniquePtr<ffi::Conformer>) -> Conformer {
        assert!(!ptr.is_null(), "RDKit returned a null conformer");
        Conformer { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::Conformer {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::Conformer> {
        self.ptr.pin_mut()
    }

    fn props(&self) -> &RDProps {
        expect_ok(conformer_as_props(self.raw()))
    }

    /// A conformer with `num_atoms` atoms at the origin. Mirrors
    /// `Conformer(numAtoms)`.
    pub fn new(num_atoms: u32) -> Conformer {
        Conformer::from_ffi(expect_ok(ffi::conformer_with_num_atoms(num_atoms)))
    }
}

/// A conformer with no atoms. Mirrors `Conformer()`.
impl Default for Conformer {
    fn default() -> Conformer {
        Conformer::from_ffi(expect_ok(ffi::conformer_new()))
    }
}

/// Mirrors `Conformer(other)`.
impl Clone for Conformer {
    fn clone(&self) -> Conformer {
        Conformer::from_ffi(expect_ok(ffi::conformer_clone(self.raw())))
    }
}

impl fmt::Debug for Conformer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Conformer")
            .field("id", &self.get_id())
            .field("num_atoms", &self.get_num_atoms())
            .finish()
    }
}

conformer_getters!(Conformer);
conformer_setters!(Conformer);
prop_getters!(Conformer);
prop_setters!(Conformer);

/// A conformer in a molecule, borrowed. Mirrors `rdchem.Conformer` reached
/// through a `Mol`.
#[derive(Clone, Copy)]
pub struct ConformerRef<'m> {
    mol: &'m Mol,
    conf: &'m ffi::Conformer,
}

impl<'m> ConformerRef<'m> {
    pub(crate) fn new(mol: &'m Mol, conf: &'m ffi::Conformer) -> ConformerRef<'m> {
        ConformerRef { mol, conf }
    }

    fn raw(&self) -> &'m ffi::Conformer {
        self.conf
    }

    fn props(&self) -> &'m RDProps {
        expect_ok(conformer_as_props(self.conf))
    }

    /// Mirrors `Conformer.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &'m Mol {
        self.mol
    }
}

conformer_getters!(ConformerRef<'m>);
prop_getters!(ConformerRef<'m>);

impl fmt::Debug for ConformerRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ConformerRef")
            .field("id", &self.get_id())
            .finish()
    }
}

/// A conformer in a molecule, borrowed exclusively for editing. It holds the
/// conformer's id and finds it on each call (W2).
pub struct ConformerMut<'m> {
    mol: &'m mut Mol,
    id: i32,
}

impl<'m> ConformerMut<'m> {
    /// A conformer with `id` exists in `mol`.
    pub(crate) fn new(mol: &'m mut Mol, id: i32) -> ConformerMut<'m> {
        ConformerMut { mol, id }
    }

    fn raw(&self) -> &ffi::Conformer {
        expect_ok_index(self.mol.get_conformer_with_id(self.id)).raw()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::Conformer> {
        expect_ok_index(self.mol.conformer_ffi_mut(self.id))
    }

    fn props(&self) -> &RDProps {
        expect_ok(conformer_as_props(self.raw()))
    }

    /// Mirrors `Conformer.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &Mol {
        self.mol
    }
}

conformer_getters!(ConformerMut<'m>);
conformer_setters!(ConformerMut<'m>);
prop_getters!(ConformerMut<'m>);
prop_setters!(ConformerMut<'m>);

impl fmt::Debug for ConformerMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ConformerMut")
            .field("id", &self.get_id())
            .finish()
    }
}
