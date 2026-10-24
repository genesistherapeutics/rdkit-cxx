use std::fmt;

use crate::Result;
use crate::chem::rdchem::{BondRef, ChiralType, HybridizationType, Mol};
use crate::error::{expect_ok, to_string};
use crate::ffi::atom::ffi;

/// Stamps the getters that need no owning molecule onto a type with
/// `fn raw(&self) -> &ffi::Atom` (CONVENTIONS.md W3).
macro_rules! atom_getters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// The atom's index in its molecule. Mirrors `Atom.GetIdx`.
            pub fn get_idx(&self) -> u32 {
                expect_ok(ffi::atom_get_idx(self.raw()))
            }

            /// Mirrors `Atom.GetAtomicNum`.
            pub fn get_atomic_num(&self) -> i32 {
                expect_ok(ffi::atom_get_atomic_num(self.raw()))
            }

            /// The element symbol. Mirrors `Atom.GetSymbol`.
            pub fn get_symbol(&self) -> String {
                to_string(expect_ok(ffi::atom_get_symbol(self.raw()))).expect("element symbols are ASCII")
            }

            /// Mirrors `Atom.GetFormalCharge`.
            pub fn get_formal_charge(&self) -> i32 {
                expect_ok(ffi::atom_get_formal_charge(self.raw()))
            }

            /// Mirrors `Atom.GetIsAromatic`.
            pub fn get_is_aromatic(&self) -> bool {
                expect_ok(ffi::atom_get_is_aromatic(self.raw()))
            }

            /// Mirrors `Atom.GetIsotope`.
            pub fn get_isotope(&self) -> u32 {
                expect_ok(ffi::atom_get_isotope(self.raw()))
            }

            /// Mirrors `Atom.GetNumExplicitHs`.
            pub fn get_num_explicit_hs(&self) -> u32 {
                expect_ok(ffi::atom_get_num_explicit_hs(self.raw()))
            }

            /// Mirrors `Atom.GetHybridization`.
            pub fn get_hybridization(&self) -> HybridizationType {
                HybridizationType::from_ffi(expect_ok(ffi::atom_get_hybridization(self.raw())))
            }

            /// Mirrors `Atom.GetChiralTag`.
            pub fn get_chiral_tag(&self) -> ChiralType {
                ChiralType::from_ffi(expect_ok(ffi::atom_get_chiral_tag(self.raw())))
            }
        }
    };
}

/// Stamps the getters that need the owning molecule onto a view with
/// `fn raw(&self) -> &ffi::Atom` and `fn mol(&self) -> &Mol` (W4).
macro_rules! atom_view_getters {
    ($type:ident) => {
        impl<'m> $type<'m> {
            /// The number of explicit neighbors. Mirrors `Atom.GetDegree`.
            pub fn get_degree(&self) -> u32 {
                expect_ok(ffi::atom_get_degree(self.raw()))
            }

            /// The degree counting hydrogens. Needs the property cache.
            /// Mirrors `Atom.GetTotalDegree`.
            pub fn get_total_degree(&self) -> Result<u32> {
                Ok(ffi::atom_get_total_degree(self.raw())?)
            }

            /// Explicit and implicit hydrogens. Needs the property cache.
            /// Mirrors `Atom.GetTotalNumHs`.
            pub fn get_total_num_hs(&self) -> Result<u32> {
                Ok(ffi::atom_get_total_num_hs(self.raw(), false)?)
            }

            /// Explicit and implicit hydrogens, with neighboring hydrogen atoms
            /// when `include_neighbors`. Needs the property cache. Mirrors
            /// `Atom.GetTotalNumHs(includeNeighbors)`.
            pub fn get_total_num_hs_with_include_neighbors(
                &self,
                include_neighbors: bool,
            ) -> Result<u32> {
                Ok(ffi::atom_get_total_num_hs(self.raw(), include_neighbors)?)
            }

            /// Needs the property cache. Mirrors `Atom.GetNumImplicitHs`.
            pub fn get_num_implicit_hs(&self) -> Result<u32> {
                Ok(ffi::atom_get_num_implicit_hs(self.raw())?)
            }

            /// Needs the property cache. Mirrors `Atom.GetTotalValence`.
            pub fn get_total_valence(&self) -> Result<u32> {
                Ok(ffi::atom_get_total_valence(self.raw())?)
            }

            /// Whether the atom is in a ring, finding the molecule's rings if
            /// needed. Mirrors `Atom.IsInRing`.
            pub fn is_in_ring(&self) -> bool {
                expect_ok(ffi::atom_is_in_ring(self.raw()))
            }
        }
    };
}

/// An atom in a molecule, borrowed. Mirrors `rdchem.Atom` reached through a
/// `Mol`.
///
/// A view cannot outlive its molecule:
///
/// ```compile_fail
/// use rdkit_cxx::chem::mol_from_smiles;
/// let atom = {
///     let mol = mol_from_smiles("C").unwrap();
///     mol.get_atom_with_idx(0).unwrap()
/// };
/// atom.get_idx();
/// ```
#[derive(Clone, Copy)]
pub struct AtomRef<'m> {
    mol: &'m Mol,
    atom: &'m ffi::Atom,
}

impl<'m> AtomRef<'m> {
    pub(crate) fn new(mol: &'m Mol, atom: &'m ffi::Atom) -> AtomRef<'m> {
        AtomRef { mol, atom }
    }

    fn raw(&self) -> &'m ffi::Atom {
        self.atom
    }

    /// The atom's neighbors. Mirrors `Atom.GetNeighbors`.
    pub fn get_neighbors(&self) -> Vec<AtomRef<'m>> {
        let mol = self.mol;
        expect_ok(ffi::atom_get_neighbors(self.atom))
            .into_iter()
            .map(|idx| expect_ok_index(mol.get_atom_with_idx(idx)))
            .collect()
    }

    /// The atom's bonds. Mirrors `Atom.GetBonds`.
    pub fn get_bonds(&self) -> Vec<BondRef<'m>> {
        let mol = self.mol;
        expect_ok(ffi::atom_get_bonds(self.atom))
            .into_iter()
            .map(|idx| expect_ok_index(mol.get_bond_with_idx(idx)))
            .collect()
    }

    /// Mirrors `Atom.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &'m Mol {
        self.mol
    }
}

atom_getters!(AtomRef<'m>);
atom_view_getters!(AtomRef);

impl fmt::Debug for AtomRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AtomRef")
            .field("idx", &self.get_idx())
            .finish()
    }
}

/// Unwraps a lookup by an index RDKit itself returned.
pub(crate) fn expect_ok_index<T>(result: Result<T>) -> T {
    result.unwrap_or_else(|err| panic!("RDKit returned a bad index: {err}"))
}
