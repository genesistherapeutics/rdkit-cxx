// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::props::{prop_getters, prop_setters};
use crate::chem::rdchem::{BondRef, ChiralType, HybridizationType, Mol};
use crate::error::{expect_ok, to_string};
use crate::ffi::atom::ffi;
use crate::ffi::props::ffi::{RDProps, atom_as_props};

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

            /// The element symbol. Fails for an atomic number RDKit's
            /// periodic table lacks. Mirrors `Atom.GetSymbol`.
            pub fn get_symbol(&self) -> Result<String> {
                to_string(ffi::atom_get_symbol(self.raw())?)
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

            /// The atom map number, or 0 for none. Mirrors `Atom.GetAtomMapNum`.
            pub fn get_atom_map_num(&self) -> i32 {
                expect_ok(ffi::atom_get_atom_map_num(self.raw()))
            }

            /// The isotope's mass, or the element's average atomic weight.
            /// Fails for an atomic number RDKit's periodic table lacks.
            /// Mirrors `Atom.GetMass`.
            pub fn get_mass(&self) -> Result<f64> {
                Ok(ffi::atom_get_mass(self.raw())?)
            }

            /// Whether implicit hydrogens are disallowed. Mirrors
            /// `Atom.GetNoImplicit`.
            pub fn get_no_implicit(&self) -> bool {
                expect_ok(ffi::atom_get_no_implicit(self.raw()))
            }

            /// Mirrors `Atom.GetNumRadicalElectrons`.
            pub fn get_num_radical_electrons(&self) -> u32 {
                expect_ok(ffi::atom_get_num_radical_electrons(self.raw()))
            }

            /// Mirrors `Atom.HasOwningMol`.
            pub fn has_owning_mol(&self) -> bool {
                expect_ok(ffi::atom_has_owning_mol(self.raw()))
            }

            /// A copy of the atom's monomer information, if any. Mirrors
            /// `Atom.GetMonomerInfo`, which returns the atom's own object.
            pub fn get_monomer_info(&self) -> Option<$crate::chem::rdchem::AtomMonomerInfo> {
                let info = expect_ok($crate::ffi::monomer_info::ffi::atom_get_monomer_info(self.raw()));
                (!info.is_null()).then(|| $crate::chem::rdchem::AtomMonomerInfo::from_ffi(info))
            }

            /// A copy of the atom's PDB residue information, if any; other
            /// kinds of monomer information are `Error::Value`. Mirrors
            /// `Atom.GetPDBResidueInfo`.
            pub fn get_pdb_residue_info(
                &self,
            ) -> Result<Option<$crate::chem::rdchem::AtomPDBResidueInfo>> {
                let info = $crate::ffi::monomer_info::ffi::atom_get_pdb_residue_info(self.raw())?;
                Ok((!info.is_null()).then(|| $crate::chem::rdchem::AtomPDBResidueInfo::from_ffi(info)))
            }

            /// Mirrors `Atom.HasQuery`.
            pub fn has_query(&self) -> bool {
                expect_ok(ffi::atom_has_query(self.raw()))
            }

            /// Whether the valences need recomputing. Mirrors
            /// `Atom.NeedsUpdatePropertyCache`.
            pub fn needs_update_property_cache(&self) -> bool {
                expect_ok(ffi::atom_needs_update_property_cache(self.raw()))
            }

            /// The explicit valence: the bonds' and explicit hydrogens'.
            /// `Error::Runtime` outside a molecule or before the property
            /// cache is computed. Mirrors `Atom.GetExplicitValence`.
            pub fn get_explicit_valence(&self) -> Result<u32> {
                Ok(ffi::atom_get_explicit_valence(self.raw())?)
            }

            /// The implicit valence: the implicit hydrogens'. `Error::Runtime`
            /// outside a molecule; −1 before the property cache is computed,
            /// as RDKit 2024.03 returns it. Mirrors `Atom.GetImplicitValence`.
            pub fn get_implicit_valence(&self) -> Result<i32> {
                Ok(ffi::atom_get_implicit_valence(self.raw())?)
            }

            /// The atom as SMARTS, or as SMILES when it has no query. Mirrors
            /// `Atom.GetSmarts`.
            pub fn get_smarts(&self) -> Result<String> {
                self.get_smarts_with_params(&$crate::chem::rdchem::GetSmartsParams::default())
            }

            /// Mirrors `Atom.GetSmarts(doKekule, allHsExplicit, isomericSmiles)`.
            pub fn get_smarts_with_params(
                &self,
                params: &$crate::chem::rdchem::GetSmartsParams,
            ) -> Result<String> {
                to_string(ffi::atom_get_smarts(
                    self.raw(),
                    params.do_kekule,
                    params.all_hs_explicit,
                    params.isomeric_smiles,
                )?)
            }
        }
    };
}

/// Stamps the getters that need the owning molecule onto a view with
/// `fn raw(&self) -> &ffi::Atom` (W4).
macro_rules! atom_view_getters {
    ($type:ident) => {
        impl $type<'_> {
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

            /// Whether the atom is in a ring of `size` atoms, finding the
            /// molecule's rings if needed. Mirrors `Atom.IsInRingSize`.
            pub fn is_in_ring_size(&self, size: i32) -> bool {
                expect_ok(ffi::atom_is_in_ring_size(self.raw(), size))
            }

            /// Whether the atom's valence is not allowed. Needs the property
            /// cache. Mirrors `Atom.HasValenceViolation`.
            pub fn has_valence_violation(&self) -> Result<bool> {
                Ok(ffi::atom_has_valence_violation(self.raw())?)
            }
        }
    };
}

/// Stamps the setters onto a type with
/// `fn raw_mut(&mut self) -> Pin<&mut ffi::Atom>` (W3).
macro_rules! atom_setters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `Atom.SetAtomicNum`.
            pub fn set_atomic_num(&mut self, new_num: i32) {
                expect_ok(ffi::atom_set_atomic_num(self.raw_mut(), new_num))
            }

            /// Mirrors `Atom.SetFormalCharge`.
            pub fn set_formal_charge(&mut self, what: i32) {
                expect_ok(ffi::atom_set_formal_charge(self.raw_mut(), what))
            }

            /// Mirrors `Atom.SetIsAromatic`.
            pub fn set_is_aromatic(&mut self, what: bool) {
                expect_ok(ffi::atom_set_is_aromatic(self.raw_mut(), what))
            }

            /// Mirrors `Atom.SetIsotope`.
            pub fn set_isotope(&mut self, what: u32) {
                expect_ok(ffi::atom_set_isotope(self.raw_mut(), what))
            }

            /// Mirrors `Atom.SetNumExplicitHs`.
            pub fn set_num_explicit_hs(&mut self, what: u32) {
                expect_ok(ffi::atom_set_num_explicit_hs(self.raw_mut(), what))
            }

            /// Disallows implicit hydrogens when `what`. Mirrors
            /// `Atom.SetNoImplicit`.
            pub fn set_no_implicit(&mut self, what: bool) {
                expect_ok(ffi::atom_set_no_implicit(self.raw_mut(), what))
            }

            /// Mirrors `Atom.SetHybridization`.
            pub fn set_hybridization(&mut self, what: HybridizationType) {
                expect_ok(ffi::atom_set_hybridization(self.raw_mut(), what.to_ffi()))
            }

            /// Mirrors `Atom.SetChiralTag`.
            pub fn set_chiral_tag(&mut self, what: ChiralType) {
                expect_ok(ffi::atom_set_chiral_tag(self.raw_mut(), what.to_ffi()))
            }

            /// Sets the atom map number; 0 clears it. Mirrors
            /// `Atom.SetAtomMapNum`.
            pub fn set_atom_map_num(&mut self, mapno: i32) {
                expect_ok(ffi::atom_set_atom_map_num(self.raw_mut(), mapno, false))
            }

            /// With `strict`, fails unless `mapno` is in 0..1000. Mirrors
            /// `Atom.SetAtomMapNum(mapno, strict)`.
            pub fn set_atom_map_num_with_strict(&mut self, mapno: i32, strict: bool) -> Result<()> {
                Ok(ffi::atom_set_atom_map_num(self.raw_mut(), mapno, strict)?)
            }

            /// Mirrors `Atom.SetNumRadicalElectrons`.
            pub fn set_num_radical_electrons(&mut self, num: u32) {
                expect_ok(ffi::atom_set_num_radical_electrons(self.raw_mut(), num))
            }

            /// Gives the atom a copy of `info`. Mirrors `Atom.SetMonomerInfo`.
            pub fn set_monomer_info(&mut self, info: &$crate::chem::rdchem::AtomMonomerInfo) {
                expect_ok($crate::ffi::monomer_info::ffi::atom_set_monomer_info(
                    self.raw_mut(),
                    info.monomer(),
                ))
            }

            /// Gives the atom a copy of `info`. Mirrors
            /// `Atom.SetPDBResidueInfo`.
            pub fn set_pdb_residue_info(&mut self, info: &$crate::chem::rdchem::AtomPDBResidueInfo) {
                expect_ok($crate::ffi::monomer_info::ffi::atom_set_monomer_info(
                    self.raw_mut(),
                    info.monomer(),
                ))
            }

            /// Inverts tetrahedral chirality, returning whether there was any.
            /// Mirrors `Atom.InvertChirality`.
            pub fn invert_chirality(&mut self) -> bool {
                expect_ok(ffi::atom_invert_chirality(self.raw_mut()))
            }
        }
    };
}

/// Parameters for `get_smarts_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetSmartsParams {
    /// Writes the Kekulé form.
    pub do_kekule: bool,
    /// Writes hydrogens in brackets.
    pub all_hs_explicit: bool,
    /// Includes stereochemistry and isotopes.
    pub isomeric_smiles: bool,
}

impl Default for GetSmartsParams {
    fn default() -> GetSmartsParams {
        GetSmartsParams {
            do_kekule: false,
            all_hs_explicit: false,
            isomeric_smiles: true,
        }
    }
}

/// `StereoInfo::NOATOM`, read from C++ for a test that pins the constant.
#[cfg(test)]
pub(crate) fn noatom_for_tests() -> u32 {
    expect_ok(ffi::atom_noatom())
}

/// An atom outside any molecule. Mirrors `rdchem.Atom()`; `RWMol::add_atom`
/// copies it into a molecule.
pub struct Atom {
    /// Never null.
    ptr: UniquePtr<ffi::Atom>,
}

impl Atom {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::Atom>) -> Atom {
        assert!(!ptr.is_null(), "RDKit returned a null atom");
        Atom { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::Atom {
        self.ptr.as_ref().unwrap()
    }

    pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::Atom> {
        self.ptr.pin_mut()
    }

    /// An atom of the element with symbol `what`. Mirrors `Atom(what)`.
    pub fn new(what: &str) -> Result<Atom> {
        Ok(Atom::from_ffi(ffi::atom_new(what)?))
    }

    /// An atom with atomic number `num`. Mirrors `Atom(num)`.
    pub fn from_num(num: u32) -> Atom {
        Atom::from_ffi(expect_ok(ffi::atom_from_num(num)))
    }
}

/// Mirrors `Atom(other)`.
impl Clone for Atom {
    fn clone(&self) -> Atom {
        Atom::from_ffi(expect_ok(ffi::atom_clone(self.raw())))
    }
}

impl fmt::Debug for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Atom")
            .field("atomic_num", &self.get_atomic_num())
            .finish()
    }
}

impl Atom {
    fn props(&self) -> &RDProps {
        expect_ok(atom_as_props(self.raw()))
    }
}

atom_getters!(Atom);
atom_setters!(Atom);
prop_getters!(Atom);
prop_setters!(Atom);

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

    pub(crate) fn raw(&self) -> &'m ffi::Atom {
        self.atom
    }

    /// The atom's neighbors. Mirrors `Atom.GetNeighbors`.
    pub fn get_neighbors(&self) -> Vec<AtomRef<'m>> {
        neighbors(self.mol, self.atom)
    }

    /// The atom's bonds. Mirrors `Atom.GetBonds`.
    pub fn get_bonds(&self) -> Vec<BondRef<'m>> {
        bonds(self.mol, self.atom)
    }

    /// Mirrors `Atom.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &'m Mol {
        self.mol
    }
}

impl<'m> AtomRef<'m> {
    fn props(&self) -> &'m RDProps {
        expect_ok(atom_as_props(self.atom))
    }
}

atom_getters!(AtomRef<'m>);
atom_view_getters!(AtomRef);
prop_getters!(AtomRef<'m>);

impl fmt::Debug for AtomRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AtomRef")
            .field("idx", &self.get_idx())
            .finish()
    }
}

/// An atom in a molecule, borrowed exclusively for editing. Mirrors
/// `rdchem.Atom` reached through a `Mol`.
///
/// It holds an index and finds the atom on each call, so it needs no
/// `unsafe` (CONVENTIONS.md W2).
pub struct AtomMut<'m> {
    mol: &'m mut Mol,
    idx: u32,
}

impl<'m> AtomMut<'m> {
    /// `idx` has been checked against the molecule.
    pub(crate) fn new(mol: &'m mut Mol, idx: u32) -> AtomMut<'m> {
        AtomMut { mol, idx }
    }

    pub(crate) fn raw(&self) -> &ffi::Atom {
        expect_ok_index(self.mol.get_atom_with_idx(self.idx)).raw()
    }

    pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::Atom> {
        expect_ok_index(self.mol.atom_ffi_mut(self.idx))
    }

    /// The atom's neighbors. Mirrors `Atom.GetNeighbors`.
    pub fn get_neighbors(&self) -> Vec<AtomRef<'_>> {
        neighbors(self.mol, self.raw())
    }

    /// The atom's bonds. Mirrors `Atom.GetBonds`.
    pub fn get_bonds(&self) -> Vec<BondRef<'_>> {
        bonds(self.mol, self.raw())
    }

    /// Mirrors `Atom.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &Mol {
        self.mol
    }

    /// Recomputes implicit valence and related properties, failing on a bad
    /// valence. Mirrors `Atom.UpdatePropertyCache`.
    pub fn update_property_cache(&mut self) -> Result<()> {
        self.update_property_cache_with_strict(true)
    }

    /// Recomputes implicit valence and related properties, failing on a bad
    /// valence only when `strict`. Mirrors `Atom.UpdatePropertyCache(strict)`.
    pub fn update_property_cache_with_strict(&mut self, strict: bool) -> Result<()> {
        Ok(ffi::atom_update_property_cache(self.raw_mut(), strict)?)
    }
}

impl AtomMut<'_> {
    fn props(&self) -> &RDProps {
        expect_ok(atom_as_props(self.raw()))
    }
}

atom_getters!(AtomMut<'m>);
atom_view_getters!(AtomMut);
atom_setters!(AtomMut<'m>);
prop_getters!(AtomMut<'m>);
prop_setters!(AtomMut<'m>);

impl fmt::Debug for AtomMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AtomMut").field("idx", &self.idx).finish()
    }
}

fn neighbors<'m>(mol: &'m Mol, atom: &ffi::Atom) -> Vec<AtomRef<'m>> {
    expect_ok(ffi::atom_get_neighbors(atom))
        .into_iter()
        .map(|idx| expect_ok_index(mol.get_atom_with_idx(idx)))
        .collect()
}

fn bonds<'m>(mol: &'m Mol, atom: &ffi::Atom) -> Vec<BondRef<'m>> {
    expect_ok(ffi::atom_get_bonds(atom))
        .into_iter()
        .map(|idx| expect_ok_index(mol.get_bond_with_idx(idx)))
        .collect()
}

/// Unwraps a lookup by an index that RDKit returned or the crate checked.
pub(crate) fn expect_ok_index<T>(result: Result<T>) -> T {
    result.unwrap_or_else(|err| panic!("RDKit returned a bad index: {err}"))
}

/// Stamps the `ExplicitBitVect` property methods, which Python defines on
/// `Atom` only, onto a type with `fn props(&self) -> &RDProps` (P1).
macro_rules! atom_bit_vect_props {
    ($type:ident $(<$lifetime:lifetime>)?, setter) => {
        atom_bit_vect_props!($type $(<$lifetime>)?);

        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `Atom.SetExplicitBitVectProp`.
            pub fn set_explicit_bit_vect_prop(
                &mut self,
                key: &str,
                val: &crate::data_structs::ExplicitBitVect,
            ) {
                expect_ok(crate::ffi::props::ffi::props_set_explicit_bit_vect_prop(
                    self.props(),
                    key,
                    val.raw(),
                ))
            }
        }
    };
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// A missing key is `Error::Key`, and another type `Error::Value`.
            /// Mirrors `Atom.GetExplicitBitVectProp`.
            pub fn get_explicit_bit_vect_prop(
                &self,
                key: &str,
            ) -> Result<crate::data_structs::ExplicitBitVect> {
                Ok(crate::data_structs::ExplicitBitVect::from_ffi(
                    crate::ffi::props::ffi::props_get_explicit_bit_vect_prop(self.props(), key)?,
                ))
            }
        }
    };
}

atom_bit_vect_props!(Atom, setter);
atom_bit_vect_props!(AtomRef<'m>);
atom_bit_vect_props!(AtomMut<'m>, setter);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noatom_matches_rdkit() {
        assert_eq!(crate::chem::rdchem::StereoInfo::NOATOM, noatom_for_tests());
    }
}
