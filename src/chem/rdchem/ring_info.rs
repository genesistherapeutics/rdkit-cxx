// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;

use crate::Result;
use crate::chem::rdchem::Mol;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::ring_info::ffi;

/// Stamps the ring queries onto a type with `fn mol(&self) -> &Mol`
/// (CONVENTIONS.md W3).
macro_rules! ring_info_getters {
    ($type:ident) => {
        impl $type<'_> {
            /// Mirrors `RingInfo.NumRings`.
            pub fn num_rings(&self) -> Result<u32> {
                Ok(ffi::ring_info_num_rings(self.mol().as_ffi())?)
            }

            /// Each ring's atom indices, empty when rings were never found.
            /// Mirrors `RingInfo.AtomRings`.
            pub fn atom_rings(&self) -> Vec<Vec<u32>> {
                IndexList::into_nested(expect_ok(ffi::ring_info_atom_rings(self.mol().as_ffi())))
            }

            /// Each ring's bond indices, empty when rings were never found.
            /// Mirrors `RingInfo.BondRings`.
            pub fn bond_rings(&self) -> Vec<Vec<u32>> {
                IndexList::into_nested(expect_ok(ffi::ring_info_bond_rings(self.mol().as_ffi())))
            }

            /// The number of rings an atom is in. Mirrors `RingInfo.NumAtomRings`.
            pub fn num_atom_rings(&self, idx: u32) -> Result<u32> {
                Ok(ffi::ring_info_num_atom_rings(self.mol().as_ffi(), idx)?)
            }

            /// The number of rings a bond is in. Mirrors `RingInfo.NumBondRings`.
            pub fn num_bond_rings(&self, idx: u32) -> Result<u32> {
                Ok(ffi::ring_info_num_bond_rings(self.mol().as_ffi(), idx)?)
            }

            /// Mirrors `RingInfo.IsAtomInRingOfSize`.
            pub fn is_atom_in_ring_of_size(&self, idx: u32, size: u32) -> Result<bool> {
                Ok(ffi::ring_info_is_atom_in_ring_of_size(
                    self.mol().as_ffi(),
                    idx,
                    size,
                )?)
            }

            /// Mirrors `RingInfo.IsBondInRingOfSize`.
            pub fn is_bond_in_ring_of_size(&self, idx: u32, size: u32) -> Result<bool> {
                Ok(ffi::ring_info_is_bond_in_ring_of_size(
                    self.mol().as_ffi(),
                    idx,
                    size,
                )?)
            }

            /// The indices of the rings an atom is in. Mirrors
            /// `RingInfo.AtomMembers`.
            pub fn atom_members(&self, idx: u32) -> Result<Vec<i32>> {
                Ok(ffi::ring_info_atom_members(self.mol().as_ffi(), idx)?)
            }

            /// The indices of the rings a bond is in. Mirrors
            /// `RingInfo.BondMembers`.
            pub fn bond_members(&self, idx: u32) -> Result<Vec<i32>> {
                Ok(ffi::ring_info_bond_members(self.mol().as_ffi(), idx)?)
            }

            /// The sizes of the rings an atom is in. Mirrors
            /// `RingInfo.AtomRingSizes`.
            pub fn atom_ring_sizes(&self, idx: u32) -> Result<Vec<i32>> {
                Ok(ffi::ring_info_atom_ring_sizes(self.mol().as_ffi(), idx)?)
            }

            /// The sizes of the rings a bond is in. Mirrors
            /// `RingInfo.BondRingSizes`.
            pub fn bond_ring_sizes(&self, idx: u32) -> Result<Vec<i32>> {
                Ok(ffi::ring_info_bond_ring_sizes(self.mol().as_ffi(), idx)?)
            }

            /// The smallest ring an atom is in, or 0. Mirrors
            /// `RingInfo.MinAtomRingSize`.
            pub fn min_atom_ring_size(&self, idx: u32) -> Result<u32> {
                Ok(ffi::ring_info_min_atom_ring_size(self.mol().as_ffi(), idx)?)
            }

            /// The smallest ring a bond is in, or 0. Mirrors
            /// `RingInfo.MinBondRingSize`.
            pub fn min_bond_ring_size(&self, idx: u32) -> Result<u32> {
                Ok(ffi::ring_info_min_bond_ring_size(self.mol().as_ffi(), idx)?)
            }

            /// Mirrors `RingInfo.AreAtomsInSameRing`.
            pub fn are_atoms_in_same_ring(&self, idx1: u32, idx2: u32) -> Result<bool> {
                Ok(ffi::ring_info_are_atoms_in_same_ring(
                    self.mol().as_ffi(),
                    idx1,
                    idx2,
                )?)
            }

            /// Mirrors `RingInfo.AreAtomsInSameRingOfSize`.
            pub fn are_atoms_in_same_ring_of_size(
                &self,
                idx1: u32,
                idx2: u32,
                size: u32,
            ) -> Result<bool> {
                Ok(ffi::ring_info_are_atoms_in_same_ring_of_size(
                    self.mol().as_ffi(),
                    idx1,
                    idx2,
                    size,
                )?)
            }

            /// Mirrors `RingInfo.AreBondsInSameRing`.
            pub fn are_bonds_in_same_ring(&self, idx1: u32, idx2: u32) -> Result<bool> {
                Ok(ffi::ring_info_are_bonds_in_same_ring(
                    self.mol().as_ffi(),
                    idx1,
                    idx2,
                )?)
            }

            /// Mirrors `RingInfo.AreBondsInSameRingOfSize`.
            pub fn are_bonds_in_same_ring_of_size(
                &self,
                idx1: u32,
                idx2: u32,
                size: u32,
            ) -> Result<bool> {
                Ok(ffi::ring_info_are_bonds_in_same_ring_of_size(
                    self.mol().as_ffi(),
                    idx1,
                    idx2,
                    size,
                )?)
            }

            /// Whether a ring shares a bond with another. Mirrors
            /// `RingInfo.IsRingFused`.
            pub fn is_ring_fused(&self, ring_idx: u32) -> Result<bool> {
                Ok(ffi::ring_info_is_ring_fused(self.mol().as_ffi(), ring_idx)?)
            }

            /// Whether two rings share a bond. Mirrors `RingInfo.AreRingsFused`.
            pub fn are_rings_fused(&self, ring1_idx: u32, ring2_idx: u32) -> Result<bool> {
                Ok(ffi::ring_info_are_rings_fused(
                    self.mol().as_ffi(),
                    ring1_idx,
                    ring2_idx,
                )?)
            }

            /// The number of a ring's bonds shared with other rings. Mirrors
            /// `RingInfo.NumFusedBonds`.
            pub fn num_fused_bonds(&self, ring_idx: u32) -> Result<u32> {
                Ok(ffi::ring_info_num_fused_bonds(
                    self.mol().as_ffi(),
                    ring_idx,
                )?)
            }

            /// Whether ring families have been found, as `find_ring_families`
            /// finds them; parsing does not. Mirrors
            /// `RingInfo.AreRingFamiliesInitialized`.
            pub fn are_ring_families_initialized(&self) -> bool {
                expect_ok(ffi::ring_info_are_ring_families_initialized(
                    self.mol().as_ffi(),
                ))
            }

            /// Needs the ring families. Mirrors `RingInfo.NumRingFamilies`.
            pub fn num_ring_families(&self) -> Result<u32> {
                Ok(ffi::ring_info_num_ring_families(self.mol().as_ffi())?)
            }

            /// Needs the ring families: `Error::Runtime` before they are found,
            /// where RDKit 2025.09 converts RingDecomposerLib's failure value
            /// out of `unsigned int`'s range, and `Error::Value` for more
            /// relevant cycles than `u32` holds. Mirrors
            /// `RingInfo.NumRelevantCycles`.
            pub fn num_relevant_cycles(&self) -> Result<u32> {
                Ok(ffi::ring_info_num_relevant_cycles(self.mol().as_ffi())?)
            }

            /// Each ring family's atom indices, empty when they were never
            /// found. Mirrors `RingInfo.AtomRingFamilies`.
            pub fn atom_ring_families(&self) -> Vec<Vec<u32>> {
                IndexList::into_nested(expect_ok(ffi::ring_info_atom_ring_families(
                    self.mol().as_ffi(),
                )))
            }

            /// Each ring family's bond indices, empty when they were never
            /// found. Mirrors `RingInfo.BondRingFamilies`.
            pub fn bond_ring_families(&self) -> Vec<Vec<u32>> {
                IndexList::into_nested(expect_ok(ffi::ring_info_bond_ring_families(
                    self.mol().as_ffi(),
                )))
            }
        }
    };
}

/// A molecule's ring information, borrowed. Mirrors `rdchem.RingInfo`.
///
/// Parsing finds rings. On a molecule built atom by atom, the counting queries
/// return `Error::Runtime` ("RingInfo not initialized") until something finds
/// them, and the ring lists are empty.
///
/// ```compile_fail
/// use rdkit_cxx::chem::mol_from_smiles;
/// let info = {
///     let mol = mol_from_smiles("C1CC1").unwrap();
///     mol.get_ring_info()
/// };
/// info.num_rings();
/// ```
#[derive(Clone, Copy)]
pub struct RingInfoRef<'m> {
    mol: &'m Mol,
}

impl<'m> RingInfoRef<'m> {
    pub(crate) fn new(mol: &'m Mol) -> RingInfoRef<'m> {
        RingInfoRef { mol }
    }

    fn mol(&self) -> &'m Mol {
        self.mol
    }
}

ring_info_getters!(RingInfoRef);

impl fmt::Debug for RingInfoRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RingInfoRef")
            .field("num_rings", &self.num_rings().ok())
            .finish()
    }
}

/// A molecule's ring information, borrowed exclusively to add rings. Mirrors
/// `rdchem.RingInfo` (W2).
pub struct RingInfoMut<'m> {
    mol: &'m mut Mol,
}

impl<'m> RingInfoMut<'m> {
    pub(crate) fn new(mol: &'m mut Mol) -> RingInfoMut<'m> {
        RingInfoMut { mol }
    }

    fn mol(&self) -> &Mol {
        self.mol
    }

    /// Records a ring, given its atoms and bonds in order, initializing the
    /// ring information if needed. Mirrors `RingInfo.AddRing`.
    ///
    /// `Error::Value` unless the ring is a cycle: 3 or more distinct atoms,
    /// with bond `i` joining atoms `i` and `i + 1`, and the last bond joining
    /// the last atom to the first. RDKit's ring algorithms and its pickle
    /// reader assume one, and read out of bounds otherwise. `Error::Index`
    /// for an atom or bond index out of range.
    pub fn add_ring(&mut self, atom_ids: &[u32], bond_ids: &[u32]) -> Result<()> {
        Ok(ffi::ring_info_add_ring(
            self.mol.pin_mut(),
            atom_ids,
            bond_ids,
        )?)
    }
}

ring_info_getters!(RingInfoMut);

impl fmt::Debug for RingInfoMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RingInfoMut")
            .field("num_rings", &self.num_rings().ok())
            .finish()
    }
}
