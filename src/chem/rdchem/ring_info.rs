use std::fmt;

use crate::Result;
use crate::chem::rdchem::Mol;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::ring_info::ffi;

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

    /// Mirrors `RingInfo.NumRings`.
    pub fn num_rings(&self) -> Result<u32> {
        Ok(ffi::ring_info_num_rings(self.mol.as_ffi())?)
    }

    /// Each ring's atom indices, empty when rings were never found. Mirrors
    /// `RingInfo.AtomRings`.
    pub fn atom_rings(&self) -> Vec<Vec<u32>> {
        IndexList::into_nested(expect_ok(ffi::ring_info_atom_rings(self.mol.as_ffi())))
    }

    /// Each ring's bond indices, empty when rings were never found. Mirrors
    /// `RingInfo.BondRings`.
    pub fn bond_rings(&self) -> Vec<Vec<u32>> {
        IndexList::into_nested(expect_ok(ffi::ring_info_bond_rings(self.mol.as_ffi())))
    }

    /// The number of rings an atom is in. Mirrors `RingInfo.NumAtomRings`.
    pub fn num_atom_rings(&self, idx: u32) -> Result<u32> {
        Ok(ffi::ring_info_num_atom_rings(self.mol.as_ffi(), idx)?)
    }

    /// The number of rings a bond is in. Mirrors `RingInfo.NumBondRings`.
    pub fn num_bond_rings(&self, idx: u32) -> Result<u32> {
        Ok(ffi::ring_info_num_bond_rings(self.mol.as_ffi(), idx)?)
    }

    /// Mirrors `RingInfo.IsAtomInRingOfSize`.
    pub fn is_atom_in_ring_of_size(&self, idx: u32, size: u32) -> Result<bool> {
        Ok(ffi::ring_info_is_atom_in_ring_of_size(
            self.mol.as_ffi(),
            idx,
            size,
        )?)
    }

    /// Mirrors `RingInfo.IsBondInRingOfSize`.
    pub fn is_bond_in_ring_of_size(&self, idx: u32, size: u32) -> Result<bool> {
        Ok(ffi::ring_info_is_bond_in_ring_of_size(
            self.mol.as_ffi(),
            idx,
            size,
        )?)
    }
}

impl fmt::Debug for RingInfoRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RingInfoRef")
            .field("num_rings", &self.num_rings().ok())
            .finish()
    }
}
