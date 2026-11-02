// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolCatalog`: hierarchical catalogs of molecules.
//!
//! A catalog's entries form a directed graph. Entries have indices, their
//! order in the catalog, and bit ids, which `add_entry` assigns in order.
//!
//! RDKit reads the pickles' lengths and indices unchecked, so the wrappers
//! check them first (CONVENTIONS.md U2): a description's length must fit in
//! the pickle, and an edge must name an entry. Molecules in pickles are
//! checked as `Mol::from_pkl_string` checks them (U10).

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::mol_catalog::ffi;

/// A hierarchical catalog of molecules. Mirrors `rdMolCatalog.MolCatalog`.
///
/// It is `Send` but not `Sync`, as the molecules it holds are not
/// (CONVENTIONS.md S1).
pub struct MolCatalog {
    /// Never null.
    ptr: UniquePtr<ffi::MolCatalog>,
}

impl MolCatalog {
    fn raw(&self) -> &ffi::MolCatalog {
        self.ptr.as_ref().unwrap()
    }

    /// A catalog from `serialize`'s bytes. Mirrors `MolCatalog(pickle)`.
    /// `Error::Value` for a pickle whose lengths or indices RDKit would
    /// read out of bounds; Python's reads them. A truncated pickle is
    /// `Error::Runtime`.
    ///
    /// # Safety
    ///
    /// `pickle` must be bytes `serialize` (this crate's or RDKit's)
    /// returned, unmodified: each entry holds a molecule pickle, which RDKit
    /// reads as `Mol::from_pkl_string` does.
    pub unsafe fn from_pickle(pickle: &[u8]) -> Result<MolCatalog> {
        Ok(MolCatalog {
            ptr: ffi::mol_catalog_from_pickle(pickle)?,
        })
    }

    /// The number of entries. Mirrors `MolCatalog.GetNumEntries`.
    pub fn get_num_entries(&self) -> u32 {
        expect_ok(ffi::mol_catalog_get_num_entries(self.raw()))
    }

    /// The number of bit ids assigned. Mirrors `MolCatalog.GetFPLength`.
    pub fn get_fp_length(&self) -> u32 {
        expect_ok(ffi::mol_catalog_get_fp_length(self.raw()))
    }

    /// The catalog as bytes. Mirrors `MolCatalog.Serialize`.
    /// `Error::Runtime` when an entry has no molecule, as in Python, and
    /// `Error::Value` for a description of 2^31 - 1 bytes or more, whose
    /// length RDKit writes in 32 bits, or a molecule `Mol::to_binary`
    /// refuses.
    pub fn serialize(&self) -> Result<Vec<u8>> {
        Ok(ffi::mol_catalog_serialize(self.raw())?)
    }

    /// The description of the entry with a bit id. Mirrors
    /// `MolCatalog.GetBitDescription`. `Error::Index` for a bit past
    /// `get_fp_length`; `Error::Value` when no entry at or after the bit's
    /// index has it, as a pickle can arrange, where Python's dereferences
    /// null.
    pub fn get_bit_description(&self, idx: u32) -> Result<String> {
        to_string(ffi::mol_catalog_get_bit_description(self.raw(), idx)?)
    }

    /// The index of the entry with a bit id, or `None` when no entry at or
    /// after the bit's index has it (where Python returns 2^32 - 1). Mirrors
    /// `MolCatalog.GetBitEntryId`. `Error::Index` for a bit past
    /// `get_fp_length`.
    pub fn get_bit_entry_id(&self, idx: u32) -> Result<Option<u32>> {
        let id = ffi::mol_catalog_get_bit_entry_id(self.raw(), idx)?;
        Ok(u32::try_from(id).ok())
    }

    /// An entry's bit id, or `None` for an entry without one, whose id is
    /// negative (Python returns it modulo 2^32). Mirrors
    /// `MolCatalog.GetEntryBitId`.
    /// `Error::Index` for an index past the entries.
    pub fn get_entry_bit_id(&self, idx: u32) -> Result<Option<u32>> {
        let id = ffi::mol_catalog_get_entry_bit_id(self.raw(), idx)?;
        Ok(u32::try_from(id).ok())
    }

    /// An entry's description. Mirrors `MolCatalog.GetEntryDescription`.
    /// `Error::Index` for an index past the entries.
    pub fn get_entry_description(&self, idx: u32) -> Result<String> {
        to_string(ffi::mol_catalog_get_entry_description(self.raw(), idx)?)
    }

    /// The indices of the entries an entry has edges to. Mirrors
    /// `MolCatalog.GetEntryDownIds`. `Error::Index` for an index past the
    /// entries; Python lets one past the end through, and RDKit reads past
    /// its graph.
    pub fn get_entry_down_ids(&self, idx: u32) -> Result<Vec<u32>> {
        Ok(ffi::mol_catalog_get_entry_down_ids(self.raw(), idx)?)
    }

    /// Adds a copy of an entry with the next bit id, returning its index.
    /// Mirrors `MolCatalog.AddEntry`.
    pub fn add_entry(&mut self, entry: &MolCatalogEntry) -> u32 {
        expect_ok(ffi::mol_catalog_add_entry(self.ptr.pin_mut(), entry.raw()))
    }

    /// Adds an edge from one entry to another, unless it exists. Mirrors
    /// `MolCatalog.AddEdge`. `Error::Index` for an index past the entries.
    pub fn add_edge(&mut self, id1: u32, id2: u32) -> Result<()> {
        Ok(ffi::mol_catalog_add_edge(self.ptr.pin_mut(), id1, id2)?)
    }
}

impl fmt::Debug for MolCatalog {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MolCatalog")
            .field("num_entries", &self.get_num_entries())
            .field("fp_length", &self.get_fp_length())
            .finish()
    }
}

/// An empty catalog. Mirrors `rdMolCatalog.CreateMolCatalog`.
pub fn create_mol_catalog() -> MolCatalog {
    MolCatalog {
        ptr: expect_ok(ffi::create_mol_catalog()),
    }
}

/// A molecule with a description and an order, for a `MolCatalog`. Mirrors
/// `rdMolCatalog.MolCatalogEntry`.
///
/// It is `Send` but not `Sync`, as the molecule it holds is not
/// (CONVENTIONS.md S1).
pub struct MolCatalogEntry {
    /// Never null.
    ptr: UniquePtr<ffi::MolCatalogEntry>,
}

impl MolCatalogEntry {
    fn raw(&self) -> &ffi::MolCatalogEntry {
        self.ptr.as_ref().unwrap()
    }

    /// An entry from `serialize`'s bytes. Mirrors
    /// `MolCatalogEntry(pickle)`. `Error::Value` for a description whose
    /// length RDKit would read out of bounds; Python's reads it. A
    /// truncated pickle is `Error::Runtime`.
    ///
    /// # Safety
    ///
    /// `pickle` must be bytes `serialize` (this crate's or RDKit's)
    /// returned, unmodified: it holds a molecule pickle, which RDKit reads
    /// as `Mol::from_pkl_string` does.
    pub unsafe fn from_pickle(pickle: &[u8]) -> Result<MolCatalogEntry> {
        Ok(MolCatalogEntry {
            ptr: ffi::mol_catalog_entry_from_pickle(pickle)?,
        })
    }

    /// Mirrors `MolCatalogEntry.GetDescription`. `Error::Value` for one
    /// from a pickle that is not UTF-8.
    pub fn get_description(&self) -> Result<String> {
        to_string(expect_ok(ffi::mol_catalog_entry_get_description(
            self.raw(),
        )))
    }

    /// Mirrors `MolCatalogEntry.SetDescription`.
    pub fn set_description(&mut self, val: &str) {
        expect_ok(ffi::mol_catalog_entry_set_description(
            self.ptr.pin_mut(),
            val,
        ))
    }

    /// A copy of the entry's molecule, or `None` before `set_mol` (where
    /// Python's dereferences null). Mirrors `MolCatalogEntry.GetMol`, which
    /// returns the molecule itself.
    pub fn get_mol(&self) -> Option<Mol> {
        let mol = expect_ok(ffi::mol_catalog_entry_get_mol(self.raw()));
        (!mol.is_null()).then(|| Mol::from_ffi(mol))
    }

    /// Gives the entry a copy of a molecule. Mirrors `MolCatalogEntry.SetMol`.
    pub fn set_mol(&mut self, mol: &Mol) {
        expect_ok(ffi::mol_catalog_entry_set_mol(
            self.ptr.pin_mut(),
            mol.as_ffi(),
        ))
    }

    /// Mirrors `MolCatalogEntry.GetOrder`.
    pub fn get_order(&self) -> u32 {
        expect_ok(ffi::mol_catalog_entry_get_order(self.raw()))
    }

    /// Mirrors `MolCatalogEntry.SetOrder`.
    pub fn set_order(&mut self, order: u32) {
        expect_ok(ffi::mol_catalog_entry_set_order(self.ptr.pin_mut(), order))
    }

    /// The entry as bytes. Mirrors `MolCatalogEntry.Serialize`.
    /// `Error::Runtime` without a molecule, as in Python, and `Error::Value`
    /// for a description of 2^31 - 1 bytes or more, whose length RDKit
    /// writes in 32 bits, or a molecule `Mol::to_binary` refuses.
    pub fn serialize(&self) -> Result<Vec<u8>> {
        Ok(ffi::mol_catalog_entry_serialize(self.raw())?)
    }
}

/// An entry without a molecule. Mirrors `MolCatalogEntry()`.
impl Default for MolCatalogEntry {
    fn default() -> MolCatalogEntry {
        MolCatalogEntry {
            ptr: expect_ok(ffi::mol_catalog_entry_new()),
        }
    }
}

impl fmt::Debug for MolCatalogEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MolCatalogEntry")
            .field("description", &self.get_description())
            .field("order", &self.get_order())
            .finish()
    }
}
