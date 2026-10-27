// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdfragcatalog`: catalogs of the fragments of
//! molecules, and fingerprints of which fragments a molecule has.
//!
//! A `FragCatParams` reads functional groups from a file: lines of a name, a
//! tab, and a SMARTS whose first atom is where the group attaches, with
//! comment lines starting `//`. A `FragCatGenerator` adds the fragments of
//! molecules, with their functional groups removed and noted, to a
//! `FragCatalog`; fragments within the parameters' lengths get bits, which a
//! `FragFPGenerator` sets for the fragments a molecule has.
//!
//! RDKit reads the pickles' lengths and indices unchecked, so the wrappers
//! check them first (CONVENTIONS.md U2), as `rd_mol_catalog` does.

use std::fmt;
use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::data_structs::ExplicitBitVect;
use crate::error::{expect_ok, to_string};
use crate::ffi::frag_catalog::ffi;
use crate::util::path_bytes;

/// The fragment lengths, tolerance, and functional groups of a fragment
/// catalog. Mirrors `rdfragcatalog.FragCatParams`.
///
/// It is `Send` but not `Sync`, as the molecules it holds are not
/// (CONVENTIONS.md S1).
pub struct FragCatParams {
    /// Never null.
    ptr: UniquePtr<ffi::FragCatParams>,
}

impl FragCatParams {
    fn raw(&self) -> &ffi::FragCatParams {
        self.ptr.as_ref().unwrap()
    }

    /// Parameters for fragments of `l_len` to `u_len` bonds, with the
    /// functional groups in a file and the default tolerance, 1e-8. Mirrors
    /// `FragCatParams(lLen, uLen, fgroupFilename)`. `Error::Os` for a file
    /// that cannot be read; `Error::Runtime` for `l_len > u_len`, a line
    /// without a SMARTS, or a SMARTS that does not parse, as in Python.
    /// The generators set aside a list for every length up to `u_len`, so a
    /// very large one exhausts memory, as in Python.
    pub fn new<P: AsRef<Path>>(
        l_len: u32,
        u_len: u32,
        fgroup_filename: P,
    ) -> Result<FragCatParams> {
        FragCatParams::new_with_tol(l_len, u_len, fgroup_filename, 1e-8)
    }

    /// `tol` is the tolerance within which fragments' discriminators are
    /// equal. Mirrors `FragCatParams(lLen, uLen, fgroupFilename, tol)`.
    pub fn new_with_tol<P: AsRef<Path>>(
        l_len: u32,
        u_len: u32,
        fgroup_filename: P,
        tol: f64,
    ) -> Result<FragCatParams> {
        Ok(FragCatParams {
            ptr: ffi::frag_cat_params_new(l_len, u_len, path_bytes(fgroup_filename.as_ref()), tol)?,
        })
    }

    /// Mirrors `FragCatParams.GetTypeString`.
    pub fn get_type_string(&self) -> String {
        to_string(expect_ok(ffi::frag_cat_params_get_type_string(self.raw())))
            .expect("the type string is RDKit's")
    }

    /// The most bonds of a fragment with a bit. Mirrors
    /// `FragCatParams.GetUpperFragLength`.
    pub fn get_upper_frag_length(&self) -> u32 {
        expect_ok(ffi::frag_cat_params_get_upper_frag_length(self.raw()))
    }

    /// The fewest bonds of a fragment with a bit. Mirrors
    /// `FragCatParams.GetLowerFragLength`.
    pub fn get_lower_frag_length(&self) -> u32 {
        expect_ok(ffi::frag_cat_params_get_lower_frag_length(self.raw()))
    }

    /// Mirrors `FragCatParams.GetTolerance`.
    pub fn get_tolerance(&self) -> f64 {
        expect_ok(ffi::frag_cat_params_get_tolerance(self.raw()))
    }

    /// Mirrors `FragCatParams.GetNumFuncGroups`.
    pub fn get_num_func_groups(&self) -> u32 {
        expect_ok(ffi::frag_cat_params_get_num_func_groups(self.raw()))
    }

    /// A copy of a functional group's query molecule, with its name in
    /// `_Name` and its SMARTS in `_fragSMARTS`. Mirrors
    /// `FragCatParams.GetFuncGroup`, which returns the molecule itself.
    /// `Error::Index` for an index past the groups.
    pub fn get_func_group(&self, fid: u32) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::frag_cat_params_get_func_group(
            self.raw(),
            fid,
        )?))
    }

    /// The parameters as text: the lengths and tolerance (to six
    /// significant digits), then each group's name and SMARTS. Mirrors
    /// `FragCatParams.Serialize`.
    pub fn serialize(&self) -> Vec<u8> {
        expect_ok(ffi::frag_cat_params_serialize(self.raw()))
    }
}

impl fmt::Debug for FragCatParams {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FragCatParams")
            .field("lower_frag_length", &self.get_lower_frag_length())
            .field("upper_frag_length", &self.get_upper_frag_length())
            .field("tolerance", &self.get_tolerance())
            .field("num_func_groups", &self.get_num_func_groups())
            .finish()
    }
}

/// A hierarchical catalog of molecular fragments. Mirrors
/// `rdfragcatalog.FragCatalog`.
///
/// It is `Send` but not `Sync`: RDKit caches each fragment's discriminators
/// in its properties when it first reads them, through `&self` methods.
pub struct FragCatalog {
    /// Never null.
    ptr: UniquePtr<ffi::FragCatalog>,
}

impl FragCatalog {
    fn raw(&self) -> &ffi::FragCatalog {
        self.ptr.as_ref().unwrap()
    }

    /// An empty catalog with a copy of the parameters. Mirrors
    /// `FragCatalog(params)`.
    pub fn new(params: &FragCatParams) -> FragCatalog {
        FragCatalog {
            ptr: expect_ok(ffi::frag_catalog_from_params(params.raw())),
        }
    }

    /// A catalog from `serialize`'s bytes. Mirrors `FragCatalog(pickle)`.
    /// `Error::Value` for a pickle whose lengths or indices RDKit would
    /// read out of bounds; Python's reads them. A truncated pickle is
    /// `Error::Runtime`, as is a catalog with no functional groups, whose
    /// pickles RDKit cannot read back.
    ///
    /// # Safety
    ///
    /// `pickle` must be bytes `serialize` (this crate's or RDKit's)
    /// returned, unmodified: each entry holds a molecule pickle, which RDKit
    /// reads as `Mol::from_pkl_string` does.
    pub unsafe fn from_pickle(pickle: &[u8]) -> Result<FragCatalog> {
        Ok(FragCatalog {
            ptr: ffi::frag_catalog_from_pickle(pickle)?,
        })
    }

    /// The number of fragments. Mirrors `FragCatalog.GetNumEntries`.
    pub fn get_num_entries(&self) -> u32 {
        expect_ok(ffi::frag_catalog_get_num_entries(self.raw()))
    }

    /// The number of bits fragments have. Mirrors
    /// `FragCatalog.GetFPLength`.
    pub fn get_fp_length(&self) -> u32 {
        expect_ok(ffi::frag_catalog_get_fp_length(self.raw()))
    }

    /// A copy of the parameters, or `None` for a catalog made by `default`.
    /// Mirrors `FragCatalog.GetCatalogParams`, which returns them
    /// themselves.
    pub fn get_catalog_params(&self) -> Option<FragCatParams> {
        let ptr = expect_ok(ffi::frag_catalog_get_catalog_params(self.raw()));
        (!ptr.is_null()).then_some(FragCatParams { ptr })
    }

    /// The catalog as bytes. Mirrors `FragCatalog.Serialize`.
    /// `Error::Runtime` for a catalog without parameters, as in Python.
    pub fn serialize(&self) -> Result<Vec<u8>> {
        Ok(ffi::frag_catalog_serialize(self.raw())?)
    }

    /// The SMILES of the fragment with a bit, its functional groups in angle
    /// brackets. Mirrors `FragCatalog.GetBitDescription`. `Error::Index` for
    /// a bit past `get_fp_length`; `Error::Value` when no fragment at or
    /// after the bit's index has it, as a pickle can arrange, where Python's
    /// dereferences null (and likewise for the other bit lookups).
    pub fn get_bit_description(&self, idx: u32) -> Result<String> {
        to_string(ffi::frag_catalog_get_bit_description(self.raw(), idx)?)
    }

    /// The number of bonds of the fragment with a bit. Mirrors
    /// `FragCatalog.GetBitOrder`.
    pub fn get_bit_order(&self, idx: u32) -> Result<u32> {
        Ok(ffi::frag_catalog_get_bit_order(self.raw(), idx)?)
    }

    /// The functional groups of the fragment with a bit, by atom. Mirrors
    /// `FragCatalog.GetBitFuncGroupIds`.
    pub fn get_bit_func_group_ids(&self, idx: u32) -> Result<Vec<i32>> {
        Ok(ffi::frag_catalog_get_bit_func_group_ids(self.raw(), idx)?)
    }

    /// The index of the fragment with a bit, or `None` when no fragment at
    /// or after the bit's index has it (where Python returns 2^32 - 1).
    /// Mirrors `FragCatalog.GetBitEntryId`.
    pub fn get_bit_entry_id(&self, idx: u32) -> Result<Option<u32>> {
        let id = ffi::frag_catalog_get_bit_entry_id(self.raw(), idx)?;
        Ok(u32::try_from(id).ok())
    }

    /// The discriminators of the fragment with a bit, which tell fragments
    /// apart. Mirrors `FragCatalog.GetBitDiscrims`.
    pub fn get_bit_discrims(&self, idx: u32) -> Result<[f64; 3]> {
        let values = ffi::frag_catalog_get_bit_discrims(self.raw(), idx)?;
        Ok([values[0], values[1], values[2]])
    }

    /// A fragment's bit, or `None` for a fragment without one, whose id is
    /// negative (Python returns it modulo 2^32). Mirrors
    /// `FragCatalog.GetEntryBitId`. `Error::Index` for an index past the
    /// fragments (and likewise for the other entry lookups).
    pub fn get_entry_bit_id(&self, idx: u32) -> Result<Option<u32>> {
        let id = ffi::frag_catalog_get_entry_bit_id(self.raw(), idx)?;
        Ok(u32::try_from(id).ok())
    }

    /// Mirrors `FragCatalog.GetEntryDescription`.
    pub fn get_entry_description(&self, idx: u32) -> Result<String> {
        to_string(ffi::frag_catalog_get_entry_description(self.raw(), idx)?)
    }

    /// Mirrors `FragCatalog.GetEntryOrder`.
    pub fn get_entry_order(&self, idx: u32) -> Result<u32> {
        Ok(ffi::frag_catalog_get_entry_order(self.raw(), idx)?)
    }

    /// Mirrors `FragCatalog.GetEntryFuncGroupIds`.
    pub fn get_entry_func_group_ids(&self, idx: u32) -> Result<Vec<i32>> {
        Ok(ffi::frag_catalog_get_entry_func_group_ids(self.raw(), idx)?)
    }

    /// The indices of the fragments one bond larger that contain a
    /// fragment. Mirrors `FragCatalog.GetEntryDownIds`. Python lets an index
    /// one past the end through, and RDKit reads past its graph.
    pub fn get_entry_down_ids(&self, idx: u32) -> Result<Vec<u32>> {
        Ok(ffi::frag_catalog_get_entry_down_ids(self.raw(), idx)?)
    }
}

/// A catalog without parameters, to which no fragments can be added; a
/// Rust addition (Python's constructors take parameters or a pickle).
impl Default for FragCatalog {
    fn default() -> FragCatalog {
        FragCatalog {
            ptr: expect_ok(ffi::frag_catalog_new()),
        }
    }
}

impl fmt::Debug for FragCatalog {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FragCatalog")
            .field("num_entries", &self.get_num_entries())
            .field("fp_length", &self.get_fp_length())
            .finish()
    }
}

/// Adds molecules' fragments to catalogs. Mirrors
/// `rdfragcatalog.FragCatGenerator`.
#[derive(Clone, Copy, Debug, Default)]
pub struct FragCatGenerator;

impl FragCatGenerator {
    /// Adds the fragments of a molecule, up to the catalog's upper length,
    /// that the catalog lacks, returning how many it added. Mirrors
    /// `FragCatGenerator.AddFragsFromMol`. `Error::Value` for a catalog
    /// without parameters, where Python's dereferences null, and for a
    /// molecule of more than 1000 bonds when the upper length is 2 or more:
    /// RDKit then indexes its table of 1000 primes by bond index.
    /// `Error::Value` for a molecule whose fragments a pickle may not hold: a
    /// query of more than 255 children, or an atom or bond of more than 65535
    /// properties (`Mol::to_binary`). RDKit 2024.03's fragments keep no
    /// substance or stereo groups, so those are not checked.
    pub fn add_frags_from_mol(&self, mol: &Mol, fcat: &mut FragCatalog) -> Result<u32> {
        Ok(ffi::frag_cat_generator_add_frags_from_mol(
            mol.as_ffi(),
            fcat.ptr.pin_mut(),
        )?)
    }
}

/// Fingerprints molecules by the fragments of a catalog they have. Mirrors
/// `rdfragcatalog.FragFPGenerator`.
#[derive(Clone, Copy, Debug, Default)]
pub struct FragFPGenerator;

impl FragFPGenerator {
    /// A fingerprint of `fcat.get_fp_length()` bits, with the bit of each of
    /// the catalog's fragments the molecule has. Mirrors
    /// `FragFPGenerator.GetFPForMol`. `Error::Value` as for
    /// `FragCatGenerator::add_frags_from_mol`, and for a catalog from a
    /// pickle with more than 2^31 - 1 bits (U3) or with a fragment whose bit
    /// is past them; `Error::Runtime` for a catalog without fragments of one
    /// bond, as in Python.
    pub fn get_fp_for_mol(&self, mol: &Mol, fcat: &FragCatalog) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(
            ffi::frag_fp_generator_get_fp_for_mol(mol.as_ffi(), fcat.raw())?,
        ))
    }
}
