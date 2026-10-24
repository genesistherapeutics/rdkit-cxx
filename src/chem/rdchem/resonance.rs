// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::Mol;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::resonance::ffi;
use crate::ffi::substruct::ffi as substruct;

/// Which resonance structures to allow, as bits. Mirrors
/// `rdchem.ResonanceFlags` (CONVENTIONS.md F1).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct ResonanceFlags(u32);

macro_rules! resonance_flags {
    ($($name:ident),* $(,)?) => {
        impl ResonanceFlags {
            $(
                #[doc = concat!("Mirrors `ResonanceFlags.", stringify!($name), "`.")]
                pub const $name: ResonanceFlags = ResonanceFlags(ffi::ResonanceFlags::$name.repr);
            )*

            const NAMES: &[(&str, ResonanceFlags)] =
                &[$((stringify!($name), ResonanceFlags::$name)),*];
        }
    };
}

resonance_flags!(
    ALLOW_INCOMPLETE_OCTETS,
    ALLOW_CHARGE_SEPARATION,
    KEKULE_ALL,
    UNCONSTRAINED_CATIONS,
    UNCONSTRAINED_ANIONS,
);

crate::flags::flag_ops!(ResonanceFlags);

/// Parameters for `ResonanceMolSupplier::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ResonanceMolSupplierParams {
    pub flags: ResonanceFlags,
    /// The most structures to enumerate; RDKit caps it at 1,000,000.
    pub max_structs: u32,
}

impl Default for ResonanceMolSupplierParams {
    fn default() -> ResonanceMolSupplierParams {
        ResonanceMolSupplierParams {
            flags: ResonanceFlags::default(),
            max_structs: 1000,
        }
    }
}

/// Parameters for `ResonanceMolSupplier::get_substruct_match_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetSubstructMatchParams {
    pub use_chirality: bool,
    pub use_query_query_matches: bool,
}

/// Parameters for `ResonanceMolSupplier::get_substruct_matches_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetSubstructMatchesParams {
    pub uniquify: bool,
    pub use_chirality: bool,
    pub use_query_query_matches: bool,
    pub max_matches: u32,
    pub num_threads: i32,
}

impl Default for GetSubstructMatchesParams {
    fn default() -> GetSubstructMatchesParams {
        GetSubstructMatchesParams {
            uniquify: false,
            use_chirality: false,
            use_query_query_matches: false,
            max_matches: 1000,
            num_threads: 1,
        }
    }
}

/// A molecule's resonance structures. Mirrors `rdchem.ResonanceMolSupplier`.
///
/// `structures` iterates over them from the first, as Python's `for` loop
/// does. Enumeration happens on first use, as in RDKit. The progress
/// callbacks are not bound.
///
/// The supplier is not itself an iterator: `Iterator::enumerate` would take
/// the place of `enumerate`.
pub struct ResonanceMolSupplier {
    /// Never null.
    ptr: UniquePtr<ffi::ResonanceMolSupplier>,
}

impl ResonanceMolSupplier {
    fn raw(&self) -> &ffi::ResonanceMolSupplier {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> std::pin::Pin<&mut ffi::ResonanceMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier for a copy of `mol`, which it kekulizes; a molecule that
    /// cannot be kekulized is `Error::Kekulize`. Mirrors
    /// `ResonanceMolSupplier(mol)`.
    pub fn new(mol: &Mol) -> Result<ResonanceMolSupplier> {
        ResonanceMolSupplier::new_with_params(mol, &ResonanceMolSupplierParams::default())
    }

    /// Mirrors `ResonanceMolSupplier(mol, flags, maxStructs)`.
    pub fn new_with_params(
        mol: &Mol,
        params: &ResonanceMolSupplierParams,
    ) -> Result<ResonanceMolSupplier> {
        let ptr =
            ffi::resonance_mol_supplier_new(mol.as_ffi(), params.flags.bits(), params.max_structs)?;
        assert!(!ptr.is_null(), "RDKit returned a null supplier");
        Ok(ResonanceMolSupplier { ptr })
    }

    /// The number of conjugated groups. Mirrors
    /// `ResonanceMolSupplier.GetNumConjGrps`.
    pub fn get_num_conj_grps(&self) -> u32 {
        expect_ok(ffi::resonance_mol_supplier_get_num_conj_grps(self.raw()))
    }

    /// The bond's conjugated group, or -1. Mirrors
    /// `ResonanceMolSupplier.GetBondConjGrpIdx`.
    pub fn get_bond_conj_grp_idx(&self, bi: u32) -> Result<i32> {
        Ok(ffi::resonance_mol_supplier_get_bond_conj_grp_idx(
            self.raw(),
            bi,
        )?)
    }

    /// The atom's conjugated group, or -1. Mirrors
    /// `ResonanceMolSupplier.GetAtomConjGrpIdx`.
    pub fn get_atom_conj_grp_idx(&self, ai: u32) -> Result<i32> {
        Ok(ffi::resonance_mol_supplier_get_atom_conj_grp_idx(
            self.raw(),
            ai,
        )?)
    }

    /// Threads for the enumeration: 0 for every core, and negative numbers
    /// count back from it. Mirrors `ResonanceMolSupplier.SetNumThreads`.
    pub fn set_num_threads(&mut self, num_threads: i32) {
        expect_ok(ffi::resonance_mol_supplier_set_num_threads(
            self.raw_mut(),
            num_threads,
        ))
    }

    /// Whether a progress callback stopped the enumeration; always false
    /// here, since the callbacks are not bound. Mirrors
    /// `ResonanceMolSupplier.WasCanceled`.
    pub fn was_canceled(&self) -> bool {
        expect_ok(ffi::resonance_mol_supplier_was_canceled(self.raw()))
    }

    /// Enumerates the structures now, rather than on first use. Mirrors
    /// `ResonanceMolSupplier.Enumerate`.
    pub fn enumerate(&mut self) {
        expect_ok(ffi::resonance_mol_supplier_enumerate(self.raw_mut()))
    }

    /// Mirrors `ResonanceMolSupplier.GetIsEnumerated`.
    pub fn get_is_enumerated(&self) -> bool {
        expect_ok(ffi::resonance_mol_supplier_get_is_enumerated(self.raw()))
    }

    /// The number of structures. Mirrors `len(ResonanceMolSupplier)`.
    pub fn len(&mut self) -> u32 {
        expect_ok(ffi::resonance_mol_supplier_length(self.raw_mut()))
    }

    /// Whether there are no structures; there is always at least one, for a
    /// molecule with atoms.
    pub fn is_empty(&mut self) -> bool {
        self.len() == 0
    }

    /// The structures, from the first. Mirrors `iter(ResonanceMolSupplier)`,
    /// which resets the supplier.
    pub fn structures(&mut self) -> ResonanceStructures<'_> {
        self.reset();
        self.structures_from_here()
    }

    /// The structures from where the supplier stands, as Python's `next`
    /// continues.
    pub fn structures_from_here(&mut self) -> ResonanceStructures<'_> {
        ResonanceStructures { suppl: self }
    }

    /// Returns to the first structure. Mirrors `ResonanceMolSupplier.reset`.
    pub fn reset(&mut self) {
        expect_ok(ffi::resonance_mol_supplier_reset(self.raw_mut()))
    }

    /// Whether the iterator is done. Mirrors `ResonanceMolSupplier.atEnd`.
    pub fn at_end(&mut self) -> bool {
        expect_ok(ffi::resonance_mol_supplier_at_end(self.raw_mut()))
    }

    /// The structure at `idx`. Mirrors `ResonanceMolSupplier[idx]`; negative
    /// indices from the end are not bound.
    pub fn get(&mut self, idx: u32) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::resonance_mol_supplier_get_item(
            self.raw_mut(),
            idx,
        )?))
    }

    /// The first match of `query` in any structure, or empty. Mirrors
    /// `ResonanceMolSupplier.GetSubstructMatch`.
    pub fn get_substruct_match(&mut self, query: &Mol) -> Result<Vec<u32>> {
        self.get_substruct_match_with_params(query, &GetSubstructMatchParams::default())
    }

    /// Mirrors `ResonanceMolSupplier.GetSubstructMatch(query, useChirality,
    /// useQueryQueryMatches)`.
    pub fn get_substruct_match_with_params(
        &mut self,
        query: &Mol,
        params: &GetSubstructMatchParams,
    ) -> Result<Vec<u32>> {
        Ok(substruct::resonance_mol_supplier_get_substruct_match(
            self.raw_mut(),
            query.as_ffi(),
            params.use_chirality,
            params.use_query_query_matches,
        )?)
    }

    /// The matches of `query` across the structures. Mirrors
    /// `ResonanceMolSupplier.GetSubstructMatches`.
    pub fn get_substruct_matches(&mut self, query: &Mol) -> Result<Vec<Vec<u32>>> {
        self.get_substruct_matches_with_params(query, &GetSubstructMatchesParams::default())
    }

    /// Mirrors `ResonanceMolSupplier.GetSubstructMatches(query, uniquify,
    /// useChirality, useQueryQueryMatches, maxMatches, numThreads)`.
    pub fn get_substruct_matches_with_params(
        &mut self,
        query: &Mol,
        params: &GetSubstructMatchesParams,
    ) -> Result<Vec<Vec<u32>>> {
        Ok(IndexList::into_nested(
            substruct::resonance_mol_supplier_get_substruct_matches(
                self.raw_mut(),
                query.as_ffi(),
                params.uniquify,
                params.use_chirality,
                params.use_query_query_matches,
                params.max_matches,
                params.num_threads,
            )?,
        ))
    }
}

/// The structures of a `ResonanceMolSupplier`, from where it stands. A
/// structure RDKit fails to build is an `Err`, where Python yields `None`, and
/// the iteration goes on.
#[derive(Debug)]
pub struct ResonanceStructures<'s> {
    suppl: &'s mut ResonanceMolSupplier,
}

/// The next structure. Mirrors `next(ResonanceMolSupplier)`.
impl Iterator for ResonanceStructures<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.suppl.at_end() {
            return None;
        }
        Some(
            ffi::resonance_mol_supplier_next(self.suppl.raw_mut())
                .map(Mol::from_ffi)
                .map_err(Into::into),
        )
    }
}

impl fmt::Debug for ResonanceMolSupplier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ResonanceMolSupplier")
            .field("num_conj_grps", &self.get_num_conj_grps())
            .field("is_enumerated", &self.get_is_enumerated())
            .finish()
    }
}
