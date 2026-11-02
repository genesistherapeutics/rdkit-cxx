// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdSubstructLibrary`: libraries of molecules searched
//! by substructure, screened by pattern fingerprints.
//!
//! A library owns its holders, which Python's shares with the holder
//! objects it was made from; `get_mol_holder_mut` and `downcast_mut` reach
//! them.
//!
//! Every search runs on the calling thread, whatever its `num_threads`,
//! which is accepted for parity. RDKit's search threads each construct a
//! `ControlCHandler`, which the wrappers' lock cannot order
//! (CONVENTIONS.md U4), and an exception in one of them (a stored SMILES
//! that does not sanitize, a library with fewer fingerprints than
//! molecules) unwinds into freed memory while the others still run; Python
//! aborts with heap corruption, or segfaults. The results are RDKit's for
//! any number of threads, except that a `MolBundle`'s matches come member
//! by member, as RDKit's single thread finds them, where its threads sort
//! them and can cut earlier members' matches short.
//!
//! RDKit's tautomer matcher reads the matches of a recursive query in a
//! tautomer query's template as the template's own, past their end, so a
//! search for a `TautomerQuery`, or an `ExtendedQueryMol` of tautomers,
//! whose template has a recursive SMARTS (`$(...)`) is `Error::Value`.
//!
//! ```no_run
//! use rdkit_cxx::chem::rd_substruct_library::{
//!     CachedSmilesMolHolder, PatternHolder, SubstructLibrary,
//! };
//! use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles};
//!
//! let mut library = SubstructLibrary::new(
//!     CachedSmilesMolHolder::default().into(),
//!     Some(PatternHolder::default().into()),
//!     None,
//! );
//! for smiles in ["CCO", "c1ccccc1O", "CCN"] {
//!     library.add_mol(&mol_from_smiles(smiles)?)?;
//! }
//! let hydroxyls = library.get_matches(&mol_from_smarts("[OX2H]")?)?;
//! assert_eq!(hydroxyls, [0, 1]);
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

mod holders;

use std::fmt;
use std::io::{Read, Write};

use cxx::UniquePtr;

pub use holders::{
    CachedMolHolder, CachedSmilesMolHolder, CachedTrustedSmilesMolHolder, FPHolderBase,
    FPHolderType, KeyFromPropHolder, KeyHolderBase, KeyHolderType, MolHolder, MolHolderBase,
    MolHolderType, PatternHolder, TautomerPatternHolder,
};

use crate::chem::rd_generalized_substruct::ExtendedQueryMol;
use crate::chem::rd_tautomer_query::TautomerQuery;
use crate::chem::rdchem::{Mol, MolBundle, SubstructMatchParameters};
use crate::error::expect_ok;
use crate::ffi::substruct_library::ffi;
use crate::{Error, Result};

/// The keyword arguments of `SubstructLibrary.GetMatches`, whose
/// `parameters` overload takes `parameters`' fields (A3). `Default` is the
/// keyword overload's: `SubstructMatchParameters`' defaults with
/// `use_chirality` set.
#[derive(Clone, Debug, PartialEq)]
pub struct GetMatchesParams {
    /// How molecules are matched.
    pub parameters: SubstructMatchParameters,
    /// Accepted for parity: searches run on the calling thread (see the
    /// module documentation).
    pub num_threads: i32,
    /// The most matches to return; a negative number returns all.
    pub max_results: i32,
}

impl Default for GetMatchesParams {
    fn default() -> GetMatchesParams {
        GetMatchesParams {
            parameters: keyword_parameters(),
            num_threads: -1,
            max_results: 1000,
        }
    }
}

/// The keyword arguments of `SubstructLibrary.CountMatches`, whose
/// `parameters` overload takes `parameters`' fields (A3). `Default` is the
/// keyword overload's: `SubstructMatchParameters`' defaults with
/// `use_chirality` set.
#[derive(Clone, Debug, PartialEq)]
pub struct CountMatchesParams {
    /// How molecules are matched.
    pub parameters: SubstructMatchParameters,
    /// Accepted for parity: searches run on the calling thread (see the
    /// module documentation).
    pub num_threads: i32,
}

impl Default for CountMatchesParams {
    fn default() -> CountMatchesParams {
        CountMatchesParams {
            parameters: keyword_parameters(),
            num_threads: -1,
        }
    }
}

/// The keyword arguments of `SubstructLibrary.HasMatch`: those of
/// `CountMatches`.
pub type HasMatchParams = CountMatchesParams;

// The keyword overloads set recursionPossible, useChirality, and
// useQueryQueryMatches on default parameters, with Python's defaults of
// true, true, and false.
fn keyword_parameters() -> SubstructMatchParameters {
    SubstructMatchParameters {
        recursion_possible: true,
        use_chirality: true,
        use_query_query_matches: false,
        ..SubstructMatchParameters::default()
    }
}

/// A library of molecules to search by substructure: a molecule holder,
/// with pattern fingerprints to screen them and keys to name them if it
/// has those holders, and the order to search them in. Mirrors
/// `rdSubstructLibrary.SubstructLibrary`.
///
/// `SubstructLibrary` is `Send` but not `Sync`: a search writes the ring
/// information of a `MolHolder`'s molecules (CONVENTIONS.md S1).
pub struct SubstructLibrary {
    molecules: MolHolderBase,
    fingerprints: Option<FPHolderBase>,
    keys: Option<KeyHolderBase>,
    search_order: Vec<u32>,
}

/// A query, borrowed as C++ for one call.
#[derive(Clone, Copy)]
enum RawQuery<'a> {
    Mol(&'a crate::ffi::mol::ffi::RWMol),
    Bundle(&'a crate::ffi::mol_bundle::ffi::MolBundle),
    Tautomer(&'a crate::ffi::tautomer_query::ffi::TautomerQuery),
    Extended(&'a crate::ffi::generalized_substruct::ffi::ExtendedQueryMol),
}

/// Whether a search has a range, and the range.
#[derive(Clone, Copy)]
struct Range {
    ranged: bool,
    start_idx: u32,
    end_idx: u32,
}

type Data = crate::ffi::substruct::ffi::SubstructMatchParametersData;
type Raw<T> = std::result::Result<T, cxx::Exception>;

impl RawQuery<'_> {
    fn get_matches(
        self,
        view: &ffi::SubstructLibrary,
        r: Range,
        p: &Data,
        max: i32,
    ) -> Raw<Vec<u32>> {
        let (ranged, start, end) = (r.ranged, r.start_idx, r.end_idx);
        match self {
            RawQuery::Mol(q) => {
                ffi::substruct_library_get_matches(view, q, ranged, start, end, p, max)
            }
            RawQuery::Bundle(q) => {
                ffi::substruct_library_get_matches_mol_bundle(view, q, ranged, start, end, p, max)
            }
            RawQuery::Tautomer(q) => ffi::substruct_library_get_matches_tautomer_query(
                view, q, ranged, start, end, p, max,
            ),
            RawQuery::Extended(q) => ffi::substruct_library_get_matches_extended_query_mol(
                view, q, ranged, start, end, p, max,
            ),
        }
    }

    fn count_matches(self, view: &ffi::SubstructLibrary, r: Range, p: &Data) -> Raw<u32> {
        let (ranged, start, end) = (r.ranged, r.start_idx, r.end_idx);
        match self {
            RawQuery::Mol(q) => {
                ffi::substruct_library_count_matches(view, q, ranged, start, end, p)
            }
            RawQuery::Bundle(q) => {
                ffi::substruct_library_count_matches_mol_bundle(view, q, ranged, start, end, p)
            }
            RawQuery::Tautomer(q) => {
                ffi::substruct_library_count_matches_tautomer_query(view, q, ranged, start, end, p)
            }
            RawQuery::Extended(q) => ffi::substruct_library_count_matches_extended_query_mol(
                view, q, ranged, start, end, p,
            ),
        }
    }

    fn has_match(self, view: &ffi::SubstructLibrary, r: Range, p: &Data) -> Raw<bool> {
        let (ranged, start, end) = (r.ranged, r.start_idx, r.end_idx);
        match self {
            RawQuery::Mol(q) => ffi::substruct_library_has_match(view, q, ranged, start, end, p),
            RawQuery::Bundle(q) => {
                ffi::substruct_library_has_match_mol_bundle(view, q, ranged, start, end, p)
            }
            RawQuery::Tautomer(q) => {
                ffi::substruct_library_has_match_tautomer_query(view, q, ranged, start, end, p)
            }
            RawQuery::Extended(q) => {
                ffi::substruct_library_has_match_extended_query_mol(view, q, ranged, start, end, p)
            }
        }
    }
}

/// A query a library searches for (Python's overloads by query type, A5).
trait LibraryQuery {
    fn with_raw<T>(&self, f: impl FnOnce(RawQuery<'_>) -> T) -> T;
}

impl LibraryQuery for Mol {
    fn with_raw<T>(&self, f: impl FnOnce(RawQuery<'_>) -> T) -> T {
        f(RawQuery::Mol(self.as_ffi()))
    }
}

impl LibraryQuery for MolBundle {
    fn with_raw<T>(&self, f: impl FnOnce(RawQuery<'_>) -> T) -> T {
        self.with_ffi(|bundle| f(RawQuery::Bundle(bundle)))
    }
}

impl LibraryQuery for TautomerQuery {
    fn with_raw<T>(&self, f: impl FnOnce(RawQuery<'_>) -> T) -> T {
        f(RawQuery::Tautomer(self.as_ffi()))
    }
}

impl LibraryQuery for ExtendedQueryMol {
    fn with_raw<T>(&self, f: impl FnOnce(RawQuery<'_>) -> T) -> T {
        f(RawQuery::Extended(self.as_ffi()))
    }
}

/// Defines the twelve searches for one query type: `GetMatches`,
/// `CountMatches`, and `HasMatch`, each over the whole library or a range,
/// with Python's defaults or parameters.
macro_rules! searches {
    (
        $type:ty, $what:literal;
        $gm:ident, $gmp:ident, $gmr:ident, $gmrp:ident;
        $cm:ident, $cmp:ident, $cmr:ident, $cmrp:ident;
        $hm:ident, $hmp:ident, $hmr:ident, $hmrp:ident
    ) => {
        #[doc = concat!("The molecules that match ", $what, ", at most 1000. Mirrors `SubstructLibrary.GetMatches(query)`.")]
        pub fn $gm(&self, query: &$type) -> Result<Vec<u32>> {
            self.$gmp(query, &GetMatchesParams::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.GetMatches(query, parameters, numThreads, maxResults)`, for ", $what, ", and its keyword form.")]
        pub fn $gmp(&self, query: &$type, params: &GetMatchesParams) -> Result<Vec<u32>> {
            self.search(None, query, |view, range, raw| {
                raw.get_matches(view, range, &params.parameters.to_ffi(), params.max_results)
            })
        }

        #[doc = concat!(
            "The molecules that match ", $what, " from `start_idx` up to `end_idx` (positions in the search ",
            "order, if it is set), which is clipped to the library. A `start_idx` past the end is `Error::Index`, ",
            "and an `end_idx` at or before it `Error::Runtime`. Mirrors `SubstructLibrary.GetMatches(query, startIdx, endIdx)`."
        )]
        pub fn $gmr(&self, query: &$type, start_idx: u32, end_idx: u32) -> Result<Vec<u32>> {
            self.$gmrp(query, start_idx, end_idx, &GetMatchesParams::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.GetMatches(query, startIdx, endIdx, parameters, numThreads, maxResults)`, for ", $what, ", and its keyword form.")]
        pub fn $gmrp(
            &self,
            query: &$type,
            start_idx: u32,
            end_idx: u32,
            params: &GetMatchesParams,
        ) -> Result<Vec<u32>> {
            self.search(Some((start_idx, end_idx)), query, |view, range, raw| {
                raw.get_matches(view, range, &params.parameters.to_ffi(), params.max_results)
            })
        }

        #[doc = concat!("The number of molecules that match ", $what, ". Mirrors `SubstructLibrary.CountMatches(query)`.")]
        pub fn $cm(&self, query: &$type) -> Result<u32> {
            self.$cmp(query, &CountMatchesParams::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.CountMatches(query, parameters, numThreads)`, for ", $what, ", and its keyword form.")]
        pub fn $cmp(&self, query: &$type, params: &CountMatchesParams) -> Result<u32> {
            self.search(None, query, |view, range, raw| {
                raw.count_matches(view, range, &params.parameters.to_ffi())
            })
        }

        #[doc = concat!("As `", stringify!($gmr), "`, counted. Mirrors `SubstructLibrary.CountMatches(query, startIdx, endIdx)`.")]
        pub fn $cmr(&self, query: &$type, start_idx: u32, end_idx: u32) -> Result<u32> {
            self.$cmrp(query, start_idx, end_idx, &CountMatchesParams::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.CountMatches(query, startIdx, endIdx, parameters, numThreads)`, for ", $what, ", and its keyword form.")]
        pub fn $cmrp(
            &self,
            query: &$type,
            start_idx: u32,
            end_idx: u32,
            params: &CountMatchesParams,
        ) -> Result<u32> {
            self.search(Some((start_idx, end_idx)), query, |view, range, raw| {
                raw.count_matches(view, range, &params.parameters.to_ffi())
            })
        }

        #[doc = concat!("Whether any molecule matches ", $what, ". Mirrors `SubstructLibrary.HasMatch(query)`.")]
        pub fn $hm(&self, query: &$type) -> Result<bool> {
            self.$hmp(query, &HasMatchParams::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.HasMatch(query, parameters, numThreads)`, for ", $what, ", and its keyword form.")]
        pub fn $hmp(&self, query: &$type, params: &HasMatchParams) -> Result<bool> {
            self.search(None, query, |view, range, raw| {
                raw.has_match(view, range, &params.parameters.to_ffi())
            })
        }

        #[doc = concat!("As `", stringify!($gmr), "`, for any match. Mirrors `SubstructLibrary.HasMatch(query, startIdx, endIdx)`.")]
        pub fn $hmr(&self, query: &$type, start_idx: u32, end_idx: u32) -> Result<bool> {
            self.$hmrp(query, start_idx, end_idx, &HasMatchParams::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.HasMatch(query, startIdx, endIdx, parameters, numThreads)`, for ", $what, ", and its keyword form.")]
        pub fn $hmrp(
            &self,
            query: &$type,
            start_idx: u32,
            end_idx: u32,
            params: &HasMatchParams,
        ) -> Result<bool> {
            self.search(Some((start_idx, end_idx)), query, |view, range, raw| {
                raw.has_match(view, range, &params.parameters.to_ffi())
            })
        }
    };
}

impl SubstructLibrary {
    /// A library of the molecules, with fingerprints and keys if given:
    /// Python's `(molecules)`, `(molecules, fingerprints_or_keys)`, and
    /// `(molecules, fingerprints, keys)` overloads. The holders should hold
    /// the same number of entries, in the same order. Mirrors
    /// `SubstructLibrary(molecules, ...)`.
    pub fn new(
        molecules: MolHolderBase,
        fingerprints: Option<FPHolderBase>,
        keys: Option<KeyHolderBase>,
    ) -> SubstructLibrary {
        SubstructLibrary {
            molecules,
            fingerprints,
            keys,
            search_order: Vec::new(),
        }
    }

    /// The library a pickle holds (`serialize`'s bytes, here or in
    /// Python). A pickle whose molecules' plain monomer infos have the PDB
    /// residue type (CONVENTIONS.md U10), whose search order or
    /// fingerprints do not fit its holders, or that has no molecule holder,
    /// is `Error::Value`. Mirrors `SubstructLibrary(pickle)`.
    ///
    /// # Safety
    ///
    /// `pickle` must be bytes `serialize` (this crate's or RDKit's)
    /// returned. RDKit reads it as a Boost.Serialization archive of
    /// polymorphic holders, which trusts its class and object ids
    /// (CONVENTIONS.md U2).
    pub unsafe fn from_pickle(pickle: &[u8]) -> Result<SubstructLibrary> {
        let lib = ffi::substruct_library_from_pickle(pickle)?;
        let lib = lib.as_ref().expect("RDKit returned a null library");
        let fingerprints = expect_ok(ffi::substruct_library_fp_holder(lib));
        let keys = expect_ok(ffi::substruct_library_key_holder(lib));
        Ok(SubstructLibrary {
            molecules: MolHolderBase::from_ffi(expect_ok(ffi::substruct_library_mol_holder(lib))),
            fingerprints: (!fingerprints.is_null()).then(|| FPHolderBase::from_ffi(fingerprints)),
            keys: (!keys.is_null()).then(|| KeyHolderBase::from_ffi(keys)),
            search_order: expect_ok(ffi::substruct_library_search_order(lib)),
        })
    }

    /// An RDKit library sharing the holders, for one call, with the search
    /// order if `with_order`. A search order the molecule holder no longer
    /// fits (one swapped in through `get_mol_holder_mut`) is `Error::Index`.
    fn view(&self, with_order: bool) -> Result<UniquePtr<ffi::SubstructLibrary>> {
        let mut view = expect_ok(ffi::substruct_library_view(self.molecules.raw()));
        if let Some(fingerprints) = &self.fingerprints {
            expect_ok(ffi::substruct_library_view_set_fp_holder(
                view.pin_mut(),
                fingerprints.raw(),
            ));
        }
        if let Some(keys) = &self.keys {
            expect_ok(ffi::substruct_library_view_set_key_holder(
                view.pin_mut(),
                keys.raw(),
            ));
        }
        if with_order && !self.search_order.is_empty() {
            ffi::substruct_library_view_set_search_order(view.pin_mut(), &self.search_order)?;
        }
        Ok(view)
    }

    // Runs one search of the query on a view of the library, over the
    // range if there is one.
    fn search<T>(
        &self,
        range: Option<(u32, u32)>,
        query: &impl LibraryQuery,
        f: impl FnOnce(&ffi::SubstructLibrary, Range, RawQuery<'_>) -> Raw<T>,
    ) -> Result<T> {
        let view = self.view(true)?;
        let (start_idx, end_idx) = range.unwrap_or((0, 0));
        let range = Range {
            ranged: range.is_some(),
            start_idx,
            end_idx,
        };
        Ok(query.with_raw(|raw| f(&view, range, raw))?)
    }

    /// The molecule holder. Mirrors `SubstructLibrary.GetMolHolder`.
    pub fn get_mol_holder(&self) -> &MolHolderBase {
        &self.molecules
    }

    /// The molecule holder, to add to it alone or to replace it. Mirrors
    /// `SubstructLibrary.GetMolHolder`.
    pub fn get_mol_holder_mut(&mut self) -> &mut MolHolderBase {
        &mut self.molecules
    }

    /// The fingerprint holder, if the library has one. Mirrors
    /// `SubstructLibrary.GetFpHolder`.
    pub fn get_fp_holder(&self) -> Option<&FPHolderBase> {
        self.fingerprints.as_ref()
    }

    /// The fingerprint holder, if the library has one. Mirrors
    /// `SubstructLibrary.GetFpHolder`.
    pub fn get_fp_holder_mut(&mut self) -> Option<&mut FPHolderBase> {
        self.fingerprints.as_mut()
    }

    /// The key holder, if the library has one. Mirrors
    /// `SubstructLibrary.GetKeyHolder`.
    pub fn get_key_holder(&self) -> Option<&KeyHolderBase> {
        self.keys.as_ref()
    }

    /// The key holder, if the library has one. Mirrors
    /// `SubstructLibrary.GetKeyHolder`.
    pub fn get_key_holder_mut(&mut self) -> Option<&mut KeyHolderBase> {
        self.keys.as_mut()
    }

    /// Adds a molecule to each holder and returns its index. Holders whose
    /// sizes differ are `Error::Runtime` after the molecule holder has the
    /// molecule, as in RDKit, and full ones `Error::Value`, as is a molecule
    /// a `CachedMolHolder` refuses (`MolHolderBase::add_mol`). Mirrors
    /// `SubstructLibrary.AddMol`.
    pub fn add_mol(&mut self, mol: &Mol) -> Result<u32> {
        let mut view = self.view(false)?;
        Ok(ffi::substruct_library_add_mol(
            view.pin_mut(),
            mol.as_ffi(),
        )?)
    }

    searches!(
        Mol, "the query";
        get_matches, get_matches_with_params, get_matches_range, get_matches_range_with_params;
        count_matches, count_matches_with_params, count_matches_range,
        count_matches_range_with_params;
        has_match, has_match_with_params, has_match_range, has_match_range_with_params
    );

    searches!(
        MolBundle, "any of the bundle's members (member by member, each skipping the molecules earlier members matched)";
        get_matches_mol_bundle, get_matches_mol_bundle_with_params, get_matches_mol_bundle_range,
        get_matches_mol_bundle_range_with_params;
        count_matches_mol_bundle, count_matches_mol_bundle_with_params,
        count_matches_mol_bundle_range, count_matches_mol_bundle_range_with_params;
        has_match_mol_bundle, has_match_mol_bundle_with_params, has_match_mol_bundle_range,
        has_match_mol_bundle_range_with_params
    );

    searches!(
        TautomerQuery, "any of the tautomer query's tautomers (a template with a recursive query is `Error::Value`; see the module docs)";
        get_matches_tautomer_query, get_matches_tautomer_query_with_params,
        get_matches_tautomer_query_range, get_matches_tautomer_query_range_with_params;
        count_matches_tautomer_query, count_matches_tautomer_query_with_params,
        count_matches_tautomer_query_range, count_matches_tautomer_query_range_with_params;
        has_match_tautomer_query, has_match_tautomer_query_with_params,
        has_match_tautomer_query_range, has_match_tautomer_query_range_with_params
    );

    searches!(
        ExtendedQueryMol, "the extended query (one whose tautomer templates have a recursive query is `Error::Value`; see the module docs)";
        get_matches_extended_query_mol, get_matches_extended_query_mol_with_params,
        get_matches_extended_query_mol_range, get_matches_extended_query_mol_range_with_params;
        count_matches_extended_query_mol, count_matches_extended_query_mol_with_params,
        count_matches_extended_query_mol_range,
        count_matches_extended_query_mol_range_with_params;
        has_match_extended_query_mol, has_match_extended_query_mol_with_params,
        has_match_extended_query_mol_range, has_match_extended_query_mol_range_with_params
    );

    /// A copy of the molecule at `idx`, as `MolHolderBase::get_mol`.
    /// Mirrors `SubstructLibrary.GetMol`.
    pub fn get_mol(&self, idx: u32) -> Result<Mol> {
        self.molecules.get_mol(idx)
    }

    /// Searches the molecules in this order, which may leave some out or
    /// repeat them; an empty order (Python's `None`) searches them in
    /// turn. An index past the molecules is `Error::Index`. Mirrors
    /// `SubstructLibrary.SetSearchOrder`.
    pub fn set_search_order(&mut self, seq: &[u32]) -> Result<()> {
        let len = self.len();
        if let Some(&index) = seq.iter().find(|&&idx| idx >= len) {
            return Err(Error::index(index, len as usize));
        }
        self.search_order = seq.to_vec();
        Ok(())
    }

    /// The search order, empty when unset. Mirrors
    /// `SubstructLibrary.GetSearchOrder`.
    pub fn get_search_order(&self) -> &[u32] {
        &self.search_order
    }

    /// The number of molecules. Mirrors `SubstructLibrary.__len__`.
    pub fn len(&self) -> u32 {
        self.molecules.len()
    }

    /// Whether the library has no molecules.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The library as RDKit pickles it, a Boost.Serialization text archive;
    /// `Error::Value` if RDKit was built without Boost.Serialization
    /// (`substruct_library_can_serialize`), or when a `MolHolder` holds a
    /// molecule `Mol::to_binary` refuses. Mirrors
    /// `SubstructLibrary.Serialize`.
    pub fn serialize(&self) -> Result<Vec<u8>> {
        let view = self.view(true)?;
        Ok(ffi::substruct_library_serialize(&view)?)
    }

    /// Writes `serialize`'s bytes; a failed write is `Error::Os`. Mirrors
    /// `SubstructLibrary.ToStream`.
    pub fn to_stream<W: Write + ?Sized>(&self, stream: &mut W) -> Result<()> {
        let pickle = self.serialize()?;
        stream
            .write_all(&pickle)
            .and_then(|()| stream.flush())
            .map_err(|err| Error::Os {
                message: err.to_string(),
            })
    }

    /// Replaces the library with the one a stream holds, read to its end
    /// and checked as `from_pickle` checks it; on an error the library is
    /// unchanged. A failed read is `Error::Os`. Mirrors
    /// `SubstructLibrary.InitFromStream`.
    ///
    /// # Safety
    ///
    /// As `from_pickle`: the stream must hold bytes `serialize` returned.
    pub unsafe fn init_from_stream<R: Read + ?Sized>(&mut self, stream: &mut R) -> Result<()> {
        let mut pickle = Vec::new();
        stream.read_to_end(&mut pickle).map_err(|err| Error::Os {
            message: err.to_string(),
        })?;
        // SAFETY: the caller's contract.
        *self = unsafe { SubstructLibrary::from_pickle(&pickle)? };
        Ok(())
    }
}

/// An empty library of molecules held as copies, without fingerprints or
/// keys. Mirrors `SubstructLibrary()`.
impl Default for SubstructLibrary {
    fn default() -> SubstructLibrary {
        SubstructLibrary::new(MolHolder::default().into(), None, None)
    }
}

impl fmt::Debug for SubstructLibrary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SubstructLibrary")
            .field("molecules", &self.molecules)
            .field("fingerprints", &self.fingerprints)
            .field("keys", &self.keys)
            .field("search_order", &self.search_order)
            .finish()
    }
}

/// Gives the library pattern fingerprints of its molecules (an empty
/// molecule's for one whose SMILES do not parse); a library with
/// fingerprints is `Error::Runtime`. Mirrors `rdSubstructLibrary.AddPatterns`.
pub fn add_patterns(sslib: &mut SubstructLibrary) -> Result<()> {
    add_patterns_holder(sslib, PatternHolder::default().into())
}

/// As `add_patterns`; `num_threads` is accepted for parity, and the
/// fingerprints are made on the calling thread (see the module
/// documentation). Mirrors `rdSubstructLibrary.AddPatterns(sslib,
/// numThreads)`.
pub fn add_patterns_with_num_threads(sslib: &mut SubstructLibrary, num_threads: i32) -> Result<()> {
    let _ = num_threads;
    add_patterns(sslib)
}

/// Gives the library `patterns`, whose fingerprints are replaced by those
/// of the library's molecules. On an error (a stored SMILES that does not
/// sanitize, where Python leaves null fingerprints in the holder that later
/// reads crash on), the library is unchanged. Mirrors
/// `rdSubstructLibrary.AddPatterns(sslib, patterns)`.
pub fn add_patterns_holder(sslib: &mut SubstructLibrary, patterns: FPHolderBase) -> Result<()> {
    if sslib.fingerprints.is_some() {
        // As RDKit's precondition.
        return Err(Error::Runtime {
            message: "Substruct Library already has fingerprints".to_owned(),
        });
    }
    let mut patterns = patterns;
    let view = sslib.view(false)?;
    ffi::add_patterns(&view, patterns.raw_mut())?;
    sslib.fingerprints = Some(patterns);
    Ok(())
}

/// As `add_patterns_holder`, with `num_threads` accepted for parity.
/// Mirrors `rdSubstructLibrary.AddPatterns(sslib, patterns, numThreads)`.
pub fn add_patterns_holder_with_num_threads(
    sslib: &mut SubstructLibrary,
    patterns: FPHolderBase,
    num_threads: i32,
) -> Result<()> {
    let _ = num_threads;
    add_patterns_holder(sslib, patterns)
}

/// Whether libraries can be serialized (RDKit built with Boost
/// serialization). Mirrors `rdSubstructLibrary.SubstructLibraryCanSerialize`.
pub fn substruct_library_can_serialize() -> bool {
    expect_ok(ffi::substruct_library_can_serialize())
}
