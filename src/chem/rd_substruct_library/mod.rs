// Copyright (c) 2026 Steven Kearnes
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

/// Defines a search's eight forms: a molecule or a bundle, the whole
/// library or a range, Python's defaults or parameters.
macro_rules! searches {
    (
        $python:literal, $out:ty, $params:ty;
        $name:ident, $with_params:ident, $range:ident, $range_with_params:ident,
        $bundle:ident, $bundle_with_params:ident, $bundle_range:ident,
        $bundle_range_with_params:ident;
        |$view:ident, $query:ident, $ranged:ident, $start:ident, $end:ident, $p:ident| $mol_call:expr,
        $bundle_call:expr
    ) => {
        #[doc = concat!("Mirrors `SubstructLibrary.", $python, "(query)`.")]
        pub fn $name(&self, query: &Mol) -> Result<$out> {
            self.$with_params(query, &<$params>::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.", $python, "(query, parameters, ...)` and its keyword form.")]
        pub fn $with_params(&self, query: &Mol, params: &$params) -> Result<$out> {
            self.search(None, |$view, $ranged, $start, $end| {
                let ($query, $p) = (query.as_ffi(), params);
                $mol_call
            })
        }

        #[doc = concat!(
            "Searches the molecules (positions in the search order, if it is set) from `start_idx` up to `end_idx`, ",
            "which is clipped to the library. A `start_idx` past the end is `Error::Index`, and an `end_idx` ",
            "at or before it `Error::Runtime`. Mirrors `SubstructLibrary.", $python, "(query, startIdx, endIdx)`."
        )]
        pub fn $range(&self, query: &Mol, start_idx: u32, end_idx: u32) -> Result<$out> {
            self.$range_with_params(query, start_idx, end_idx, &<$params>::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.", $python, "(query, startIdx, endIdx, parameters, ...)` and its keyword form.")]
        pub fn $range_with_params(
            &self,
            query: &Mol,
            start_idx: u32,
            end_idx: u32,
            params: &$params,
        ) -> Result<$out> {
            self.search(Some((start_idx, end_idx)), |$view, $ranged, $start, $end| {
                let ($query, $p) = (query.as_ffi(), params);
                $mol_call
            })
        }

        #[doc = concat!(
            "The molecules any of the bundle's members match. A bundle of more than one member needs ",
            "(members - 1) times the molecules searched below 2^31, or RDKit's count overflows `int` ",
            "(`Error::Value`). Mirrors `SubstructLibrary.", $python, "(query: MolBundle)`."
        )]
        pub fn $bundle(&self, query: &MolBundle) -> Result<$out> {
            self.$bundle_with_params(query, &<$params>::default())
        }

        #[doc = concat!("Mirrors `SubstructLibrary.", $python, "(query: MolBundle, parameters, ...)` and its keyword form.")]
        pub fn $bundle_with_params(&self, query: &MolBundle, params: &$params) -> Result<$out> {
            query.with_ffi(|$query| {
                self.search(None, |$view, $ranged, $start, $end| {
                    let $p = params;
                    $bundle_call
                })
            })
        }

        #[doc = concat!("Mirrors `SubstructLibrary.", $python, "(query: MolBundle, startIdx, endIdx)`.")]
        pub fn $bundle_range(&self, query: &MolBundle, start_idx: u32, end_idx: u32) -> Result<$out> {
            self.$bundle_range_with_params(query, start_idx, end_idx, &<$params>::default())
        }

        #[doc = concat!(
            "Mirrors `SubstructLibrary.", $python,
            "(query: MolBundle, startIdx, endIdx, parameters, ...)` and its keyword form."
        )]
        pub fn $bundle_range_with_params(
            &self,
            query: &MolBundle,
            start_idx: u32,
            end_idx: u32,
            params: &$params,
        ) -> Result<$out> {
            query.with_ffi(|$query| {
                self.search(Some((start_idx, end_idx)), |$view, $ranged, $start, $end| {
                    let $p = params;
                    $bundle_call
                })
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

    fn search<T>(
        &self,
        range: Option<(u32, u32)>,
        f: impl FnOnce(&ffi::SubstructLibrary, bool, u32, u32) -> std::result::Result<T, cxx::Exception>,
    ) -> Result<T> {
        let view = self.view(true)?;
        let (start_idx, end_idx) = range.unwrap_or((0, 0));
        Ok(f(&view, range.is_some(), start_idx, end_idx)?)
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
    /// molecule, as in RDKit, and full ones `Error::Value`. Mirrors
    /// `SubstructLibrary.AddMol`.
    pub fn add_mol(&mut self, mol: &Mol) -> Result<u32> {
        let mut view = self.view(false)?;
        Ok(ffi::substruct_library_add_mol(
            view.pin_mut(),
            mol.as_ffi(),
        )?)
    }

    searches!(
        "GetMatches", Vec<u32>, GetMatchesParams;
        get_matches, get_matches_with_params, get_matches_range, get_matches_range_with_params,
        get_matches_mol_bundle, get_matches_mol_bundle_with_params, get_matches_mol_bundle_range,
        get_matches_mol_bundle_range_with_params;
        |view, query, ranged, start, end, p| ffi::substruct_library_get_matches(
            view, query, ranged, start, end, &p.parameters.to_ffi(), p.max_results,
        ),
        ffi::substruct_library_get_matches_mol_bundle(
            view, query, ranged, start, end, &p.parameters.to_ffi(), p.max_results,
        )
    );

    searches!(
        "CountMatches", u32, CountMatchesParams;
        count_matches, count_matches_with_params, count_matches_range,
        count_matches_range_with_params, count_matches_mol_bundle,
        count_matches_mol_bundle_with_params, count_matches_mol_bundle_range,
        count_matches_mol_bundle_range_with_params;
        |view, query, ranged, start, end, p| ffi::substruct_library_count_matches(
            view, query, ranged, start, end, &p.parameters.to_ffi(),
        ),
        ffi::substruct_library_count_matches_mol_bundle(
            view, query, ranged, start, end, &p.parameters.to_ffi(),
        )
    );

    searches!(
        "HasMatch", bool, HasMatchParams;
        has_match, has_match_with_params, has_match_range, has_match_range_with_params,
        has_match_mol_bundle, has_match_mol_bundle_with_params, has_match_mol_bundle_range,
        has_match_mol_bundle_range_with_params;
        |view, query, ranged, start, end, p| ffi::substruct_library_has_match(
            view, query, ranged, start, end, &p.parameters.to_ffi(),
        ),
        ffi::substruct_library_has_match_mol_bundle(
            view, query, ranged, start, end, &p.parameters.to_ffi(),
        )
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
            return Err(Error::Index { index, len });
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
    /// (`substruct_library_can_serialize`). Mirrors
    /// `SubstructLibrary.Serialize`.
    pub fn serialize(&self) -> Result<Vec<u8>> {
        let view = self.view(true)?;
        Ok(ffi::substruct_library_serialize(&view)?)
    }

    /// Writes `serialize`'s bytes; a failed write is `Error::Os`. Mirrors
    /// `SubstructLibrary.ToStream`.
    pub fn to_stream(&self, stream: &mut dyn Write) -> Result<()> {
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
    pub unsafe fn init_from_stream(&mut self, stream: &mut dyn Read) -> Result<()> {
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
