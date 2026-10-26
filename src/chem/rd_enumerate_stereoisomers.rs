// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdEnumerateStereoisomers`: the stereoisomers of a
//! molecule, one at a time. Chiral centers and double bonds are enumerated
//! when they are unassigned (or, with `only_unassigned` off, also when they
//! are assigned); atropisomers only when they are assigned.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smiles;
//! use rdkit_cxx::chem::rd_enumerate_stereoisomers::StereoisomerEnumerator;
//!
//! let mol = mol_from_smiles("CC(F)Cl")?;
//! for isomer in StereoisomerEnumerator::new(&mol)? {
//!     let isomer = isomer?;
//! }
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

use std::fmt;
use std::iter::FusedIterator;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::enumerate_stereoisomers::ffi;

/// Options for `StereoisomerEnumerator::new_with_params`, with Python's
/// defaults, read from RDKit's (A3). Mirrors
/// `rdEnumerateStereoisomers.StereoEnumerationOptions`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoEnumerationOptions {
    /// Whether each isomer is embedded (with hydrogens, and its heavy atoms'
    /// coordinates kept as a conformer), and dropped if that fails: slow,
    /// and a heuristic that can lose isomers. Mirrors
    /// `StereoEnumerationOptions.tryEmbedding`.
    pub try_embedding: bool,
    /// Whether assigned stereocenters stay as they are, unless they belong to
    /// a relative stereo group. Mirrors `StereoEnumerationOptions.onlyUnassigned`.
    pub only_unassigned: bool,
    /// Whether only the stereo groups vary. Mirrors
    /// `StereoEnumerationOptions.onlyStereoGroups`.
    pub only_stereo_groups: bool,
    /// Whether isomers with the same canonical CXSMILES are made once.
    /// Mirrors `StereoEnumerationOptions.unique`.
    pub unique: bool,
    /// The most isomers to make, a random subset when there are more; 0 for
    /// no limit. Mirrors `StereoEnumerationOptions.maxIsomers`.
    pub max_isomers: u64,
    /// The seed of the enumerator's own random generator; -1 for a random
    /// seed. Mirrors `StereoEnumerationOptions.randomSeed`.
    pub random_seed: i32,
    /// Whether isomers that fail to embed are logged, to RDKit's info log;
    /// Python's `StereoisomerEnumerator(mol, options, verbose)` argument.
    pub verbose: bool,
}

impl Default for StereoEnumerationOptions {
    fn default() -> StereoEnumerationOptions {
        let data = expect_ok(ffi::stereo_enumeration_options_default());
        StereoEnumerationOptions {
            try_embedding: data.try_embedding,
            only_unassigned: data.only_unassigned,
            only_stereo_groups: data.only_stereo_groups,
            unique: data.unique,
            max_isomers: data.max_isomers,
            random_seed: data.random_seed,
            verbose: false,
        }
    }
}

impl StereoEnumerationOptions {
    fn to_ffi(&self) -> ffi::StereoEnumerationOptionsData {
        ffi::StereoEnumerationOptionsData {
            try_embedding: self.try_embedding,
            only_unassigned: self.only_unassigned,
            only_stereo_groups: self.only_stereo_groups,
            unique: self.unique,
            max_isomers: self.max_isomers,
            random_seed: self.random_seed,
        }
    }
}

/// The stereoisomers of a molecule, in random order, as an iterator: each
/// item is an isomer, or the error RDKit raised making it. The enumerator
/// owns a copy of the molecule. Mirrors
/// `rdEnumerateStereoisomers.StereoisomerEnumerator`.
///
/// With `try_embedding`, the embedder's limits apply to the molecule with
/// hydrogens (`Error::Value` from the first isomer, as `embed_multiple_confs`
/// gives), and making each isomer holds the locks the embedder takes, so
/// threads that embed run one at a time.
pub struct StereoisomerEnumerator {
    /// Never null.
    ptr: UniquePtr<ffi::StereoisomerEnumerator>,
}

impl fmt::Debug for StereoisomerEnumerator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("StereoisomerEnumerator")
            .field("stereoisomer_count", &self.get_stereoisomer_count())
            .finish()
    }
}

impl StereoisomerEnumerator {
    /// An enumerator with the default options. Fails where RDKit's stereo
    /// perception raises for the molecule. Mirrors
    /// `StereoisomerEnumerator(mol)`.
    pub fn new(mol: &Mol) -> Result<StereoisomerEnumerator> {
        StereoisomerEnumerator::new_with_params(mol, &StereoEnumerationOptions::default())
    }

    /// Mirrors `StereoisomerEnumerator(mol, verbose)`.
    pub fn new_with_verbose(mol: &Mol, verbose: bool) -> Result<StereoisomerEnumerator> {
        StereoisomerEnumerator::new_with_params(
            mol,
            &StereoEnumerationOptions {
                verbose,
                ..Default::default()
            },
        )
    }

    /// Mirrors `StereoisomerEnumerator(mol, options, verbose)`; see
    /// `StereoEnumerationOptions`.
    pub fn new_with_params(
        mol: &Mol,
        options: &StereoEnumerationOptions,
    ) -> Result<StereoisomerEnumerator> {
        let ptr =
            ffi::stereoisomer_enumerator_new(mol.as_ffi(), &options.to_ffi(), options.verbose)?;
        assert!(!ptr.is_null(), "RDKit returned a null enumerator");
        Ok(StereoisomerEnumerator { ptr })
    }

    /// The number of possible stereoisomers, 2 to the number of stereo
    /// elements to vary (`u64::MAX` past 63), which is at least the number the
    /// enumerator makes. Mirrors `StereoisomerEnumerator.GetStereoisomerCount`.
    pub fn get_stereoisomer_count(&self) -> u64 {
        expect_ok(ffi::stereoisomer_enumerator_get_stereoisomer_count(
            &self.ptr,
        ))
    }
}

/// The next isomer; `None` once there are no more. Mirrors
/// `StereoisomerEnumerator.next`, which returns `None` then.
impl Iterator for StereoisomerEnumerator {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        match ffi::stereoisomer_enumerator_next(self.ptr.pin_mut()) {
            Ok(ptr) if ptr.is_null() => None,
            Ok(ptr) => Some(Ok(Mol::from_ffi(ptr))),
            Err(err) => Some(Err(err.into())),
        }
    }
}

// RDKit returns no isomer again once it has returned none.
impl FusedIterator for StereoisomerEnumerator {}
