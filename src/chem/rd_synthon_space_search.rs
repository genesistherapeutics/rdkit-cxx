// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdSynthonSpaceSearch`: searches of synthon spaces,
//! combinatorial libraries such as Enamine REAL held as reactions and the
//! synthons each combines, without enumerating their products. RDKit calls
//! this functionality experimental.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smarts;
//! use rdkit_cxx::chem::rd_synthon_space_search::SynthonSpace;
//!
//! let mut space = SynthonSpace::default();
//! space.read_text_file("amide_space.txt")?;
//! let query = mol_from_smarts("c1ccccc1C(=O)N")?;
//! let hits = space.substructure_search(&query)?.get_hit_molecules();
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```
//!
//! Each search, read, and fingerprint build runs on the calling thread,
//! whatever `num_threads` asks for, and holds the crate's Ctrl-C lock:
//! RDKit's searches swap the process's SIGINT handler through static
//! variables, as every substructure match does, and its own threads would do
//! so at once and call `std::terminate` on an exception. The hits do not
//! depend on the thread count. RDKit's handler stops a search at a Ctrl-C,
//! which then returns `Error::Runtime`, as Python raises
//! `KeyboardInterrupt`.
//!
//! The incremental searches hand each batch of hits to a closure, on the
//! calling thread, which returns `Ok(true)` to stop the search. Its error,
//! or a panic in it (as `Error::Runtime`), ends the search, which returns
//! it. The substructure and fingerprint searches call the closure with
//! RDKit's Ctrl-C handler installed, so the closure cannot call functions
//! that match substructures, depict, or take the shape modules' locks: they
//! return `Error::Runtime`. The RASCAL search's closure can.
//!
//! The reads, the searches, and shape building remove atoms and bonds from
//! the synthons, the products they combine, the query's fragments, and the
//! molecules a conformer generator returns, so what RDKit would misread then
//! is refused with `Error::Value`: a bond's `_MolFileBondEndPts`, or
//! substance groups whose `PARENT`s form a cycle (RDKit follows it forever).
//! A reaction's synthons' chains are checked together, as a product carries
//! them all.
//!
//! Shapes are built (`build_synthon_shapes`) for the database files, but
//! the shape search and the progress bars are not bound (the manifest says
//! why).

use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rd_enumerate_stereoisomers::StereoEnumerationOptions;
use crate::chem::rd_fingerprint_generator::FingerprintGenerator64;
use crate::chem::rd_generalized_substruct::ExtendedQueryMol;
use crate::chem::rd_rascal_mces::RascalOptions;
use crate::chem::rdchem::SubstructMatchParameters;
use crate::chem::rdmolops::mol_list_into_vec;
use crate::error::{expect_ok, to_string};
use crate::ffi::synthon_space_search::{ConformerCallback, ConformerFn, SearchCallback, ffi};
use crate::util::{path_bytes, path_from_bytes};

/// Parameters of the searches. Mirrors
/// `rdSynthonSpaceSearch.SynthonSpaceSearchParams` (A3); `Default` reads
/// RDKit's defaults. The shape search's fields are not bound (the manifest
/// says why).
#[derive(Clone, Debug, PartialEq)]
pub struct SynthonSpaceSearchParams {
    /// The most hits to return; -1 for no limit, and below -1 is
    /// `Error::Value`.
    pub max_hits: i64,
    /// The most fragment sets the query is split into; a search that
    /// reaches it is likely incomplete.
    pub max_num_frags: u64,
    /// The number of the first hit to return, to page through hits; negative
    /// is `Error::Value`, as is a sum with `max_hits` past `i64::MAX`.
    pub hit_start: i64,
    /// Returns a random sample of the hits, up to `max_hits` (which may not
    /// then be -1).
    pub random_sample: bool,
    /// The seed of the random sample's generator, which is the search's own;
    /// -1 for a random seed.
    pub random_seed: i32,
    /// Builds the hits; without it the search only counts them.
    pub build_hits: bool,
    /// The similarity a fingerprint search's hit needs.
    pub similarity_cutoff: f64,
    /// How much lower the similarity cutoff is for the fragments.
    pub frag_similarity_adjuster: f64,
    /// How much lower the cutoff is for the approximate similarity the
    /// fingerprint search screens possible hits with.
    pub approx_similarity_adjuster: f64,
    /// The fewest heavy atoms in a hit.
    pub min_hit_heavy_atoms: u32,
    /// The most heavy atoms in a hit; -1 for no limit.
    pub max_hit_heavy_atoms: i32,
    /// The lowest molecular weight of a hit.
    pub min_hit_mol_wt: f64,
    /// The highest molecular weight of a hit; 0 or less for no limit.
    pub max_hit_mol_wt: f64,
    /// The fewest chiral atoms in a hit.
    pub min_hit_chiral_atoms: u32,
    /// The most chiral atoms in a hit; -1 for no limit.
    pub max_hit_chiral_atoms: i32,
    /// The most seconds a search takes; 0 for no limit, and a limit past the
    /// clock's range is `Error::Value`.
    pub time_out: u64,
    /// How many possible hits are built and checked at a time, and the
    /// largest batch an incremental search hands its closure.
    pub to_try_chunk_size: u64,
    /// Not read: the search runs on the calling thread (see the module
    /// documentation).
    pub num_threads: i32,
    /// A file the search writes its possible hits to, one line of synthon
    /// SMILES and a name each, emptied first, for a later `*_lines` search;
    /// `None` for none. RDKit opens the file itself and ignores a failed
    /// write: a file that cannot be opened is `Error::Os` before the search,
    /// and one that does not end with a complete line after it.
    pub possible_hits_file: Option<PathBuf>,
    /// The most lines written to `possible_hits_file`.
    pub max_possible_hits_to_write: u64,
    /// Writes `possible_hits_file` and stops, without building the hits.
    pub write_possible_hits_and_stop: bool,
}

impl Default for SynthonSpaceSearchParams {
    fn default() -> SynthonSpaceSearchParams {
        let d = expect_ok(ffi::synthon_space_search_params_default());
        SynthonSpaceSearchParams {
            max_hits: d.max_hits,
            max_num_frags: d.max_num_frags,
            hit_start: d.hit_start,
            random_sample: d.random_sample,
            random_seed: d.random_seed,
            build_hits: d.build_hits,
            similarity_cutoff: d.similarity_cutoff,
            frag_similarity_adjuster: d.frag_similarity_adjuster,
            approx_similarity_adjuster: d.approx_similarity_adjuster,
            min_hit_heavy_atoms: d.min_hit_heavy_atoms,
            max_hit_heavy_atoms: d.max_hit_heavy_atoms,
            min_hit_mol_wt: d.min_hit_mol_wt,
            max_hit_mol_wt: d.max_hit_mol_wt,
            min_hit_chiral_atoms: d.min_hit_chiral_atoms,
            max_hit_chiral_atoms: d.max_hit_chiral_atoms,
            time_out: d.time_out,
            to_try_chunk_size: d.to_try_chunk_size,
            num_threads: d.num_threads,
            possible_hits_file: (!d.possible_hits_file.is_empty())
                .then(|| path_from_bytes(d.possible_hits_file)),
            max_possible_hits_to_write: d.max_possible_hits_to_write,
            write_possible_hits_and_stop: d.write_possible_hits_and_stop,
        }
    }
}

impl SynthonSpaceSearchParams {
    fn to_ffi(&self) -> ffi::SynthonSpaceSearchParamsData {
        ffi::SynthonSpaceSearchParamsData {
            max_hits: self.max_hits,
            max_num_frags: self.max_num_frags,
            hit_start: self.hit_start,
            random_sample: self.random_sample,
            random_seed: self.random_seed,
            build_hits: self.build_hits,
            similarity_cutoff: self.similarity_cutoff,
            frag_similarity_adjuster: self.frag_similarity_adjuster,
            approx_similarity_adjuster: self.approx_similarity_adjuster,
            min_hit_heavy_atoms: self.min_hit_heavy_atoms,
            max_hit_heavy_atoms: self.max_hit_heavy_atoms,
            min_hit_mol_wt: self.min_hit_mol_wt,
            max_hit_mol_wt: self.max_hit_mol_wt,
            min_hit_chiral_atoms: self.min_hit_chiral_atoms,
            max_hit_chiral_atoms: self.max_hit_chiral_atoms,
            time_out: self.time_out,
            to_try_chunk_size: self.to_try_chunk_size,
            num_threads: self.num_threads,
            possible_hits_file: self
                .possible_hits_file
                .as_deref()
                .map(|path| path_bytes(path).to_vec())
                .unwrap_or_default(),
            max_possible_hits_to_write: self.max_possible_hits_to_write,
            write_possible_hits_and_stop: self.write_possible_hits_and_stop,
        }
    }
}

/// Parameters for building the synthons' shapes, Gaussian shapes of their
/// conformers. Mirrors `rdSynthonSpaceSearch.ShapeBuildParams` (A3);
/// `Default` reads RDKit's defaults.
#[derive(Clone)]
pub struct ShapeBuildParams {
    /// The most conformers per synthon.
    pub num_confs: u32,
    /// The RMS threshold that prunes the conformers.
    pub rms_threshold: f64,
    /// The similarity below which a synthon's shapes must stay apart; NaN is
    /// `Error::Value`.
    pub shape_sim_threshold: f64,
    /// Not read: shapes are built on the calling thread (see the module
    /// documentation).
    pub num_threads: i32,
    /// The embedding's seed; -1 for a random one.
    pub random_seed: i32,
    /// How the stereoisomers whose shapes are built are enumerated. Its
    /// `verbose` is not read: RDKit's enumerator logs nothing here.
    pub stereo_enum_opts: StereoEnumerationOptions,
    /// Builds shapes only for synthons of at most this many heavy atoms; 0
    /// for all.
    pub max_synthon_atoms: u32,
    /// The most attempts at embedding a synthon.
    pub max_embed_attempts: u32,
    /// The most seconds an embedding takes.
    pub time_out: u32,
    /// A database file the build writes the space to as it goes, to restart
    /// from; `None` for none. RDKit writes it itself and ignores a failed
    /// write, so the last one is checked: `Error::Os` if it does not hold a
    /// space.
    pub interim_file: Option<PathBuf>,
    /// How many shapes the build makes between writes of the interim file; 0
    /// for no writes.
    pub interim_writes: u64,
    /// Generates the synthons' conformers instead of RDKit's embedder;
    /// `set_user_conformer_generator` sets it.
    pub user_conformer_generator: Option<UserConformerGenerator>,
}

/// A conformer generator for `ShapeBuildParams`, which clones share: it
/// gets an isomer's SMILES and the most conformers to make, and returns a
/// molecule with conformers, or `None` when it fails.
/// `ShapeBuildParams::set_user_conformer_generator` says how the build calls
/// it. The C++ counterpart is `SynthonSpaceSearch::UserConfGenerator`.
/// A Rust addition, for the callable Python passes.
#[derive(Clone)]
pub struct UserConformerGenerator(Arc<ConformerFn>);

impl UserConformerGenerator {
    /// A generator that runs `func`.
    pub fn new<F>(func: F) -> UserConformerGenerator
    where
        F: Fn(&str, u32) -> Result<Option<Mol>> + Send + Sync + 'static,
    {
        UserConformerGenerator(Arc::new(func))
    }

    /// Runs the generator.
    pub fn call(&self, smiles: &str, num_conformers: u32) -> Result<Option<Mol>> {
        (self.0)(smiles, num_conformers)
    }
}

impl fmt::Debug for UserConformerGenerator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("UserConformerGenerator")
            .finish_non_exhaustive()
    }
}

impl fmt::Debug for ShapeBuildParams {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShapeBuildParams")
            .field("num_confs", &self.num_confs)
            .field("rms_threshold", &self.rms_threshold)
            .field("shape_sim_threshold", &self.shape_sim_threshold)
            .field("random_seed", &self.random_seed)
            .field("stereo_enum_opts", &self.stereo_enum_opts)
            .field("interim_file", &self.interim_file)
            .field("user_conformer_generator", &self.user_conformer_generator)
            .finish_non_exhaustive()
    }
}

impl Default for ShapeBuildParams {
    fn default() -> ShapeBuildParams {
        let d = expect_ok(ffi::shape_build_params_default());
        ShapeBuildParams {
            num_confs: d.num_confs,
            rms_threshold: d.rms_threshold,
            shape_sim_threshold: d.shape_sim_threshold,
            num_threads: d.num_threads,
            random_seed: d.random_seed,
            stereo_enum_opts: StereoEnumerationOptions {
                try_embedding: d.stereo_try_embedding,
                only_unassigned: d.stereo_only_unassigned,
                only_stereo_groups: d.stereo_only_stereo_groups,
                unique: d.stereo_unique,
                max_isomers: d.stereo_max_isomers,
                random_seed: d.stereo_random_seed,
                verbose: false,
            },
            max_synthon_atoms: d.max_synthon_atoms,
            max_embed_attempts: d.max_embed_attempts,
            time_out: d.time_out,
            interim_file: (!d.interim_file.is_empty()).then(|| path_from_bytes(d.interim_file)),
            interim_writes: d.interim_writes,
            user_conformer_generator: None,
        }
    }
}

impl ShapeBuildParams {
    fn to_ffi(&self) -> ffi::ShapeBuildParamsData {
        ffi::ShapeBuildParamsData {
            num_confs: self.num_confs,
            rms_threshold: self.rms_threshold,
            shape_sim_threshold: self.shape_sim_threshold,
            num_threads: self.num_threads,
            random_seed: self.random_seed,
            stereo_try_embedding: self.stereo_enum_opts.try_embedding,
            stereo_only_unassigned: self.stereo_enum_opts.only_unassigned,
            stereo_only_stereo_groups: self.stereo_enum_opts.only_stereo_groups,
            stereo_unique: self.stereo_enum_opts.unique,
            stereo_max_isomers: self.stereo_enum_opts.max_isomers,
            stereo_random_seed: self.stereo_enum_opts.random_seed,
            max_synthon_atoms: self.max_synthon_atoms,
            max_embed_attempts: self.max_embed_attempts,
            time_out: self.time_out,
            interim_file: self
                .interim_file
                .as_deref()
                .map(|path| path_bytes(path).to_vec())
                .unwrap_or_default(),
            interim_writes: self.interim_writes,
            has_user_conformer_generator: self.user_conformer_generator.is_some(),
        }
    }

    /// Generates the synthons' conformers with `func` instead of RDKit's
    /// embedder: it gets an isomer's SMILES and the most conformers to make,
    /// and returns a molecule with conformers, or `None` when it fails, which
    /// drops the isomer as a failed embedding does (Python's returns null,
    /// which RDKit dereferences). RDKit maps the molecule onto the isomer by
    /// a substructure match. `func` runs on the calling thread, with the
    /// build's locks released, so it may embed, match, or build shapes
    /// itself; its error, or a panic in it (as `Error::Runtime`), ends the
    /// build, which returns it, and coordinates that are not finite or lie
    /// beyond 1e8 are `Error::Value`. Clones of the parameters share it.
    /// Mirrors `ShapeBuildParams.setUserConformerGenerator`.
    pub fn set_user_conformer_generator<F>(&mut self, func: F)
    where
        F: Fn(&str, u32) -> Result<Option<Mol>> + Send + Sync + 'static,
    {
        self.user_conformer_generator = Some(UserConformerGenerator::new(func));
    }
}

/// Parameters for `convert_text_to_db_file_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default)]
pub struct ConvertTextToDBFileParams<'a> {
    /// Fingerprints the synthons with this generator.
    pub fp_gen: Option<&'a FingerprintGenerator64>,
    /// Builds the synthons' shapes with these parameters.
    pub py_shape_params: Option<ShapeBuildParams>,
}

/// Parameters for `SynthonSpace::substructure_search_with_params` and the
/// other substructure searches, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SubstructureSearchParams {
    /// How the query is matched to the synthons and the hits.
    pub substruct_match_params: SubstructMatchParameters,
    /// The search's other parameters.
    pub params: SynthonSpaceSearchParams,
}

/// Parameters for `SynthonSpace::substructure_search_incremental_with_params`,
/// with Python's defaults (A2), which are `SubstructureSearch`'s.
pub type SubstructureSearchIncrementalParams = SubstructureSearchParams;

/// The results of a search: the hits, as copies of RDKit's, and how the
/// search went. Python's searches return it for every kind of search.
/// Mirrors `rdSynthonSpaceSearch.SubstructureResult`.
#[non_exhaustive]
pub struct SubstructureResult {
    /// Never null.
    ptr: UniquePtr<ffi::SearchResults>,
}

impl fmt::Debug for SubstructureResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SubstructureResult")
            .field("max_num_results", &self.get_max_num_results())
            .field("timed_out", &self.get_timed_out())
            .finish_non_exhaustive()
    }
}

impl SubstructureResult {
    fn from_ffi(ptr: UniquePtr<ffi::SearchResults>) -> SubstructureResult {
        assert!(!ptr.is_null(), "RDKit returned null search results");
        SubstructureResult { ptr }
    }

    /// The hits, named by their synthons and reaction (`_Name`), with a
    /// `Similarity` property from a similarity search, most similar first.
    /// Mirrors `SubstructureResult.GetHitMolecules`.
    pub fn get_hit_molecules(&self) -> Vec<Mol> {
        mol_list_into_vec(expect_ok(ffi::substructure_result_get_hit_molecules(
            &self.ptr,
        )))
    }

    /// The most hits the search could have found, an upper bound: duplicate
    /// synthon sets and products that do not match fall away. Mirrors
    /// `SubstructureResult.GetMaxNumResults`.
    pub fn get_max_num_results(&self) -> u64 {
        expect_ok(ffi::substructure_result_get_max_num_results(&self.ptr))
    }

    /// Whether the search ran out of time. Mirrors
    /// `SubstructureResult.GetTimedOut`.
    pub fn get_timed_out(&self) -> bool {
        expect_ok(ffi::substructure_result_get_timed_out(&self.ptr))
    }

    /// Whether a Ctrl-C canceled the search; a search that returns
    /// canceled is an error instead, so this is false. Mirrors
    /// `SubstructureResult.GetCancelled`.
    pub fn get_cancelled(&self) -> bool {
        expect_ok(ffi::substructure_result_get_cancelled(&self.ptr))
    }

    /// The most similar molecule a similarity search built, even below the
    /// cutoff; `None` for a substructure search or when nothing came close.
    /// Mirrors `SubstructureResult.GetBestHit`.
    pub fn get_best_hit(&self) -> Option<Mol> {
        let ptr = expect_ok(ffi::substructure_result_get_best_hit(&self.ptr));
        (!ptr.is_null()).then(|| Mol::from_ffi(ptr))
    }
}

/// A synthon space: reactions, each a list of synthon sets, whose products
/// combine one synthon from each set. Mirrors
/// `rdSynthonSpaceSearch.SynthonSpace`.
///
/// It is `Send` but not `Sync`: it holds RDKit molecules (CONVENTIONS.md
/// S1).
pub struct SynthonSpace {
    /// Never null.
    ptr: UniquePtr<ffi::SynthonSpace>,
}

impl Default for SynthonSpace {
    /// An empty space. Mirrors `SynthonSpace()`.
    fn default() -> SynthonSpace {
        let ptr = expect_ok(ffi::synthon_space_new());
        assert!(!ptr.is_null(), "RDKit returned a null synthon space");
        SynthonSpace { ptr }
    }
}

impl fmt::Debug for SynthonSpace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SynthonSpace")
            .field("num_reactions", &self.get_num_reactions())
            .field("num_synthons", &self.get_num_synthons())
            .field("num_products", &self.get_num_products())
            .finish()
    }
}

/// Runs an incremental search with `callback` lent to C++, and returns the
/// error the callback returned, if it returned one, in place of the
/// exception that ended the search.
fn incremental<F>(
    mut callback: F,
    search: impl FnOnce(&mut SearchCallback<'_>) -> std::result::Result<(), cxx::Exception>,
) -> Result<()>
where
    F: FnMut(Vec<Mol>) -> Result<bool>,
{
    let mut lent = SearchCallback {
        callback: &mut callback,
        error: None,
    };
    let result = search(&mut lent);
    if let Some(error) = lent.error.take() {
        return Err(error);
    }
    Ok(result?)
}

impl SynthonSpace {
    /// Replaces the space with the one in a text file of the kind ChemSpace
    /// and Enamine publish: a header line (`SMILES synton_id synton#
    /// reaction_id`, with an optional `release` column, or `synton_role` for
    /// `synton#`), separated by tabs, other whitespace, or commas, and then
    /// one synthon a line, whose connection points are `[U]`, `[Np]`, `[Pu]`,
    /// and `[Am]` or `[1*]` to `[4*]`. A file RDKit cannot parse is
    /// `Error::Runtime` (or `Error::Value` for a synthon SMILES that does not
    /// parse), and one that cannot be read `Error::Os`; the space is then
    /// unchanged. Python's adds to what the space held; here a read replaces
    /// it, so a failed read leaves no half-read space, whose synthons RDKit
    /// would search without their molecules. RDKit trusts a few things
    /// about the file, which are `Error::Value` here: a synthon number (the
    /// third column) from 0 to 15, at most 4 synthon sets a reaction, a
    /// synthon with an atom other than a dummy, at most 2^31 - 2^20 lines,
    /// and no run of more than 4000 separators. `summarise` names no file
    /// for a space read from text: RDKit is handed the bytes that were
    /// checked, not the file. Mirrors `SynthonSpace.ReadTextFile`.
    pub fn read_text_file<P: AsRef<Path>>(&mut self, in_file: P) -> Result<()> {
        Ok(ffi::synthon_space_read_text_file(
            self.ptr.pin_mut(),
            path_bytes(in_file.as_ref()),
        )?)
    }

    /// Replaces the space with the one in a database file `write_db_file`
    /// wrote. The header must give at least one reaction, a format version
    /// from 3.0 to 3.2 (those this RDKit reads), and record offsets that
    /// ascend through the file, and the space read must hold reactions of 1
    /// to 4 synthon sets with their molecules; anything else is
    /// `Error::Value`, and a file that cannot be opened (RDKit opens it for
    /// writing too) `Error::Os` or `Error::Runtime`. The space is then
    /// unchanged. The file is read on the calling thread; Python's
    /// `numThreads` is not bound. Mirrors `SynthonSpace.ReadDBFile`.
    ///
    /// # Safety
    ///
    /// The file must hold bytes RDKit's `writeDBFile` wrote (this crate's
    /// `write_db_file`, or Python's or C++'s `WriteDBFile` or
    /// `ConvertTextToDBFile` with this RDKit or an earlier one), and must not
    /// change while it is read: RDKit maps it into memory and reads its
    /// records' molecule and bit-vector pickles, Boost archives of shapes,
    /// counts, and bools without validating them.
    pub unsafe fn read_db_file<P: AsRef<Path>>(&mut self, in_file: P) -> Result<()> {
        Ok(ffi::synthon_space_read_db_file(
            self.ptr.pin_mut(),
            path_bytes(in_file.as_ref()),
        )?)
    }

    /// Writes the space to a database file in RDKit's format, which
    /// `read_db_file` reads. RDKit opens the file itself and ignores a failed
    /// write, so the file is read back and checked: one that cannot be
    /// written, or that does not hold the space, is `Error::Os`. Mirrors
    /// `SynthonSpace.WriteDBFile`.
    pub fn write_db_file<P: AsRef<Path>>(&self, out_file: P) -> Result<()> {
        Ok(ffi::synthon_space_write_db_file(
            &self.ptr,
            path_bytes(out_file.as_ref()),
        )?)
    }

    /// Writes every product of the space, one SMILES and name a line;
    /// `Error::Os` if the file cannot be written. Mirrors
    /// `SynthonSpace.WriteEnumeratedFile`.
    pub fn write_enumerated_file<P: AsRef<Path>>(&self, out_file: P) -> Result<()> {
        Ok(ffi::synthon_space_write_enumerated_file(
            &self.ptr,
            path_bytes(out_file.as_ref()),
        )?)
    }

    /// Mirrors `SynthonSpace.GetNumReactions`.
    pub fn get_num_reactions(&self) -> u64 {
        expect_ok(ffi::synthon_space_get_num_reactions(&self.ptr))
    }

    /// The number of products the space can make, counting duplicates.
    /// Mirrors `SynthonSpace.GetNumProducts`.
    pub fn get_num_products(&self) -> u64 {
        expect_ok(ffi::synthon_space_get_num_products(&self.ptr))
    }

    /// The number of distinct synthons, some of which several reactions use.
    /// Mirrors `SynthonSpace.GetNumSynthons`.
    pub fn get_num_synthons(&self) -> u64 {
        expect_ok(ffi::synthon_space_get_num_synthons(&self.ptr))
    }

    /// The number of synthons with shapes. Mirrors
    /// `SynthonSpace.GetNumSynthonsWithShapes`.
    pub fn get_num_synthons_with_shapes(&self) -> u64 {
        expect_ok(ffi::synthon_space_get_num_synthons_with_shapes(&self.ptr))
    }

    /// Writes a summary of the space to standard output (C++'s
    /// `std::cout`). Mirrors `SynthonSpace.Summarise`.
    pub fn summarise(&self) {
        expect_ok(ffi::synthon_space_summarise(&self.ptr));
    }

    /// Writes how many reactions the synthons are used in to standard
    /// output (C++'s `std::cout`). Mirrors `SynthonSpace.ReportSynthonUsage`.
    pub fn report_synthon_usage(&self) {
        expect_ok(ffi::synthon_space_report_synthon_usage(&self.ptr));
    }

    /// The info string of the generator the synthons' fingerprints were made
    /// with; empty when they have none. Mirrors
    /// `SynthonSpace.GetSynthonFingerprintType`.
    pub fn get_synthon_fingerprint_type(&self) -> String {
        to_string(expect_ok(ffi::synthon_space_get_synthon_fingerprint_type(
            &self.ptr,
        )))
        .expect("a fingerprint info string is not UTF-8")
    }

    /// Fingerprints the synthons with the generator, as the first
    /// fingerprint search does, unless they already have its fingerprints.
    /// For atom pairs, a space whose products can pass 46340 atoms is
    /// `Error::Value` (U3). Python's `progressBarWidth` is not bound (the
    /// manifest says why). Mirrors `SynthonSpace.BuildSynthonFingerprints`.
    pub fn build_synthon_fingerprints(
        &mut self,
        fingerprint_generator: &FingerprintGenerator64,
    ) -> Result<()> {
        Ok(ffi::synthon_space_build_synthon_fingerprints(
            self.ptr.pin_mut(),
            fingerprint_generator.as_ffi(),
        )?)
    }

    /// Builds the shapes of the synthons that have none, with RDKit's
    /// defaults: each synthon's conformers in a sample product, for each of
    /// up to 3 stereoisomers, by ETKDGv3. Mirrors
    /// `SynthonSpace.BuildSynthonShapes`.
    pub fn build_synthon_shapes(&mut self) -> Result<()> {
        self.build_synthon_shapes_with_py_params(&ShapeBuildParams::default())
    }

    /// Mirrors `SynthonSpace.BuildSynthonShapes(py_params)`. Every synthon
    /// needs a heavy atom and a connector, and the products at most 613
    /// atoms (so that, with their hydrogens, the embedder takes them), or the
    /// build is `Error::Value`: RDKit sorts samples by their first, and its
    /// shapes need volume. A space already built with `num_confs` conformers
    /// is left as it is. Optimizers may land in other local optima on other
    /// platforms, so the shapes may differ there.
    pub fn build_synthon_shapes_with_py_params(
        &mut self,
        py_params: &ShapeBuildParams,
    ) -> Result<()> {
        let data = py_params.to_ffi();
        let mut callback = ConformerCallback {
            generator: py_params
                .user_conformer_generator
                .as_ref()
                .map(|generator| &*generator.0),
            error: None,
        };
        // SAFETY: the generator is lent for this call only.
        let result = unsafe {
            ffi::synthon_space_build_synthon_shapes(self.ptr.pin_mut(), &data, &mut callback)
        };
        if let Some(error) = callback.error.take() {
            return Err(error);
        }
        Ok(result?)
    }

    /// The products that have the query as a substructure. A query without
    /// atoms is `Error::Runtime`. Mirrors `SynthonSpace.SubstructureSearch`.
    pub fn substructure_search(&mut self, query: &Mol) -> Result<SubstructureResult> {
        self.substructure_search_with_params(query, &SubstructureSearchParams::default())
    }

    /// Mirrors `SynthonSpace.SubstructureSearch(query, substructMatchParams,
    /// params)`; see `SynthonSpaceSearchParams` for the errors.
    pub fn substructure_search_with_params(
        &mut self,
        query: &Mol,
        params: &SubstructureSearchParams,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_substructure_search(
                self.ptr.pin_mut(),
                query.as_ffi(),
                &params.substruct_match_params.to_ffi(),
                &params.params.to_ffi(),
            )?,
        ))
    }

    /// As `substructure_search`, with each molecule of the extended query
    /// (its molecule, bundle, or tautomers) searched in turn, and the hits
    /// merged. Mirrors `SynthonSpace.SubstructureSearch(query:
    /// ExtendedQueryMol)`.
    pub fn substructure_search_extended_query_mol(
        &mut self,
        query: &ExtendedQueryMol,
    ) -> Result<SubstructureResult> {
        self.substructure_search_extended_query_mol_with_params(
            query,
            &SubstructureSearchParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.SubstructureSearch(query: ExtendedQueryMol,
    /// substructMatchParams, params)`.
    pub fn substructure_search_extended_query_mol_with_params(
        &mut self,
        query: &ExtendedQueryMol,
        params: &SubstructureSearchParams,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_substructure_search_extended(
                self.ptr.pin_mut(),
                query.as_ffi(),
                &params.substruct_match_params.to_ffi(),
                &params.params.to_ffi(),
            )?,
        ))
    }

    /// Builds the possible hits from lines `start_line` (from 0) up to
    /// `finish_line` of the params' `possible_hits_file`, which an earlier
    /// search wrote, and returns those that match. RDKit dereferences null
    /// for a line whose SMILES does not parse, which is `Error::Value` here,
    /// so the lines are checked and RDKit reads a copy of them; with no
    /// file the result is empty, and with a file that cannot be opened it
    /// is empty too, after RDKit's error message. `start_line` not below
    /// `finish_line` is `Error::Runtime`. Mirrors
    /// `SynthonSpace.SubstructureSearch(query, substructMatchParams, params,
    /// startLine, finishLine)`.
    pub fn substructure_search_lines(
        &mut self,
        query: &Mol,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SubstructureResult> {
        self.substructure_search_lines_with_params(
            query,
            start_line,
            finish_line,
            &SubstructureSearchParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.SubstructureSearch(query, substructMatchParams,
    /// params, startLine, finishLine)`.
    pub fn substructure_search_lines_with_params(
        &mut self,
        query: &Mol,
        start_line: u64,
        finish_line: u64,
        params: &SubstructureSearchParams,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_substructure_search_lines(
                self.ptr.pin_mut(),
                query.as_ffi(),
                &params.substruct_match_params.to_ffi(),
                &params.params.to_ffi(),
                start_line,
                finish_line,
            )?,
        ))
    }

    /// As `substructure_search_lines`, for an extended query. Mirrors
    /// `SynthonSpace.SubstructureSearch(query: ExtendedQueryMol,
    /// substructMatchParams, params, startLine, finishLine)`.
    pub fn substructure_search_extended_query_mol_lines(
        &mut self,
        query: &ExtendedQueryMol,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SubstructureResult> {
        self.substructure_search_extended_query_mol_lines_with_params(
            query,
            start_line,
            finish_line,
            &SubstructureSearchParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.SubstructureSearch(query: ExtendedQueryMol,
    /// substructMatchParams, params, startLine, finishLine)`.
    pub fn substructure_search_extended_query_mol_lines_with_params(
        &mut self,
        query: &ExtendedQueryMol,
        start_line: u64,
        finish_line: u64,
        params: &SubstructureSearchParams,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_substructure_search_extended_lines(
                self.ptr.pin_mut(),
                query.as_ffi(),
                &params.substruct_match_params.to_ffi(),
                &params.params.to_ffi(),
                start_line,
                finish_line,
            )?,
        ))
    }

    /// As `substructure_search`, handing the hits to `callback` in batches
    /// of at most `to_try_chunk_size` (see the module documentation).
    /// Mirrors `SynthonSpace.SubstructureSearchIncremental`.
    pub fn substructure_search_incremental<F>(&mut self, query: &Mol, callback: F) -> Result<()>
    where
        F: FnMut(Vec<Mol>) -> Result<bool>,
    {
        self.substructure_search_incremental_with_params(
            query,
            callback,
            &SubstructureSearchIncrementalParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.SubstructureSearchIncremental(query, callback,
    /// substructMatchParams, params)`.
    pub fn substructure_search_incremental_with_params<F>(
        &mut self,
        query: &Mol,
        callback: F,
        params: &SubstructureSearchIncrementalParams,
    ) -> Result<()>
    where
        F: FnMut(Vec<Mol>) -> Result<bool>,
    {
        let match_params = params.substruct_match_params.to_ffi();
        let search_params = params.params.to_ffi();
        let space = self.ptr.pin_mut();
        incremental(callback, |lent| {
            // SAFETY: the callback is lent for this call only.
            unsafe {
                ffi::synthon_space_substructure_search_incremental(
                    space,
                    query.as_ffi(),
                    lent,
                    &match_params,
                    &search_params,
                )
            }
        })
    }

    /// The products whose fingerprints, from the generator, are at least
    /// `similarity_cutoff` similar to the query's, by Tanimoto; the synthons
    /// are fingerprinted first if they have no fingerprints from it, and a
    /// space with another generator's is `Error::Runtime`. The generator's
    /// errors apply to the query. Mirrors `SynthonSpace.FingerprintSearch`.
    pub fn fingerprint_search(
        &mut self,
        query: &Mol,
        fingerprint_generator: &FingerprintGenerator64,
    ) -> Result<SubstructureResult> {
        self.fingerprint_search_with_params(
            query,
            fingerprint_generator,
            &SynthonSpaceSearchParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.FingerprintSearch(query, fingerprintGenerator,
    /// params)`. RDKit converts `similarity_cutoff -
    /// frag_similarity_adjuster` times the fingerprint size, and the size
    /// over it, to `unsigned int`, so it must be 0, or above 0 with both
    /// below 2^32, or the search is `Error::Value`; atom pairs need a query
    /// below 46333 atoms and products below 46341 (U3).
    pub fn fingerprint_search_with_params(
        &mut self,
        query: &Mol,
        fingerprint_generator: &FingerprintGenerator64,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_fingerprint_search(
                self.ptr.pin_mut(),
                query.as_ffi(),
                fingerprint_generator.as_ffi(),
                &params.to_ffi(),
            )?,
        ))
    }

    /// As `substructure_search_lines`, by fingerprint similarity; for atom
    /// pairs, a line of more than 46340 atoms is `Error::Value`. Mirrors
    /// `SynthonSpace.FingerprintSearch(query, fingerprintGenerator, params,
    /// startLine, finishLine)`.
    pub fn fingerprint_search_lines(
        &mut self,
        query: &Mol,
        fingerprint_generator: &FingerprintGenerator64,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_fingerprint_search_lines(
                self.ptr.pin_mut(),
                query.as_ffi(),
                fingerprint_generator.as_ffi(),
                &params.to_ffi(),
                start_line,
                finish_line,
            )?,
        ))
    }

    /// As `fingerprint_search`, handing the hits to `callback` (see the
    /// module documentation). Mirrors
    /// `SynthonSpace.FingerprintSearchIncremental`.
    pub fn fingerprint_search_incremental<F>(
        &mut self,
        query: &Mol,
        fingerprint_generator: &FingerprintGenerator64,
        callback: F,
    ) -> Result<()>
    where
        F: FnMut(Vec<Mol>) -> Result<bool>,
    {
        self.fingerprint_search_incremental_with_params(
            query,
            fingerprint_generator,
            callback,
            &SynthonSpaceSearchParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.FingerprintSearchIncremental(query,
    /// fingerprintGenerator, callback, params)`.
    pub fn fingerprint_search_incremental_with_params<F>(
        &mut self,
        query: &Mol,
        fingerprint_generator: &FingerprintGenerator64,
        callback: F,
        params: &SynthonSpaceSearchParams,
    ) -> Result<()>
    where
        F: FnMut(Vec<Mol>) -> Result<bool>,
    {
        let search_params = params.to_ffi();
        let space = self.ptr.pin_mut();
        incremental(callback, |lent| {
            // SAFETY: the callback is lent for this call only.
            unsafe {
                ffi::synthon_space_fingerprint_search_incremental(
                    space,
                    query.as_ffi(),
                    fingerprint_generator.as_ffi(),
                    lent,
                    &search_params,
                )
            }
        })
    }

    /// The products whose RASCAL similarity to the query reaches the
    /// options' `similarity_threshold` (the params' `similarity_cutoff` is
    /// not read). RDKit converts `similarity_threshold -
    /// frag_similarity_adjuster` times a fragment's atoms and bonds, and
    /// those over it, to `unsigned int`, so it must be above 0 with both
    /// below 2^32, or the search is `Error::Value`, as are a query of more
    /// than 46324 atoms and bonds and products of more than 46340 atoms
    /// (U3), and the options' own errors. Mirrors
    /// `SynthonSpace.RascalSearch`.
    pub fn rascal_search(
        &mut self,
        query: &Mol,
        rascal_options: &RascalOptions,
    ) -> Result<SubstructureResult> {
        self.rascal_search_with_params(query, rascal_options, &SynthonSpaceSearchParams::default())
    }

    /// Mirrors `SynthonSpace.RascalSearch(query, rascalOptions, params)`.
    pub fn rascal_search_with_params(
        &mut self,
        query: &Mol,
        rascal_options: &RascalOptions,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_rascal_search(
                self.ptr.pin_mut(),
                query.as_ffi(),
                &rascal_options.to_ffi(),
                &params.to_ffi(),
            )?,
        ))
    }

    /// As `substructure_search_lines`, by RASCAL similarity; a line of more
    /// than 46340 atoms is `Error::Value`, and so is `start_line` above
    /// `finish_line`, which RDKit refuses with `Error::Runtime`. Mirrors
    /// `SynthonSpace.RascalSearch(query, rascalOptions, params, startLine,
    /// finishLine)`.
    pub fn rascal_search_lines(
        &mut self,
        query: &Mol,
        rascal_options: &RascalOptions,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SubstructureResult> {
        Ok(SubstructureResult::from_ffi(
            ffi::synthon_space_rascal_search_lines(
                self.ptr.pin_mut(),
                query.as_ffi(),
                &rascal_options.to_ffi(),
                &params.to_ffi(),
                start_line,
                finish_line,
            )?,
        ))
    }

    /// As `rascal_search`, with RDKit's default options for `None`, handing
    /// the hits to `callback` (see the module documentation). Mirrors
    /// `SynthonSpace.RascalSearchIncremental`.
    pub fn rascal_search_incremental<F>(
        &mut self,
        query: &Mol,
        rascal_options: Option<&RascalOptions>,
        callback: F,
    ) -> Result<()>
    where
        F: FnMut(Vec<Mol>) -> Result<bool>,
    {
        self.rascal_search_incremental_with_params(
            query,
            rascal_options,
            callback,
            &SynthonSpaceSearchParams::default(),
        )
    }

    /// Mirrors `SynthonSpace.RascalSearchIncremental(query, rascalOptions,
    /// callback, params)`.
    pub fn rascal_search_incremental_with_params<F>(
        &mut self,
        query: &Mol,
        rascal_options: Option<&RascalOptions>,
        callback: F,
        params: &SynthonSpaceSearchParams,
    ) -> Result<()>
    where
        F: FnMut(Vec<Mol>) -> Result<bool>,
    {
        let options = rascal_options.cloned().unwrap_or_default().to_ffi();
        let search_params = params.to_ffi();
        let space = self.ptr.pin_mut();
        incremental(callback, |lent| {
            // SAFETY: the callback is lent for this call only.
            unsafe {
                ffi::synthon_space_rascal_search_incremental(
                    space,
                    query.as_ffi(),
                    &options,
                    lent,
                    &search_params,
                )
            }
        })
    }
}

/// Reads a synthon text file, as `SynthonSpace::read_text_file` does, and
/// writes it as a database file, as `SynthonSpace::write_db_file` does, with
/// their errors; Python's says it uses less memory than the two, but RDKit's
/// is the two as well. Mirrors `rdSynthonSpaceSearch.ConvertTextToDBFile`.
pub fn convert_text_to_db_file(
    in_filename: impl AsRef<Path>,
    out_filename: impl AsRef<Path>,
) -> Result<()> {
    convert_text_to_db_file_with_params(
        in_filename,
        out_filename,
        &ConvertTextToDBFileParams::default(),
    )
}

/// As `convert_text_to_db_file`, with the synthons fingerprinted, as
/// `SynthonSpace::build_synthon_fingerprints` does, and their shapes built,
/// as `SynthonSpace::build_synthon_shapes_with_py_params` does, before the
/// file is written. Mirrors `rdSynthonSpaceSearch.ConvertTextToDBFile(
/// inFilename, outFilename, fpGen, py_shapeParams)`.
pub fn convert_text_to_db_file_with_params(
    in_filename: impl AsRef<Path>,
    out_filename: impl AsRef<Path>,
    params: &ConvertTextToDBFileParams<'_>,
) -> Result<()> {
    let mut space = SynthonSpace::default();
    space.read_text_file(in_filename)?;
    if let Some(fp_gen) = params.fp_gen {
        space.build_synthon_fingerprints(fp_gen)?;
    }
    if let Some(shape_params) = &params.py_shape_params {
        space.build_synthon_shapes_with_py_params(shape_params)?;
    }
    space.write_db_file(out_filename)
}

/// The integer with a space between each group of 3 digits (`"1 234
/// 567"`). Mirrors `rdSynthonSpaceSearch.FormattedIntegerString`.
pub fn formatted_integer_string(value: i64) -> String {
    to_string(expect_ok(ffi::formatted_integer_string(value)))
        .expect("a formatted integer is not UTF-8")
}
