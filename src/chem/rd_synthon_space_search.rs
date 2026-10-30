// Copyright (c) 2026 Genesis Molecular AI, Inc.
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
//! depend on the thread count. RDKit's handler stops a substructure or
//! fingerprint search, or a text read, at a Ctrl-C, which then returns
//! `Error::Runtime`, as Python raises `KeyboardInterrupt`; the RASCAL search
//! installs no handler.
//!
//! The reads and the searches remove atoms and bonds from the synthons,
//! the products they combine, and the query's fragments, so what RDKit
//! would misread then is refused with `Error::Value`: a bond's
//! `_MolFileBondEndPts`, or substance groups whose `PARENT`s form a cycle
//! (RDKit follows it forever). A reaction's synthons' chains are checked
//! together, as a product carries them all.
//!
//! RDKit's fingerprint build stops between reactions at a Ctrl-C, leaving a
//! space whose fingerprint search would dereference the missing
//! fingerprints. Here the synthons are fingerprinted without RDKit's Ctrl-C
//! handler, so a Ctrl-C does not cut the build short, and a space whose
//! fingerprints are incomplete anyway (a build that failed partway) is
//! `Error::Value` for `build_synthon_fingerprints` and `fingerprint_search`.
//!
//! RDKit 2025.03 has no shape search, incremental searches, or searches of
//! a possible-hits file; the progress bar of `BuildSynthonFingerprints` is
//! not bound (the manifest says why).

use std::fmt;
use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rd_fingerprint_generator::FingerprintGenerator64;
use crate::chem::rd_generalized_substruct::ExtendedQueryMol;
use crate::chem::rd_rascal_mces::RascalOptions;
use crate::chem::rdchem::SubstructMatchParameters;
use crate::chem::rdmolops::mol_list_into_vec;
use crate::error::{expect_ok, to_string};
use crate::ffi::synthon_space_search::ffi;
use crate::util::path_bytes;

/// Parameters of the searches. Mirrors
/// `rdSynthonSpaceSearch.SynthonSpaceSearchParams` (A3); `Default` reads
/// RDKit's defaults.
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
    /// Not read by RDKit 2025.03's searches.
    pub num_random_sweeps: i32,
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
    /// Not read: the search runs on the calling thread (see the module
    /// documentation).
    pub num_threads: i32,
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
            num_random_sweeps: d.num_random_sweeps,
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
            num_threads: d.num_threads,
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
            num_random_sweeps: self.num_random_sweeps,
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
            num_threads: self.num_threads,
        }
    }
}

/// Parameters for `convert_text_to_db_file_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default)]
pub struct ConvertTextToDBFileParams<'a> {
    /// Fingerprints the synthons with this generator.
    pub fp_gen: Option<&'a FingerprintGenerator64>,
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
            .field("num_products", &self.get_num_products())
            .finish()
    }
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
    /// from 3.0 to 3.1 (those this RDKit reads), and record offsets that
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
    /// records' molecule and bit-vector pickles, counts, and bools without
    /// validating them.
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
    /// `Error::Value` for a synthon a pickle cannot hold (`Mol::to_binary`).
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

    /// Writes a summary of the space to standard output (C++'s
    /// `std::cout`). Mirrors `SynthonSpace.Summarise`.
    pub fn summarise(&self) {
        expect_ok(ffi::synthon_space_summarise(&self.ptr));
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

    /// Fingerprints the synthons with the generator, as a fingerprint search
    /// with it does, unless they already have its fingerprints. For atom
    /// pairs, a space whose products can pass 46340 atoms is `Error::Value`
    /// (U3), as is a space whose fingerprints are incomplete (see the module
    /// documentation). Mirrors `SynthonSpace.BuildSynthonFingerprints`.
    pub fn build_synthon_fingerprints(
        &mut self,
        fingerprint_generator: &FingerprintGenerator64,
    ) -> Result<()> {
        Ok(ffi::synthon_space_build_synthon_fingerprints(
            self.ptr.pin_mut(),
            fingerprint_generator.as_ffi(),
        )?)
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

    /// The products whose fingerprints, from the generator, are at least
    /// `similarity_cutoff` similar to the query's, by Tanimoto; the synthons
    /// are fingerprinted first, as `build_synthon_fingerprints` does, if they
    /// have no fingerprints, and a space with another generator's is
    /// `Error::Runtime`. The generator's errors apply to the query. Mirrors
    /// `SynthonSpace.FingerprintSearch`.
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
    /// params)`. Atom pairs need a query below 46333 atoms and products
    /// below 46341 (U3), or the search is `Error::Value`.
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

    /// The products whose RASCAL similarity to the query reaches the
    /// options' `similarity_threshold` (the params' `similarity_cutoff` is
    /// not read); its fragments' threshold is `similarity_threshold -
    /// frag_similarity_adjuster`. A NaN threshold (RASCAL casts its product
    /// to `unsigned int`) is `Error::Value`, as are a query of more than 46324
    /// atoms and bonds and products of more than 46340 atoms (U3), and the
    /// options' own errors. Mirrors `SynthonSpace.RascalSearch`.
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
/// `SynthonSpace::build_synthon_fingerprints` does, before the file is
/// written. Mirrors `rdSynthonSpaceSearch.ConvertTextToDBFile(inFilename,
/// outFilename, fpGen)`.
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
    space.write_db_file(out_filename)
}

/// The integer with a space between each group of 3 digits (`"1 234
/// 567"`). Mirrors `rdSynthonSpaceSearch.FormattedIntegerString`.
pub fn formatted_integer_string(value: i64) -> String {
    to_string(expect_ok(ffi::formatted_integer_string(value)))
        .expect("a formatted integer is not UTF-8")
}
