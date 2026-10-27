// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Rust bindings for RDKit's C++ API, with Rust ports for much of the
//! functionality RDKit currently offers only in Python. The API mirrors
//! RDKit's Python interface in Rust spelling, so Python code ports line by
//! line.
//!
//! ```no_run
//! use rdkit_cxx::chem::rd_fingerprint_generator::get_morgan_generator;
//! use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles, mol_to_smiles};
//! use rdkit_cxx::data_structs::tanimoto_similarity;
//!
//! let mol = mol_from_smiles("Oc1ccccc1O")?;
//! let hydroxyl = mol_from_smarts("[OX2H]")?;
//! println!("{} {:?}", mol_to_smiles(&mol)?, mol.get_substruct_matches(&hydroxyl)?);
//!
//! let morgan = get_morgan_generator();
//! let phenol = mol_from_smiles("Oc1ccccc1")?;
//! let similarity = tanimoto_similarity(
//!     &morgan.get_fingerprint(&mol)?,
//!     &morgan.get_fingerprint(&phenol)?,
//! )?;
//! println!("{similarity:.3}");
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```
//!
//! The crate compiles a C++ wrapper layer against RDKit from conda-forge
//! and links RDKit's shared libraries; it reimplements none of RDKit's C++.
//! What RDKit implements only in Python (QED, BRICS, `Draw`, Butina
//! clustering, and more) is ported to Rust on top of the bindings, and
//! computes what Python computes. Every name in RDKit 2025.09's Python
//! interface has a Rust counterpart or a recorded reason it has none.
//!
//! # From Python to Rust
//!
//! - **Paths** follow Python's, each component in snake_case:
//!   `rdkit.Chem.rdmolfiles.MolFromSmiles` is
//!   [`chem::rdmolfiles::mol_from_smiles`], also re-exported from [`chem`]
//!   as `rdkit.Chem` re-exports it, and `rdkit.Chem.QED` is [`chem::qed`].
//!   Classes keep their Python names ([`chem::Mol`], not `ROMol`).
//! - **Optional arguments** become `_with_` functions.
//!   `GetNumAtoms(onlyExplicit=False)` is
//!   `get_num_atoms_with_only_explicit(false)`, and a call with several
//!   optional arguments takes a params struct with Python's defaults:
//!
//!   ```no_run
//!   # use rdkit_cxx::chem::{SmilesWriteParams, mol_from_smiles, mol_to_smiles_with_params};
//!   # let mol = mol_from_smiles("c1ccccc1")?;
//!   let kekule = mol_to_smiles_with_params(
//!       &mol,
//!       &SmilesWriteParams { do_kekule: true, ..Default::default() },
//!   )?;
//!   # Ok::<(), rdkit_cxx::Error>(())
//!   ```
//!
//!   Params structs gain fields only with a new RDKit series, never in a
//!   patch release, so that idiom keeps compiling. Enums, and the structs
//!   functions return, are `#[non_exhaustive]`.
//! - **Errors:** exceptions become [`Result<T, Error>`](Result). A parser
//!   that returns `None` in Python returns [`Error::Parse`], and
//!   sanitization errors keep their atom indices
//!   ([`Error::AtomValence`]`{ atom_idx, .. }`).
//! - **Views and `RWMol`:** atoms, bonds, and conformers reached through a
//!   molecule are views that borrow it ([`chem::AtomRef`],
//!   [`chem::AtomMut`]). Editing goes through `&mut Mol`, one atom at a
//!   time by index, or through [`chem::RWMol`], which dereferences to
//!   `Mol` and converts to and from it without copying:
//!
//!   ```no_run
//!   use rdkit_cxx::chem::{Atom, RWMol, mol_from_smiles};
//!
//!   let mut mol = mol_from_smiles("CCO")?;
//!   mol.get_atom_with_idx_mut(2)?.set_formal_charge(-1);
//!   let mut rw = RWMol::from(mol);
//!   let n = rw.add_atom(&Atom::new("N")?);
//!   rw.add_bond(0, n)?;
//!   for atom in rw.get_atoms() {
//!       println!("{} {}", atom.get_symbol()?, atom.get_degree());
//!   }
//!   # Ok::<(), rdkit_cxx::Error>(())
//!   ```
//! - **Threads:** [`chem::Mol`] is `Send` but not `Sync`. RDKit updates
//!   cached data inside molecules it treats as read-only, so parallel code
//!   gives each thread its own molecule.
//! - **Callbacks** Python takes as any callable take closures: `impl Fn`
//!   arguments, or `&dyn Fn` fields in params structs.
//!
//! # Building and linking
//!
//! Depend on the crate with a tilde requirement, `rdkit-cxx = "~2025.9"`,
//! which stays within one RDKit series. Each crate version supports one
//! RDKit release series, and within it the releases it has checked: the
//! build finds RDKit through `CONDA_PREFIX` (or `RDKIT_PREFIX`), in an
//! environment with conda-forge's `librdkit-dev`, `libboost-devel`,
//! `cxx-compiler`, and `cmake`, and refuses any other release.
//!
//! RDKit's shared libraries stay in that environment, so a binary that
//! depends on this crate finds them through an rpath its own build script
//! adds; Cargo gives `DEP_RDKIT_CXX_LIBDIR` only to packages that depend on
//! rdkit-cxx directly:
//!
//! ```ignore
//! // build.rs
//! fn main() {
//!     // Absent where nothing links, as when docs.rs builds the documentation.
//!     if let Ok(dir) = std::env::var("DEP_RDKIT_CXX_LIBDIR") {
//!         println!("cargo::rustc-link-arg=-Wl,-rpath,{dir}");
//!     }
//! }
//! ```
//!
//! The [README](https://github.com/skearnes/rdkit-cxx#readme) has the rest:
//! the environment, how the wrappers work, the known limitations, and how
//! coverage is measured.
//!
//! # Modules
//!
//! The modules follow RDKit's Python packages. The ones most code starts
//! from:
//!
//! - [`chem::rdchem`]: molecules, atoms, bonds, conformers, and their
//!   views; properties; substructure matching; the periodic table.
//! - [`chem::rdmolfiles`]: reading and writing SMILES, SMARTS, mol blocks,
//!   and other formats, and the SD, SMILES, and TDT suppliers and writers.
//! - [`chem::rdmolops`]: hydrogens, sanitization, kekulization,
//!   aromaticity, fragments, stereochemistry, and the classic fingerprints.
//! - [`data_structs`]: bit and count vectors, and the similarity metrics
//!   on them.
//! - [`ml`]: ports of `rdkit.ML`: clustering, data sets, information
//!   theory, and scoring.
//! - [`chem::draw`]: drawing molecules and reactions to SVG and PNG, and
//!   similarity maps.
//!
//! Others mirror RDKit's modules of the same names, such as
//! [`chem::rd_fingerprint_generator`], [`chem::rd_dist_geom`],
//! [`chem::rd_chem_reactions`], and [`chem::mol_standardize`].

#![warn(missing_docs)]

pub mod avalon;
pub mod chem;
pub mod data_manip;
pub mod data_structs;
pub mod distance_geometry;
pub mod error;
mod ffi;
mod flags;
pub mod force_field;
pub mod geometry;
pub mod ml;
pub mod numerics;
mod py;
pub mod rd_base;
pub mod rd_logger;
pub mod rd_paths;
pub mod rd_random;
pub mod sim_div_filters;
mod util;
pub mod utils;

pub use error::{Error, Result};
