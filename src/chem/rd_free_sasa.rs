// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdFreeSASA`: solvent accessible surface areas from
//! the FreeSASA library.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_pdb_file;
//! use rdkit_cxx::chem::rd_free_sasa::{calc_sasa, classify_atoms};
//!
//! let mut mol = mol_from_pdb_file("protein.pdb").unwrap();
//! let radii = classify_atoms(&mut mol);
//! println!("{}", calc_sasa(&mut mol, &radii).unwrap());
//! ```
//!
//! FreeSASA trusts what it is given, so `calc_sasa` refuses with
//! `Error::Value` what would make it read or write out of bounds: radii
//! that are not one per atom, coordinates that are not finite, a largest
//! radius plus probe radius that is not positive and finite, and
//! coordinates spread over more of FreeSASA's cells (each twice that
//! largest radius wide) than `i32` counts, or too large for cells that
//! small. `classify_atoms` classifies an atom whose PDB residue or atom
//! name is blank, which FreeSASA would compare uninitialized memory with,
//! as an unknown name.

use crate::Result;
use crate::chem::rdchem::{AsAtom, Atom, Mol};
use crate::error::expect_ok;
use crate::ffi::free_sasa::ffi;

/// The algorithm FreeSASA computes areas with. Mirrors
/// `rdFreeSASA.SASAAlgorithm`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SASAAlgorithm {
    /// Lee and Richards' slices. Mirrors `LeeRichards`.
    LeeRichards,
    /// Shrake and Rupley's test points. Mirrors `ShrakeRupley`.
    ShrakeRupley,
}

impl SASAAlgorithm {
    fn to_ffi(self) -> i32 {
        match self {
            SASAAlgorithm::LeeRichards => 0,
            SASAAlgorithm::ShrakeRupley => 1,
        }
    }

    fn from_ffi(value: i32) -> SASAAlgorithm {
        match value {
            0 => SASAAlgorithm::LeeRichards,
            1 => SASAAlgorithm::ShrakeRupley,
            _ => panic!("unknown SASAAlgorithm value {value}"),
        }
    }
}

/// The classifier that gives atoms of PDB residues their radii and
/// classes. Mirrors `rdFreeSASA.SASAClassifier`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SASAClassifier {
    /// Mirrors `Protor`.
    Protor,
    /// Mirrors `NACCESS`.
    Naccess,
    /// Mirrors `OONS`.
    Oons,
}

impl SASAClassifier {
    fn to_ffi(self) -> i32 {
        match self {
            SASAClassifier::Protor => 0,
            SASAClassifier::Naccess => 1,
            SASAClassifier::Oons => 2,
        }
    }

    fn from_ffi(value: i32) -> SASAClassifier {
        match value {
            0 => SASAClassifier::Protor,
            1 => SASAClassifier::Naccess,
            2 => SASAClassifier::Oons,
            _ => panic!("unknown SASAClassifier value {value}"),
        }
    }
}

/// An atom's class. Mirrors `rdFreeSASA.SASAClass`.
///
/// `classify_atoms` records FreeSASA's own class numbers in the atoms'
/// `SASAClass` property (0 apolar, 1 polar, 2 unknown), not these; the
/// `SASAClassName` property holds the name (`Apolar`, `Polar`, or
/// `Unclassified`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SASAClass {
    /// Mirrors `Unclassified` (0).
    Unclassified,
    /// Mirrors `APolar` (1).
    APolar,
    /// Mirrors `Polar` (2).
    Polar,
}

/// How FreeSASA classifies atoms and computes areas. Mirrors
/// `rdFreeSASA.SASAOpts` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct SASAOpts {
    pub algorithm: SASAAlgorithm,
    pub classifier: SASAClassifier,
    /// The radius of the solvent probe, in Å.
    pub probe_radius: f64,
}

impl Default for SASAOpts {
    /// Lee and Richards' algorithm, the ProtOr classifier, and a probe of
    /// radius 1.4. Mirrors `SASAOpts()`.
    fn default() -> SASAOpts {
        let data = expect_ok(ffi::sasa_opts_default());
        SASAOpts {
            algorithm: SASAAlgorithm::from_ffi(data.algorithm),
            classifier: SASAClassifier::from_ffi(data.classifier),
            probe_radius: data.probe_radius,
        }
    }
}

impl SASAOpts {
    /// Options with the default probe radius. Mirrors `SASAOpts(alg, cls)`.
    pub fn new(alg: SASAAlgorithm, cls: SASAClassifier) -> SASAOpts {
        SASAOpts {
            algorithm: alg,
            classifier: cls,
            ..SASAOpts::default()
        }
    }

    /// Mirrors `SASAOpts(alg, cls, pr)`, with probe radius `pr`.
    pub fn new_with_pr(alg: SASAAlgorithm, cls: SASAClassifier, pr: f64) -> SASAOpts {
        SASAOpts {
            algorithm: alg,
            classifier: cls,
            probe_radius: pr,
        }
    }

    fn to_ffi(&self) -> ffi::SASAOptsData {
        ffi::SASAOptsData {
            algorithm: self.algorithm.to_ffi(),
            classifier: self.classifier.to_ffi(),
            probe_radius: self.probe_radius,
        }
    }
}

/// Classifies the atoms with the ProtOr classifier and returns their radii,
/// indexed by atom, or none when some atom could not be classified. Mirrors
/// `rdFreeSASA.classifyAtoms`.
///
/// Only atoms with PDB residue information are classified; the others get
/// radius 0. Each atom's `SASAClass` and `SASAClassName` properties record
/// its class.
pub fn classify_atoms(mol: &mut Mol) -> Vec<f64> {
    classify_atoms_with_options(mol, &SASAOpts::default())
}

/// Mirrors `rdFreeSASA.classifyAtoms(mol, options)`, which uses the
/// options' classifier.
pub fn classify_atoms_with_options(mol: &mut Mol, options: &SASAOpts) -> Vec<f64> {
    expect_ok(ffi::classify_atoms(mol.pin_mut(), &options.to_ffi()))
}

/// Parameters for `calc_sasa_with_params`, with Python's defaults (A2).
#[derive(Clone)]
pub struct CalcSASAParams<'a> {
    /// The conformer; -1 for the first.
    pub conf_idx: i32,
    /// A query atom that limits the result to the atoms it matches, such as
    /// `make_free_sasa_polar_atom_query`'s; `None` for every atom.
    pub query: Option<&'a dyn AsAtom>,
    pub opts: SASAOpts,
}

impl Default for CalcSASAParams<'_> {
    fn default() -> Self {
        CalcSASAParams {
            conf_idx: -1,
            query: None,
            opts: SASAOpts::default(),
        }
    }
}

/// The solvent accessible surface area of the molecule's first conformer,
/// with one radius per atom (from `classify_atoms`, for example), by Lee
/// and Richards' algorithm with a probe of radius 1.4. Mirrors
/// `rdFreeSASA.CalcSASA`.
///
/// Sets each atom's `SASA` property to its area, and the molecule's to the
/// total, so it takes the molecule mutably (S4). A molecule without atoms or conformers is `Error::Runtime`, a
/// missing conformer `Error::Value`, and so are the inputs the module
/// documentation lists. FreeSASA failing for want of memory gives 0.
pub fn calc_sasa(mol: &mut Mol, radii: &[f64]) -> Result<f64> {
    calc_sasa_with_params(mol, radii, &CalcSASAParams::default())
}

/// Mirrors `rdFreeSASA.CalcSASA(mol, radii, confIdx, query, opts)`. With a
/// query, the result is the sum of the areas of the atoms it matches; an
/// atom that is not a query atom is `Error::Value`.
pub fn calc_sasa_with_params(
    mol: &mut Mol,
    radii: &[f64],
    params: &CalcSASAParams<'_>,
) -> Result<f64> {
    let query = params
        .query
        .map_or(std::ptr::null(), |atom| atom.atom() as *const _);
    // SAFETY: the query is null or borrowed for the call.
    Ok(unsafe {
        ffi::calc_sasa(
            mol.pin_mut(),
            radii,
            params.conf_idx,
            query,
            &params.opts.to_ffi(),
        )
    }?)
}

/// A query atom that matches the atoms `classify_atoms` found apolar.
/// Mirrors `rdFreeSASA.MakeFreeSasaAPolarAtomQuery`.
pub fn make_free_sasa_a_polar_atom_query() -> Atom {
    Atom::from_ffi(expect_ok(ffi::make_free_sasa_a_polar_atom_query()))
}

/// A query atom that matches the atoms `classify_atoms` found polar.
/// Mirrors `rdFreeSASA.MakeFreeSasaPolarAtomQuery`.
pub fn make_free_sasa_polar_atom_query() -> Atom {
    Atom::from_ffi(expect_ok(ffi::make_free_sasa_polar_atom_query()))
}
