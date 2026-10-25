// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! `CleanupParameters` and the module functions that take them.

use std::path::PathBuf;

use crate::chem::Mol;
use crate::chem::rdmolops::with_borrowed_mols;
use crate::error::{expect_ok, to_string};
use crate::ffi::standardize::ffi;
use crate::util::{path_bytes, path_from_bytes};
use crate::{Error, Result};

/// Parameters controlling molecular standardization. Mirrors
/// `rdMolStandardize.CleanupParameters`.
///
/// An empty file name chooses RDKit's built-in definitions. The four data
/// lists are what `update_params_from_json` reads from `normalizationData`,
/// `fragmentData`, `acidbaseData`, and `tautomerTransformData`: Python keeps
/// them on the C++ object, with no attributes. A list that is not empty
/// replaces the matching file.
#[derive(Clone, Debug, PartialEq)]
pub struct CleanupParameters {
    /// The file of normalization transformations.
    pub normalizations_file: PathBuf,
    /// The file of acid and base definitions.
    pub acidbase_file: PathBuf,
    /// The file of fragments to remove.
    pub fragment_file: PathBuf,
    /// The file of tautomer transformations.
    pub tautomer_transforms_file: PathBuf,
    /// The most times to restart the series of normalizations.
    pub max_restarts: i32,
    /// Prefers organic fragments to inorganic ones when choosing which to
    /// keep.
    pub prefer_organic: bool,
    /// Applies atom-order dependent normalizations (like uncharging) in a
    /// canonical order.
    pub do_canonical: bool,
    /// The most tautomers to generate.
    pub max_tautomers: i32,
    /// The most transforms to apply while enumerating tautomers.
    pub max_transforms: i32,
    /// Removes stereochemistry from sp3 centers involved in tautomerism.
    pub tautomer_remove_sp3_stereo: bool,
    /// Removes stereochemistry from double bonds involved in tautomerism.
    pub tautomer_remove_bond_stereo: bool,
    /// Removes isotopic hydrogens from centers involved in tautomerism.
    pub tautomer_remove_isotopic_hs: bool,
    /// Assigns stereochemistry again on every generated tautomer.
    pub tautomer_reassign_stereo: bool,
    /// Has `LargestFragmentChooser` compare atom counts before masses.
    pub largest_fragment_chooser_use_atom_count: bool,
    /// Has `LargestFragmentChooser` count only heavy atoms.
    pub largest_fragment_chooser_count_heavy_atoms_only: bool,
    /// RDKit's installation directory, from the `RDBASE` environment
    /// variable; no standardization reads it.
    pub rdbase: PathBuf,
    /// Normalizations, as (name, reaction SMARTS).
    pub normalization_data: Vec<(String, String)>,
    /// Fragments to remove, as (name, SMARTS).
    pub fragment_data: Vec<(String, String)>,
    /// Acid-base pairs, as (name, acid SMARTS, base SMARTS).
    pub acidbase_data: Vec<(String, String, String)>,
    /// Tautomer transformations, as (name, SMARTS, bonds, charges): the
    /// bonds (`-`, `=`, `#`, `:`) to give the SMARTS's bonds, in order, and
    /// the charges (`+`, `0`, `-`) to add to its atoms, either one empty.
    pub tautomer_transform_data: Vec<(String, String, String, String)>,
}

/// Mirrors `CleanupParameters()`.
impl Default for CleanupParameters {
    fn default() -> CleanupParameters {
        CleanupParameters::from_ffi(expect_ok(ffi::cleanup_parameters_default()))
            .expect("RDKit's default CleanupParameters are UTF-8")
    }
}

fn utf8(bytes: Vec<u8>) -> Result<String> {
    String::from_utf8(bytes).map_err(|err| Error::Value {
        message: format!("invalid UTF-8 in CleanupParameters: {err}"),
    })
}

fn entry(first: &str, second: &str, third: &str, fourth: &str) -> ffi::CleanupDataEntry {
    ffi::CleanupDataEntry {
        name: first.as_bytes().to_vec(),
        first: second.as_bytes().to_vec(),
        second: third.as_bytes().to_vec(),
        third: fourth.as_bytes().to_vec(),
    }
}

impl CleanupParameters {
    pub(crate) fn from_ffi(data: ffi::CleanupParametersData) -> Result<CleanupParameters> {
        let pairs = |entries: Vec<ffi::CleanupDataEntry>| -> Result<Vec<(String, String)>> {
            entries
                .into_iter()
                .map(|e| Ok((utf8(e.name)?, utf8(e.first)?)))
                .collect()
        };
        Ok(CleanupParameters {
            normalizations_file: path_from_bytes(data.normalizations),
            acidbase_file: path_from_bytes(data.acidbase_file),
            fragment_file: path_from_bytes(data.fragment_file),
            tautomer_transforms_file: path_from_bytes(data.tautomer_transforms),
            max_restarts: data.max_restarts,
            prefer_organic: data.prefer_organic,
            do_canonical: data.do_canonical,
            max_tautomers: data.max_tautomers,
            max_transforms: data.max_transforms,
            tautomer_remove_sp3_stereo: data.tautomer_remove_sp3_stereo,
            tautomer_remove_bond_stereo: data.tautomer_remove_bond_stereo,
            tautomer_remove_isotopic_hs: data.tautomer_remove_isotopic_hs,
            tautomer_reassign_stereo: data.tautomer_reassign_stereo,
            largest_fragment_chooser_use_atom_count: data.largest_fragment_chooser_use_atom_count,
            largest_fragment_chooser_count_heavy_atoms_only: data
                .largest_fragment_chooser_count_heavy_atoms_only,
            rdbase: path_from_bytes(data.rdbase),
            normalization_data: pairs(data.normalization_data)?,
            fragment_data: pairs(data.fragment_data)?,
            acidbase_data: data
                .acidbase_data
                .into_iter()
                .map(|e| Ok((utf8(e.name)?, utf8(e.first)?, utf8(e.second)?)))
                .collect::<Result<_>>()?,
            tautomer_transform_data: data
                .tautomer_transform_data
                .into_iter()
                .map(|e| {
                    Ok((
                        utf8(e.name)?,
                        utf8(e.first)?,
                        utf8(e.second)?,
                        utf8(e.third)?,
                    ))
                })
                .collect::<Result<_>>()?,
        })
    }

    pub(crate) fn to_ffi(&self) -> ffi::CleanupParametersData {
        let bytes = |path: &PathBuf| path_bytes(path).to_vec();
        ffi::CleanupParametersData {
            rdbase: bytes(&self.rdbase),
            normalizations: bytes(&self.normalizations_file),
            acidbase_file: bytes(&self.acidbase_file),
            fragment_file: bytes(&self.fragment_file),
            tautomer_transforms: bytes(&self.tautomer_transforms_file),
            max_restarts: self.max_restarts,
            prefer_organic: self.prefer_organic,
            do_canonical: self.do_canonical,
            max_tautomers: self.max_tautomers,
            max_transforms: self.max_transforms,
            tautomer_remove_sp3_stereo: self.tautomer_remove_sp3_stereo,
            tautomer_remove_bond_stereo: self.tautomer_remove_bond_stereo,
            tautomer_remove_isotopic_hs: self.tautomer_remove_isotopic_hs,
            tautomer_reassign_stereo: self.tautomer_reassign_stereo,
            largest_fragment_chooser_use_atom_count: self.largest_fragment_chooser_use_atom_count,
            largest_fragment_chooser_count_heavy_atoms_only: self
                .largest_fragment_chooser_count_heavy_atoms_only,
            normalization_data: self
                .normalization_data
                .iter()
                .map(|(name, smarts)| entry(name, smarts, "", ""))
                .collect(),
            fragment_data: self
                .fragment_data
                .iter()
                .map(|(name, smarts)| entry(name, smarts, "", ""))
                .collect(),
            acidbase_data: self
                .acidbase_data
                .iter()
                .map(|(name, acid, base)| entry(name, acid, base, ""))
                .collect(),
            tautomer_transform_data: self
                .tautomer_transform_data
                .iter()
                .map(|(name, smarts, bonds, charges)| entry(name, smarts, bonds, charges))
                .collect(),
        }
    }
}

/// Updates the parameters from a JSON object, as RDKit reads it: the keys
/// are the C++ field names (`normalizations`, `acidbaseFile`, `maxRestarts`,
/// ...; not the `largestFragmentChooser` options), a value that does not
/// convert keeps the current one, and the entries of the data lists
/// (`normalizationData`, `fragmentData`, `acidbaseData`, and
/// `tautomerTransformData`, each a list of objects) are appended to the
/// current ones, skipping any without a name or SMARTS. Empty JSON changes
/// nothing. `Error::Runtime` for JSON that does not parse, and
/// `Error::Value` for strings that are not UTF-8 (the parser decodes `\u`
/// escapes, lone surrogates among them). Mirrors
/// `rdMolStandardize.UpdateParamsFromJSON`.
pub fn update_params_from_json(params: &mut CleanupParameters, json: &str) -> Result<()> {
    let data = ffi::update_params_from_json(&params.to_ffi(), json)?;
    *params = CleanupParameters::from_ffi(data)?;
    Ok(())
}

/// Parameters for the parents, with Python's defaults (A2).
macro_rules! parent_params {
    ($($python:literal $name:ident),* $(,)?) => {
        $(
            #[doc = concat!("Parameters for `", stringify!($name), "`'s function, with Python's defaults (A2). Mirrors `rdMolStandardize.", $python, "`'s `params` and `skipStandardize`.")]
            #[derive(Clone, Debug, Default, PartialEq)]
            pub struct $name {
                /// The parameters; Python's `None` is their default.
                pub params: CleanupParameters,
                /// Skips the cleanup that comes first.
                pub skip_standardize: bool,
            }
        )*
    };
}

parent_params!(
    "TautomerParent" TautomerParentParams,
    "TautomerParentInPlace" TautomerParentInPlaceParams,
    "FragmentParent" FragmentParentParams,
    "FragmentParentInPlace" FragmentParentInPlaceParams,
    "StereoParent" StereoParentParams,
    "StereoParentInPlace" StereoParentInPlaceParams,
    "IsotopeParent" IsotopeParentParams,
    "IsotopeParentInPlace" IsotopeParentInPlaceParams,
    "ChargeParent" ChargeParentParams,
    "ChargeParentInPlace" ChargeParentInPlaceParams,
    "SuperParent" SuperParentParams,
    "SuperParentInPlace" SuperParentInPlaceParams,
);

fn copy(
    op: ffi::StandardizeOp,
    mol: &Mol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::standardize_copy(
        op,
        mol.as_ffi(),
        &params.to_ffi(),
        skip_standardize,
    )?))
}

fn in_place(
    op: ffi::StandardizeOp,
    mol: &mut Mol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    Ok(ffi::standardize_in_place(
        op,
        mol.pin_mut(),
        &params.to_ffi(),
        skip_standardize,
    )?)
}

fn in_place_mols(
    op: ffi::StandardizeOp,
    mols: &mut [Mol],
    num_threads: i32,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    let data = params.to_ffi();
    Ok(with_borrowed_mols(
        mols.iter_mut().map(|mol| ("", mol)),
        |list| ffi::standardize_in_place_mols(op, list, num_threads, &data, skip_standardize),
    )?)
}

/// The functions that take only `CleanupParameters`: a copy, in place, and
/// for several molecules.
macro_rules! cleanup_functions {
    (
        $op:ident, $python:literal, $doc:literal, $errors:literal:
        $copy:ident, $copy_params:ident, $in_place:ident, $in_place_params:ident,
        $mols:ident, $mols_params:ident
    ) => {
        #[doc = concat!($doc, " Mirrors `rdMolStandardize.", $python, "`.")]
        #[doc = ""]
        #[doc = $errors]
        pub fn $copy(mol: &Mol) -> Result<Mol> {
            $copy_params(mol, &CleanupParameters::default())
        }

        #[doc = concat!("Mirrors `rdMolStandardize.", $python, "(mol, params)`.")]
        pub fn $copy_params(mol: &Mol, params: &CleanupParameters) -> Result<Mol> {
            copy(ffi::StandardizeOp::$op, mol, params, false)
        }

        #[doc = concat!("As `", stringify!($copy), "`, in place. Mirrors `rdMolStandardize.", $python, "InPlace`.")]
        pub fn $in_place(mol: &mut Mol) -> Result<()> {
            $in_place_params(mol, &CleanupParameters::default())
        }

        #[doc = concat!("Mirrors `rdMolStandardize.", $python, "InPlace(mol, params)`.")]
        pub fn $in_place_params(mol: &mut Mol, params: &CleanupParameters) -> Result<()> {
            in_place(ffi::StandardizeOp::$op, mol, params, false)
        }

        #[doc = concat!("As `", stringify!($in_place), "`, for several molecules. Mirrors `rdMolStandardize.", $python, "InPlace(mols, numThreads)`.")]
        #[doc = ""]
        #[doc = "The molecules are split among `num_threads` threads (0 or less uses all but \
                 that many of the cores), each standardized by one. When some fail, the others \
                 may still have been standardized, and one of the errors is returned."]
        pub fn $mols(mols: &mut [Mol], num_threads: i32) -> Result<()> {
            $mols_params(mols, num_threads, &CleanupParameters::default())
        }

        #[doc = concat!("Mirrors `rdMolStandardize.", $python, "InPlace(mols, numThreads, params)`.")]
        pub fn $mols_params(
            mols: &mut [Mol],
            num_threads: i32,
            params: &CleanupParameters,
        ) -> Result<()> {
            in_place_mols(ffi::StandardizeOp::$op, mols, num_threads, params, false)
        }
    };
}

cleanup_functions!(
    Cleanup, "Cleanup",
    "Standardizes a molecule: removes hydrogens, disconnects metals, normalizes, reionizes, and \
     assigns stereochemistry, as MolVS's `Standardizer.standardize` does.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     `Error::Value` for a molecule of more than `i32::MAX / 255` atoms, whose charge sums could \
     overflow RDKit's `int`s (CONVENTIONS.md U9).":
    cleanup, cleanup_with_params, cleanup_in_place, cleanup_in_place_with_params,
    cleanup_in_place_mols, cleanup_in_place_mols_with_params
);

cleanup_functions!(
    Normalize, "Normalize",
    "Applies a series of standard transformations to correct functional groups and recombine \
     charges.",
    "The copy normalizes each fragment, keeping each transformation's first product by \
     SMILES. In place, the transformations must have one reactant and no more product atoms \
     than reactant atoms (`Error::Value`). RDKit's errors (a file that cannot be read, a \
     reaction that does not parse) are returned.":
    normalize, normalize_with_params, normalize_in_place, normalize_in_place_with_params,
    normalize_in_place_mols, normalize_in_place_mols_with_params
);

cleanup_functions!(
    Reionize, "Reionize",
    "Ensures the strongest acid groups are charged first.",
    "RDKit's errors (a file that cannot be read, SMARTS that do not parse) are returned, and \
     `Error::Value` for a molecule of more than `i32::MAX / 255` atoms, whose charge sums could \
     overflow RDKit's `int`s (CONVENTIONS.md U9).":
    reionize, reionize_with_params, reionize_in_place, reionize_in_place_with_params,
    reionize_in_place_mols, reionize_in_place_mols_with_params
);

cleanup_functions!(
    RemoveFragments, "RemoveFragments",
    "Removes the fragments the fragment definitions match, leaving the last one when all \
     match.",
    "RDKit's errors (a file that cannot be read, SMARTS that do not parse) are returned.":
    remove_fragments, remove_fragments_with_params, remove_fragments_in_place,
    remove_fragments_in_place_with_params, remove_fragments_in_place_mols,
    remove_fragments_in_place_mols_with_params
);

/// The canonical tautomer of a molecule. Mirrors
/// `rdMolStandardize.CanonicalTautomer`.
///
/// `Error::Value` for a tautomer transformation with fewer bonds or charges
/// than its SMARTS has bonds or atoms, which RDKit would read past their
/// ends, and for a molecule whose tautomer scores could overflow RDKit's
/// `int`s (one of millions of atoms or rings).
pub fn canonical_tautomer(mol: &Mol) -> Result<Mol> {
    canonical_tautomer_with_params(mol, &CleanupParameters::default())
}

/// Mirrors `rdMolStandardize.CanonicalTautomer(mol, params)`.
pub fn canonical_tautomer_with_params(mol: &Mol, params: &CleanupParameters) -> Result<Mol> {
    copy(ffi::StandardizeOp::CanonicalTautomer, mol, params, false)
}

/// The functions that take a parent's parameters.
macro_rules! parent_functions {
    (
        $op:ident, $python:literal, $doc:literal, $errors:literal,
        $params:ident, $in_place_params:ident:
        $copy:ident, $copy_params:ident, $in_place:ident, $in_place_with:ident,
        $mols:ident, $mols_params:ident
    ) => {
        #[doc = concat!($doc, " Mirrors `rdMolStandardize.", $python, "`.")]
        #[doc = ""]
        #[doc = $errors]
        pub fn $copy(mol: &Mol) -> Result<Mol> {
            $copy_params(mol, &$params::default())
        }

        #[doc = concat!("Mirrors `rdMolStandardize.", $python, "(mol, params, skipStandardize)`.")]
        pub fn $copy_params(mol: &Mol, params: &$params) -> Result<Mol> {
            copy(ffi::StandardizeOp::$op, mol, &params.params, params.skip_standardize)
        }

        #[doc = concat!("As `", stringify!($copy), "`, in place. Mirrors `rdMolStandardize.", $python, "InPlace`.")]
        pub fn $in_place(mol: &mut Mol) -> Result<()> {
            $in_place_with(mol, &$in_place_params::default())
        }

        #[doc = concat!("Mirrors `rdMolStandardize.", $python, "InPlace(mol, params, skipStandardize)`.")]
        pub fn $in_place_with(mol: &mut Mol, params: &$in_place_params) -> Result<()> {
            in_place(ffi::StandardizeOp::$op, mol, &params.params, params.skip_standardize)
        }

        #[doc = concat!("As `", stringify!($in_place), "`, for several molecules. Mirrors `rdMolStandardize.", $python, "InPlace(mols, numThreads)`.")]
        #[doc = ""]
        #[doc = "The molecules are split among `num_threads` threads (0 or less uses all but \
                 that many of the cores), each standardized by one. When some fail, the others \
                 may still have been standardized, and one of the errors is returned."]
        pub fn $mols(mols: &mut [Mol], num_threads: i32) -> Result<()> {
            $mols_params(mols, num_threads, &$in_place_params::default())
        }

        #[doc = concat!("Mirrors `rdMolStandardize.", $python, "InPlace(mols, numThreads, params, skipStandardize)`.")]
        pub fn $mols_params(
            mols: &mut [Mol],
            num_threads: i32,
            params: &$in_place_params,
        ) -> Result<()> {
            in_place_mols(
                ffi::StandardizeOp::$op,
                mols,
                num_threads,
                &params.params,
                params.skip_standardize,
            )
        }
    };
}

parent_functions!(
    TautomerParent, "TautomerParent",
    "The tautomer parent: the standardized canonical tautomer of the molecule.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     `Error::Value` for a molecule of more than `i32::MAX / 255` atoms, whose charge sums could \
     overflow RDKit's `int`s (CONVENTIONS.md U9); for a tautomer transformation with fewer bonds \
     or charges than its SMARTS has bonds or atoms, which RDKit would read past their ends; and \
     for a molecule whose tautomer scores could overflow RDKit's `int`s.",
    TautomerParentParams, TautomerParentInPlaceParams:
    tautomer_parent, tautomer_parent_with_params, tautomer_parent_in_place,
    tautomer_parent_in_place_with_params, tautomer_parent_in_place_mols,
    tautomer_parent_in_place_mols_with_params
);

parent_functions!(
    FragmentParent, "FragmentParent",
    "The fragment parent: the largest fragment after a cleanup.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     unless the cleanup is skipped, `Error::Value` for a molecule of more than `i32::MAX / 255` \
     atoms, whose charge sums could overflow RDKit's `int`s (CONVENTIONS.md U9).",
    FragmentParentParams, FragmentParentInPlaceParams:
    fragment_parent, fragment_parent_with_params, fragment_parent_in_place,
    fragment_parent_in_place_with_params, fragment_parent_in_place_mols,
    fragment_parent_in_place_mols_with_params
);

parent_functions!(
    StereoParent, "StereoParent",
    "The stereo parent: the molecule after a cleanup, without stereochemistry.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     unless the cleanup is skipped, `Error::Value` for a molecule of more than `i32::MAX / 255` \
     atoms, whose charge sums could overflow RDKit's `int`s (CONVENTIONS.md U9).",
    StereoParentParams, StereoParentInPlaceParams:
    stereo_parent, stereo_parent_with_params, stereo_parent_in_place,
    stereo_parent_in_place_with_params, stereo_parent_in_place_mols,
    stereo_parent_in_place_mols_with_params
);

parent_functions!(
    IsotopeParent, "IsotopeParent",
    "The isotope parent: the molecule after a cleanup, without isotopes.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     unless the cleanup is skipped, `Error::Value` for a molecule of more than `i32::MAX / 255` \
     atoms, whose charge sums could overflow RDKit's `int`s (CONVENTIONS.md U9).",
    IsotopeParentParams, IsotopeParentInPlaceParams:
    isotope_parent, isotope_parent_with_params, isotope_parent_in_place,
    isotope_parent_in_place_with_params, isotope_parent_in_place_mols,
    isotope_parent_in_place_mols_with_params
);

parent_functions!(
    ChargeParent, "ChargeParent",
    "The charge parent: the uncharged fragment parent, cleaned up again.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     `Error::Value` for a molecule of more than `i32::MAX / 255` atoms, whose charge sums could \
     overflow RDKit's `int`s (CONVENTIONS.md U9).",
    ChargeParentParams, ChargeParentInPlaceParams:
    charge_parent, charge_parent_with_params, charge_parent_in_place,
    charge_parent_in_place_with_params, charge_parent_in_place_mols,
    charge_parent_in_place_mols_with_params
);

parent_functions!(
    SuperParent, "SuperParent",
    "The super parent: the fragment, charge, isotope, stereo, and tautomer parent of the \
     molecule.",
    "RDKit's errors (a failed sanitization, a file that cannot be read) are returned, and \
     `Error::Value` for a molecule of more than `i32::MAX / 255` atoms, whose charge sums could \
     overflow RDKit's `int`s (CONVENTIONS.md U9); for a tautomer transformation with fewer bonds \
     or charges than its SMARTS has bonds or atoms, which RDKit would read past their ends; and \
     for a molecule whose tautomer scores could overflow RDKit's `int`s.",
    SuperParentParams, SuperParentInPlaceParams:
    super_parent, super_parent_with_params, super_parent_in_place,
    super_parent_in_place_with_params, super_parent_in_place_mols,
    super_parent_in_place_mols_with_params
);

/// Standardizes a SMILES string with `cleanup`'s defaults, returning the
/// canonical SMILES. `Error::Value` for SMILES that do not parse, as in
/// Python, and for more than `i32::MAX / 255` atoms. Mirrors
/// `rdMolStandardize.StandardizeSmiles`.
pub fn standardize_smiles(smiles: &str) -> Result<String> {
    to_string(ffi::standardize_smiles(smiles)?)
}
