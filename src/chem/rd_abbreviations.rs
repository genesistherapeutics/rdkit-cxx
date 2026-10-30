// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdAbbreviations`: abbreviations (and linkers) found in
//! molecules, condensed to labeled dummy atoms or labeled with `SUP`
//! substance groups.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smiles;
//! use rdkit_cxx::chem::rd_abbreviations::{
//!     condense_mol_abbreviations, get_default_abbreviations,
//! };
//!
//! let mol = mol_from_smiles("c1ccccc1C(F)(F)F")?;
//! let abbrevs = get_default_abbreviations();
//! let condensed = condense_mol_abbreviations(&mol, &abbrevs)?;
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::abbreviations::ffi;

/// An abbreviation: its label, the labels a drawing shows for it, and the
/// query molecule that finds it. Mirrors
/// `rdAbbreviations.AbbreviationDefinition`.
///
/// The functions that find abbreviations refuse one without a molecule with
/// `Error::Runtime`, as Python's raise `RuntimeError`, and one whose molecule
/// lacks the `_numDummies` property that the parsers set with `Error::Key`,
/// unless `max_coverage` is zero or less.
#[derive(Clone, Debug)]
pub struct AbbreviationDefinition {
    /// The label. Mirrors `AbbreviationDefinition.label`.
    pub label: String,
    /// The label in a drawing when the bond comes from the right. Mirrors
    /// `AbbreviationDefinition.displayLabel`.
    pub display_label: String,
    /// The label in a drawing when the bond comes from the west. Mirrors
    /// `AbbreviationDefinition.displayLabelW`.
    pub display_label_w: String,
    /// The query molecule, whose first atom is the dummy for the bond to the
    /// rest of the molecule. Mirrors `AbbreviationDefinition.mol`.
    pub mol: Option<Mol>,
}

/// RDKit's default-constructed definition: empty labels and no molecule.
/// Mirrors `AbbreviationDefinition()`.
impl Default for AbbreviationDefinition {
    fn default() -> AbbreviationDefinition {
        AbbreviationDefinition::from_ffi(expect_ok(ffi::abbreviation_definition_default()))
    }
}

// Every string is UTF-8: the labels come from RDKit's ASCII tables, or from a
// &str split at ASCII separators.
fn text(value: UniquePtr<cxx::CxxString>) -> String {
    to_string(value).expect("an abbreviation's label is not UTF-8")
}

impl AbbreviationDefinition {
    fn from_ffi(data: ffi::AbbreviationDefinitionData) -> AbbreviationDefinition {
        AbbreviationDefinition {
            label: text(data.label),
            display_label: text(data.display_label),
            display_label_w: text(data.display_label_w),
            mol: (!data.mol.is_null()).then(|| Mol::from_ffi(data.mol)),
        }
    }
}

fn from_ffi_vec(data: Vec<ffi::AbbreviationDefinitionData>) -> Vec<AbbreviationDefinition> {
    data.into_iter()
        .map(AbbreviationDefinition::from_ffi)
        .collect()
}

/// Calls `f` with a C++ list of the definitions that borrows their molecules
/// (W8).
fn with_list<T>(
    abbrevs: &[AbbreviationDefinition],
    f: impl FnOnce(&ffi::AbbreviationList) -> T,
) -> T {
    let mut list = expect_ok(ffi::abbreviation_list_new());
    for def in abbrevs {
        expect_ok(match &def.mol {
            Some(mol) => ffi::abbreviation_list_push_with_mol(
                list.pin_mut(),
                &def.label,
                &def.display_label,
                &def.display_label_w,
                mol.as_ffi(),
            ),
            None => ffi::abbreviation_list_push(
                list.pin_mut(),
                &def.label,
                &def.display_label,
                &def.display_label_w,
            ),
        });
    }
    // The list drops here, before the borrow of the molecules ends.
    f(&list)
}

/// RDKit's default abbreviations (`CO2Et`, `tBu`, `CF3`, ...), in the order
/// they are tried. Mirrors `rdAbbreviations.GetDefaultAbbreviations`.
pub fn get_default_abbreviations() -> Vec<AbbreviationDefinition> {
    from_ffi_vec(expect_ok(ffi::get_default_abbreviations()))
}

/// RDKit's default linkers (`PEG6`, ..., `Hept`), whose definitions have a
/// dummy atom at each end. Mirrors `rdAbbreviations.GetDefaultLinkers`.
pub fn get_default_linkers() -> Vec<AbbreviationDefinition> {
    from_ffi_vec(expect_ok(ffi::get_default_linkers()))
}

/// Parameters for `parse_abbreviations_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ParseAbbreviationsParams {
    /// Whether dummy atoms after the first are removed, as for linkers.
    pub remove_extra_dummies: bool,
    /// Whether an abbreviation may attach to a dummy atom, as for linkers.
    pub allow_connection_to_dummies: bool,
}

/// Definitions from text, one per line: a label, a SMARTS, and optionally
/// the display labels, separated by spaces or tabs. A SMARTS that does not
/// start with a dummy atom gets one, and one RDKit cannot parse, or with
/// fewer than two atoms, is skipped. So is a line with only a label, which
/// Python's binding reads again as its SMARTS (RDKit reads past the line's
/// fields). Mirrors `rdAbbreviations.ParseAbbreviations`.
pub fn parse_abbreviations(text: &str) -> Result<Vec<AbbreviationDefinition>> {
    parse_abbreviations_with_params(text, &ParseAbbreviationsParams::default())
}

/// Mirrors `rdAbbreviations.ParseAbbreviations(text, removeExtraDummies,
/// allowConnectionToDummies)`; see `ParseAbbreviationsParams`.
pub fn parse_abbreviations_with_params(
    text: &str,
    params: &ParseAbbreviationsParams,
) -> Result<Vec<AbbreviationDefinition>> {
    Ok(from_ffi_vec(ffi::parse_abbreviations(
        text,
        params.remove_extra_dummies,
        params.allow_connection_to_dummies,
    )?))
}

/// Linker definitions from text, as `parse_abbreviations` with both of
/// `ParseAbbreviationsParams`' options set. Mirrors
/// `rdAbbreviations.ParseLinkers`.
pub fn parse_linkers(text: &str) -> Result<Vec<AbbreviationDefinition>> {
    parse_abbreviations_with_params(
        text,
        &ParseAbbreviationsParams {
            remove_extra_dummies: true,
            allow_connection_to_dummies: true,
        },
    )
}

/// Parameters for `condense_mol_abbreviations_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CondenseMolAbbreviationsParams {
    /// Abbreviations that would cover this fraction of the molecule's atoms
    /// (dummies aside) or more are skipped; zero or less for no limit.
    pub max_coverage: f64,
    /// Whether the result's rings are found again (symmetrized SSSR); the
    /// result is not sanitized either way.
    pub sanitize: bool,
}

impl Default for CondenseMolAbbreviationsParams {
    fn default() -> CondenseMolAbbreviationsParams {
        CondenseMolAbbreviationsParams {
            max_coverage: 0.4,
            sanitize: true,
        }
    }
}

/// A copy of the molecule with the abbreviations it matches, tried in order,
/// condensed to dummy atoms labeled with them (`atomLabel`,
/// `_displayLabel`, and `_displayLabelW`). The result is not sanitized.
/// Mirrors `rdAbbreviations.CondenseMolAbbreviations`.
///
/// RDKit marks the atoms and bonds it removes in sets sized for the
/// molecule as it was, but adds bonds as it goes, so abbreviations whose
/// matches it would bond into each other are refused with `Error::Value`, as
/// is what RDKit would misread when it removes the atoms: a bond's
/// `_MolFileBondEndPts`, or substance groups whose `PARENT`s form a cycle.
pub fn condense_mol_abbreviations(mol: &Mol, abbrevs: &[AbbreviationDefinition]) -> Result<Mol> {
    condense_mol_abbreviations_with_params(mol, abbrevs, &CondenseMolAbbreviationsParams::default())
}

/// Mirrors `rdAbbreviations.CondenseMolAbbreviations(mol, abbrevs,
/// maxCoverage, sanitize)`; see `CondenseMolAbbreviationsParams`. With no
/// abbreviations, Python returns its copy without the ring perception
/// `sanitize` runs.
pub fn condense_mol_abbreviations_with_params(
    mol: &Mol,
    abbrevs: &[AbbreviationDefinition],
    params: &CondenseMolAbbreviationsParams,
) -> Result<Mol> {
    let ptr = with_list(abbrevs, |list| {
        ffi::condense_mol_abbreviations(mol.as_ffi(), list, params.max_coverage, params.sanitize)
    })?;
    Ok(Mol::from_ffi(ptr))
}

/// A copy of the molecule with a `SUP` substance group, labeled with its
/// abbreviation (`LABEL`), for each abbreviation it matches, tried in order;
/// abbreviations covering 40% of the atoms or more are skipped. Mirrors
/// `rdAbbreviations.LabelMolAbbreviations`.
pub fn label_mol_abbreviations(mol: &Mol, abbrevs: &[AbbreviationDefinition]) -> Result<Mol> {
    label_mol_abbreviations_with_max_coverage(mol, abbrevs, 0.4)
}

/// As `label_mol_abbreviations`, skipping abbreviations that cover
/// `max_coverage` of the atoms or more (none for zero or less). Mirrors
/// `rdAbbreviations.LabelMolAbbreviations(mol, abbrevs, maxCoverage)`.
pub fn label_mol_abbreviations_with_max_coverage(
    mol: &Mol,
    abbrevs: &[AbbreviationDefinition],
    max_coverage: f64,
) -> Result<Mol> {
    let ptr = with_list(abbrevs, |list| {
        ffi::label_mol_abbreviations(mol.as_ffi(), list, max_coverage)
    })?;
    Ok(Mol::from_ffi(ptr))
}

/// A copy of the molecule with each `SUP` substance group condensed to a
/// dummy atom labeled with the group's `LABEL` (`"abbrev"` without one). The
/// result is not sanitized. Mirrors
/// `rdAbbreviations.CondenseAbbreviationSubstanceGroups`.
///
/// As for `condense_mol_abbreviations`, groups RDKit would bond into each
/// other as it condenses them, groups with atoms out of range, and what RDKit
/// would misread when it removes the atoms are refused with `Error::Value`.
pub fn condense_abbreviation_substance_groups(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::condense_abbreviation_substance_groups(
        mol.as_ffi(),
    )?))
}
