// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdDeprotect`: removing protecting groups with
//! reactions.

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::deprotect::ffi;

/// One deprotection: a reaction that removes a protecting group. Mirrors
/// `rdDeprotect.DeprotectData`.
pub struct DeprotectData {
    /// Never null.
    ptr: UniquePtr<ffi::DeprotectData>,
}

impl Clone for DeprotectData {
    fn clone(&self) -> DeprotectData {
        DeprotectData {
            ptr: expect_ok(ffi::deprotect_data_clone(&self.ptr)),
        }
    }
}

impl fmt::Debug for DeprotectData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("DeprotectData")
            .field("deprotection_class", &self.deprotection_class())
            .field("reaction_smarts", &self.reaction_smarts())
            .field("abbreviation", &self.abbreviation())
            .field("full_name", &self.full_name())
            .finish()
    }
}

// Every string is UTF-8: Rust gave it, or RDKit's built-in list did.
fn text(value: std::result::Result<UniquePtr<cxx::CxxString>, cxx::Exception>) -> String {
    to_string(expect_ok(value)).expect("a DeprotectData string is not UTF-8")
}

impl DeprotectData {
    /// A deprotection from reaction SMARTS with one product. Mirrors
    /// `DeprotectData(deprotection_class, reaction_smarts, abbreviation,
    /// full_name)`. `Error::Parse` for reaction SMARTS that do not parse;
    /// a reaction without one product parses but is not valid.
    pub fn new(
        deprotection_class: &str,
        reaction_smarts: &str,
        abbreviation: &str,
        full_name: &str,
    ) -> Result<DeprotectData> {
        Ok(DeprotectData {
            ptr: ffi::deprotect_data_new(
                deprotection_class,
                reaction_smarts,
                abbreviation,
                full_name,
            )?,
        })
    }

    /// The functional group protected. Mirrors
    /// `DeprotectData.deprotection_class`.
    pub fn deprotection_class(&self) -> String {
        text(ffi::deprotect_data_deprotection_class(&self.ptr))
    }

    /// The protecting group's name. Mirrors `DeprotectData.full_name`.
    pub fn full_name(&self) -> String {
        text(ffi::deprotect_data_full_name(&self.ptr))
    }

    /// The protecting group's abbreviation. Mirrors
    /// `DeprotectData.abbreviation`.
    pub fn abbreviation(&self) -> String {
        text(ffi::deprotect_data_abbreviation(&self.ptr))
    }

    /// The deprotection's reaction SMARTS. Mirrors
    /// `DeprotectData.reaction_smarts`.
    pub fn reaction_smarts(&self) -> String {
        text(ffi::deprotect_data_reaction_smarts(&self.ptr))
    }

    /// An example reaction, empty for deprotections made here. Mirrors
    /// `DeprotectData.example`.
    pub fn example(&self) -> String {
        text(ffi::deprotect_data_example(&self.ptr))
    }

    /// Whether the reaction parsed with exactly one product. Mirrors
    /// `DeprotectData.isValid`.
    pub fn is_valid(&self) -> bool {
        expect_ok(ffi::deprotect_data_is_valid(&self.ptr))
    }
}

/// RDKit's built-in deprotections. Mirrors `rdDeprotect.GetDeprotections`.
pub fn get_deprotections() -> Vec<DeprotectData> {
    (0..expect_ok(ffi::num_default_deprotections()))
        .map(|i| DeprotectData {
            ptr: expect_ok(ffi::default_deprotection(i)),
        })
        .collect()
}

fn to_list(deprotections: &[DeprotectData]) -> UniquePtr<ffi::DeprotectList> {
    let mut list = expect_ok(ffi::deprotect_list_new());
    for data in deprotections {
        expect_ok(ffi::deprotect_list_push(list.pin_mut(), &data.ptr));
    }
    list
}

/// The molecule with every protecting group RDKit's built-in deprotections
/// remove, applied until none applies (at most 100 times), and the
/// abbreviations of those used in its `DEPROTECTIONS` property and their
/// count in `DEPROTECTION_COUNT`. Mirrors `rdDeprotect.Deprotect`.
pub fn deprotect(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::deprotect_default(mol.as_ffi())?))
}

/// Mirrors `rdDeprotect.Deprotect(mol, deprotections)`. Invalid
/// deprotections are skipped.
pub fn deprotect_with_deprotections(mol: &Mol, deprotections: &[DeprotectData]) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::deprotect_with(
        mol.as_ffi(),
        &to_list(deprotections),
    )?))
}

/// Deprotects the molecule in place, returning whether it changed. Mirrors
/// `rdDeprotect.DeprotectInPlace`. `Error::Value` for a deprotection
/// whose reaction has more than one reactant.
pub fn deprotect_in_place(mol: &mut Mol) -> Result<bool> {
    Ok(ffi::deprotect_in_place_default(mol.pin_mut())?)
}

/// Mirrors `rdDeprotect.DeprotectInPlace(mol, deprotections)`.
pub fn deprotect_in_place_with_deprotections(
    mol: &mut Mol,
    deprotections: &[DeprotectData],
) -> Result<bool> {
    Ok(ffi::deprotect_in_place_with(
        mol.pin_mut(),
        &to_list(deprotections),
    )?)
}
