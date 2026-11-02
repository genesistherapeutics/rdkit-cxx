// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FragmentMatcher`, a Python module: a matcher for a
//! fragment, given as SMARTS, that a molecule must not also contain any of a
//! list of excluded fragments.

use crate::Result;
use crate::chem::{BondRef, Mol, SubstructMatchParameters, mol_from_smarts};

/// A SMARTS pattern to match, and patterns that veto a match. Mirrors
/// `FragmentMatcher.FragmentMatcher`; `Default` mirrors `FragmentMatcher()`,
/// which matches nothing until `init` gives it a pattern.
#[derive(Clone, Debug, Default)]
pub struct FragmentMatcher {
    on_patt: Option<Mol>,
    off_patts: Vec<Mol>,
}

impl FragmentMatcher {
    /// Adds a pattern that vetoes a match. SMARTS RDKit cannot parse is
    /// `Error::Parse`, and is not added; Python adds `None`, and its later
    /// matches fail. Mirrors `FragmentMatcher.AddExclusion`.
    pub fn add_exclusion(&mut self, sma: &str) -> Result<()> {
        self.off_patts.push(mol_from_smarts(sma)?);
        Ok(())
    }

    /// Sets the pattern to match. SMARTS RDKit cannot parse is
    /// `Error::Parse` and leaves the matcher without a pattern, as Python's
    /// `None` does. Mirrors `FragmentMatcher.Init`.
    pub fn init(&mut self, sma: &str) -> Result<()> {
        self.on_patt = None;
        self.on_patt = Some(mol_from_smarts(sma)?);
        Ok(())
    }

    /// Whether the molecule matches the pattern and none of the exclusions;
    /// false without a pattern. Python returns 1 or 0. Mirrors
    /// `FragmentMatcher.HasMatch`.
    pub fn has_match(&self, mol: &Mol) -> Result<bool> {
        let Some(on_patt) = &self.on_patt else {
            return Ok(false);
        };
        if !mol.has_substruct_match(on_patt)? {
            return Ok(false);
        }
        for patt in &self.off_patts {
            if mol.has_substruct_match(patt)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// The pattern's first match in the molecule, empty for none, whatever
    /// the exclusions; `None` without a pattern. Mirrors
    /// `FragmentMatcher.GetMatch`.
    pub fn get_match(&self, mol: &Mol) -> Result<Option<Vec<u32>>> {
        self.on_patt
            .as_ref()
            .map(|patt| mol.get_substruct_match(patt))
            .transpose()
    }

    /// The pattern's matches in the molecule, each once, whatever the
    /// exclusions; `None` without a pattern. Mirrors
    /// `FragmentMatcher.GetMatches(mol)`.
    pub fn get_matches(&self, mol: &Mol) -> Result<Option<Vec<Vec<u32>>>> {
        self.get_matches_with_uniquify(mol, true)
    }

    /// As [`FragmentMatcher::get_matches`]; `uniquify` false returns every
    /// match, in each atom order. Mirrors
    /// `FragmentMatcher.GetMatches(mol, uniquify)`.
    pub fn get_matches_with_uniquify(
        &self,
        mol: &Mol,
        uniquify: bool,
    ) -> Result<Option<Vec<Vec<u32>>>> {
        let params = SubstructMatchParameters {
            uniquify,
            ..Default::default()
        };
        self.on_patt
            .as_ref()
            .map(|patt| mol.get_substruct_matches_with_params(patt, &params))
            .transpose()
    }

    /// The pattern's bond `idx`; `None` without a pattern, and
    /// `Error::Index` for an index out of range. Mirrors
    /// `FragmentMatcher.GetBond`.
    pub fn get_bond(&self, idx: u32) -> Result<Option<BondRef<'_>>> {
        self.on_patt
            .as_ref()
            .map(|patt| patt.get_bond_with_idx(idx))
            .transpose()
    }
}
