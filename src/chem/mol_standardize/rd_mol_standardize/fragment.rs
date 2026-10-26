// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Removing fragments and choosing the largest.

use std::path::PathBuf;

use cxx::UniquePtr;

use super::CleanupParameters;
use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::standardize_fragment::ffi;
use crate::util::path_bytes;

/// Removes the fragments that match its definitions (counter-ions,
/// solvents). Mirrors `rdMolStandardize.FragmentRemover`.
pub struct FragmentRemover {
    /// Never null.
    ptr: UniquePtr<ffi::FragmentRemover>,
}

impl std::fmt::Debug for FragmentRemover {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FragmentRemover").finish_non_exhaustive()
    }
}

/// RDKit's built-in fragments, leaving the last fragment when all match.
/// Mirrors `FragmentRemover()`.
impl Default for FragmentRemover {
    fn default() -> FragmentRemover {
        FragmentRemover::from_ffi(expect_ok(ffi::fragment_remover_default()))
    }
}

/// Parameters for `FragmentRemover::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentRemoverParams {
    /// The file of fragment definitions (name and SMARTS, tab-separated, one
    /// per line); empty for RDKit's built-in ones.
    pub fragment_filename: PathBuf,
    /// Leaves the last fragments when every fragment left matches.
    pub leave_last: bool,
    /// Returns the molecule unchanged when every fragment matches.
    pub skip_if_all_match: bool,
}

impl Default for FragmentRemoverParams {
    fn default() -> FragmentRemoverParams {
        FragmentRemoverParams {
            fragment_filename: PathBuf::new(),
            leave_last: true,
            skip_if_all_match: false,
        }
    }
}

/// Parameters for `fragment_remover_from_data_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentRemoverFromDataParams {
    /// Leaves the last fragments when every fragment left matches.
    pub leave_last: bool,
    /// Returns the molecule unchanged when every fragment matches.
    pub skip_if_all_match: bool,
}

impl Default for FragmentRemoverFromDataParams {
    fn default() -> FragmentRemoverFromDataParams {
        FragmentRemoverFromDataParams {
            leave_last: true,
            skip_if_all_match: false,
        }
    }
}

impl FragmentRemover {
    fn from_ffi(ptr: UniquePtr<ffi::FragmentRemover>) -> FragmentRemover {
        assert!(!ptr.is_null(), "RDKit returned a null FragmentRemover");
        FragmentRemover { ptr }
    }

    fn raw(&self) -> &ffi::FragmentRemover {
        self.ptr.as_ref().unwrap()
    }

    /// Mirrors `FragmentRemover(fragmentFilename, leave_last,
    /// skip_if_all_match)`. `Error::Os` for a file that cannot be read and
    /// `Error::Value` for SMARTS that do not parse.
    pub fn new_with_params(params: &FragmentRemoverParams) -> Result<FragmentRemover> {
        Ok(FragmentRemover::from_ffi(ffi::fragment_remover_new(
            path_bytes(&params.fragment_filename),
            params.leave_last,
            params.skip_if_all_match,
        )?))
    }

    /// A copy of the molecule without the matching fragments. Mirrors
    /// `FragmentRemover.remove`.
    pub fn remove(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::fragment_remover_remove(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// Removes the matching fragments. Mirrors
    /// `FragmentRemover.removeInPlace`.
    pub fn remove_in_place(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::fragment_remover_remove_in_place(
            self.raw(),
            mol.pin_mut(),
        )?)
    }
}

/// A fragment remover from fragment definitions, one per line: a name and
/// SMARTS, tab-separated (lines starting `//` are comments). `Error::Value`
/// for a line without SMARTS or SMARTS that do not parse. Mirrors
/// `rdMolStandardize.FragmentRemoverFromData`.
pub fn fragment_remover_from_data(fragment_data: &str) -> Result<FragmentRemover> {
    fragment_remover_from_data_with_params(fragment_data, &Default::default())
}

/// Mirrors `rdMolStandardize.FragmentRemoverFromData(fragmentData,
/// leave_last, skip_if_all_match)`.
pub fn fragment_remover_from_data_with_params(
    fragment_data: &str,
    params: &FragmentRemoverFromDataParams,
) -> Result<FragmentRemover> {
    Ok(FragmentRemover::from_ffi(ffi::fragment_remover_from_data(
        fragment_data,
        params.leave_last,
        params.skip_if_all_match,
    )?))
}

/// Keeps a molecule's largest fragment: the one with the most atoms
/// (counting hydrogens, unless only heavy atoms count), then the heaviest,
/// then the first by SMILES. Mirrors `rdMolStandardize.LargestFragmentChooser`,
/// whose options Python does not expose.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LargestFragmentChooser {
    prefer_organic: bool,
    use_atom_count: bool,
    count_heavy_atoms_only: bool,
}

/// Mirrors `LargestFragmentChooser()`.
impl Default for LargestFragmentChooser {
    fn default() -> LargestFragmentChooser {
        LargestFragmentChooser::new_with_prefer_organic(false)
    }
}

impl LargestFragmentChooser {
    /// With `prefer_organic`, a fragment with carbon beats any without.
    /// Mirrors `LargestFragmentChooser(preferOrganic)`.
    pub fn new_with_prefer_organic(prefer_organic: bool) -> LargestFragmentChooser {
        LargestFragmentChooser {
            prefer_organic,
            use_atom_count: true,
            count_heavy_atoms_only: false,
        }
    }

    /// The chooser `prefer_organic` and the `largest_fragment_chooser_`
    /// options describe. Mirrors `LargestFragmentChooser(params)`.
    pub fn from_params(params: &CleanupParameters) -> LargestFragmentChooser {
        LargestFragmentChooser {
            prefer_organic: params.prefer_organic,
            use_atom_count: params.largest_fragment_chooser_use_atom_count,
            count_heavy_atoms_only: params.largest_fragment_chooser_count_heavy_atoms_only,
        }
    }

    /// The largest fragment, sanitized (`Error::MolSanitize` and the like
    /// when it cannot be). Mirrors `LargestFragmentChooser.choose`.
    pub fn choose(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::largest_fragment_chooser_choose(
            mol.as_ffi(),
            self.prefer_organic,
            self.use_atom_count,
            self.count_heavy_atoms_only,
        )?))
    }

    /// Removes every fragment but the largest, without sanitizing. Mirrors
    /// `LargestFragmentChooser.chooseInPlace`.
    pub fn choose_in_place(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::largest_fragment_chooser_choose_in_place(
            mol.pin_mut(),
            self.prefer_organic,
            self.use_atom_count,
            self.count_heavy_atoms_only,
        )?)
    }
}
