// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.SaltRemover`, a pure-Python module: removes salt
//! fragments, as definitions of salts match them, from molecules.
//!
//! Python reads its default salts from `Salts.txt` in RDKit's data
//! directory; the crate embeds that file
//! (`wrapper/third_party/rdkit/Data/Salts.txt`), so the default needs no data
//! directory.

use std::path::{Path, PathBuf};

use crate::chem::{
    DeleteSubstructsParams, Mol, SDMolSupplier, SmilesMolSupplier, delete_substructs_with_params,
    get_mol_frags, mol_from_smarts, mol_from_smiles, sanitize_mol,
};
use crate::{Error, Result};

/// The salts RDKit ships, as SMARTS, one per line.
const DEFAULT_SALTS: &str = include_str!("../../wrapper/third_party/rdkit/Data/Salts.txt");

/// The format of salt definitions. Python's `InputFormat` is a class of
/// string constants. Mirrors `SaltRemover.InputFormat`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum InputFormat {
    /// SMARTS, one per line, with `//` comments and anything after the
    /// first space or tab ignored. Mirrors `InputFormat.SMARTS` (`'smarts'`).
    #[default]
    Smarts,
    /// An SD file. Mirrors `InputFormat.MOL` (`'mol'`).
    Mol,
    /// SMILES, one per line. Mirrors `InputFormat.SMILES` (`'smiles'`).
    Smiles,
}

/// A line's SMARTS: its text up to the first space or tab, after removing
/// `//` comments; `None` for a line without one.
fn smarts_from_smarts_line(line: &str) -> Result<Option<Mol>> {
    // Python strips the line, cuts it at "//", and splits at runs of spaces
    // and tabs.
    let line = line.trim().split("//").next().unwrap_or_default();
    if line.is_empty() {
        return Ok(None);
    }
    let smarts = line.split([' ', '\t']).next().unwrap_or_default();
    match mol_from_smarts(smarts) {
        Ok(salt) => Ok(Some(salt)),
        Err(_) => Err(Error::Value {
            message: line.to_owned(),
        }),
    }
}

fn smarts_salts(text: &str) -> Result<Vec<Mol>> {
    let mut salts = Vec::new();
    for line in text.lines() {
        salts.extend(smarts_from_smarts_line(line)?);
    }
    Ok(salts)
}

/// Parameters for `SaltRemover::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SaltRemoverParams {
    /// The file of salt definitions; `None` for RDKit's own salts.
    pub defn_filename: Option<PathBuf>,
    /// Salt definitions as text, one per line, used instead of any file.
    pub defn_data: Option<String>,
    /// The definitions' format.
    pub defn_format: InputFormat,
    /// Matches the salts with chirality.
    pub use_chirality: bool,
}

/// A molecule with its salts removed, and the salts that were. Mirrors the
/// `StrippedMol` named tuple that `SaltRemover.StripMolWithDeleted` returns.
#[derive(Clone, Debug)]
pub struct StrippedMol {
    /// The molecule without its salts.
    pub mol: Mol,
    /// The salt definitions that removed atoms, in their order.
    pub deleted: Vec<Mol>,
}

/// Parameters for `SaltRemover::strip_mol_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct StripMolParams {
    /// Keeps the last fragment when removing salts would remove every one;
    /// when the last fragments are the same salt, it keeps them all.
    pub dont_remove_everything: bool,
    /// Sanitizes the molecule after removing salts from it.
    pub sanitize: bool,
}

impl Default for StripMolParams {
    fn default() -> StripMolParams {
        StripMolParams {
            dont_remove_everything: false,
            sanitize: true,
        }
    }
}

/// Removes fragments that salt definitions match whole. Mirrors
/// `SaltRemover.SaltRemover`.
#[derive(Clone, Debug)]
pub struct SaltRemover {
    defn_filename: Option<PathBuf>,
    salts: Vec<Mol>,
    use_chirality: bool,
}

impl SaltRemover {
    /// A remover for RDKit's own salts. Mirrors `SaltRemover()`.
    pub fn new() -> SaltRemover {
        SaltRemover::new_with_params(&SaltRemoverParams::default()).expect("RDKit's salts parse")
    }

    /// A remover for the salts `params` defines. A definition RDKit cannot
    /// parse is `Error::Value` with its line, as in Python, and so is a
    /// line of `defn_data` without one (a blank line or a comment, which
    /// Python rejects in data but not in files), and `defn_data` in the
    /// `Mol` format. A file that cannot be read is `Error::Os`, and a record
    /// of an SD or SMILES file that does not parse is its error, where
    /// Python keeps `None` and fails when it strips. RDKit's own salts are
    /// SMARTS: in another format, without `defn_filename` or `defn_data`,
    /// they are `Error::Value`, where Python misreads the file. Mirrors
    /// `SaltRemover(defnFilename, defnData, defnFormat, useChirality)`.
    pub fn new_with_params(params: &SaltRemoverParams) -> Result<SaltRemover> {
        let salts = match (&params.defn_data, &params.defn_filename) {
            (Some(data), _) if !data.is_empty() => {
                // Python reads the data line by line, each with its newline.
                let mut salts = Vec::new();
                for line in data.split_inclusive('\n') {
                    let salt = match params.defn_format {
                        InputFormat::Smarts => smarts_from_smarts_line(line)?,
                        InputFormat::Smiles => mol_from_smiles(line).ok(),
                        InputFormat::Mol => {
                            return Err(Error::Value {
                                message: "Unsupported format for supplier.".to_owned(),
                            });
                        }
                    };
                    salts.push(salt.ok_or_else(|| Error::Value {
                        message: line.to_owned(),
                    })?);
                }
                salts
            }
            (_, Some(file)) => match params.defn_format {
                InputFormat::Smarts => {
                    smarts_salts(&std::fs::read_to_string(file).map_err(|e| Error::Os {
                        message: format!("{}: {e}", file.display()),
                    })?)?
                }
                InputFormat::Mol => SDMolSupplier::new(file)?.mols().collect::<Result<_>>()?,
                InputFormat::Smiles => SmilesMolSupplier::new(file)?
                    .mols()
                    .collect::<Result<_>>()?,
            },
            (_, None) => match params.defn_format {
                InputFormat::Smarts => smarts_salts(DEFAULT_SALTS)?,
                _ => {
                    return Err(Error::Value {
                        message: "RDKit's own salts are SMARTS".to_owned(),
                    });
                }
            },
        };
        Ok(SaltRemover {
            defn_filename: params.defn_filename.clone(),
            salts,
            use_chirality: params.use_chirality,
        })
    }

    /// The file the definitions were read from; `None` for RDKit's own
    /// salts, which Python reads from `Salts.txt` in its data directory and
    /// the crate embeds. Mirrors `SaltRemover.defnFilename`.
    pub fn defn_filename(&self) -> Option<&Path> {
        self.defn_filename.as_deref()
    }

    /// The salt definitions, in the order they are tried. Mirrors the
    /// `salts` attribute.
    pub fn salts(&self) -> &[Mol] {
        &self.salts
    }

    /// The molecule without the fragments the salts match whole. Mirrors
    /// `SaltRemover.StripMol(mol)` and calling the remover.
    pub fn strip_mol(&self, mol: &Mol) -> Result<Mol> {
        self.strip_mol_with_params(mol, &StripMolParams::default())
    }

    /// As [`SaltRemover::strip_mol`]. Mirrors
    /// `SaltRemover.StripMol(mol, dontRemoveEverything, sanitize)`.
    pub fn strip_mol_with_params(&self, mol: &Mol, params: &StripMolParams) -> Result<Mol> {
        Ok(self.strip(mol, params)?.mol)
    }

    /// The molecule without its salts, and the salts removed. Mirrors
    /// `SaltRemover.StripMolWithDeleted(mol)`.
    pub fn strip_mol_with_deleted(&self, mol: &Mol) -> Result<StrippedMol> {
        self.strip_mol_with_deleted_with_dont_remove_everything(mol, false)
    }

    /// As [`SaltRemover::strip_mol_with_deleted`], keeping the last fragment
    /// when `dont_remove_everything` is set. Mirrors
    /// `SaltRemover.StripMolWithDeleted(mol, dontRemoveEverything)`.
    pub fn strip_mol_with_deleted_with_dont_remove_everything(
        &self,
        mol: &Mol,
        dont_remove_everything: bool,
    ) -> Result<StrippedMol> {
        self.strip(
            mol,
            &StripMolParams {
                dont_remove_everything,
                sanitize: true,
            },
        )
    }

    /// Removes one salt's fragments, again while that removes atoms.
    fn apply_pattern(&self, m: Mol, salt: &Mol, not_everything: bool) -> Result<Mol> {
        let mut n_ats = m.get_num_atoms();
        if n_ats == 0 {
            return Ok(m);
        }
        let params = DeleteSubstructsParams {
            only_frags: true,
            use_chirality: self.use_chirality,
        };
        let t = delete_substructs_with_params(&m, salt, &params)?;
        if not_everything && t.get_num_atoms() == 0 {
            return Ok(m);
        }
        let mut res = t;
        while res.get_num_atoms() > 0 && n_ats > res.get_num_atoms() {
            n_ats = res.get_num_atoms();
            let t = delete_substructs_with_params(&res, salt, &params)?;
            if not_everything && t.get_num_atoms() == 0 {
                break;
            }
            res = t;
        }
        Ok(res)
    }

    fn strip(&self, mol: &Mol, params: &StripMolParams) -> Result<StrippedMol> {
        let mut deleted = Vec::new();
        if params.dont_remove_everything && get_mol_frags(mol).len() <= 1 {
            return Ok(StrippedMol {
                mol: mol.clone(),
                deleted,
            });
        }
        let mut modified = false;
        let mut mol = mol.clone();
        let mut natoms = mol.get_num_atoms();
        for salt in &self.salts {
            mol = self.apply_pattern(mol, salt, params.dont_remove_everything)?;
            if natoms != mol.get_num_atoms() {
                natoms = mol.get_num_atoms();
                modified = true;
                deleted.push(salt.clone());
                if params.dont_remove_everything && get_mol_frags(&mol).len() <= 1 {
                    break;
                }
            }
        }
        if params.sanitize && modified && mol.get_num_atoms() > 0 {
            sanitize_mol(&mut mol)?;
        }
        Ok(StrippedMol { mol, deleted })
    }
}

impl Default for SaltRemover {
    /// Mirrors `SaltRemover()`.
    fn default() -> SaltRemover {
        SaltRemover::new()
    }
}
