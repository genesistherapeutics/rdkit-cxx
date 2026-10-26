// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdinchi`: InChI strings and InChI keys, through the
//! InChI library RDKit links.
//!
//! `MolBlockToInchi` is not bound: it parses the mol block with the InChI
//! library's own mol-file reader, which crashes on malformed blocks. Parse
//! the block with [`mol_from_mol_block`](crate::chem::mol_from_mol_block)
//! and call [`mol_to_inchi`] instead.

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::inchi::ffi;

/// What `InchiToMol` returns: the molecule, and the InChI library's return
/// code, message, and log. Mirrors the tuple `rdinchi.InchiToMol` returns.
#[derive(Debug)]
pub struct InchiToMolResult {
    /// The molecule; `None` when InChI cannot read the string (Python's
    /// `None`).
    pub mol: Option<Mol>,
    /// InChI's return code: 0 for success, 1 for a warning, 2 or more for
    /// an error, and negative when the string holds no structure.
    pub retcode: i32,
    /// InChI's error or warning message.
    pub message: String,
    /// InChI's log.
    pub log: String,
}

/// What `MolToInchi` returns: the InChI, and the InChI library's return
/// code, message, log, and auxiliary information. Mirrors the tuple
/// `rdinchi.MolToInchi` returns.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MolToInchiResult {
    /// The InChI; empty when InChI failed.
    pub inchi: String,
    /// InChI's return code: 0 for success, 1 for a warning, 2 or more for
    /// an error.
    pub retcode: i32,
    /// InChI's error or warning message.
    pub message: String,
    /// InChI's log.
    pub log: String,
    /// InChI's AuxInfo string; empty with `/AuxNone`.
    pub aux_info: String,
}

/// Parameters for [`inchi_to_mol_with_params`], with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InchiToMolParams {
    /// Sanitizes the molecule.
    pub sanitize: bool,
    /// Removes hydrogens (which sanitizes too) when `sanitize` is set.
    pub remove_hs: bool,
}

impl Default for InchiToMolParams {
    fn default() -> InchiToMolParams {
        InchiToMolParams {
            sanitize: true,
            remove_hs: true,
        }
    }
}

/// The molecule an InChI describes, with InChI's return code, message, and
/// log. A string InChI cannot read gives no molecule (`mol` is `None`), not
/// an error; a molecule that fails sanitization is the sanitization error.
/// Mirrors `rdinchi.InchiToMol`.
///
/// ```no_run
/// use rdkit_cxx::chem::rdinchi::inchi_to_mol;
///
/// let result = inchi_to_mol("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)")?;
/// assert!(result.mol.is_some());
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub fn inchi_to_mol(inchi: &str) -> Result<InchiToMolResult> {
    inchi_to_mol_with_params(inchi, &InchiToMolParams::default())
}

/// As [`inchi_to_mol`], with `sanitize` and `remove_hs`. Mirrors
/// `rdinchi.InchiToMol(inchi, sanitize, removeHs)`.
pub fn inchi_to_mol_with_params(
    inchi: &str,
    params: &InchiToMolParams,
) -> Result<InchiToMolResult> {
    let data = ffi::inchi_to_mol(inchi, params.sanitize, params.remove_hs)?;
    Ok(InchiToMolResult {
        mol: (!data.mol.is_null()).then(|| Mol::from_ffi(data.mol)),
        retcode: data.retcode,
        message: to_string(data.message)?,
        log: to_string(data.log)?,
    })
}

/// The molecule's standard InChI, with InChI's return code, message, log,
/// and auxiliary information. When InChI fails, `inchi` is empty and
/// `retcode` and `message` say why. `Error::Value` for a molecule InChI
/// cannot represent, which RDKit would pass on garbled: more than 32766
/// atoms; more than 32767 stereo centers and stereo bonds; more than 127
/// hydrogens on an atom whose count RDKit passes (any atom but C, N, O, F,
/// Cl, Br, and I, whose hydrogens InChI counts itself, and any radical); an
/// isotope more than 100 from its element's average mass; a coordinate that
/// is not finite or exceeds 1e18 in magnitude; a tetrahedral stereo center
/// with three or four neighbors, counting hydrogens, but fewer than three
/// neighbor atoms, whose other neighbors RDKit leaves unset; or a stereo
/// bond whose stereo atoms are not in the molecule. `Error::Kekulize` (or
/// another sanitization variant) for aromatic systems RDKit cannot kekulize.
/// Mirrors `rdinchi.MolToInchi`.
///
/// ```no_run
/// use rdkit_cxx::chem::{mol_from_smiles, rdinchi::mol_to_inchi};
///
/// let result = mol_to_inchi(&mol_from_smiles("CC(=O)O")?)?;
/// assert_eq!(result.inchi, "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)");
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub fn mol_to_inchi(mol: &Mol) -> Result<MolToInchiResult> {
    mol_to_inchi_with_options(mol, "")
}

/// As [`mol_to_inchi`], with InChI options separated by spaces or tabs, each
/// prefixed with `/` or `-`, in any case. Only the structure-perception
/// options are accepted: `SNon`, `SRel`, `SRac`, `SUCF`, `ChiralFlagON`,
/// `ChiralFlagOFF`, `NEWPSOFF`, `DoNotAddH`, `SUU`, `SLUUD`, `FixedH`,
/// `RecMet`, `KET`, `15T`, `LargeMolecules`, `Polymers`, `Polymers105`,
/// `NPZz`, `SAtZZ`, `AuxNone`, `WarnOnEmptyStructure`, `NoWarnings`, and
/// the timeouts `W<seconds>` (a decimal number) and `WM<milliseconds>` (an
/// integer), at most 10^12 ms. Any other token is `Error::Value`: InChI's
/// other options steer its file input and output, and some corrupt memory
/// through the API RDKit calls (`/OutputSDF /Key`), and InChI reads some
/// timeouts into undefined conversions. Mirrors
/// `rdinchi.MolToInchi(mol, options)`.
///
/// ```no_run
/// use rdkit_cxx::chem::{mol_from_smiles, rdinchi::mol_to_inchi_with_options};
///
/// let mol = mol_from_smiles("OC1=NC=CC=C1")?;
/// let result = mol_to_inchi_with_options(&mol, "/FixedH /SUU")?;
/// assert!(result.inchi.ends_with("/f/h7H"));
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub fn mol_to_inchi_with_options(mol: &Mol, options: &str) -> Result<MolToInchiResult> {
    let data = ffi::mol_to_inchi(mol.as_ffi(), options)?;
    Ok(MolToInchiResult {
        inchi: to_string(data.inchi)?,
        retcode: data.retcode,
        message: to_string(data.message)?,
        log: to_string(data.log)?,
        aux_info: to_string(data.aux_info)?,
    })
}

/// The molecule's standard InChI key; empty when InChI fails. Errors as
/// [`mol_to_inchi`]. Mirrors `rdinchi.MolToInchiKey`.
pub fn mol_to_inchi_key(mol: &Mol) -> Result<String> {
    mol_to_inchi_key_with_options(mol, "")
}

/// As [`mol_to_inchi_key`], with the options
/// [`mol_to_inchi_with_options`] accepts. Mirrors
/// `rdinchi.MolToInchiKey(mol, options)`.
pub fn mol_to_inchi_key_with_options(mol: &Mol, options: &str) -> Result<String> {
    to_string(ffi::mol_to_inchi_key(mol.as_ffi(), options)?)
}

/// The InChI key of an InChI; empty when InChI rejects the string (RDKit
/// logs why). InChI reads the string only up to a NUL, as from Python.
/// Mirrors `rdinchi.InchiToInchiKey`.
pub fn inchi_to_inchi_key(inchi: &str) -> String {
    to_string(expect_ok(ffi::inchi_to_inchi_key(inchi))).expect("an InChI key is ASCII")
}

/// The version of the InChI library, such as `"1.07.3"`. Mirrors
/// `rdinchi.GetInchiVersion`.
pub fn get_inchi_version() -> String {
    to_string(expect_ok(ffi::get_inchi_version())).expect("the InChI version is ASCII")
}
