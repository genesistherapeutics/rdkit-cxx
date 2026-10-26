// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.ChemUtils.SDFToCSV`, a Python script: SD files as
//! CSV, a SMILES column and one column per property. Its command-line entry
//! point is not ported.

use std::io::Write;

use crate::chem::{
    Mol, SDMolSupplier, SmilesWriteParams, mol_from_smiles, mol_to_smiles_with_params,
};
use crate::{Error, Result};

/// Parameters for `convert_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ConvertParams {
    /// A property to write first, as the key; `None` for none.
    pub key_col: Option<String>,
    /// The number of records after which to stop; -1 for all.
    pub stop_after: i64,
    /// Writes stereochemistry in the SMILES.
    pub include_chirality: bool,
    /// A property whose SMILES to write, canonicalized, instead of the
    /// molecule's; empty, or a molecule without it, for the molecule's.
    pub smiles_from: String,
}

impl Default for ConvertParams {
    fn default() -> ConvertParams {
        ConvertParams {
            key_col: None,
            stop_after: -1,
            include_chirality: false,
            smiles_from: String::new(),
        }
    }
}

/// A row as Python's `csv.writer` writes it: fields with a comma, quote,
/// or line break quoted, quotes doubled, and `\r\n` after.
fn write_row(out: &mut impl Write, row: &[String]) -> std::io::Result<()> {
    let fields: Vec<String> = row
        .iter()
        .map(|f| {
            if f.contains([',', '"', '\r', '\n']) || (row.len() == 1 && f.is_empty()) {
                format!("\"{}\"", f.replace('"', "\"\""))
            } else {
                f.clone()
            }
        })
        .collect();
    write!(out, "{}\r\n", fields.join(","))
}

fn os_error(e: std::io::Error) -> Error {
    Error::Os {
        message: e.to_string(),
    }
}

/// Writes the supplier's molecules as CSV: a header of `SMILES` and the
/// first molecule's properties, then a row for each molecule that parses,
/// with an empty field for a property it lacks. A supplier without records
/// is `Error::Index`, where Python's indexing raises, and a write that
/// fails is `Error::Os`. Mirrors `SDFToCSV.Convert(suppl, outFile)`.
pub fn convert(suppl: &mut SDMolSupplier, out_file: &mut impl Write) -> Result<()> {
    convert_with_params(suppl, out_file, &ConvertParams::default())
}

/// As [`convert`]. Mirrors `SDFToCSV.Convert(suppl, outFile, keyCol,
/// stopAfter, includeChirality, smilesFrom)`.
pub fn convert_with_params(
    suppl: &mut SDMolSupplier,
    out_file: &mut impl Write,
    params: &ConvertParams,
) -> Result<()> {
    let mol = suppl.get(0)?;
    let mut prop_names = mol.get_prop_names()?;
    let key_col = params.key_col.as_deref().filter(|k| !k.is_empty());
    if let Some(key) = key_col {
        if let Some(i) = prop_names.iter().position(|p| p == key) {
            prop_names.remove(i);
        }
    }
    let mut out_l = Vec::new();
    if let Some(key) = key_col {
        out_l.push(key.to_owned());
    }
    out_l.push("SMILES".to_owned());
    out_l.extend(prop_names.iter().cloned());
    write_row(out_file, &out_l).map_err(os_error)?;
    let smiles_params = SmilesWriteParams {
        do_isomeric_smiles: params.include_chirality,
        ..Default::default()
    };
    for (i, mol) in suppl.mols().enumerate() {
        let n_done = i as i64 + 1;
        // Python's supplier yields None for a record that does not parse.
        let Ok(mol) = mol else {
            continue;
        };
        let smi = row_smiles(&mol, params, &smiles_params)?;
        let mut out_l = Vec::new();
        if let Some(key) = key_col {
            out_l.push(mol.get_prop(key)?);
        }
        out_l.push(smi);
        for prop in &prop_names {
            out_l.push(if mol.has_prop(prop) {
                mol.get_prop(prop)?
            } else {
                String::new()
            });
        }
        write_row(out_file, &out_l).map_err(os_error)?;
        if n_done == params.stop_after {
            break;
        }
    }
    Ok(())
}

fn row_smiles(
    mol: &Mol,
    params: &ConvertParams,
    smiles_params: &SmilesWriteParams,
) -> Result<String> {
    if params.smiles_from.is_empty() || !mol.has_prop(&params.smiles_from) {
        return mol_to_smiles_with_params(mol, smiles_params);
    }
    let smi = mol.get_prop(&params.smiles_from)?;
    mol_to_smiles_with_params(&mol_from_smiles(&smi)?, smiles_params)
}
