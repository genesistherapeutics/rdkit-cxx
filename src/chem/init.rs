// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/Chem/__init__.py`, re-exported from
//! `rdkit_cxx::chem` as Python's is from `rdkit.Chem`.
//!
//! Python's package also patches faster iterators onto Boost.Python's
//! molecules; Rust's `get_atoms` and `get_bonds` are iterators already.

use std::path::Path;

use crate::chem::{
    AssignCIPLabelsParams, AssignStereochemistryParams, Mol, SDMolSupplier, SDMolSupplierParams,
    SmilesMolSupplier, SmilesMolSupplierParams, SmilesWriteParams, StereoDescriptor,
    StereoSpecified, StereoType, SubstructMatchParameters, TDTMolSupplier, TDTMolSupplierParams,
    assign_cip_labels_with_params, assign_stereochemistry_with_params, find_potential_stereo,
    get_use_legacy_stereo_perception, mol_from_smarts, mol_from_smiles, mol_to_smiles_with_params,
};
use crate::{Error, Result};

/// The matches of a SMARTS pattern in a molecule given as SMILES, as atom
/// indices, each match once. `Error::Parse` for SMILES or SMARTS RDKit
/// cannot parse, where Python fails on `None`. Mirrors
/// `Chem.QuickSmartsMatch(smi, sma)`.
pub fn quick_smarts_match(smi: &str, sma: &str) -> Result<Vec<Vec<u32>>> {
    quick_smarts_match_with_unique(smi, sma, true)
}

/// As [`quick_smarts_match`]; `unique` false returns every match, in each
/// atom order. Python's `display` argument does nothing and is not ported.
/// Mirrors `Chem.QuickSmartsMatch(smi, sma, unique)`.
pub fn quick_smarts_match_with_unique(smi: &str, sma: &str, unique: bool) -> Result<Vec<Vec<u32>>> {
    let m = mol_from_smiles(smi)?;
    let p = mol_from_smarts(sma)?;
    m.get_substruct_matches_with_params(
        &p,
        &SubstructMatchParameters {
            uniquify: unique,
            ..Default::default()
        },
    )
}

/// The canonical SMILES for a SMILES, with stereochemistry. `Error::Parse`
/// for SMILES RDKit cannot parse. Mirrors `Chem.CanonSmiles(smi)`.
pub fn canon_smiles(smi: &str) -> Result<String> {
    canon_smiles_with_use_chiral(smi, true)
}

/// As [`canon_smiles`]; `use_chiral` false leaves out stereochemistry and
/// isotopes. Mirrors `Chem.CanonSmiles(smi, useChiral)`.
pub fn canon_smiles_with_use_chiral(smi: &str, use_chiral: bool) -> Result<String> {
    let m = mol_from_smiles(smi)?;
    mol_to_smiles_with_params(
        &m,
        &SmilesWriteParams {
            do_isomeric_smiles: use_chiral,
            ..Default::default()
        },
    )
}

/// The supplier [`supplier_from_filename`] makes, by the file's extension.
/// Python returns the supplier itself.
#[derive(Debug)]
#[non_exhaustive]
pub enum FileMolSupplier {
    /// For `.sdf` files.
    Sd(SDMolSupplier),
    /// For `.csv` and `.txt` files.
    Smiles(SmilesMolSupplier),
    /// For `.tdt` files.
    Tdt(TDTMolSupplier),
}

impl FileMolSupplier {
    /// The supplier's molecules from its start, each `Error::Parse` (or a
    /// sanitization variant) where Python's supplier yields `None`.
    pub fn mols(&mut self) -> Box<dyn Iterator<Item = Result<Mol>> + '_> {
        match self {
            FileMolSupplier::Sd(suppl) => Box::new(suppl.mols()),
            FileMolSupplier::Smiles(suppl) => Box::new(suppl.mols()),
            FileMolSupplier::Tdt(suppl) => Box::new(suppl.mols()),
        }
    }
}

/// Parameters for `supplier_from_filename_with_params`, with Python's
/// defaults (A2). Python passes its other keyword arguments to the
/// supplier's constructor; here each kind of supplier has its own.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SupplierFromFilenameParams {
    /// The delimiter of a SMILES file; empty for `,` in `.csv` files and a
    /// tab in `.txt` files. It replaces `smiles.delimiter`.
    pub delim: String,
    /// The parameters of an SD supplier.
    pub sd: SDMolSupplierParams,
    /// The parameters of a SMILES supplier but its delimiter.
    pub smiles: SmilesMolSupplierParams,
    /// The parameters of a TDT supplier.
    pub tdt: TDTMolSupplierParams,
}

/// A supplier for a file, picked by its extension, in any case: `.sdf`,
/// `.csv`, `.txt`, or `.tdt`. Another extension is `Error::Value`, as in
/// Python, and a file the supplier cannot open is `Error::Os`. Mirrors
/// `Chem.SupplierFromFilename(fileN)`.
pub fn supplier_from_filename(file_n: impl AsRef<Path>) -> Result<FileMolSupplier> {
    supplier_from_filename_with_params(file_n, &SupplierFromFilenameParams::default())
}

/// As [`supplier_from_filename`]. Python passes its delimiter to
/// `TDTMolSupplier`, which takes none, so it fails for every `.tdt` file;
/// here a TDT supplier ignores `delim`. Mirrors
/// `Chem.SupplierFromFilename(fileN, delim, **kwargs)`.
pub fn supplier_from_filename_with_params(
    file_n: impl AsRef<Path>,
    params: &SupplierFromFilenameParams,
) -> Result<FileMolSupplier> {
    let file_n = file_n.as_ref();
    // Python splits the name at its last dot, so a name without one is its
    // own extension.
    let name = file_n.to_string_lossy();
    let ext = name.rsplit('.').next().unwrap_or_default().to_lowercase();
    let smiles = |default: &str| -> Result<FileMolSupplier> {
        let delimiter = if params.delim.is_empty() {
            default.to_owned()
        } else {
            params.delim.clone()
        };
        let smiles = SmilesMolSupplierParams {
            delimiter,
            ..params.smiles.clone()
        };
        Ok(FileMolSupplier::Smiles(SmilesMolSupplier::new_with_params(
            file_n, &smiles,
        )?))
    };
    match ext.as_str() {
        "sdf" => Ok(FileMolSupplier::Sd(SDMolSupplier::new_with_params(
            file_n, &params.sd,
        )?)),
        "csv" => smiles(","),
        "txt" => smiles("\t"),
        "tdt" => Ok(FileMolSupplier::Tdt(TDTMolSupplier::new_with_params(
            file_n,
            &params.tdt,
        )?)),
        _ => Err(Error::Value {
            message: format!("unrecognized extension: {ext}"),
        }),
    }
}

/// Parameters for `find_mol_chiral_centers_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindMolChiralCentersParams {
    /// Assigns stereochemistry even if it has been already.
    pub force: bool,
    /// Includes unassigned stereocenters, labeled `?`.
    pub include_unassigned: bool,
    /// Labels each center with its CIP code; with the new implementation,
    /// `false` labels specified centers with their `StereoDescriptor`
    /// instead (`Tet_CW`, `Tet_CCW`).
    pub include_cip: bool,
    /// Uses the legacy stereo perception code; `None` for the process-wide
    /// setting (`get_use_legacy_stereo_perception`).
    pub use_legacy_implementation: Option<bool>,
}

impl Default for FindMolChiralCentersParams {
    fn default() -> FindMolChiralCentersParams {
        FindMolChiralCentersParams {
            force: true,
            include_unassigned: false,
            include_cip: true,
            use_legacy_implementation: None,
        }
    }
}

/// The molecule's chiral centers as `(atom index, CIP code)` pairs. Mirrors
/// `Chem.FindMolChiralCenters`.
pub fn find_mol_chiral_centers(mol: &mut Mol) -> Result<Vec<(u32, String)>> {
    find_mol_chiral_centers_with_params(mol, &FindMolChiralCentersParams::default())
}

/// Mirrors `Chem.FindMolChiralCenters(mol, force, includeUnassigned,
/// includeCIP, useLegacyImplementation)`.
///
/// Python switches the process-wide legacy-perception setting for the call,
/// which safe Rust cannot do (CONVENTIONS.md U5). The new implementation
/// does not read it. The legacy one runs `assign_stereochemistry`, which
/// does, so asking for it while the setting is off is `Error::Value`; turn
/// it on with `set_use_legacy_stereo_perception` first.
pub fn find_mol_chiral_centers_with_params(
    mol: &mut Mol,
    params: &FindMolChiralCentersParams,
) -> Result<Vec<(u32, String)>> {
    let use_legacy_setting = get_use_legacy_stereo_perception();
    if params
        .use_legacy_implementation
        .unwrap_or(use_legacy_setting)
    {
        if !use_legacy_setting {
            return Err(Error::Value {
                message: "the legacy implementation needs set_use_legacy_stereo_perception(true)"
                    .to_owned(),
            });
        }
        legacy_chiral_centers(mol, params)
    } else {
        chiral_centers(mol, params)
    }
}

fn legacy_chiral_centers(
    mol: &mut Mol,
    params: &FindMolChiralCentersParams,
) -> Result<Vec<(u32, String)>> {
    assign_stereochemistry_with_params(
        mol,
        &AssignStereochemistryParams {
            clean_it: false,
            force: params.force,
            flag_possible_stereo_centers: params.include_unassigned,
        },
    )?;
    let mut centers = Vec::new();
    for atom in mol.get_atoms() {
        if atom.has_prop("_CIPCode") {
            centers.push((atom.get_idx(), atom.get_prop("_CIPCode")?));
        } else if params.include_unassigned && atom.has_prop("_ChiralityPossible") {
            centers.push((atom.get_idx(), "?".to_owned()));
        }
    }
    Ok(centers)
}

fn chiral_centers(
    mol: &mut Mol,
    params: &FindMolChiralCentersParams,
) -> Result<Vec<(u32, String)>> {
    let itms = find_potential_stereo(mol)?;
    if params.include_cip {
        let centered_on = |type_| {
            itms.iter()
                .filter(|si| si.type_ == type_)
                .map(|si| si.centered_on)
                .collect()
        };
        let label = AssignCIPLabelsParams {
            atoms_to_label: Some(centered_on(StereoType::AtomTetrahedral)),
            bonds_to_label: Some(centered_on(StereoType::BondDouble)),
            ..Default::default()
        };
        assign_cip_labels_with_params(mol, &label)?;
    }
    let mut centers = Vec::new();
    for si in &itms {
        if si.type_ != StereoType::AtomTetrahedral
            || !(params.include_unassigned || si.specified == StereoSpecified::Specified)
        {
            continue;
        }
        let idx = si.centered_on;
        let mut atm = mol.get_atom_with_idx_mut(idx)?;
        let code = if params.include_cip && atm.has_prop("_CIPCode") {
            atm.get_prop("_CIPCode")?
        } else if si.specified != StereoSpecified::Unspecified {
            python_descriptor_name(si.descriptor).to_owned()
        } else {
            atm.set_int_prop("_ChiralityPossible", 1);
            "?".to_owned()
        };
        centers.push((idx, code));
    }
    Ok(centers)
}

/// The descriptor's name in Python's `rdchem.StereoDescriptor`.
fn python_descriptor_name(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::NoValue => "NoValue",
        StereoDescriptor::TetCw => "Tet_CW",
        StereoDescriptor::TetCcw => "Tet_CCW",
        StereoDescriptor::BondCis => "Bond_Cis",
        StereoDescriptor::BondTrans => "Bond_Trans",
        StereoDescriptor::BondAtropCw => "Bond_AtropCW",
        StereoDescriptor::BondAtropCcw => "Bond_AtropCCW",
    }
}
