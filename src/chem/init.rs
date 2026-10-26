// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/Chem/__init__.py`, re-exported from
//! `rdkit_cxx::chem` as Python's is from `rdkit.Chem`.

use crate::chem::{
    AssignCIPLabelsParams, AssignStereochemistryParams, Mol, StereoDescriptor, StereoSpecified,
    StereoType, assign_cip_labels_with_params, assign_stereochemistry_with_params,
    find_potential_stereo, get_use_legacy_stereo_perception,
};
use crate::{Error, Result};

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
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)?
            .clear_prop("_ChiralityPossible");
    }
    assign_stereochemistry_with_params(
        mol,
        &AssignStereochemistryParams {
            clean_it: true,
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
