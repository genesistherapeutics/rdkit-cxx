// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.TemplateAlign`, a Python module: 2D depictions
//! that place a molecule's atoms where a template's matching atoms are.

use std::collections::BTreeMap;

use crate::chem::Mol;
use crate::chem::rd_depictor::{Compute2DCoordsParams, compute_2d_coords_with_params};
use crate::geometry::rd_geometry::Point2D;
use crate::{Error, Result};

/// Parameters for `align_mol_to_template_2d_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolToTemplate2DParams {
    /// The molecule's atoms that match the template's, in the template's
    /// order; `None` (or empty, as in Python) for the molecule's first
    /// substructure match of the template.
    pub match_: Option<Vec<u32>>,
    /// Removes the molecule's conformers first, rather than adding one.
    pub clear_confs: bool,
    /// The template's conformer; -1 for its first.
    pub template_conf_id: i32,
}

impl Default for AlignMolToTemplate2DParams {
    fn default() -> AlignMolToTemplate2DParams {
        AlignMolToTemplate2DParams {
            match_: None,
            clear_confs: false,
            template_conf_id: -1,
        }
    }
}

/// Gives the molecule 2D coordinates in which the atoms that match the
/// template sit at the template's positions, and returns the new
/// conformer's id. A molecule that does not match is `Error::Value`, as in
/// Python. Mirrors `TemplateAlign.AlignMolToTemplate2D(mol, template)`.
pub fn align_mol_to_template_2d(mol: &mut Mol, template: &Mol) -> Result<u32> {
    align_mol_to_template_2d_with_params(mol, template, &AlignMolToTemplate2DParams::default())
}

/// As [`align_mol_to_template_2d`]. A match longer than the template is
/// `Error::Index`, where Python's raises from RDKit. Mirrors
/// `TemplateAlign.AlignMolToTemplate2D(mol, template, match, clearConfs,
/// templateConfId)`.
pub fn align_mol_to_template_2d_with_params(
    mol: &mut Mol,
    template: &Mol,
    params: &AlignMolToTemplate2DParams,
) -> Result<u32> {
    let match_ = match &params.match_ {
        Some(m) if !m.is_empty() => m.clone(),
        _ => {
            let m = mol.get_substruct_match(template)?;
            if m.is_empty() {
                return Err(Error::Value {
                    message: "no match between mol and template".to_owned(),
                });
            }
            m
        }
    };
    let template_conf = template.get_conformer_with_id(params.template_conf_id)?;
    let mut atom_map = BTreeMap::new();
    for (i, &idx) in match_.iter().enumerate() {
        let p = template_conf.get_atom_position(i as u32)?;
        atom_map.insert(idx, Point2D::new(p.x, p.y));
    }
    compute_2d_coords_with_params(
        mol,
        &Compute2DCoordsParams {
            clear_confs: params.clear_confs,
            coord_map: atom_map,
            ..Default::default()
        },
    )
}
