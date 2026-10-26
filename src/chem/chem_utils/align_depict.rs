// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.ChemUtils.AlignDepict`, a Python script: 2D
//! depictions with a core's atoms where the core's depiction has them. Its
//! command-line entry point is not ported.

use std::collections::BTreeMap;

use crate::chem::Mol;
use crate::chem::rd_depictor::{Compute2DCoordsParams, compute_2d_coords_with_params};
use crate::geometry::rd_geometry::Point2D;
use crate::{Error, Result};

/// Parameters for `align_depict_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default)]
pub struct AlignDepictParams {
    /// A pattern that maps the molecule's atoms to the core's, matching
    /// both; it must have as many atoms as the core. `None` matches the
    /// core itself.
    pub core_pattern: Option<Mol>,
    /// Depicts a molecule that does not match freely instead of failing.
    pub accept_failure: bool,
}

/// Replaces the molecule's conformers with one 2D depiction in which its
/// atoms that match the core sit where the core's depiction has them. A
/// molecule that does not match is `Error::Value`, as in Python. Mirrors
/// `AlignDepict.AlignDepict(mol, core)`.
pub fn align_depict(mol: &mut Mol, core: &Mol) -> Result<()> {
    align_depict_with_params(mol, core, &AlignDepictParams::default())
}

/// As [`align_depict`]. A pattern with another number of atoms than the
/// core, or that the core does not match, is `Error::Value`, as in Python.
/// Mirrors `AlignDepict.AlignDepict(mol, core, corePattern, acceptFailure)`.
pub fn align_depict_with_params(
    mol: &mut Mol,
    core: &Mol,
    params: &AlignDepictParams,
) -> Result<()> {
    let core_match: Vec<u32> = match &params.core_pattern {
        Some(pattern) => {
            if core.get_num_atoms_with_only_explicit(true)?
                != pattern.get_num_atoms_with_only_explicit(true)?
            {
                return Err(Error::Value {
                    message: "When a pattern is provided, it must have the same number of atoms as the core"
                        .to_owned(),
                });
            }
            let m = core.get_substruct_match(pattern)?;
            if m.is_empty() {
                return Err(Error::Value {
                    message: "Core does not map to itself".to_owned(),
                });
            }
            m
        }
        None => (0..core.get_num_atoms_with_only_explicit(true)?).collect(),
    };
    let match_ = mol.get_substruct_match(params.core_pattern.as_ref().unwrap_or(core))?;
    let mut coord_map = BTreeMap::new();
    if match_.is_empty() {
        if !params.accept_failure {
            return Err(Error::Value {
                message: "Substructure match with core not found.".to_owned(),
            });
        }
    } else {
        let conf = core.get_conformer()?;
        for (i, &idx) in match_.iter().enumerate() {
            let core_idx = *core_match.get(i).ok_or(Error::index(i, core_match.len()))?;
            let pt3 = conf.get_atom_position(core_idx)?;
            coord_map.insert(idx, Point2D::new(pt3.x, pt3.y));
        }
    }
    compute_2d_coords_with_params(
        mol,
        &Compute2DCoordsParams {
            clear_confs: true,
            coord_map,
            canon_orient: false,
            ..Default::default()
        },
    )?;
    Ok(())
}
