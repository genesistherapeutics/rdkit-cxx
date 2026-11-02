// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdConformerParser`: conformers from molecular
//! dynamics trajectories.

use std::path::Path;

use crate::Result;
use crate::chem::Mol;
use crate::ffi::conformer_parser::ffi;

/// Parameters for `add_conformers_from_amber_trajectory_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AddConformersFromAmberTrajectoryParams {
    /// How many conformers to read; negative for all of them.
    pub num_confs: i32,
    /// Removes the molecule's conformers first.
    pub clear_confs: bool,
}

impl Default for AddConformersFromAmberTrajectoryParams {
    fn default() -> AddConformersFromAmberTrajectoryParams {
        AddConformersFromAmberTrajectoryParams {
            num_confs: -1,
            clear_confs: true,
        }
    }
}

/// Adds a conformer for each frame of an Amber trajectory (a title line,
/// then coordinates), returning their ids. Mirrors
/// `rdConformerParser.AddConformersFromAmberTrajectory`.
///
/// `Error::Os` for a file that cannot be opened. `Error::Value` for a
/// molecule with no atoms (RDKit divides by its atom count), a coordinate
/// count that is not a multiple of three times the atoms, text that is not
/// a number, or more frames asked for than the file holds. On an error the
/// molecule is unchanged; Python clears its conformers before reading.
pub fn add_conformers_from_amber_trajectory(
    mol: &mut Mol,
    traj: impl AsRef<Path>,
) -> Result<Vec<u32>> {
    add_conformers_from_amber_trajectory_with_params(
        mol,
        traj,
        &AddConformersFromAmberTrajectoryParams::default(),
    )
}

/// Mirrors `rdConformerParser.AddConformersFromAmberTrajectory(mol, traj,
/// numConfs, clearConfs)`.
pub fn add_conformers_from_amber_trajectory_with_params(
    mol: &mut Mol,
    traj: impl AsRef<Path>,
    params: &AddConformersFromAmberTrajectoryParams,
) -> Result<Vec<u32>> {
    Ok(ffi::add_conformers_from_amber_trajectory(
        mol.pin_mut(),
        traj.as_ref().as_os_str().as_encoded_bytes(),
        params.num_confs,
        params.clear_confs,
    )?)
}
