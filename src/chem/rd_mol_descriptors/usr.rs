// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's ultrafast shape recognition (USR) descriptors: the
//! moments of the atoms' distances from four reference points.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::descriptors_3d::ffi;
use crate::ffi::geometry::ffi::IndexList;
use crate::geometry::Point3D;

fn to_nested(lists: Vec<ffi::DoubleList>) -> Vec<Vec<f64>> {
    lists.into_iter().map(|list| list.values).collect()
}

/// The 12 USR descriptors of the default conformer. `Error::Value` for a
/// molecule without conformers or with fewer than three atoms. Mirrors
/// `rdMolDescriptors.GetUSR`.
pub fn get_usr(mol: &Mol) -> Result<Vec<f64>> {
    get_usr_with_conf_id(mol, -1)
}

/// Mirrors `rdMolDescriptors.GetUSR(mol, confId)`.
pub fn get_usr_with_conf_id(mol: &Mol, conf_id: i32) -> Result<Vec<f64>> {
    Ok(ffi::get_usr(mol.as_ffi(), conf_id)?)
}

/// Parameters for `get_usrcat_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetUSRCATParams {
    /// The atom classes, each a list of atom numbers that start at 1, as in
    /// Python (its wrapper subtracts 1); `None` for RDKit's four classes
    /// (hydrophobic, aromatic, acceptor, donor). A number outside 1 to the
    /// number of atoms, which RDKit would read coordinates for unchecked, is
    /// `Error::Index` for its atom index, the number − 1 (wrapping for 0).
    /// An empty list is `Error::Value`.
    pub atom_selections: Option<Vec<Vec<u32>>>,
    /// The conformer; −1 is the default conformer.
    pub conf_id: i32,
}

impl Default for GetUSRCATParams {
    fn default() -> GetUSRCATParams {
        GetUSRCATParams {
            atom_selections: None,
            conf_id: -1,
        }
    }
}

/// The 60 USRCAT descriptors: USR for all atoms and for each of RDKit's
/// four atom classes. Mirrors `rdMolDescriptors.GetUSRCAT`.
pub fn get_usrcat(mol: &Mol) -> Result<Vec<f64>> {
    get_usrcat_with_params(mol, &GetUSRCATParams::default())
}

/// With the given classes, 12 values for each and 12 for all atoms. Mirrors
/// `rdMolDescriptors.GetUSRCAT(mol, atomSelections, confId)`.
pub fn get_usrcat_with_params(mol: &Mol, params: &GetUSRCATParams) -> Result<Vec<f64>> {
    let selections: Vec<IndexList> = params
        .atom_selections
        .iter()
        .flatten()
        .map(|indices| IndexList {
            indices: indices.clone(),
        })
        .collect();
    Ok(ffi::get_usrcat(
        mol.as_ffi(),
        params.atom_selections.is_some(),
        &selections,
        params.conf_id,
    )?)
}

/// What `get_usr_distributions_with_outputs` fills in (A9).
/// A Rust addition, for Python's output argument (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct GetUSRDistributionsOutputs {
    /// The reference points: the centroid, the atom nearest it, the atom
    /// farthest from it, and the atom farthest from that. Python's `points`.
    pub points: Vec<Point3D>,
}

/// The distances of the points from the four USR reference points.
/// `Error::Value` for no points. Mirrors
/// `rdMolDescriptors.GetUSRDistributions`.
pub fn get_usr_distributions(coords: &[Point3D]) -> Result<Vec<Vec<f64>>> {
    Ok(get_usr_distributions_with_outputs(coords)?.0)
}

/// With the reference points. Mirrors
/// `rdMolDescriptors.GetUSRDistributions(coords, points)`.
pub fn get_usr_distributions_with_outputs(
    coords: &[Point3D],
) -> Result<(Vec<Vec<f64>>, GetUSRDistributionsOutputs)> {
    let mut points = Vec::new();
    let distances = ffi::get_usr_distributions(coords, &mut points)?;
    Ok((to_nested(distances), GetUSRDistributionsOutputs { points }))
}

/// The distances of the points from each of `points`. `Error::Value` when
/// either is empty. Mirrors `rdMolDescriptors.GetUSRDistributionsFromPoints`.
pub fn get_usr_distributions_from_points(
    coords: &[Point3D],
    points: &[Point3D],
) -> Result<Vec<Vec<f64>>> {
    Ok(to_nested(ffi::get_usr_distributions_from_points(
        coords, points,
    )?))
}

/// The 12 USR descriptors from four distance distributions: each one's
/// mean, standard deviation, and cube root of skewness. `Error::Value` for
/// no distributions or an empty one, and `Error::Runtime` for other than
/// four. Mirrors `rdMolDescriptors.GetUSRFromDistributions`.
pub fn get_usr_from_distributions(distances: &[Vec<f64>]) -> Result<Vec<f64>> {
    let lists: Vec<ffi::DoubleList> = distances
        .iter()
        .map(|values| ffi::DoubleList {
            values: values.clone(),
        })
        .collect();
    Ok(ffi::get_usr_from_distributions(&lists)?)
}

/// The similarity of two USR or USRCAT descriptors, 1 / (1 + mean absolute
/// difference), with each set of 12 values weighted equally.
/// `Error::Value` for descriptors of different lengths. Mirrors
/// `rdMolDescriptors.GetUSRScore`.
pub fn get_usr_score(descriptor1: &[f64], descriptor2: &[f64]) -> Result<f64> {
    get_usr_score_with_weights(descriptor1, descriptor2, &[])
}

/// With a weight for each set of 12 values; `Error::Value` for another
/// number of weights. Empty means 1 for each, as Python's empty list does.
/// Mirrors `rdMolDescriptors.GetUSRScore(descriptor1, descriptor2,
/// weights)`.
pub fn get_usr_score_with_weights(
    descriptor1: &[f64],
    descriptor2: &[f64],
    weights: &[f64],
) -> Result<f64> {
    Ok(ffi::get_usr_score(descriptor1, descriptor2, weights)?)
}
