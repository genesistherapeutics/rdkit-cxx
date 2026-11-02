// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Descriptors3D`, a Python module: descriptors of a
//! conformer's shape, from its principal moments of inertia and plane of
//! best fit.
//!
//! The functions call `rdkit_cxx::chem::rd_mol_descriptors`, as Python's
//! do, and take its parameter structs. Python's `force` arguments only
//! steer RDKit's descriptor caches, which the crate never keeps
//! (CONVENTIONS.md U8). Each fails with `Error::Value` when the molecule
//! has no such conformer.

use crate::chem::Mol;
use crate::chem::rd_mol_descriptors as descriptors;
use crate::{Error, Result};

/// The smallest principal moment of inertia, mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.PMI1`, which calls `rdMolDescriptors.CalcPMI1`.
pub fn pmi1(mol: &Mol) -> Result<f64> {
    descriptors::calc_pmi1(mol)
}

/// Mirrors `Descriptors3D.PMI1(mol, confId, useAtomicMasses)`.
pub fn pmi1_with_params(mol: &Mol, params: &descriptors::CalcPMI1Params) -> Result<f64> {
    descriptors::calc_pmi1_with_params(mol, params)
}

/// The middle principal moment of inertia, mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.PMI2`, which calls `rdMolDescriptors.CalcPMI2`.
pub fn pmi2(mol: &Mol) -> Result<f64> {
    descriptors::calc_pmi2(mol)
}

/// Mirrors `Descriptors3D.PMI2(mol, confId, useAtomicMasses)`.
pub fn pmi2_with_params(mol: &Mol, params: &descriptors::CalcPMI2Params) -> Result<f64> {
    descriptors::calc_pmi2_with_params(mol, params)
}

/// The largest principal moment of inertia, mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.PMI3`, which calls `rdMolDescriptors.CalcPMI3`.
pub fn pmi3(mol: &Mol) -> Result<f64> {
    descriptors::calc_pmi3(mol)
}

/// Mirrors `Descriptors3D.PMI3(mol, confId, useAtomicMasses)`.
pub fn pmi3_with_params(mol: &Mol, params: &descriptors::CalcPMI3Params) -> Result<f64> {
    descriptors::calc_pmi3_with_params(mol, params)
}

/// The normalized principal moments ratio PMI1 / PMI3 (Sauer and Schwarz, JCIM 43, 987-1003 (2003)), mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.NPR1`, which calls `rdMolDescriptors.CalcNPR1`.
pub fn npr1(mol: &Mol) -> Result<f64> {
    descriptors::calc_npr1(mol)
}

/// Mirrors `Descriptors3D.NPR1(mol, confId, useAtomicMasses)`.
pub fn npr1_with_params(mol: &Mol, params: &descriptors::CalcNPR1Params) -> Result<f64> {
    descriptors::calc_npr1_with_params(mol, params)
}

/// The normalized principal moments ratio PMI2 / PMI3, mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.NPR2`, which calls `rdMolDescriptors.CalcNPR2`.
pub fn npr2(mol: &Mol) -> Result<f64> {
    descriptors::calc_npr2(mol)
}

/// Mirrors `Descriptors3D.NPR2(mol, confId, useAtomicMasses)`.
pub fn npr2_with_params(mol: &Mol, params: &descriptors::CalcNPR2Params) -> Result<f64> {
    descriptors::calc_npr2_with_params(mol, params)
}

/// The radius of gyration (Todeschini and Consonni, _Descriptors from Molecular Geometry_), mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.RadiusOfGyration`, which calls `rdMolDescriptors.CalcRadiusOfGyration`.
pub fn radius_of_gyration(mol: &Mol) -> Result<f64> {
    descriptors::calc_radius_of_gyration(mol)
}

/// Mirrors `Descriptors3D.RadiusOfGyration(mol, confId, useAtomicMasses)`.
pub fn radius_of_gyration_with_params(
    mol: &Mol,
    params: &descriptors::CalcRadiusOfGyrationParams,
) -> Result<f64> {
    descriptors::calc_radius_of_gyration_with_params(mol, params)
}

/// The inertial shape factor, PMI2 / (PMI1 PMI3), mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.InertialShapeFactor`, which calls `rdMolDescriptors.CalcInertialShapeFactor`.
pub fn inertial_shape_factor(mol: &Mol) -> Result<f64> {
    descriptors::calc_inertial_shape_factor(mol)
}

/// Mirrors `Descriptors3D.InertialShapeFactor(mol, confId, useAtomicMasses)`.
pub fn inertial_shape_factor_with_params(
    mol: &Mol,
    params: &descriptors::CalcInertialShapeFactorParams,
) -> Result<f64> {
    descriptors::calc_inertial_shape_factor_with_params(mol, params)
}

/// The eccentricity, sqrt(PMI3² − PMI1²) / PMI3, mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.Eccentricity`, which calls `rdMolDescriptors.CalcEccentricity`.
pub fn eccentricity(mol: &Mol) -> Result<f64> {
    descriptors::calc_eccentricity(mol)
}

/// Mirrors `Descriptors3D.Eccentricity(mol, confId, useAtomicMasses)`.
pub fn eccentricity_with_params(
    mol: &Mol,
    params: &descriptors::CalcEccentricityParams,
) -> Result<f64> {
    descriptors::calc_eccentricity_with_params(mol, params)
}

/// The asphericity, mass-weighted, for the default conformer.
/// Mirrors `Descriptors3D.Asphericity`, which calls `rdMolDescriptors.CalcAsphericity`.
pub fn asphericity(mol: &Mol) -> Result<f64> {
    descriptors::calc_asphericity(mol)
}

/// Mirrors `Descriptors3D.Asphericity(mol, confId, useAtomicMasses)`.
pub fn asphericity_with_params(
    mol: &Mol,
    params: &descriptors::CalcAsphericityParams,
) -> Result<f64> {
    descriptors::calc_asphericity_with_params(mol, params)
}

/// The spherocity index, from the unweighted moments, for the default
/// conformer. Mirrors `Descriptors3D.SpherocityIndex`, which calls
/// `rdMolDescriptors.CalcSpherocityIndex`.
pub fn spherocity_index(mol: &Mol) -> Result<f64> {
    descriptors::calc_spherocity_index(mol)
}

/// For the conformer `conf_id`. Mirrors
/// `Descriptors3D.SpherocityIndex(mol, confId)`.
pub fn spherocity_index_with_conf_id(mol: &Mol, conf_id: i32) -> Result<f64> {
    descriptors::calc_spherocity_index_with_conf_id(mol, conf_id)
}

/// The plane of best fit: the atoms' mean distance from the plane that
/// fits them best (Firth et al., JCIM 52, 2516-2525 (2012)). Mirrors
/// `Descriptors3D.PBF`, which calls `rdMolDescriptors.CalcPBF`.
pub fn pbf(mol: &Mol) -> Result<f64> {
    descriptors::calc_pbf(mol)
}

/// For the conformer `conf_id`. Mirrors `Descriptors3D.PBF(mol, confId)`.
pub fn pbf_with_conf_id(mol: &Mol, conf_id: i32) -> Result<f64> {
    descriptors::calc_pbf_with_conf_id(mol, conf_id)
}

/// A 3D descriptor: a function of a molecule and a conformer id.
/// A Rust addition, naming the functions of `descList`.
pub type Descriptor3DFn = fn(&Mol, i32) -> Result<f64>;

/// Defines the conformer-id forms of the moment descriptors that
/// `desc_list` holds.
macro_rules! with_conf_id {
    ($($name:ident = $with_params:ident, $params:ident;)*) => {
        $(
            pub(super) fn $name(mol: &Mol, conf_id: i32) -> Result<f64> {
                $with_params(
                    mol,
                    &descriptors::$params {
                        conf_id,
                        ..descriptors::$params::default()
                    },
                )
            }
        )*
    };
}

mod conf_id {
    use super::*;

    with_conf_id! {
        pmi1 = pmi1_with_params, CalcPMI1Params;
        pmi2 = pmi2_with_params, CalcPMI2Params;
        pmi3 = pmi3_with_params, CalcPMI3Params;
        npr1 = npr1_with_params, CalcNPR1Params;
        npr2 = npr2_with_params, CalcNPR2Params;
        radius_of_gyration = radius_of_gyration_with_params, CalcRadiusOfGyrationParams;
        inertial_shape_factor = inertial_shape_factor_with_params, CalcInertialShapeFactorParams;
        eccentricity = eccentricity_with_params, CalcEccentricityParams;
        asphericity = asphericity_with_params, CalcAsphericityParams;
    }
}

/// The 3D descriptors' names and functions, in Python's order. Mirrors
/// `Descriptors3D.descList`.
pub fn desc_list() -> Vec<(&'static str, Descriptor3DFn)> {
    vec![
        ("PMI1", conf_id::pmi1),
        ("PMI2", conf_id::pmi2),
        ("PMI3", conf_id::pmi3),
        ("NPR1", conf_id::npr1),
        ("NPR2", conf_id::npr2),
        ("RadiusOfGyration", conf_id::radius_of_gyration),
        ("InertialShapeFactor", conf_id::inertial_shape_factor),
        ("Eccentricity", conf_id::eccentricity),
        ("Asphericity", conf_id::asphericity),
        ("SpherocityIndex", spherocity_index_with_conf_id),
        ("PBF", pbf_with_conf_id),
    ]
}

/// Every 3D descriptor of the default conformer, by name, in
/// [`desc_list`]'s order, as Python's dictionary has them. `Error::Value`
/// for a molecule without conformers, as Python raises. Mirrors
/// `Descriptors3D.CalcMolDescriptors3D`.
pub fn calc_mol_descriptors_3d(mol: &Mol) -> Result<Vec<(&'static str, f64)>> {
    calc_mol_descriptors_3d_with_conf_id(mol, -1)
}

/// Of the conformer `conf_id`. Mirrors
/// `Descriptors3D.CalcMolDescriptors3D(mol, confId)`.
pub fn calc_mol_descriptors_3d_with_conf_id(
    mol: &Mol,
    conf_id: i32,
) -> Result<Vec<(&'static str, f64)>> {
    if mol.get_num_conformers() == 0 {
        return Err(Error::Value {
            message: "Computing 3D Descriptors requires a structure with at least 1 conformer"
                .to_owned(),
        });
    }
    desc_list()
        .into_iter()
        .map(|(name, function)| Ok((name, function(mol, conf_id)?)))
        .collect()
}
