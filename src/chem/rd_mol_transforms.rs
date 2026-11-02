// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolTransforms`: centroids, principal axes,
//! canonical orientations, and bond lengths, angles, and dihedrals of
//! conformers.
//!
//! The functions take any conformer: `Conformer`, `ConformerRef`, or
//! `ConformerMut`. Those that look at atoms (hydrogens, bonds, rings) need
//! one in a molecule, and fail with `Error::Runtime` for a conformer outside
//! any. Transforms are 4 x 4 matrices, as rows.
//!
//! The principal axes, and the canonical transforms built on them, come from
//! Eigen's eigensolver, which leaves its output uninitialized for input that
//! is not finite, so they take coordinates and weights within 1e9; others are
//! `Error::Value`.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::{AsConformer, AsConformerMut};
use crate::error::Error;
use crate::ffi::mol_transforms::ffi;
use crate::geometry::Point3D;

/// A 4 x 4 transform's 16 values, row-major.
pub(crate) fn transform_values(trans: &[Vec<f64>]) -> Result<Vec<f64>> {
    if trans.len() != 4 || trans.iter().any(|row| row.len() != 4) {
        return Err(Error::Value {
            message: "the transform must be a 4 x 4 matrix".to_owned(),
        });
    }
    Ok(trans.concat())
}

fn rows(values: &[f64], n: usize) -> Vec<Vec<f64>> {
    values.chunks(n).map(<[f64]>::to_vec).collect()
}

/// Parameters for `compute_centroid_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputeCentroidParams {
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
    /// A weight for each atom; `None` weighs them equally.
    pub weights: Option<Vec<f64>>,
}

impl Default for ComputeCentroidParams {
    fn default() -> ComputeCentroidParams {
        ComputeCentroidParams {
            ignore_hs: true,
            weights: None,
        }
    }
}

/// The centroid of the conformer's heavy atoms. Mirrors
/// `rdMolTransforms.ComputeCentroid`.
pub fn compute_centroid(conf: &impl AsConformer) -> Result<Point3D> {
    compute_centroid_with_params(conf, &ComputeCentroidParams::default())
}

/// Mirrors `rdMolTransforms.ComputeCentroid(conf, ignoreHs, weights)`;
/// fewer weights than atoms is `Error::Runtime`.
pub fn compute_centroid_with_params(
    conf: &impl AsConformer,
    params: &ComputeCentroidParams,
) -> Result<Point3D> {
    let weights = params.weights.as_deref().unwrap_or_default();
    Ok(ffi::compute_centroid(
        conf.conformer(),
        params.ignore_hs,
        params.weights.is_some(),
        weights,
    )?)
}

/// Parameters for `compute_canonical_transform_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputeCanonicalTransformParams {
    /// The center of the principal axes; `None` for the centroid.
    pub center: Option<Point3D>,
    /// Normalizes the covariance matrix by the number of atoms.
    pub normalize_covar: bool,
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
}

impl Default for ComputeCanonicalTransformParams {
    fn default() -> ComputeCanonicalTransformParams {
        ComputeCanonicalTransformParams {
            center: None,
            normalize_covar: false,
            ignore_hs: true,
        }
    }
}

/// The transform that aligns the conformer's principal axes with x, y, and
/// z, as rows; the conformer is unchanged. Mirrors
/// `rdMolTransforms.ComputeCanonicalTransform`.
pub fn compute_canonical_transform(conf: &impl AsConformer) -> Result<Vec<Vec<f64>>> {
    compute_canonical_transform_with_params(conf, &ComputeCanonicalTransformParams::default())
}

/// Mirrors `rdMolTransforms.ComputeCanonicalTransform(conf, center,
/// normalizeCovar, ignoreHs)`.
pub fn compute_canonical_transform_with_params(
    conf: &impl AsConformer,
    params: &ComputeCanonicalTransformParams,
) -> Result<Vec<Vec<f64>>> {
    let values = ffi::compute_canonical_transform(
        conf.conformer(),
        params.center.is_some(),
        &params.center.unwrap_or_default(),
        params.normalize_covar,
        params.ignore_hs,
    )?;
    Ok(rows(&values, 4))
}

/// Parameters for `compute_principal_axes_and_moments_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputePrincipalAxesAndMomentsParams {
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
    /// A weight for each atom; `None` weighs them equally.
    pub weights: Option<Vec<f64>>,
}

impl Default for ComputePrincipalAxesAndMomentsParams {
    fn default() -> ComputePrincipalAxesAndMomentsParams {
        ComputePrincipalAxesAndMomentsParams {
            ignore_hs: true,
            weights: None,
        }
    }
}

/// Parameters for
/// `compute_principal_axes_and_moments_from_gyration_matrix_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputePrincipalAxesAndMomentsFromGyrationMatrixParams {
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
    /// A weight for each atom; `None` weighs them equally.
    pub weights: Option<Vec<f64>>,
}

impl Default for ComputePrincipalAxesAndMomentsFromGyrationMatrixParams {
    fn default() -> ComputePrincipalAxesAndMomentsFromGyrationMatrixParams {
        ComputePrincipalAxesAndMomentsFromGyrationMatrixParams {
            ignore_hs: true,
            weights: None,
        }
    }
}

type AxesAndMoments = Option<(Vec<Vec<f64>>, Vec<f64>)>;

fn principal_axes(
    conf: &impl AsConformer,
    from_gyration_matrix: bool,
    ignore_hs: bool,
    weights: Option<&[f64]>,
) -> Result<AxesAndMoments> {
    let mut axes = Vec::new();
    let mut moments = Vec::new();
    let found = ffi::compute_principal_axes_and_moments(
        conf.conformer(),
        from_gyration_matrix,
        ignore_hs,
        weights.is_some(),
        weights.unwrap_or_default(),
        &mut axes,
        &mut moments,
    )?;
    Ok(found.then(|| (rows(&axes, 3), moments)))
}

/// The principal axes, as rows, and moments of the inertia tensor, or
/// `None` when the eigensolver fails (Python's `(None, None)`). Mirrors
/// `rdMolTransforms.ComputePrincipalAxesAndMoments`. Python caches the
/// unweighted result in the molecule's computed properties, where the PMI
/// descriptors read it; this leaves none (U8).
pub fn compute_principal_axes_and_moments(conf: &impl AsConformer) -> Result<AxesAndMoments> {
    compute_principal_axes_and_moments_with_params(
        conf,
        &ComputePrincipalAxesAndMomentsParams::default(),
    )
}

/// Mirrors `rdMolTransforms.ComputePrincipalAxesAndMoments(conf, ignoreHs,
/// weights)`; a weight count other than the number of atoms is
/// `Error::Value`.
pub fn compute_principal_axes_and_moments_with_params(
    conf: &impl AsConformer,
    params: &ComputePrincipalAxesAndMomentsParams,
) -> Result<AxesAndMoments> {
    principal_axes(conf, false, params.ignore_hs, params.weights.as_deref())
}

/// The principal axes and moments of the gyration tensor. Mirrors
/// `rdMolTransforms.ComputePrincipalAxesAndMomentsFromGyrationMatrix`.
pub fn compute_principal_axes_and_moments_from_gyration_matrix(
    conf: &impl AsConformer,
) -> Result<AxesAndMoments> {
    compute_principal_axes_and_moments_from_gyration_matrix_with_params(
        conf,
        &ComputePrincipalAxesAndMomentsFromGyrationMatrixParams::default(),
    )
}

/// Mirrors
/// `rdMolTransforms.ComputePrincipalAxesAndMomentsFromGyrationMatrix(conf,
/// ignoreHs, weights)`.
pub fn compute_principal_axes_and_moments_from_gyration_matrix_with_params(
    conf: &impl AsConformer,
    params: &ComputePrincipalAxesAndMomentsFromGyrationMatrixParams,
) -> Result<AxesAndMoments> {
    principal_axes(conf, true, params.ignore_hs, params.weights.as_deref())
}

/// Applies a 4 x 4 transform, as rows, to every position. Mirrors
/// `rdMolTransforms.TransformConformer`, which copies rows * rows values into
/// its 4 x 4 transform whatever the shape; another shape is `Error::Value`.
pub fn transform_conformer(conf: &mut impl AsConformerMut, trans: &[Vec<f64>]) -> Result<()> {
    let values = transform_values(trans)?;
    Ok(ffi::transform_conformer(conf.conformer_mut(), &values)?)
}

/// Parameters for `canonicalize_conformer_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CanonicalizeConformerParams {
    /// The center of the principal axes; `None` for the centroid.
    pub center: Option<Point3D>,
    /// Normalizes the covariance matrix by the number of atoms.
    pub normalize_covar: bool,
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
}

impl Default for CanonicalizeConformerParams {
    fn default() -> CanonicalizeConformerParams {
        CanonicalizeConformerParams {
            center: None,
            normalize_covar: false,
            ignore_hs: true,
        }
    }
}

/// Moves the conformer so that its principal axes are x, y, and z. Mirrors
/// `rdMolTransforms.CanonicalizeConformer`.
pub fn canonicalize_conformer(conf: &mut impl AsConformerMut) -> Result<()> {
    canonicalize_conformer_with_params(conf, &CanonicalizeConformerParams::default())
}

/// Mirrors `rdMolTransforms.CanonicalizeConformer(conf, center,
/// normalizeCovar, ignoreHs)`.
pub fn canonicalize_conformer_with_params(
    conf: &mut impl AsConformerMut,
    params: &CanonicalizeConformerParams,
) -> Result<()> {
    Ok(ffi::canonicalize_conformer(
        conf.conformer_mut(),
        params.center.is_some(),
        &params.center.unwrap_or_default(),
        params.normalize_covar,
        params.ignore_hs,
    )?)
}

/// Parameters for `canonicalize_mol_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CanonicalizeMolParams {
    /// Normalizes the covariance matrix by the number of atoms.
    pub normalize_covar: bool,
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
}

impl Default for CanonicalizeMolParams {
    fn default() -> CanonicalizeMolParams {
        CanonicalizeMolParams {
            normalize_covar: false,
            ignore_hs: true,
        }
    }
}

/// Canonicalizes each of the molecule's conformers. Mirrors
/// `rdMolTransforms.CanonicalizeMol`.
pub fn canonicalize_mol(mol: &mut Mol) -> Result<()> {
    canonicalize_mol_with_params(mol, &CanonicalizeMolParams::default())
}

/// Mirrors `rdMolTransforms.CanonicalizeMol(mol, normalizeCovar,
/// ignoreHs)`.
pub fn canonicalize_mol_with_params(mol: &mut Mol, params: &CanonicalizeMolParams) -> Result<()> {
    Ok(ffi::canonicalize_mol(
        mol.pin_mut(),
        params.normalize_covar,
        params.ignore_hs,
    )?)
}

/// The distance between two atoms. Mirrors `rdMolTransforms.GetBondLength`.
pub fn get_bond_length(conf: &impl AsConformer, i_atom_id: u32, j_atom_id: u32) -> Result<f64> {
    Ok(ffi::get_bond_length(
        conf.conformer(),
        i_atom_id,
        j_atom_id,
    )?)
}

/// Sets the length of the bond between two atoms, moving the atoms on `j`'s
/// side. Mirrors `rdMolTransforms.SetBondLength`. Atoms that are not bonded,
/// a ring bond, or atoms at the same place are `Error::Value`.
pub fn set_bond_length(
    conf: &mut impl AsConformerMut,
    i_atom_id: u32,
    j_atom_id: u32,
    value: f64,
) -> Result<()> {
    Ok(ffi::set_bond_length(
        conf.conformer_mut(),
        i_atom_id,
        j_atom_id,
        value,
    )?)
}

/// The angle i-j-k, in radians. Mirrors `rdMolTransforms.GetAngleRad`.
pub fn get_angle_rad(
    conf: &impl AsConformer,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
) -> Result<f64> {
    Ok(ffi::get_angle(
        conf.conformer(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        false,
    )?)
}

/// The angle i-j-k, in degrees. Mirrors `rdMolTransforms.GetAngleDeg`.
pub fn get_angle_deg(
    conf: &impl AsConformer,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
) -> Result<f64> {
    Ok(ffi::get_angle(
        conf.conformer(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        true,
    )?)
}

/// Sets the angle i-j-k, in radians, moving the atoms on `k`'s side.
/// Mirrors `rdMolTransforms.SetAngleRad`.
pub fn set_angle_rad(
    conf: &mut impl AsConformerMut,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
    value: f64,
) -> Result<()> {
    Ok(ffi::set_angle(
        conf.conformer_mut(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        value,
        false,
    )?)
}

/// Sets the angle i-j-k, in degrees. Mirrors `rdMolTransforms.SetAngleDeg`.
pub fn set_angle_deg(
    conf: &mut impl AsConformerMut,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
    value: f64,
) -> Result<()> {
    Ok(ffi::set_angle(
        conf.conformer_mut(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        value,
        true,
    )?)
}

/// The dihedral i-j-k-l, in radians. Mirrors
/// `rdMolTransforms.GetDihedralRad`.
pub fn get_dihedral_rad(
    conf: &impl AsConformer,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
    l_atom_id: u32,
) -> Result<f64> {
    Ok(ffi::get_dihedral(
        conf.conformer(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        l_atom_id,
        false,
    )?)
}

/// The dihedral i-j-k-l, in degrees. Mirrors
/// `rdMolTransforms.GetDihedralDeg`.
pub fn get_dihedral_deg(
    conf: &impl AsConformer,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
    l_atom_id: u32,
) -> Result<f64> {
    Ok(ffi::get_dihedral(
        conf.conformer(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        l_atom_id,
        true,
    )?)
}

/// Sets the dihedral i-j-k-l, in radians, moving the atoms on `l`'s side.
/// Mirrors `rdMolTransforms.SetDihedralRad`.
pub fn set_dihedral_rad(
    conf: &mut impl AsConformerMut,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
    l_atom_id: u32,
    value: f64,
) -> Result<()> {
    Ok(ffi::set_dihedral(
        conf.conformer_mut(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        l_atom_id,
        value,
        false,
    )?)
}

/// Sets the dihedral i-j-k-l, in degrees. Mirrors
/// `rdMolTransforms.SetDihedralDeg`.
pub fn set_dihedral_deg(
    conf: &mut impl AsConformerMut,
    i_atom_id: u32,
    j_atom_id: u32,
    k_atom_id: u32,
    l_atom_id: u32,
    value: f64,
) -> Result<()> {
    Ok(ffi::set_dihedral(
        conf.conformer_mut(),
        i_atom_id,
        j_atom_id,
        k_atom_id,
        l_atom_id,
        value,
        true,
    )?)
}
