// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdShapeHelpers`: shapes of conformers encoded on
//! grids, compared, and boxed.
//!
//! ```no_run
//! use rdkit_cxx::chem::rd_dist_geom::embed_molecule;
//! use rdkit_cxx::chem::rd_shape_helpers::shape_tanimoto_dist;
//! use rdkit_cxx::chem::{add_hs, mol_from_smiles};
//!
//! let mut mol = add_hs(&mol_from_smiles("CCO").unwrap()).unwrap();
//! embed_molecule(&mut mol).unwrap();
//! assert_eq!(shape_tanimoto_dist(&mol, &mol).unwrap(), 0.0);
//! ```
//!
//! Transforms are 4 x 4 matrices, as rows. RDKit casts grid coordinates
//! and the extents of the atoms' spheres to `int` without checks; a sphere
//! too large for its grid, or a negative step size, is `Error::Value`, and
//! `encode_shape` skips atoms that RDKit could not place, as it skips any
//! outside the grid. The comparisons orient the shapes along the first
//! conformer's principal axes, so coordinates beyond 1e9 or not finite are
//! `Error::Value` there, as in `rd_mol_transforms`.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rd_mol_transforms::transform_values;
use crate::chem::rdchem::AsConformer;
use crate::data_structs::DiscreteValueType;
use crate::error::expect_ok;
use crate::ffi::shape_helpers::ffi;
use crate::geometry::{Point3D, UniformGrid3D};

/// The 16 values of an optional transform, and whether there is one.
fn optional_transform(trans: &Option<Vec<Vec<f64>>>) -> Result<(bool, Vec<f64>)> {
    match trans {
        Some(trans) => Ok((true, transform_values(trans)?)),
        None => Ok((false, Vec::new())),
    }
}

fn corners(points: Vec<Point3D>) -> (Point3D, Point3D) {
    match points.as_slice() {
        [first, second] => (*first, *second),
        _ => unreachable!("RDKit returns two points"),
    }
}

/// Parameters for `encode_shape_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EncodeShapeParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// A transform applied to the coordinates on the way to the grid; the
    /// molecule is unchanged.
    pub trans: Option<Vec<Vec<f64>>>,
    /// Scales each atom's van der Waals radius, within which grid points
    /// take the largest value.
    pub vdw_scale: f64,
    /// The thickness of each layer of decreasing values outside that radius.
    pub step_size: f64,
    /// The most layers; -1 for as many as the grid's values allow.
    pub max_layers: i32,
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
}

impl Default for EncodeShapeParams {
    fn default() -> EncodeShapeParams {
        EncodeShapeParams {
            conf_id: -1,
            trans: None,
            vdw_scale: 0.8,
            step_size: 0.25,
            max_layers: -1,
            ignore_hs: true,
        }
    }
}

/// Encodes the shape of the molecule's first conformer onto the grid.
/// Mirrors `rdShapeHelpers.EncodeShape`.
pub fn encode_shape(mol: &Mol, grid: &mut UniformGrid3D) -> Result<()> {
    encode_shape_with_params(mol, grid, &EncodeShapeParams::default())
}

/// Mirrors `rdShapeHelpers.EncodeShape(mol, grid, confId, trans, vdwScale,
/// stepSize, maxLayers, ignoreHs)`.
pub fn encode_shape_with_params(
    mol: &Mol,
    grid: &mut UniformGrid3D,
    params: &EncodeShapeParams,
) -> Result<()> {
    let (has_trans, trans) = optional_transform(&params.trans)?;
    Ok(ffi::encode_shape(
        mol.as_ffi(),
        grid.raw_mut(),
        params.conf_id,
        has_trans,
        &trans,
        params.vdw_scale,
        params.step_size,
        params.max_layers,
        params.ignore_hs,
    )?)
}

/// Parameters for the shape comparisons' `_with_params` forms, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeCompareParams {
    /// The first molecule's conformer; -1 for the first.
    pub conf_id1: i32,
    /// The second molecule's conformer; -1 for the first.
    pub conf_id2: i32,
    /// The spacing of the grids that encode the shapes.
    pub grid_spacing: f64,
    /// The bits per grid point.
    pub bits_per_point: DiscreteValueType,
    /// Scales each atom's van der Waals radius.
    pub vdw_scale: f64,
    /// The thickness of each layer outside the scaled radius.
    pub step_size: f64,
    /// The most layers; -1 for as many as the bits per point allow.
    pub max_layers: i32,
    /// Leaves out the hydrogens.
    pub ignore_hs: bool,
}

impl Default for ShapeCompareParams {
    fn default() -> ShapeCompareParams {
        ShapeCompareParams {
            conf_id1: -1,
            conf_id2: -1,
            grid_spacing: 0.5,
            bits_per_point: DiscreteValueType::Twobitvalue,
            vdw_scale: 0.8,
            step_size: 0.25,
            max_layers: -1,
            ignore_hs: true,
        }
    }
}

/// The Tversky index of the molecules' shapes, as they are aligned. Mirrors
/// `rdShapeHelpers.ShapeTverskyIndex`.
pub fn shape_tversky_index(mol1: &Mol, mol2: &Mol, alpha: f64, beta: f64) -> Result<f64> {
    shape_tversky_index_with_params(mol1, mol2, alpha, beta, &ShapeCompareParams::default())
}

/// Mirrors `rdShapeHelpers.ShapeTverskyIndex(mol1, mol2, alpha, beta,
/// confId1, confId2, gridSpacing, bitsPerPoint, vdwScale, stepSize,
/// maxLayers, ignoreHs)`.
pub fn shape_tversky_index_with_params(
    mol1: &Mol,
    mol2: &Mol,
    alpha: f64,
    beta: f64,
    params: &ShapeCompareParams,
) -> Result<f64> {
    compare(mol1, mol2, b't', alpha, beta, params, true)
}

/// The Tanimoto distance between the molecules' shapes, as they are
/// aligned. Mirrors `rdShapeHelpers.ShapeTanimotoDist`.
pub fn shape_tanimoto_dist(mol1: &Mol, mol2: &Mol) -> Result<f64> {
    shape_tanimoto_dist_with_params(mol1, mol2, &ShapeCompareParams::default())
}

/// Mirrors `rdShapeHelpers.ShapeTanimotoDist(mol1, mol2, confId1, confId2,
/// gridSpacing, bitsPerPoint, vdwScale, stepSize, maxLayers, ignoreHs)`.
pub fn shape_tanimoto_dist_with_params(
    mol1: &Mol,
    mol2: &Mol,
    params: &ShapeCompareParams,
) -> Result<f64> {
    compare(mol1, mol2, b'd', 0.0, 0.0, params, true)
}

/// Parameters for `shape_protrude_dist_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeProtrudeDistParams {
    /// The shared parameters. RDKit's protrude distance encodes on grids of
    /// spacing 0.5 whatever `grid_spacing` says, as in Python.
    pub compare: ShapeCompareParams,
    /// Swaps the molecules when the second's shape is the smaller.
    pub allow_reordering: bool,
}

impl Default for ShapeProtrudeDistParams {
    fn default() -> ShapeProtrudeDistParams {
        ShapeProtrudeDistParams {
            compare: ShapeCompareParams::default(),
            allow_reordering: true,
        }
    }
}

/// The fraction of the first molecule's shape that protrudes from the
/// second's, as they are aligned. Mirrors `rdShapeHelpers.ShapeProtrudeDist`.
pub fn shape_protrude_dist(mol1: &Mol, mol2: &Mol) -> Result<f64> {
    shape_protrude_dist_with_params(mol1, mol2, &ShapeProtrudeDistParams::default())
}

/// Mirrors `rdShapeHelpers.ShapeProtrudeDist(mol1, mol2, confId1, confId2,
/// gridSpacing, bitsPerPoint, vdwScale, stepSize, maxLayers, ignoreHs,
/// allowReordering)`.
pub fn shape_protrude_dist_with_params(
    mol1: &Mol,
    mol2: &Mol,
    params: &ShapeProtrudeDistParams,
) -> Result<f64> {
    compare(
        mol1,
        mol2,
        b'p',
        0.0,
        0.0,
        &params.compare,
        params.allow_reordering,
    )
}

fn compare(
    mol1: &Mol,
    mol2: &Mol,
    kind: u8,
    alpha: f64,
    beta: f64,
    params: &ShapeCompareParams,
    allow_reordering: bool,
) -> Result<f64> {
    Ok(ffi::shape_compare(
        mol1.as_ffi(),
        mol2.as_ffi(),
        kind,
        alpha,
        beta,
        params.conf_id1,
        params.conf_id2,
        params.grid_spacing,
        params.bits_per_point.to_ffi(),
        params.vdw_scale,
        params.step_size,
        params.max_layers,
        params.ignore_hs,
        allow_reordering,
    )?)
}

/// Parameters for `compute_conf_box_with_params` and
/// `compute_conf_dims_and_offset_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputeConfBoxParams {
    /// A transform applied to the coordinates first.
    pub trans: Option<Vec<Vec<f64>>>,
    /// Added on every side.
    pub padding: f64,
}

impl Default for ComputeConfBoxParams {
    fn default() -> ComputeConfBoxParams {
        ComputeConfBoxParams {
            trans: None,
            padding: 2.0,
        }
    }
}

fn conf_box(
    conf: &impl AsConformer,
    dims: bool,
    params: &ComputeConfBoxParams,
) -> Result<(Point3D, Point3D)> {
    let (has_trans, trans) = optional_transform(&params.trans)?;
    Ok(corners(ffi::compute_conf_box(
        conf.conformer(),
        dims,
        has_trans,
        &trans,
        params.padding,
    )?))
}

/// The lower and upper corners of the conformer's box. Mirrors
/// `rdShapeHelpers.ComputeConfBox`.
pub fn compute_conf_box(conf: &impl AsConformer) -> Result<(Point3D, Point3D)> {
    compute_conf_box_with_params(conf, &ComputeConfBoxParams::default())
}

/// Mirrors `rdShapeHelpers.ComputeConfBox(conf, trans, padding)`; a
/// transform that is not 4 x 4 is `Error::Value`.
pub fn compute_conf_box_with_params(
    conf: &impl AsConformer,
    params: &ComputeConfBoxParams,
) -> Result<(Point3D, Point3D)> {
    conf_box(conf, false, params)
}

/// The dimensions and the lower corner of the conformer's box. Mirrors
/// `rdShapeHelpers.ComputeConfDimsAndOffset`.
pub fn compute_conf_dims_and_offset(conf: &impl AsConformer) -> Result<(Point3D, Point3D)> {
    compute_conf_dims_and_offset_with_params(conf, &ComputeConfBoxParams::default())
}

/// Mirrors `rdShapeHelpers.ComputeConfDimsAndOffset(conf, trans,
/// padding)`.
pub fn compute_conf_dims_and_offset_with_params(
    conf: &impl AsConformer,
    params: &ComputeConfBoxParams,
) -> Result<(Point3D, Point3D)> {
    conf_box(conf, true, params)
}

/// The box that holds both boxes, each given by its lower and upper
/// corners. Mirrors `rdShapeHelpers.ComputeUnionBox`.
pub fn compute_union_box(box1: (Point3D, Point3D), box2: (Point3D, Point3D)) -> (Point3D, Point3D) {
    corners(expect_ok(ffi::compute_union_box(
        &box1.0, &box1.1, &box2.0, &box2.1,
    )))
}
