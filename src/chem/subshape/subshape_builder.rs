// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Subshape.SubshapeBuilder`, a Python module:
//! building a molecule's shape and skeleton.

use crate::chem::Mol;
use crate::chem::rd_mol_chemical_features::MolChemicalFeatureFactory;
use crate::chem::rd_shape_helpers::{EncodeShapeParams, encode_shape_with_params};
use crate::chem::rdchem::ConformerRef;
use crate::chem::subshape::builder_utils::{
    append_skeleton_points, assign_mol_feats_to_points, calculate_directions_at_point,
    cluster_terminal_pts, expand_terminal_pts, find_terminal_pts_from_conformer,
    find_terminal_pts_from_shape,
};
use crate::chem::subshape::subshape_objects::ShapeWithSkeleton;
use crate::geometry::rd_geometry::{UniformGrid3D, UniformGrid3DParams};
use crate::{Error, Result};

/// How [`SubshapeBuilder::combine_subshapes_with_operation`] combines two
/// shapes' grids. Mirrors `SubshapeBuilder.SubshapeCombineOperations`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SubshapeCombineOperations {
    /// Each point's larger value.
    #[default]
    Union,
    /// The sum of the values.
    Sum,
    /// Each point's smaller value.
    Intersect,
}

/// Parameters for [`SubshapeBuilder::generate_subshape_shape_with_params`],
/// with Python's defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GenerateSubshapeShapeParams {
    /// The conformer; −1 for the default.
    pub conf_id: i32,
    /// Also builds the skeleton.
    pub add_skeleton: bool,
    /// As in [`GenerateSubshapeSkeletonParams`].
    pub terminal_pts_only: bool,
    /// As in [`GenerateSubshapeSkeletonParams`].
    pub skel_from_conf: bool,
}

impl Default for GenerateSubshapeShapeParams {
    fn default() -> GenerateSubshapeShapeParams {
        GenerateSubshapeShapeParams {
            conf_id: -1,
            add_skeleton: true,
            terminal_pts_only: false,
            skel_from_conf: true,
        }
    }
}

/// Parameters for [`SubshapeBuilder::generate_subshape_skeleton_with_params`],
/// with Python's defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GenerateSubshapeSkeletonParams {
    /// Keeps only the terminal points.
    pub terminal_pts_only: bool,
    /// Finds the terminal points from the conformer's atoms, if there is
    /// one, rather than from the grid.
    pub skel_from_conf: bool,
}

impl Default for GenerateSubshapeSkeletonParams {
    fn default() -> GenerateSubshapeSkeletonParams {
        GenerateSubshapeSkeletonParams {
            terminal_pts_only: false,
            skel_from_conf: true,
        }
    }
}

/// Builds shapes and skeletons. Python's class attributes are its fields.
/// Mirrors `SubshapeBuilder.SubshapeBuilder`.
#[derive(Debug)]
pub struct SubshapeBuilder {
    /// The grid's size in each dimension. Mirrors `SubshapeBuilder.gridDims`.
    pub grid_dims: (f64, f64, f64),
    /// The grid's spacing. Mirrors `SubshapeBuilder.gridSpacing`.
    pub grid_spacing: f64,
    /// The radius of the windows around points. Mirrors
    /// `SubshapeBuilder.winRad`.
    pub win_rad: f64,
    /// Atoms with fewer neighbors are terminal. Mirrors
    /// `SubshapeBuilder.nbrCount`.
    pub nbr_count: usize,
    /// Terminal points within this times `win_rad` merge. Mirrors
    /// `SubshapeBuilder.terminalPtRadScale`.
    pub terminal_pt_rad_scale: f64,
    /// The window fraction of grid terminal points. Mirrors
    /// `SubshapeBuilder.fraction`.
    pub fraction: f64,
    /// The spacing of skeleton points. Mirrors `SubshapeBuilder.stepSize`.
    pub step_size: f64,
    /// Features to assign to skeleton points; `None` for none. Mirrors
    /// `SubshapeBuilder.featFactory`.
    pub feat_factory: Option<MolChemicalFeatureFactory>,
}

impl Default for SubshapeBuilder {
    fn default() -> SubshapeBuilder {
        SubshapeBuilder {
            grid_dims: (20.0, 15.0, 10.0),
            grid_spacing: 0.5,
            win_rad: 3.0,
            nbr_count: 7,
            terminal_pt_rad_scale: 0.75,
            fraction: 0.25,
            step_size: 1.0,
            feat_factory: None,
        }
    }
}

fn grid_of(shape: &ShapeWithSkeleton) -> Result<&UniformGrid3D> {
    shape.grid.as_ref().ok_or_else(|| Error::Value {
        message: "the shape has no grid".to_owned(),
    })
}

impl SubshapeBuilder {
    fn new_grid(&self, spacing: f64) -> Result<UniformGrid3D> {
        let (x, y, z) = self.grid_dims;
        UniformGrid3D::new_with_params(
            x,
            y,
            z,
            &UniformGrid3DParams {
                spacing,
                ..Default::default()
            },
        )
    }

    /// The shape resampled on a grid of `new_spacing`; a point of the new
    /// grid outside the shape's is `Error::Value`. Mirrors
    /// `SubshapeBuilder.SampleSubshape`.
    pub fn sample_subshape(
        &self,
        subshape1: &ShapeWithSkeleton,
        new_spacing: f64,
    ) -> Result<ShapeWithSkeleton> {
        let ogrid = grid_of(subshape1)?;
        let mut rgrid = self.new_grid(new_spacing)?;
        for idx in 0..rgrid.get_size() {
            let l = rgrid.get_grid_point_loc(idx)?;
            let v = ogrid.get_val_point(&l);
            // A point outside the original grid reads as -1, which Python
            // cannot store (OverflowError).
            let v = u32::try_from(v).map_err(|_| Error::Value {
                message: "can't convert negative value to unsigned int".to_owned(),
            })?;
            rgrid.set_val(idx, v)?;
        }
        Ok(ShapeWithSkeleton {
            grid: Some(rgrid),
            ..Default::default()
        })
    }

    /// The shape of a conformer on a grid of `spacing`, without skeleton.
    pub(crate) fn shape_with_spacing(
        &self,
        cmpd: &Mol,
        conf_id: i32,
        spacing: f64,
    ) -> Result<ShapeWithSkeleton> {
        let mut grid = self.new_grid(spacing)?;
        let params = EncodeShapeParams {
            conf_id,
            ignore_hs: false,
            ..Default::default()
        };
        encode_shape_with_params(cmpd, &mut grid, &params)?;
        Ok(ShapeWithSkeleton {
            grid: Some(grid),
            ..Default::default()
        })
    }

    /// The molecule's shape, encoded on a grid with its hydrogens, and its
    /// skeleton. Mirrors `SubshapeBuilder.GenerateSubshapeShape`.
    pub fn generate_subshape_shape(&self, cmpd: &Mol) -> Result<ShapeWithSkeleton> {
        self.generate_subshape_shape_with_params(cmpd, &Default::default())
    }

    /// Mirrors `SubshapeBuilder.GenerateSubshapeShape(cmpd, confId,
    /// addSkeleton, terminalPtsOnly, skelFromConf)`.
    pub fn generate_subshape_shape_with_params(
        &self,
        cmpd: &Mol,
        params: &GenerateSubshapeShapeParams,
    ) -> Result<ShapeWithSkeleton> {
        let mut shape = self.shape_with_spacing(cmpd, params.conf_id, self.grid_spacing)?;
        if params.add_skeleton {
            let conf = cmpd.get_conformer_with_id(params.conf_id)?;
            let skeleton = GenerateSubshapeSkeletonParams {
                terminal_pts_only: params.terminal_pts_only,
                skel_from_conf: params.skel_from_conf,
            };
            self.generate_subshape_skeleton_with_params(&mut shape, Some(&conf), &skeleton)?;
        }
        Ok(shape)
    }

    /// As [`generate_subshape_shape`](SubshapeBuilder::generate_subshape_shape).
    /// Mirrors `SubshapeBuilder.__call__`.
    pub fn call(&self, cmpd: &Mol) -> Result<ShapeWithSkeleton> {
        self.generate_subshape_shape(cmpd)
    }

    /// Mirrors `SubshapeBuilder.__call__(cmpd, **kwargs)`.
    pub fn call_with_params(
        &self,
        cmpd: &Mol,
        params: &GenerateSubshapeShapeParams,
    ) -> Result<ShapeWithSkeleton> {
        self.generate_subshape_shape_with_params(cmpd, params)
    }

    /// Sets the shape's skeleton: its terminal points (from the conformer
    /// or the grid, clustered and expanded to at least 5), then skeleton
    /// points between them, each with its directions and, with a feature
    /// factory, the conformer's features near it. Fewer than 3 terminal
    /// points are `Error::Value`. Mirrors
    /// `SubshapeBuilder.GenerateSubshapeSkeleton`.
    pub fn generate_subshape_skeleton(
        &self,
        shape: &mut ShapeWithSkeleton,
        conf: Option<&ConformerRef<'_>>,
    ) -> Result<()> {
        self.generate_subshape_skeleton_with_params(shape, conf, &Default::default())
    }

    /// Mirrors `SubshapeBuilder.GenerateSubshapeSkeleton(shape, conf,
    /// terminalPtsOnly, skelFromConf)`.
    pub fn generate_subshape_skeleton_with_params(
        &self,
        shape: &mut ShapeWithSkeleton,
        conf: Option<&ConformerRef<'_>>,
        params: &GenerateSubshapeSkeletonParams,
    ) -> Result<()> {
        let pts = match conf {
            Some(conf) if params.skel_from_conf => {
                find_terminal_pts_from_conformer(conf, self.win_rad, self.nbr_count)?
            }
            _ => find_terminal_pts_from_shape(shape, self.win_rad, self.fraction)?,
        };
        let mut pts = cluster_terminal_pts(&pts, self.win_rad, self.terminal_pt_rad_scale);
        expand_terminal_pts(shape, &mut pts, self.win_rad)?;
        if pts.len() < 3 {
            return Err(Error::Value {
                message: format!("only found {} terminals, need at least 3", pts.len()),
            });
        }
        let grid = grid_of(shape)?;
        if !params.terminal_pts_only {
            pts = append_skeleton_points(grid, &pts, self.win_rad, self.step_size)?;
        }
        for pt in &mut pts {
            calculate_directions_at_point(pt, grid, self.win_rad)?;
        }
        if let (Some(conf), Some(factory)) = (conf, &self.feat_factory) {
            assign_mol_feats_to_points(&mut pts, conf.get_owning_mol(), factory, self.win_rad)?;
        }
        shape.skel_pts = pts;
        Ok(())
    }

    /// A copy of the first shape whose grid is the union of the two.
    /// Mirrors `SubshapeBuilder.CombineSubshapes`.
    pub fn combine_subshapes(
        &self,
        subshape1: &ShapeWithSkeleton,
        subshape2: &ShapeWithSkeleton,
    ) -> Result<ShapeWithSkeleton> {
        self.combine_subshapes_with_operation(
            subshape1,
            subshape2,
            SubshapeCombineOperations::Union,
        )
    }

    /// Combining the grids by `operation`; grids of different parameters
    /// are `Error::Value`. Mirrors `SubshapeBuilder.CombineSubshapes(subshape1,
    /// subshape2, operation)`.
    pub fn combine_subshapes_with_operation(
        &self,
        subshape1: &ShapeWithSkeleton,
        subshape2: &ShapeWithSkeleton,
        operation: SubshapeCombineOperations,
    ) -> Result<ShapeWithSkeleton> {
        let mut cs = subshape1.clone();
        let other = grid_of(subshape2)?;
        let grid = cs.grid.as_mut().ok_or_else(|| Error::Value {
            message: "the shape has no grid".to_owned(),
        })?;
        if !grid.compare_params(other) {
            return Err(Error::Value {
                message: "the grids' parameters differ".to_owned(),
            });
        }
        match operation {
            SubshapeCombineOperations::Union => *grid |= other,
            SubshapeCombineOperations::Sum => *grid += other,
            SubshapeCombineOperations::Intersect => *grid &= other,
        }
        Ok(cs)
    }
}
