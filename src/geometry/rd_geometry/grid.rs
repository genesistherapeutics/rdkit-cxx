// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::path::Path;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::data_structs::{DiscreteValueType, DiscreteValueVect, RealValueVect};
use crate::error::{expect_ok, op};
use crate::ffi::grid::ffi;
use crate::geometry::Point3D;
use crate::util::path_bytes;

/// Defines the methods both grid types share, over their wrappers (W3).
macro_rules! grid_type {
    (
        $(#[$meta:meta])* $name:ident, $ffi:ident, $value:ty, $point_value:ty,
        $from_pkl:ident, $clone:ident, $to_binary:ident, $point_index:ident, $grid_index:ident,
        $grid_indices:ident, $get_val_point:ident, $get_val:ident, $set_val:ident,
        $set_val_point:ident, $point_loc:ident, $size:ident, $num_x:ident, $num_y:ident,
        $num_z:ident, $offset:ident, $spacing:ident, $assign_op:ident
    ) => {
        $(#[$meta])*
        pub struct $name {
            /// Never null.
            ptr: UniquePtr<ffi::$ffi>,
        }

        impl $name {
            pub(crate) fn raw(&self) -> &ffi::$ffi {
                self.ptr.as_ref().unwrap()
            }

            pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::$ffi> {
                self.ptr.pin_mut()
            }

            pub(crate) fn from_ffi(ptr: UniquePtr<ffi::$ffi>) -> $name {
                assert!(!ptr.is_null(), "RDKit returned a null grid");
                $name { ptr }
            }

            /// A grid from a pickle (`to_binary`). Mirrors the class's
            /// `pkl` constructor.
            pub fn from_pkl(pkl: &[u8]) -> Result<$name> {
                Ok($name::from_ffi(ffi::$from_pkl(pkl)?))
            }

            /// RDKit's binary form, which `from_pkl` reads; Python pickles it.
            pub fn to_binary(&self) -> Vec<u8> {
                expect_ok(ffi::$to_binary(self.raw()))
            }

            /// The index of the grid point nearest `point`, or -1 outside the
            /// grid. Mirrors `GetGridPointIndex`.
            pub fn get_grid_point_index(&self, point: &Point3D) -> i32 {
                expect_ok(ffi::$point_index(self.raw(), point))
            }

            /// The index of grid point (`xi`, `yi`, `zi`), or -1 outside the
            /// grid. Mirrors `GetGridIndex`.
            pub fn get_grid_index(&self, xi: u32, yi: u32, zi: u32) -> i32 {
                expect_ok(ffi::$grid_index(self.raw(), xi, yi, zi))
            }

            /// Grid point `idx`'s x, y, and z indices. Mirrors
            /// `GetGridIndices`.
            pub fn get_grid_indices(&self, idx: u32) -> Result<(u32, u32, u32)> {
                let ids = ffi::$grid_indices(self.raw(), idx)?;
                Ok((ids.xi, ids.yi, ids.zi))
            }

            /// The value at the grid point nearest `pt`, or -1 outside the
            /// grid. Mirrors `GetValPoint`.
            pub fn get_val_point(&self, pt: &Point3D) -> $point_value {
                expect_ok(ffi::$get_val_point(self.raw(), pt))
            }

            /// Mirrors `GetVal`.
            pub fn get_val(&self, id: u32) -> Result<$value> {
                Ok(ffi::$get_val(self.raw(), id)?)
            }

            /// Mirrors `SetVal`.
            pub fn set_val(&mut self, id: u32, val: $value) -> Result<()> {
                Ok(ffi::$set_val(self.raw_mut(), id, val)?)
            }

            /// Sets the grid point nearest `pt`; a point outside the grid
            /// changes nothing. Mirrors `SetValPoint`.
            pub fn set_val_point(&mut self, pt: &Point3D, val: $value) -> Result<()> {
                Ok(ffi::$set_val_point(self.raw_mut(), pt, val)?)
            }

            /// Mirrors `GetGridPointLoc`.
            pub fn get_grid_point_loc(&self, point_id: u32) -> Result<Point3D> {
                Ok(ffi::$point_loc(self.raw(), point_id)?)
            }

            /// The number of grid points. Mirrors `GetSize`.
            pub fn get_size(&self) -> u32 {
                expect_ok(ffi::$size(self.raw()))
            }

            /// Mirrors `GetNumX`.
            pub fn get_num_x(&self) -> u32 {
                expect_ok(ffi::$num_x(self.raw()))
            }

            /// Mirrors `GetNumY`.
            pub fn get_num_y(&self) -> u32 {
                expect_ok(ffi::$num_y(self.raw()))
            }

            /// Mirrors `GetNumZ`.
            pub fn get_num_z(&self) -> u32 {
                expect_ok(ffi::$num_z(self.raw()))
            }

            /// The location of the grid's corner. Mirrors `GetOffset`.
            pub fn get_offset(&self) -> Point3D {
                expect_ok(ffi::$offset(self.raw()))
            }

            /// Mirrors `GetSpacing`.
            pub fn get_spacing(&self) -> f64 {
                expect_ok(ffi::$spacing(self.raw()))
            }
        }

        impl Clone for $name {
            fn clone(&self) -> $name {
                $name::from_ffi(expect_ok(ffi::$clone(self.raw())))
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("num_x", &self.get_num_x())
                    .field("num_y", &self.get_num_y())
                    .field("num_z", &self.get_num_z())
                    .field("spacing", &self.get_spacing())
                    .finish()
            }
        }

        /// Adds the values of a grid with the same parameters. Mirrors `+=`.
        impl std::ops::AddAssign<&$name> for $name {
            fn add_assign(&mut self, rhs: &$name) {
                op(ffi::$assign_op(self.raw_mut(), rhs.raw(), b'+'))
            }
        }

        /// Subtracts the values of a grid with the same parameters. Mirrors
        /// `-=`.
        impl std::ops::SubAssign<&$name> for $name {
            fn sub_assign(&mut self, rhs: &$name) {
                op(ffi::$assign_op(self.raw_mut(), rhs.raw(), b'-'))
            }
        }

        /// Keeps the smaller value at each point. Mirrors `&=`.
        impl std::ops::BitAndAssign<&$name> for $name {
            fn bitand_assign(&mut self, rhs: &$name) {
                op(ffi::$assign_op(self.raw_mut(), rhs.raw(), b'&'))
            }
        }

        /// Keeps the larger value at each point. Mirrors `|=`.
        impl std::ops::BitOrAssign<&$name> for $name {
            fn bitor_assign(&mut self, rhs: &$name) {
                op(ffi::$assign_op(self.raw_mut(), rhs.raw(), b'|'))
            }
        }
    };
}

grid_type!(
    /// A grid of small unsigned integers over a box. Mirrors
    /// `rdGeometry.UniformGrid3D`.
    ///
    /// The in-place operators panic for grids with other parameters.
    UniformGrid3D, UniformGrid3D, u32, i32,
    uniform_grid3d_from_pkl, uniform_grid3d_clone, uniform_grid3d_to_binary,
    uniform_grid3d_get_grid_point_index, uniform_grid3d_get_grid_index,
    uniform_grid3d_get_grid_indices, uniform_grid3d_get_val_point,
    uniform_grid3d_get_val, uniform_grid3d_set_val, uniform_grid3d_set_val_point,
    uniform_grid3d_get_grid_point_loc, uniform_grid3d_get_size, uniform_grid3d_get_num_x,
    uniform_grid3d_get_num_y, uniform_grid3d_get_num_z, uniform_grid3d_get_offset,
    uniform_grid3d_get_spacing, uniform_grid3d_assign_op
);

grid_type!(
    /// A grid of real numbers over a box. Mirrors
    /// `rdGeometry.UniformRealValueGrid3D`.
    ///
    /// The in-place operators panic for grids with other parameters.
    UniformRealValueGrid3D, UniformRealValueGrid3D, f64, f64,
    uniform_real_value_grid3d_from_pkl, uniform_real_value_grid3d_clone,
    uniform_real_value_grid3d_to_binary, uniform_real_value_grid3d_get_grid_point_index,
    uniform_real_value_grid3d_get_grid_index, uniform_real_value_grid3d_get_grid_indices,
    uniform_real_value_grid3d_get_val_point, uniform_real_value_grid3d_get_val,
    uniform_real_value_grid3d_set_val, uniform_real_value_grid3d_set_val_point,
    uniform_real_value_grid3d_get_grid_point_loc, uniform_real_value_grid3d_get_size,
    uniform_real_value_grid3d_get_num_x, uniform_real_value_grid3d_get_num_y,
    uniform_real_value_grid3d_get_num_z, uniform_real_value_grid3d_get_offset,
    uniform_real_value_grid3d_get_spacing, uniform_real_value_grid3d_assign_op
);

/// Parameters for `UniformGrid3D::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UniformGrid3DParams {
    /// The distance between grid points.
    pub spacing: f64,
    /// How many bits each grid point's value takes.
    pub val_type: DiscreteValueType,
    /// The grid's corner; `None` centers the grid on the origin.
    pub off_set: Option<Point3D>,
}

impl Default for UniformGrid3DParams {
    fn default() -> UniformGrid3DParams {
        UniformGrid3DParams {
            spacing: 0.5,
            val_type: DiscreteValueType::Twobitvalue,
            off_set: None,
        }
    }
}

/// Parameters for `UniformGrid3D::set_sphere_occupancy_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SetSphereOccupancyParams {
    /// The most layers of decreasing values around the sphere; -1 for as many
    /// as the values allow.
    pub max_layers: i32,
    /// Skips spheres centered outside the grid rather than failing.
    pub ignore_out_of_bound: bool,
}

impl Default for SetSphereOccupancyParams {
    fn default() -> SetSphereOccupancyParams {
        SetSphereOccupancyParams {
            max_layers: -1,
            ignore_out_of_bound: true,
        }
    }
}

impl UniformGrid3D {
    /// A grid of `dim_x` by `dim_y` by `dim_z`, which must be positive; a
    /// grid of 2^31 points or more is `Error::Value`. Mirrors
    /// `UniformGrid3D(dimX, dimY, dimZ)`.
    pub fn new(dim_x: f64, dim_y: f64, dim_z: f64) -> Result<UniformGrid3D> {
        UniformGrid3D::new_with_params(dim_x, dim_y, dim_z, &Default::default())
    }

    /// Mirrors `UniformGrid3D(dimX, dimY, dimZ, spacing, valType, offSet)`.
    pub fn new_with_params(
        dim_x: f64,
        dim_y: f64,
        dim_z: f64,
        params: &UniformGrid3DParams,
    ) -> Result<UniformGrid3D> {
        Ok(UniformGrid3D::from_ffi(ffi::uniform_grid3d_new(
            dim_x,
            dim_y,
            dim_z,
            params.spacing,
            params.val_type.to_ffi(),
            params.off_set.is_some(),
            &params.off_set.unwrap_or_default(),
        )?))
    }

    /// A copy of the values. Mirrors `GetOccupancyVect`, which returns the
    /// grid's own.
    pub fn get_occupancy_vect(&self) -> DiscreteValueVect {
        DiscreteValueVect::from_ffi(expect_ok(ffi::uniform_grid3d_get_occupancy_vect(
            self.raw(),
        )))
    }

    /// Whether `other` has the same size, spacing, and offset. Mirrors
    /// `UniformGrid3D.CompareParams`.
    pub fn compare_params(&self, other: &UniformGrid3D) -> bool {
        expect_ok(ffi::uniform_grid3d_compare_params(self.raw(), other.raw()))
    }

    /// Sets the points within `radius` of `center`, with layers of decreasing
    /// values `step_size` apart around it. A center outside the grid is
    /// skipped, or `Error::Runtime` without `ignore_out_of_bound`. A negative
    /// step size, which RDKit doesn't check, or a sphere too large to index
    /// is `Error::Value`. Mirrors `UniformGrid3D.SetSphereOccupancy`.
    pub fn set_sphere_occupancy(
        &mut self,
        center: &Point3D,
        radius: f64,
        step_size: f64,
    ) -> Result<()> {
        self.set_sphere_occupancy_with_params(center, radius, step_size, &Default::default())
    }

    /// Mirrors `UniformGrid3D.SetSphereOccupancy(center, radius, stepSize,
    /// maxLayers, ignoreOutOfBound)`.
    pub fn set_sphere_occupancy_with_params(
        &mut self,
        center: &Point3D,
        radius: f64,
        step_size: f64,
        params: &SetSphereOccupancyParams,
    ) -> Result<()> {
        Ok(ffi::uniform_grid3d_set_sphere_occupancy(
            self.raw_mut(),
            center,
            radius,
            step_size,
            params.max_layers,
            params.ignore_out_of_bound,
        )?)
    }
}

/// Parameters for `UniformRealValueGrid3D::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UniformRealValueGrid3DParams {
    /// The distance between grid points.
    pub spacing: f64,
    /// The grid's corner; `None` centers the grid on the origin.
    pub off_set: Option<Point3D>,
}

impl Default for UniformRealValueGrid3DParams {
    fn default() -> UniformRealValueGrid3DParams {
        UniformRealValueGrid3DParams {
            spacing: 0.5,
            off_set: None,
        }
    }
}

impl UniformRealValueGrid3D {
    /// A grid of `dim_x` by `dim_y` by `dim_z`, which must be positive; a
    /// grid of 2^31 points or more is `Error::Value`. Mirrors
    /// `UniformRealValueGrid3D(dimX, dimY, dimZ)`.
    pub fn new(dim_x: f64, dim_y: f64, dim_z: f64) -> Result<UniformRealValueGrid3D> {
        UniformRealValueGrid3D::new_with_params(dim_x, dim_y, dim_z, &Default::default())
    }

    /// Mirrors `UniformRealValueGrid3D(dimX, dimY, dimZ, spacing, offSet)`.
    pub fn new_with_params(
        dim_x: f64,
        dim_y: f64,
        dim_z: f64,
        params: &UniformRealValueGrid3DParams,
    ) -> Result<UniformRealValueGrid3D> {
        Ok(UniformRealValueGrid3D::from_ffi(
            ffi::uniform_real_value_grid3d_new(
                dim_x,
                dim_y,
                dim_z,
                params.spacing,
                params.off_set.is_some(),
                &params.off_set.unwrap_or_default(),
            )?,
        ))
    }

    /// A copy of the values. Mirrors `GetOccupancyVect`, which returns the
    /// grid's own.
    pub fn get_occupancy_vect(&self) -> RealValueVect {
        RealValueVect::from_ffi(expect_ok(
            ffi::uniform_real_value_grid3d_get_occupancy_vect(self.raw()),
        ))
    }

    /// Whether the values match. Mirrors
    /// `UniformRealValueGrid3D.CompareVectors`.
    pub fn compare_vectors(&self, other: &UniformRealValueGrid3D) -> bool {
        expect_ok(ffi::uniform_real_value_grid3d_compare(
            self.raw(),
            other.raw(),
            b'v',
        ))
    }

    /// Whether `other` has the same size, spacing, and offset. Mirrors
    /// `UniformRealValueGrid3D.CompareParams`.
    pub fn compare_params(&self, other: &UniformRealValueGrid3D) -> bool {
        expect_ok(ffi::uniform_real_value_grid3d_compare(
            self.raw(),
            other.raw(),
            b'p',
        ))
    }

    /// Whether the parameters and values match. Mirrors
    /// `UniformRealValueGrid3D.CompareGrids`.
    pub fn compare_grids(&self, other: &UniformRealValueGrid3D) -> bool {
        expect_ok(ffi::uniform_real_value_grid3d_compare(
            self.raw(),
            other.raw(),
            b'g',
        ))
    }
}

/// A 1 Å grid of 1 by 1 by 1. Mirrors `UniformRealValueGrid3D()`.
impl Default for UniformRealValueGrid3D {
    fn default() -> UniformRealValueGrid3D {
        UniformRealValueGrid3D::from_ffi(expect_ok(ffi::uniform_real_value_grid3d_default()))
    }
}

/// Writes the grid in a format PyMOL reads. Mirrors
/// `rdGeometry.WriteGridToFile`.
pub fn write_grid_to_file(grid: &UniformGrid3D, filename: impl AsRef<Path>) -> Result<()> {
    Ok(ffi::write_grid_to_file(
        grid.raw(),
        path_bytes(filename.as_ref()),
    )?)
}

/// The Tversky index of two grids with the same parameters. Mirrors
/// `rdGeometry.TverskyIndex`.
pub fn tversky_index(
    grid1: &UniformGrid3D,
    grid2: &UniformGrid3D,
    alpha: f64,
    beta: f64,
) -> Result<f64> {
    Ok(ffi::tversky_index(grid1.raw(), grid2.raw(), alpha, beta)?)
}

/// The Tanimoto distance of two grids with the same parameters. Mirrors
/// `rdGeometry.TanimotoDistance`.
pub fn tanimoto_distance(grid1: &UniformGrid3D, grid2: &UniformGrid3D) -> Result<f64> {
    Ok(ffi::tanimoto_distance(grid1.raw(), grid2.raw())?)
}

/// The protrude distance of two grids with the same parameters. Mirrors
/// `rdGeometry.ProtrudeDistance`.
pub fn protrude_distance(grid1: &UniformGrid3D, grid2: &UniformGrid3D) -> Result<f64> {
    Ok(ffi::protrude_distance(grid1.raw(), grid2.raw())?)
}

/// The weight sum and the centroid of the grid's values within
/// `window_radius` of `pt`. A window radius of more than 26,754 grid
/// spacings is `Error::Value`. Mirrors `rdGeometry.ComputeGridCentroid`.
pub fn compute_grid_centroid(
    grid: &UniformGrid3D,
    pt: &Point3D,
    window_radius: f64,
) -> Result<(f64, Point3D)> {
    let mut weight_sum = 0.0;
    let centroid = ffi::compute_grid_centroid(grid.raw(), pt, window_radius, &mut weight_sum)?;
    Ok((weight_sum, centroid))
}

/// The centroids of the windows around the grid's full points whose fraction
/// of the window's possible value is at most `inclusion_fraction`. A window
/// radius of more than 26,754 grid spacings is `Error::Value`. Mirrors
/// `rdGeometry.FindGridTerminalPoints`.
pub fn find_grid_terminal_points(
    grid: &UniformGrid3D,
    window_radius: f64,
    inclusion_fraction: f64,
) -> Result<Vec<Point3D>> {
    Ok(ffi::find_grid_terminal_points(
        grid.raw(),
        window_radius,
        inclusion_fraction,
    )?)
}
