// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::error::expect_ok;
use crate::ffi::geometry::ffi;
use crate::{Error, Result};

pub use crate::ffi::geometry::ffi::{Point2D, Point3D};

/// Defines a point type's component-wise operators, as `RDGeom`'s points
/// have them.
macro_rules! point_ops {
    ($name:ident { $($field:ident),* }) => {
        /// Mirrors `+`.
        impl std::ops::Add for $name {
            type Output = $name;

            fn add(self, rhs: $name) -> $name {
                $name { $($field: self.$field + rhs.$field),* }
            }
        }

        /// Mirrors `-`.
        impl std::ops::Sub for $name {
            type Output = $name;

            fn sub(self, rhs: $name) -> $name {
                $name { $($field: self.$field - rhs.$field),* }
            }
        }

        /// Mirrors `+=`.
        impl std::ops::AddAssign for $name {
            fn add_assign(&mut self, rhs: $name) {
                *self = *self + rhs;
            }
        }

        /// Mirrors `-=`.
        impl std::ops::SubAssign for $name {
            fn sub_assign(&mut self, rhs: $name) {
                *self = *self - rhs;
            }
        }

        /// Mirrors `* scale`.
        impl std::ops::Mul<f64> for $name {
            type Output = $name;

            fn mul(self, scale: f64) -> $name {
                $name { $($field: self.$field * scale),* }
            }
        }

        /// Mirrors `/ scale`.
        impl std::ops::Div<f64> for $name {
            type Output = $name;

            fn div(self, scale: f64) -> $name {
                $name { $($field: self.$field / scale),* }
            }
        }

        /// Mirrors `*= scale`.
        impl std::ops::MulAssign<f64> for $name {
            fn mul_assign(&mut self, scale: f64) {
                *self = *self * scale;
            }
        }

        /// Mirrors `/= scale`.
        impl std::ops::DivAssign<f64> for $name {
            fn div_assign(&mut self, scale: f64) {
                *self = *self / scale;
            }
        }
    };
}

point_ops!(Point2D { x, y });
point_ops!(Point3D { x, y, z });

impl Point2D {
    /// Mirrors `Point2D(x, y)`.
    pub fn new(x: f64, y: f64) -> Point2D {
        Point2D { x, y }
    }

    /// Scales the point to length 1; a zero vector is `Error::Runtime`.
    /// Mirrors `Point2D.Normalize`.
    pub fn normalize(&mut self) -> Result<()> {
        *self = ffi::point2d_normalize(self)?;
        Ok(())
    }

    /// Mirrors `Point2D.Length`.
    pub fn length(&self) -> f64 {
        expect_ok(ffi::point2d_length(self))
    }

    /// Mirrors `Point2D.LengthSq`.
    pub fn length_sq(&self) -> f64 {
        expect_ok(ffi::point2d_length_sq(self))
    }

    /// Mirrors `Point2D.DotProduct`.
    pub fn dot_product(&self, other: &Point2D) -> f64 {
        expect_ok(ffi::point2d_dot_product(self, other))
    }

    /// The angle between the two vectors, from 0 to π; a zero vector is
    /// `Error::Runtime`. Mirrors `Point2D.AngleTo`.
    pub fn angle_to(&self, other: &Point2D) -> Result<f64> {
        Ok(ffi::point2d_angle_to(self, other)?)
    }

    /// The angle from this vector to the other, from 0 to 2π; a zero vector
    /// is `Error::Runtime`. Mirrors `Point2D.SignedAngleTo`.
    pub fn signed_angle_to(&self, other: &Point2D) -> Result<f64> {
        Ok(ffi::point2d_signed_angle_to(self, other)?)
    }

    /// The unit vector from this point to the other; the same point is
    /// `Error::Runtime`. Mirrors `Point2D.DirectionVector`.
    pub fn direction_vector(&self, other: &Point2D) -> Result<Point2D> {
        Ok(ffi::point2d_direction_vector(self, other)?)
    }
}

/// Drops `z`. Mirrors `Point2D(Point3D)`.
impl From<Point3D> for Point2D {
    fn from(p: Point3D) -> Point2D {
        Point2D { x: p.x, y: p.y }
    }
}

impl Point3D {
    /// Mirrors `Point3D(x, y, z)`.
    pub fn new(x: f64, y: f64, z: f64) -> Point3D {
        Point3D { x, y, z }
    }

    /// Scales the point to length 1; a zero vector is `Error::Runtime`.
    /// Mirrors `Point3D.Normalize`.
    pub fn normalize(&mut self) -> Result<()> {
        *self = ffi::point3d_normalize(self)?;
        Ok(())
    }

    /// Mirrors `Point3D.Length`.
    pub fn length(&self) -> f64 {
        expect_ok(ffi::point3d_length(self))
    }

    /// Mirrors `Point3D.LengthSq`.
    pub fn length_sq(&self) -> f64 {
        expect_ok(ffi::point3d_length_sq(self))
    }

    /// Mirrors `Point3D.DotProduct`.
    pub fn dot_product(&self, other: &Point3D) -> f64 {
        expect_ok(ffi::point3d_dot_product(self, other))
    }

    /// Mirrors `Point3D.CrossProduct`.
    pub fn cross_product(&self, other: &Point3D) -> Point3D {
        expect_ok(ffi::point3d_cross_product(self, other))
    }

    /// The angle between the two vectors, from 0 to π, or NaN for a zero
    /// vector. Mirrors `Point3D.AngleTo`.
    pub fn angle_to(&self, other: &Point3D) -> f64 {
        expect_ok(ffi::point3d_angle_to(self, other))
    }

    /// The angle from this vector to the other, from 0 to 2π, by the sign of
    /// the cross product's z. Mirrors `Point3D.SignedAngleTo`.
    pub fn signed_angle_to(&self, other: &Point3D) -> f64 {
        expect_ok(ffi::point3d_signed_angle_to(self, other))
    }

    /// The unit vector from this point to the other; the same point is
    /// `Error::Runtime`. Mirrors `Point3D.DirectionVector`.
    pub fn direction_vector(&self, other: &Point3D) -> Result<Point3D> {
        Ok(ffi::point3d_direction_vector(self, other)?)
    }

    /// Mirrors `Point3D.Distance`.
    pub fn distance(&self, pt2: &Point3D) -> f64 {
        (*self - *pt2).length()
    }
}

/// The dihedral angle of four points, from 0 to π. Mirrors
/// `rdGeometry.ComputeDihedralAngle`.
pub fn compute_dihedral_angle(pt1: &Point3D, pt2: &Point3D, pt3: &Point3D, pt4: &Point3D) -> f64 {
    expect_ok(ffi::compute_dihedral_angle(pt1, pt2, pt3, pt4))
}

/// The dihedral angle of four points, from -π to π. Mirrors
/// `rdGeometry.ComputeSignedDihedralAngle`.
pub fn compute_signed_dihedral_angle(
    pt1: &Point3D,
    pt2: &Point3D,
    pt3: &Point3D,
    pt4: &Point3D,
) -> f64 {
    expect_ok(ffi::compute_signed_dihedral_angle(pt1, pt2, pt3, pt4))
}

/// A point in any number of dimensions. Mirrors `rdGeometry.PointND`.
///
/// Operations on points of different dimensions are `Error::Value` in the
/// methods, and panic in the operators.
#[derive(Clone, Debug, PartialEq)]
pub struct PointND {
    coords: Vec<f64>,
}

impl PointND {
    /// The origin in `dim` dimensions. Mirrors `PointND(dim)`.
    pub fn new(dim: u32) -> PointND {
        PointND {
            coords: vec![0.0; dim as usize],
        }
    }

    /// A point from its coordinates.
    pub fn from_coords(coords: Vec<f64>) -> PointND {
        PointND { coords }
    }

    /// The coordinates.
    pub fn coords(&self) -> &[f64] {
        &self.coords
    }

    /// The number of dimensions. Mirrors `len(PointND)`.
    pub fn len(&self) -> u32 {
        u32::try_from(self.coords.len()).expect("fewer than 2^32 dimensions")
    }

    /// Whether the point has no dimensions.
    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }

    fn check_index(&self, idx: u32) -> Result<usize> {
        if idx < self.len() {
            Ok(idx as usize)
        } else {
            Err(Error::index(idx, self.len() as usize))
        }
    }

    /// Mirrors `PointND[idx]`; negative indices are not bound.
    pub fn get(&self, idx: u32) -> Result<f64> {
        Ok(self.coords[self.check_index(idx)?])
    }

    /// Mirrors `PointND[idx] = val`.
    pub fn set(&mut self, idx: u32, val: f64) -> Result<()> {
        let idx = self.check_index(idx)?;
        self.coords[idx] = val;
        Ok(())
    }

    /// Scales the point to length 1. Mirrors `PointND.Normalize`.
    pub fn normalize(&mut self) -> Result<()> {
        self.coords = ffi::point_nd_normalize(&self.coords)?;
        Ok(())
    }

    /// Mirrors `PointND.Length`.
    pub fn length(&self) -> f64 {
        expect_ok(ffi::point_nd_length(&self.coords))
    }

    /// Mirrors `PointND.LengthSq`.
    pub fn length_sq(&self) -> f64 {
        expect_ok(ffi::point_nd_length_sq(&self.coords))
    }

    /// Mirrors `PointND.DotProduct`.
    pub fn dot_product(&self, other: &PointND) -> Result<f64> {
        Ok(ffi::point_nd_dot_product(&self.coords, &other.coords)?)
    }

    /// The angle between the two vectors, from 0 to π. Mirrors
    /// `PointND.AngleTo`.
    pub fn angle_to(&self, other: &PointND) -> Result<f64> {
        Ok(ffi::point_nd_angle_to(&self.coords, &other.coords)?)
    }

    /// The unit vector from this point to the other. Mirrors
    /// `PointND.DirectionVector`.
    pub fn direction_vector(&self, other: &PointND) -> Result<PointND> {
        Ok(PointND {
            coords: ffi::point_nd_direction_vector(&self.coords, &other.coords)?,
        })
    }

    fn zip(&self, other: &PointND, f: impl Fn(f64, f64) -> f64) -> PointND {
        assert_eq!(self.len(), other.len(), "PointND dimension mismatch");
        PointND {
            coords: self
                .coords
                .iter()
                .zip(&other.coords)
                .map(|(a, b)| f(*a, *b))
                .collect(),
        }
    }

    fn scale(&self, f: impl Fn(f64) -> f64) -> PointND {
        PointND {
            coords: self.coords.iter().map(|a| f(*a)).collect(),
        }
    }
}

/// Mirrors `+`.
impl std::ops::Add for &PointND {
    type Output = PointND;

    fn add(self, rhs: &PointND) -> PointND {
        self.zip(rhs, |a, b| a + b)
    }
}

/// Mirrors `-`.
impl std::ops::Sub for &PointND {
    type Output = PointND;

    fn sub(self, rhs: &PointND) -> PointND {
        self.zip(rhs, |a, b| a - b)
    }
}

/// Mirrors `+=`.
impl std::ops::AddAssign<&PointND> for PointND {
    fn add_assign(&mut self, rhs: &PointND) {
        *self = &*self + rhs;
    }
}

/// Mirrors `-=`.
impl std::ops::SubAssign<&PointND> for PointND {
    fn sub_assign(&mut self, rhs: &PointND) {
        *self = &*self - rhs;
    }
}

/// Mirrors `* scale`.
impl std::ops::Mul<f64> for &PointND {
    type Output = PointND;

    fn mul(self, scale: f64) -> PointND {
        self.scale(|a| a * scale)
    }
}

/// Mirrors `/ scale`.
impl std::ops::Div<f64> for &PointND {
    type Output = PointND;

    fn div(self, scale: f64) -> PointND {
        self.scale(|a| a / scale)
    }
}

/// Mirrors `*= scale`.
impl std::ops::MulAssign<f64> for PointND {
    fn mul_assign(&mut self, scale: f64) {
        *self = &*self * scale;
    }
}

/// Mirrors `/= scale`.
impl std::ops::DivAssign<f64> for PointND {
    fn div_assign(&mut self, scale: f64) {
        *self = &*self / scale;
    }
}
