//! Mirrors `rdkit.Geometry.rdGeometry`.

pub use crate::ffi::geometry::ffi::Point3D;

impl Point3D {
    /// Mirrors `Point3D(x, y, z)`.
    pub fn new(x: f64, y: f64, z: f64) -> Point3D {
        Point3D { x, y, z }
    }
}
