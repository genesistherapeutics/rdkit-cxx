// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Geometry.rdGeometry`.

mod grid;
mod point;

pub use grid::{
    SetSphereOccupancyParams, UniformGrid3D, UniformGrid3DParams, compute_grid_centroid,
    find_grid_terminal_points, protrude_distance, tanimoto_distance, tversky_index,
    write_grid_to_file,
};
pub use point::{Point2D, Point3D, PointND, compute_dihedral_angle, compute_signed_dihedral_angle};
