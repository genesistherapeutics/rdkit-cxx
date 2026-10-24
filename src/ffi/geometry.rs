// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A point in 2D. Mirrors `rdGeometry.Point2D`.
    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    struct Point2D {
        pub x: f64,
        pub y: f64,
    }

    /// A point in 3D. Mirrors `rdGeometry.Point3D`.
    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    struct Point3D {
        pub x: f64,
        pub y: f64,
        pub z: f64,
    }

    /// Atom or bond indices; nested lists cross as `Vec<IndexList>`.
    struct IndexList {
        indices: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/geometry.h");

        fn point2d_normalize(p: &Point2D) -> Result<Point2D>;
        fn point2d_length(p: &Point2D) -> Result<f64>;
        fn point2d_length_sq(p: &Point2D) -> Result<f64>;
        fn point2d_dot_product(p: &Point2D, other: &Point2D) -> Result<f64>;
        fn point2d_angle_to(p: &Point2D, other: &Point2D) -> Result<f64>;
        fn point2d_signed_angle_to(p: &Point2D, other: &Point2D) -> Result<f64>;
        fn point2d_direction_vector(p: &Point2D, other: &Point2D) -> Result<Point2D>;

        fn point3d_normalize(p: &Point3D) -> Result<Point3D>;
        fn point3d_length(p: &Point3D) -> Result<f64>;
        fn point3d_length_sq(p: &Point3D) -> Result<f64>;
        fn point3d_dot_product(p: &Point3D, other: &Point3D) -> Result<f64>;
        fn point3d_cross_product(p: &Point3D, other: &Point3D) -> Result<Point3D>;
        fn point3d_angle_to(p: &Point3D, other: &Point3D) -> Result<f64>;
        fn point3d_signed_angle_to(p: &Point3D, other: &Point3D) -> Result<f64>;
        fn point3d_direction_vector(p: &Point3D, other: &Point3D) -> Result<Point3D>;
        fn compute_dihedral_angle(
            pt1: &Point3D,
            pt2: &Point3D,
            pt3: &Point3D,
            pt4: &Point3D,
        ) -> Result<f64>;
        fn compute_signed_dihedral_angle(
            pt1: &Point3D,
            pt2: &Point3D,
            pt3: &Point3D,
            pt4: &Point3D,
        ) -> Result<f64>;

        fn point_nd_normalize(p: &[f64]) -> Result<Vec<f64>>;
        fn point_nd_length(p: &[f64]) -> Result<f64>;
        fn point_nd_length_sq(p: &[f64]) -> Result<f64>;
        fn point_nd_dot_product(p: &[f64], other: &[f64]) -> Result<f64>;
        fn point_nd_angle_to(p: &[f64], other: &[f64]) -> Result<f64>;
        fn point_nd_direction_vector(p: &[f64], other: &[f64]) -> Result<Vec<f64>>;
    }

    // Other bridges return vectors of these.
    impl Vec<Point3D> {}
    impl Vec<IndexList> {}
}

impl ffi::IndexList {
    pub(crate) fn into_nested(lists: Vec<ffi::IndexList>) -> Vec<Vec<u32>> {
        lists.into_iter().map(|list| list.indices).collect()
    }

    pub(crate) fn from_nested(lists: &[Vec<u32>]) -> Vec<ffi::IndexList> {
        lists
            .iter()
            .map(|indices| ffi::IndexList {
                indices: indices.clone(),
            })
            .collect()
    }
}
