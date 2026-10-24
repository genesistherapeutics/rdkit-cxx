// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A grid point's x, y, and z indices.
    struct GridIndices {
        xi: u32,
        yi: u32,
        zi: u32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/grid.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type UniformGrid3D;
        type UniformRealValueGrid3D;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type DiscreteValueType = crate::ffi::data_structs::ffi::DiscreteValueType;
        type DiscreteValueVect = crate::ffi::data_structs::ffi::DiscreteValueVect;
        type RealValueVect = crate::ffi::data_structs::ffi::RealValueVect;

        #[allow(clippy::too_many_arguments)]
        fn uniform_grid3d_new(
            dim_x: f64,
            dim_y: f64,
            dim_z: f64,
            spacing: f64,
            val_type: DiscreteValueType,
            has_offset: bool,
            offset: &Point3D,
        ) -> Result<UniquePtr<UniformGrid3D>>;
        fn uniform_grid3d_from_pkl(pkl: &[u8]) -> Result<UniquePtr<UniformGrid3D>>;
        fn uniform_grid3d_clone(grid: &UniformGrid3D) -> Result<UniquePtr<UniformGrid3D>>;
        fn uniform_grid3d_to_binary(grid: &UniformGrid3D) -> Result<Vec<u8>>;
        fn uniform_grid3d_get_grid_point_index(
            grid: &UniformGrid3D,
            point: &Point3D,
        ) -> Result<i32>;
        fn uniform_grid3d_get_grid_index(
            grid: &UniformGrid3D,
            xi: u32,
            yi: u32,
            zi: u32,
        ) -> Result<i32>;
        fn uniform_grid3d_get_grid_indices(grid: &UniformGrid3D, idx: u32) -> Result<GridIndices>;
        fn uniform_grid3d_get_val_point(grid: &UniformGrid3D, pt: &Point3D) -> Result<i32>;
        fn uniform_grid3d_get_val(grid: &UniformGrid3D, id: u32) -> Result<u32>;
        fn uniform_grid3d_set_val(grid: Pin<&mut UniformGrid3D>, id: u32, val: u32) -> Result<()>;
        fn uniform_grid3d_set_val_point(
            grid: Pin<&mut UniformGrid3D>,
            pt: &Point3D,
            val: u32,
        ) -> Result<()>;
        fn uniform_grid3d_get_grid_point_loc(
            grid: &UniformGrid3D,
            point_id: u32,
        ) -> Result<Point3D>;
        fn uniform_grid3d_get_size(grid: &UniformGrid3D) -> Result<u32>;
        fn uniform_grid3d_get_num_x(grid: &UniformGrid3D) -> Result<u32>;
        fn uniform_grid3d_get_num_y(grid: &UniformGrid3D) -> Result<u32>;
        fn uniform_grid3d_get_num_z(grid: &UniformGrid3D) -> Result<u32>;
        fn uniform_grid3d_get_offset(grid: &UniformGrid3D) -> Result<Point3D>;
        fn uniform_grid3d_get_spacing(grid: &UniformGrid3D) -> Result<f64>;
        fn uniform_grid3d_get_occupancy_vect(
            grid: &UniformGrid3D,
        ) -> Result<UniquePtr<DiscreteValueVect>>;
        fn uniform_grid3d_compare_params(
            grid: &UniformGrid3D,
            other: &UniformGrid3D,
        ) -> Result<bool>;
        fn uniform_grid3d_set_sphere_occupancy(
            grid: Pin<&mut UniformGrid3D>,
            center: &Point3D,
            radius: f64,
            step_size: f64,
            max_layers: i32,
            ignore_out_of_bound: bool,
        ) -> Result<()>;
        fn uniform_grid3d_assign_op(
            grid: Pin<&mut UniformGrid3D>,
            other: &UniformGrid3D,
            op: u8,
        ) -> Result<()>;
        fn write_grid_to_file(grid: &UniformGrid3D, filename: &[u8]) -> Result<()>;
        fn tversky_index(
            grid1: &UniformGrid3D,
            grid2: &UniformGrid3D,
            alpha: f64,
            beta: f64,
        ) -> Result<f64>;
        fn tanimoto_distance(grid1: &UniformGrid3D, grid2: &UniformGrid3D) -> Result<f64>;
        fn protrude_distance(grid1: &UniformGrid3D, grid2: &UniformGrid3D) -> Result<f64>;
        fn compute_grid_centroid(
            grid: &UniformGrid3D,
            pt: &Point3D,
            window_radius: f64,
            weight_sum: &mut f64,
        ) -> Result<Point3D>;
        fn find_grid_terminal_points(
            grid: &UniformGrid3D,
            window_radius: f64,
            inclusion_fraction: f64,
        ) -> Result<Vec<Point3D>>;
        fn uniform_real_value_grid3d_default() -> Result<UniquePtr<UniformRealValueGrid3D>>;
        fn uniform_real_value_grid3d_new(
            dim_x: f64,
            dim_y: f64,
            dim_z: f64,
            spacing: f64,
            has_offset: bool,
            offset: &Point3D,
        ) -> Result<UniquePtr<UniformRealValueGrid3D>>;
        fn uniform_real_value_grid3d_from_pkl(
            pkl: &[u8],
        ) -> Result<UniquePtr<UniformRealValueGrid3D>>;
        fn uniform_real_value_grid3d_clone(
            grid: &UniformRealValueGrid3D,
        ) -> Result<UniquePtr<UniformRealValueGrid3D>>;
        fn uniform_real_value_grid3d_to_binary(grid: &UniformRealValueGrid3D) -> Result<Vec<u8>>;
        fn uniform_real_value_grid3d_get_grid_point_index(
            grid: &UniformRealValueGrid3D,
            point: &Point3D,
        ) -> Result<i32>;
        fn uniform_real_value_grid3d_get_grid_index(
            grid: &UniformRealValueGrid3D,
            xi: u32,
            yi: u32,
            zi: u32,
        ) -> Result<i32>;
        fn uniform_real_value_grid3d_get_grid_indices(
            grid: &UniformRealValueGrid3D,
            idx: u32,
        ) -> Result<GridIndices>;
        fn uniform_real_value_grid3d_get_val_point(
            grid: &UniformRealValueGrid3D,
            pt: &Point3D,
        ) -> Result<f64>;
        fn uniform_real_value_grid3d_get_val(grid: &UniformRealValueGrid3D, id: u32)
        -> Result<f64>;
        fn uniform_real_value_grid3d_set_val(
            grid: Pin<&mut UniformRealValueGrid3D>,
            id: u32,
            val: f64,
        ) -> Result<()>;
        fn uniform_real_value_grid3d_set_val_point(
            grid: Pin<&mut UniformRealValueGrid3D>,
            pt: &Point3D,
            val: f64,
        ) -> Result<()>;
        fn uniform_real_value_grid3d_get_grid_point_loc(
            grid: &UniformRealValueGrid3D,
            point_id: u32,
        ) -> Result<Point3D>;
        fn uniform_real_value_grid3d_get_size(grid: &UniformRealValueGrid3D) -> Result<u32>;
        fn uniform_real_value_grid3d_get_num_x(grid: &UniformRealValueGrid3D) -> Result<u32>;
        fn uniform_real_value_grid3d_get_num_y(grid: &UniformRealValueGrid3D) -> Result<u32>;
        fn uniform_real_value_grid3d_get_num_z(grid: &UniformRealValueGrid3D) -> Result<u32>;
        fn uniform_real_value_grid3d_get_offset(grid: &UniformRealValueGrid3D) -> Result<Point3D>;
        fn uniform_real_value_grid3d_get_spacing(grid: &UniformRealValueGrid3D) -> Result<f64>;
        fn uniform_real_value_grid3d_get_occupancy_vect(
            grid: &UniformRealValueGrid3D,
        ) -> Result<UniquePtr<RealValueVect>>;
        fn uniform_real_value_grid3d_compare(
            grid: &UniformRealValueGrid3D,
            other: &UniformRealValueGrid3D,
            what: u8,
        ) -> Result<bool>;
        fn uniform_real_value_grid3d_assign_op(
            grid: Pin<&mut UniformRealValueGrid3D>,
            other: &UniformRealValueGrid3D,
            op: u8,
        ) -> Result<()>;
    }
}

// SAFETY: grids are plain data with no thread affinity, and their const
// methods write nothing.
unsafe impl Send for ffi::UniformGrid3D {}
unsafe impl Sync for ffi::UniformGrid3D {}
unsafe impl Send for ffi::UniformRealValueGrid3D {}
unsafe impl Sync for ffi::UniformRealValueGrid3D {}
