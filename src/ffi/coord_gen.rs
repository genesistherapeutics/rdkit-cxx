// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::CoordGen::CoordGenParams` without its template molecule; the
    /// coordinate map as parallel keys and values (A7). The read-only
    /// precisions are read from the defaults and ignored on the way in.
    struct CoordGenParamsData {
        coord_map_keys: Vec<u32>,
        coord_map_values: Vec<Point2D>,
        coordgen_scaling: f64,
        template_file_dir: Vec<u8>,
        minimizer_precision: f32,
        dbg_use_constrained: bool,
        dbg_use_fixed: bool,
        treat_nonterminal_bonds_to_metal_as_zobs: bool,
        sketcher_coarse_precision: f32,
        sketcher_standard_precision: f32,
        sketcher_best_precision: f32,
        sketcher_quick_precision: f32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/coord_gen.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point2D = crate::ffi::geometry::ffi::Point2D;

        fn coord_gen_params_default() -> Result<CoordGenParamsData>;
        fn add_coords(
            mol: Pin<&mut RWMol>,
            data: &CoordGenParamsData,
            template_mol: &RWMol,
            has_template_mol: bool,
        ) -> Result<()>;
        fn add_coords_default(mol: Pin<&mut RWMol>) -> Result<()>;
        fn set_default_template_file_dir(dir: &[u8]) -> Result<()>;
    }
}
