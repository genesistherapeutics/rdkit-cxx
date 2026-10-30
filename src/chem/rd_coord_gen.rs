// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdCoordGen`: 2D coordinates from Schrödinger's
//! CoordGen library.
//!
//! CoordGen keeps its templates and template directory in process-wide
//! globals, so its calls are serialized with the rest of 2D depiction
//! (CONVENTIONS.md U4). `add_coords` refuses, with `Error::Value`, what
//! CoordGen would mishandle: a `minimizer_precision` that is negative, not
//! finite, or above 3.5e8 (CoordGen casts six times it to `int`, and a
//! negative one makes it read a solution from an empty map); template
//! coordinates that are not finite or, times `coordgen_scaling`, beyond 5e6
//! (CoordGen computes in `float` and squares sums of their products);
//! molecules whose Morgan scores overflow the `int` CoordGen sums them in,
//! such as a chain of 40 carbons; and molecules whose formal charges'
//! magnitudes sum past `i32::MAX`, which CoordGen also sums in `int`
//! (CONVENTIONS.md U9).

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::coord_gen::ffi;
use crate::geometry::Point2D;
use crate::util::path_bytes;

/// How CoordGen lays out a molecule. Mirrors `rdCoordGen.CoordGenParams`
/// (A3). The template molecule, which Python's params keep alive, is owned
/// here.
#[derive(Clone, Debug)]
pub struct CoordGenParams {
    /// Template coordinates for some atoms, which CoordGen keeps them near
    /// (or at, with `dbg_use_fixed`). Indices past the molecule's atoms are
    /// ignored.
    pub coord_map: BTreeMap<u32, Point2D>,
    /// A molecule with one conformer whose match in the molecule takes its
    /// coordinates; ignored when it does not match or has another number of
    /// conformers.
    pub template_mol: Option<Mol>,
    /// CoordGen's length of a single bond in RDKit's units: CoordGen's
    /// coordinates are divided by it, and template coordinates multiplied.
    pub coordgen_scaling: f64,
    /// The directory of CoordGen's template file; empty for `$RDBASE/Data`.
    /// This RDKit build does not read template files, so only CoordGen's
    /// record of the directory changes.
    pub template_file_dir: PathBuf,
    /// How hard CoordGen works; see the `sketcher_*_precision` presets.
    pub minimizer_precision: f32,
    /// Keeps the template atoms near their coordinates (for debugging).
    pub dbg_use_constrained: bool,
    /// Fixes the template atoms at their coordinates (for debugging).
    pub dbg_use_fixed: bool,
    /// Treats single and double bonds between a metal and an atom with
    /// other bonds as zero-order.
    pub treat_nonterminal_bonds_to_metal_as_zobs: bool,
}

fn defaults() -> ffi::CoordGenParamsData {
    expect_ok(ffi::coord_gen_params_default())
}

impl Default for CoordGenParams {
    fn default() -> CoordGenParams {
        let data = defaults();
        CoordGenParams {
            coord_map: BTreeMap::new(),
            template_mol: None,
            coordgen_scaling: data.coordgen_scaling,
            template_file_dir: PathBuf::from(
                String::from_utf8_lossy(&data.template_file_dir).into_owned(),
            ),
            minimizer_precision: data.minimizer_precision,
            dbg_use_constrained: data.dbg_use_constrained,
            dbg_use_fixed: data.dbg_use_fixed,
            treat_nonterminal_bonds_to_metal_as_zobs: data.treat_nonterminal_bonds_to_metal_as_zobs,
        }
    }
}

impl CoordGenParams {
    /// Replaces the template coordinates. Mirrors
    /// `CoordGenParams.SetCoordMap`.
    pub fn set_coord_map(&mut self, coord_map: &BTreeMap<u32, Point2D>) {
        self.coord_map = coord_map.clone();
    }

    /// Sets a copy of `templ` as the template molecule. Mirrors
    /// `CoordGenParams.SetTemplateMol`.
    pub fn set_template_mol(&mut self, templ: &Mol) {
        self.template_mol = Some(templ.clone());
    }

    /// The fastest preset for `minimizer_precision`, and its default.
    /// Mirrors `CoordGenParams.sketcherCoarsePrecision`.
    pub fn sketcher_coarse_precision(&self) -> f32 {
        defaults().sketcher_coarse_precision
    }

    /// CoordGen's own default preset for `minimizer_precision`. Mirrors
    /// `CoordGenParams.sketcherStandardPrecision`.
    pub fn sketcher_standard_precision(&self) -> f32 {
        defaults().sketcher_standard_precision
    }

    /// The slowest, best preset for `minimizer_precision`. Mirrors
    /// `CoordGenParams.sketcherBestPrecision`.
    pub fn sketcher_best_precision(&self) -> f32 {
        defaults().sketcher_best_precision
    }

    /// A fast preset for `minimizer_precision`. Mirrors
    /// `CoordGenParams.sketcherQuickPrecision`.
    pub fn sketcher_quick_precision(&self) -> f32 {
        defaults().sketcher_quick_precision
    }

    fn to_ffi(&self) -> ffi::CoordGenParamsData {
        let defaults = defaults();
        ffi::CoordGenParamsData {
            coord_map_keys: self.coord_map.keys().copied().collect(),
            coord_map_values: self.coord_map.values().copied().collect(),
            coordgen_scaling: self.coordgen_scaling,
            template_file_dir: path_bytes(&self.template_file_dir).to_vec(),
            minimizer_precision: self.minimizer_precision,
            dbg_use_constrained: self.dbg_use_constrained,
            dbg_use_fixed: self.dbg_use_fixed,
            treat_nonterminal_bonds_to_metal_as_zobs: self.treat_nonterminal_bonds_to_metal_as_zobs,
            ..defaults
        }
    }
}

/// Replaces the molecule's conformers with one of 2D coordinates from
/// CoordGen, whose single bonds are 1 long. Mirrors `rdCoordGen.AddCoords`
/// with RDKit's default parameters, whose template directory
/// `set_default_template_file_dir` sets. `Error::Value` for a molecule
/// CoordGen would mishandle (see the module documentation).
pub fn add_coords(mol: &mut Mol) -> Result<()> {
    Ok(ffi::add_coords_default(mol.pin_mut())?)
}

/// Mirrors `rdCoordGen.AddCoords(mol, params)`. `Error::Value` also for
/// parameters CoordGen would mishandle.
pub fn add_coords_with_params(mol: &mut Mol, params: &CoordGenParams) -> Result<()> {
    // The bridge takes a molecule either way; an empty one stands in.
    let none = Mol::default();
    let template = params.template_mol.as_ref();
    Ok(ffi::add_coords(
        mol.pin_mut(),
        &params.to_ffi(),
        template.unwrap_or(&none).as_ffi(),
        template.is_some(),
    )?)
}

/// Sets the template directory `add_coords` uses. Mirrors
/// `rdCoordGen.SetDefaultTemplateFileDir`. `Error::Value` for a path with a
/// NUL byte.
pub fn set_default_template_file_dir(dir: impl AsRef<Path>) -> Result<()> {
    Ok(ffi::set_default_template_file_dir(path_bytes(
        dir.as_ref(),
    ))?)
}
