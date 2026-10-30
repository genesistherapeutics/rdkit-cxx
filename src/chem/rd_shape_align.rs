// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdShapeAlign`: PubChem's alignment of molecules by
//! the overlap of their Gaussian shapes and pharmacophore features
//! ("colors").
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_mol_file;
//! use rdkit_cxx::chem::rd_shape_align::align_mol;
//!
//! let reference = mol_from_mol_file("ref.mol").unwrap();
//! let mut probe = mol_from_mol_file("probe.mol").unwrap();
//! let (shape, color) = align_mol(&reference, &mut probe).unwrap();
//! println!("{shape} {color}");
//! ```
//!
//! A `ShapeInput`'s vectors are read-only here: the alignment indexes one
//! by the others and divides by the self overlaps without checks, so Python's
//! setters, which can make them disagree, are not bound.
//!
//! PubChem's code trusts its inputs, so these refuse with `Error::Value`
//! what it would mishandle: a heavy atom after a hydrogen or dummy atom
//! (RDKit 2024.09 writes each heavy atom's coordinates at its own index in
//! an array of the heavy atoms and features, past its end); shapes without
//! atoms or features (it reads past an empty vector), shapes that overlap in
//! none of its four first starting poses, and, when both shapes have
//! features of a type, an `opt_param` outside (0, 1] (for each of these it
//! can finish without a pose, and returns a matrix of uninitialized memory);
//! coordinates beyond 1e8; and a molecule's `PUBCHEM_PHARMACOPHORE_FEATURES`
//! counts that exceed its lines or tokens (it allocates them first). Its
//! global state is guarded by a lock, so alignments in several threads take
//! turns. RDKit 2024.09 computes a shape's self overlaps under whatever
//! distance cutoff its last alignment left; these always use the cutoff the
//! alignment sets, as `AlignMol` does.

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::shape_align::ffi;

/// A molecule's shape, from `prepare_conformer`. Mirrors
/// `rdShapeAlign.ShapeInput`, read-only.
pub struct ShapeInput {
    ptr: UniquePtr<ffi::PubChemShapeInput>,
}

impl std::fmt::Debug for ShapeInput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ShapeInput").finish_non_exhaustive()
    }
}

impl Clone for ShapeInput {
    fn clone(&self) -> ShapeInput {
        ShapeInput {
            ptr: expect_ok(ffi::pubchem_shape_input_clone(&self.ptr)),
        }
    }
}

impl ShapeInput {
    /// The x, y, z coordinates of each atom, then each feature. Mirrors
    /// `ShapeInput.coord`.
    pub fn coord(&self) -> Vec<f64> {
        expect_ok(ffi::pubchem_shape_input_coord(&self.ptr))
            .into_iter()
            .map(f64::from)
            .collect()
    }

    /// The Gaussian width of each atom and feature. Mirrors
    /// `ShapeInput.alpha_vector`.
    pub fn alpha_vector(&self) -> Vec<f64> {
        expect_ok(ffi::pubchem_shape_input_alpha_vector(&self.ptr))
    }

    /// The type of each atom (0) and feature (above 0). Mirrors
    /// `ShapeInput.atom_type_vector`.
    pub fn atom_type_vector(&self) -> Vec<u32> {
        expect_ok(ffi::pubchem_shape_input_atom_type_vector(&self.ptr))
    }

    /// The positions of the atoms, rather than features, among the points.
    /// Mirrors `ShapeInput.volumeAtomIndexVector`.
    pub fn volume_atom_index_vector(&self) -> Vec<u32> {
        expect_ok(ffi::pubchem_shape_input_volume_atom_index_vector(&self.ptr))
    }

    /// The translation of the conformer's center to the origin. Mirrors
    /// `ShapeInput.shift`.
    pub fn shift(&self) -> Vec<f64> {
        expect_ok(ffi::pubchem_shape_input_shift(&self.ptr))
    }

    /// The shape's overlap with itself. Mirrors `ShapeInput.sov`.
    pub fn sov(&self) -> f64 {
        expect_ok(ffi::pubchem_shape_input_sov(&self.ptr))
    }

    /// The features' overlap with themselves. Mirrors `ShapeInput.sof`.
    pub fn sof(&self) -> f64 {
        expect_ok(ffi::pubchem_shape_input_sof(&self.ptr))
    }
}

/// Parameters for `prepare_conformer_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct PrepareConformerParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// Adds pharmacophore features (colors): the molecule's
    /// `PUBCHEM_PHARMACOPHORE_FEATURES`, or RDKit's patterns' matches.
    pub use_colors: bool,
}

impl Default for PrepareConformerParams {
    fn default() -> PrepareConformerParams {
        PrepareConformerParams {
            conf_id: -1,
            use_colors: true,
        }
    }
}

/// The shape of the molecule's first conformer, centered on its heavy
/// atoms. Mirrors `rdShapeAlign.PrepareConformer`.
///
/// A missing or 2D conformer, an element without a van der Waals radius,
/// and features that do not parse are `Error::Value`, as is what the module
/// documentation lists.
pub fn prepare_conformer(mol: &Mol) -> Result<ShapeInput> {
    prepare_conformer_with_params(mol, &PrepareConformerParams::default())
}

/// Mirrors `rdShapeAlign.PrepareConformer(mol, confId, useColors)`.
pub fn prepare_conformer_with_params(
    mol: &Mol,
    params: &PrepareConformerParams,
) -> Result<ShapeInput> {
    Ok(ShapeInput {
        ptr: ffi::pubchem_prepare_conformer(mol.as_ffi(), params.conf_id, params.use_colors)?,
    })
}

fn pair(scores: Vec<f64>) -> (f64, f64) {
    (scores[0], scores[1])
}

/// Parameters for `align_mol_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolParams {
    /// The reference's conformer; -1 for the first.
    pub ref_conf_id: i32,
    /// The probe's conformer; -1 for the first.
    pub probe_conf_id: i32,
    /// Adds features (colors) to the shapes and scores them.
    pub use_colors: bool,
    /// The weight of the shape score against the color score in choosing
    /// poses: 1 for shape only.
    pub opt_param: f64,
    /// The most optimization steps for every starting pose.
    pub max_preiters: u32,
    /// The most further steps for the best poses.
    pub max_postiters: u32,
}

impl Default for AlignMolParams {
    fn default() -> AlignMolParams {
        AlignMolParams {
            ref_conf_id: -1,
            probe_conf_id: -1,
            use_colors: true,
            opt_param: 1.0,
            max_preiters: 10,
            max_postiters: 30,
        }
    }
}

/// Aligns the probe's first conformer to the reference's and returns the
/// shape and color Tanimoto scores (the color score is 0 for an
/// `opt_param` of 1). Mirrors `rdShapeAlign.AlignMol`.
///
/// The probe's conformer moves, and its `shape_align_shape_tanimoto` and
/// `shape_align_color_tanimoto` properties record the scores. Besides
/// `prepare_conformer`'s errors, the module documentation lists what is
/// `Error::Value`.
pub fn align_mol(ref_: &Mol, probe: &mut Mol) -> Result<(f64, f64)> {
    align_mol_with_params(ref_, probe, &AlignMolParams::default())
}

/// Mirrors `rdShapeAlign.AlignMol(ref, probe, refConfId, probeConfId,
/// useColors, opt_param, max_preiters, max_postiters)`.
pub fn align_mol_with_params(
    ref_: &Mol,
    probe: &mut Mol,
    params: &AlignMolParams,
) -> Result<(f64, f64)> {
    Ok(pair(ffi::pubchem_align_mol(
        ref_.as_ffi(),
        probe.pin_mut(),
        params.ref_conf_id,
        params.probe_conf_id,
        params.use_colors,
        params.opt_param,
        params.max_preiters,
        params.max_postiters,
    )?))
}

/// Parameters for `align_mol_shape_input_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolShapeInputParams {
    /// The probe's conformer; -1 for the first.
    pub probe_conf_id: i32,
    /// Adds features (colors) to the probe's shape and scores them.
    pub use_colors: bool,
    /// The weight of the shape score against the color score.
    pub opt_param: f64,
    /// The most optimization steps for every starting pose.
    pub max_preiters: u32,
    /// The most further steps for the best poses.
    pub max_postiters: u32,
}

impl Default for AlignMolShapeInputParams {
    fn default() -> AlignMolShapeInputParams {
        AlignMolShapeInputParams {
            probe_conf_id: -1,
            use_colors: true,
            opt_param: 1.0,
            max_preiters: 10,
            max_postiters: 30,
        }
    }
}

/// Aligns the probe's first conformer to a reference shape, which is
/// centered on the origin. Mirrors `rdShapeAlign.AlignMol(refShape,
/// probe)`, with the errors of `align_mol`.
pub fn align_mol_shape_input(ref_shape: &ShapeInput, probe: &mut Mol) -> Result<(f64, f64)> {
    align_mol_shape_input_with_params(ref_shape, probe, &AlignMolShapeInputParams::default())
}

/// Mirrors `rdShapeAlign.AlignMol(refShape, probe, probeConfId, useColors,
/// opt_param, max_preiters, max_postiters)`.
pub fn align_mol_shape_input_with_params(
    ref_shape: &ShapeInput,
    probe: &mut Mol,
    params: &AlignMolShapeInputParams,
) -> Result<(f64, f64)> {
    Ok(pair(ffi::pubchem_align_mol_shape_input(
        &ref_shape.ptr,
        probe.pin_mut(),
        params.probe_conf_id,
        params.use_colors,
        params.opt_param,
        params.max_preiters,
        params.max_postiters,
    )?))
}
