// Copyright (c) 2026 Steven Kearnes
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
//! what it would mishandle: shapes without atoms or features (it reads past
//! an empty vector), shapes that overlap in none of its four first starting
//! poses, radii below about 0.16 Å, non-finite shape coordinates, and, when
//! both shapes have features of a type, an `opt_param` outside (0, 1]
//! (for each of these it can finish without a pose, and returns a matrix of
//! uninitialized memory); coordinates beyond 1e8; and a molecule's
//! `PUBCHEM_PHARMACOPHORE_FEATURES` counts that exceed its lines or tokens
//! (it allocates them first). Its global state is guarded by a lock, so
//! alignments in several threads take turns.

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::AsConformerMut;
use crate::error::expect_ok;
use crate::ffi::shape_align::ffi;
use crate::geometry::Point3D;

/// How a molecule becomes a shape. Mirrors
/// `rdShapeAlign.ShapeInputOptions` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeInputOptions {
    /// Adds features (colors) to the shape and scores them.
    pub use_colors: bool,
    /// Includes dummy atoms in the shape; hydrogens never are.
    pub include_dummies: bool,
    /// The radius of dummy atoms.
    pub dummy_radius: f64,
    /// The atoms that make up the shape; empty for all.
    pub atom_subset: Vec<u32>,
    /// Atoms no feature may use.
    pub not_color_atoms: Vec<u32>,
    /// Radii for some atoms, by atom index, instead of their elements'.
    pub atom_radii: Vec<(u32, f64)>,
    /// Features to use instead of the default ones: (type, position,
    /// radius).
    pub custom_features: Vec<(u32, Point3D, f64)>,
    /// Moves the shape to its principal axes, centered on the origin.
    pub normalize: bool,
}

impl Default for ShapeInputOptions {
    fn default() -> ShapeInputOptions {
        let data = expect_ok(ffi::pubchem_shape_input_options_default());
        ShapeInputOptions {
            use_colors: data.use_colors,
            include_dummies: data.include_dummies,
            dummy_radius: data.dummy_radius,
            atom_subset: data.atom_subset,
            not_color_atoms: data.not_color_atoms,
            atom_radii: Vec::new(),
            custom_features: Vec::new(),
            normalize: data.normalize,
        }
    }
}

impl ShapeInputOptions {
    fn to_ffi(&self) -> ffi::ShapeInputOptionsData {
        ffi::ShapeInputOptionsData {
            use_colors: self.use_colors,
            include_dummies: self.include_dummies,
            dummy_radius: self.dummy_radius,
            atom_subset: self.atom_subset.clone(),
            not_color_atoms: self.not_color_atoms.clone(),
            atom_radii_atoms: self.atom_radii.iter().map(|&(atom, _)| atom).collect(),
            atom_radii_values: self.atom_radii.iter().map(|&(_, radius)| radius).collect(),
            custom_feature_types: self.custom_features.iter().map(|f| f.0).collect(),
            custom_feature_positions: self.custom_features.iter().map(|f| f.1).collect(),
            custom_feature_radii: self.custom_features.iter().map(|f| f.2).collect(),
            normalize: self.normalize,
        }
    }
}

/// A molecule's shape, from `prepare_conformer`. Mirrors
/// `rdShapeAlign.ShapeInput`, read-only.
pub struct ShapeInput {
    ptr: UniquePtr<ffi::PubChemShapeInput>,
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

    /// The rotation, as 9 values by rows, that put the conformer in its
    /// principal axes. Mirrors `ShapeInput.inertialRot`.
    pub fn inertial_rot(&self) -> Vec<f64> {
        expect_ok(ffi::pubchem_shape_input_inertial_rot(&self.ptr))
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
    pub opts: ShapeInputOptions,
}

impl Default for PrepareConformerParams {
    fn default() -> PrepareConformerParams {
        PrepareConformerParams {
            conf_id: -1,
            opts: ShapeInputOptions::default(),
        }
    }
}

/// The shape of the molecule's first conformer, centered and in its
/// principal axes. Mirrors `rdShapeAlign.PrepareConformer`.
///
/// A missing or 2D conformer, an element without a van der Waals radius,
/// and features that do not parse are `Error::Value`.
pub fn prepare_conformer(mol: &Mol) -> Result<ShapeInput> {
    prepare_conformer_with_params(mol, &PrepareConformerParams::default())
}

/// Mirrors `rdShapeAlign.PrepareConformer(mol, confId, opts)`.
pub fn prepare_conformer_with_params(
    mol: &Mol,
    params: &PrepareConformerParams,
) -> Result<ShapeInput> {
    Ok(ShapeInput {
        ptr: ffi::pubchem_prepare_conformer(mol.as_ffi(), params.conf_id, &params.opts.to_ffi())?,
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
    let opts = ShapeInputOptions {
        use_colors: params.use_colors,
        ..Default::default()
    }
    .to_ffi();
    Ok(pair(ffi::pubchem_align_mol(
        ref_.as_ffi(),
        probe.pin_mut(),
        &opts,
        &opts,
        params.ref_conf_id,
        params.probe_conf_id,
        params.opt_param,
        params.max_preiters,
        params.max_postiters,
    )?))
}

/// Parameters for `align_mol_shape_opts_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolShapeOptsParams {
    /// The reference's conformer; -1 for the first.
    pub ref_conf_id: i32,
    /// The probe's conformer; -1 for the first.
    pub probe_conf_id: i32,
    /// The weight of the shape score against the color score.
    pub opt_param: f64,
    /// The most optimization steps for every starting pose.
    pub max_preiters: u32,
    /// The most further steps for the best poses.
    pub max_postiters: u32,
}

impl Default for AlignMolShapeOptsParams {
    fn default() -> AlignMolShapeOptsParams {
        AlignMolShapeOptsParams {
            ref_conf_id: -1,
            probe_conf_id: -1,
            opt_param: 1.0,
            max_preiters: 10,
            max_postiters: 30,
        }
    }
}

/// `align_mol` with options for each molecule's shape; colors count only
/// when both use them. Mirrors `rdShapeAlign.AlignMol(ref, probe,
/// refShapeOpts, probeShapeOpts)`.
pub fn align_mol_shape_opts(
    ref_: &Mol,
    probe: &mut Mol,
    ref_shape_opts: &ShapeInputOptions,
    probe_shape_opts: &ShapeInputOptions,
) -> Result<(f64, f64)> {
    align_mol_shape_opts_with_params(
        ref_,
        probe,
        ref_shape_opts,
        probe_shape_opts,
        &AlignMolShapeOptsParams::default(),
    )
}

/// Mirrors `rdShapeAlign.AlignMol(ref, probe, refShapeOpts,
/// probeShapeOpts, refConfId, probeConfId, opt_param, max_preiters,
/// max_postiters)`.
pub fn align_mol_shape_opts_with_params(
    ref_: &Mol,
    probe: &mut Mol,
    ref_shape_opts: &ShapeInputOptions,
    probe_shape_opts: &ShapeInputOptions,
    params: &AlignMolShapeOptsParams,
) -> Result<(f64, f64)> {
    Ok(pair(ffi::pubchem_align_mol(
        ref_.as_ffi(),
        probe.pin_mut(),
        &ref_shape_opts.to_ffi(),
        &probe_shape_opts.to_ffi(),
        params.ref_conf_id,
        params.probe_conf_id,
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
    /// Moves the probe into the reference's original frame, rather than the
    /// reference shape's.
    pub apply_ref_shift: bool,
}

impl Default for AlignMolShapeInputParams {
    fn default() -> AlignMolShapeInputParams {
        AlignMolShapeInputParams {
            probe_conf_id: -1,
            use_colors: true,
            opt_param: 1.0,
            max_preiters: 10,
            max_postiters: 30,
            apply_ref_shift: false,
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
/// opt_param, max_preiters, max_postiters, applyRefShift)`.
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
        params.apply_ref_shift,
    )?))
}

/// Parameters for `align_shapes_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignShapesParams {
    /// The weight of the shape score against the color score.
    pub opt_param: f64,
    /// The most optimization steps for every starting pose.
    pub max_preiters: u32,
    /// The most further steps for the best poses.
    pub max_postiters: u32,
}

impl Default for AlignShapesParams {
    fn default() -> AlignShapesParams {
        AlignShapesParams {
            opt_param: 1.0,
            max_preiters: 10,
            max_postiters: 30,
        }
    }
}

/// Aligns the probe shape to the reference shape, both centered on the
/// origin, and returns the shape and color scores and the 12 values of the
/// transformation (a 3 x 3 rotation by rows, then the translation). The
/// probe's coordinates move. Mirrors `rdShapeAlign.AlignShapes`, with the
/// errors of `align_mol`; on an error the probe is unchanged.
pub fn align_shapes(
    ref_shape: &ShapeInput,
    probe_shape: &mut ShapeInput,
) -> Result<(f64, f64, Vec<f64>)> {
    align_shapes_with_params(ref_shape, probe_shape, &AlignShapesParams::default())
}

/// Mirrors `rdShapeAlign.AlignShapes(refShape, probeShape, opt_param,
/// max_preiters, max_postiters)`.
pub fn align_shapes_with_params(
    ref_shape: &ShapeInput,
    probe_shape: &mut ShapeInput,
    params: &AlignShapesParams,
) -> Result<(f64, f64, Vec<f64>)> {
    let result = ffi::pubchem_align_shapes(
        &ref_shape.ptr,
        probe_shape.ptr.pin_mut(),
        params.opt_param,
        params.max_preiters,
        params.max_postiters,
    )?;
    Ok((result[0], result[1], result[2..].to_vec()))
}

/// Moves the conformer the probe shape came from as `align_shapes` moved
/// the shape (`matrix`), then by `final_rot` (9 values, by rows) and
/// `final_trans`. Mirrors `rdShapeAlign.TransformConformer`.
///
/// Vectors of other lengths, and matrix values beyond `f32`, are
/// `Error::Value`.
pub fn transform_conformer(
    final_trans: &[f64],
    final_rot: &[f64],
    matrix: &[f64],
    probe_shape: &ShapeInput,
    probe_conformer: &mut impl AsConformerMut,
) -> Result<()> {
    Ok(ffi::pubchem_transform_conformer_to_shape(
        final_trans,
        final_rot,
        matrix,
        &probe_shape.ptr,
        probe_conformer.conformer_mut(),
    )?)
}

/// Parameters for `score_mol_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ScoreMolParams {
    /// The first molecule's conformer; -1 for the first.
    pub mol1_conf_id: i32,
    /// The second molecule's conformer; -1 for the first.
    pub mol2_conf_id: i32,
}

impl Default for ScoreMolParams {
    fn default() -> ScoreMolParams {
        ScoreMolParams {
            mol1_conf_id: -1,
            mol2_conf_id: -1,
        }
    }
}

/// The shape and color scores of two molecules' first conformers, where
/// they are; colors count only when both options use them. Mirrors
/// `rdShapeAlign.ScoreMol`.
pub fn score_mol(
    mol1: &Mol,
    mol2: &Mol,
    mol1_shape_opts: &ShapeInputOptions,
    mol2_shape_opts: &ShapeInputOptions,
) -> Result<(f64, f64)> {
    score_mol_with_params(
        mol1,
        mol2,
        mol1_shape_opts,
        mol2_shape_opts,
        &ScoreMolParams::default(),
    )
}

/// Mirrors `rdShapeAlign.ScoreMol(mol1, mol2, mol1ShapeOpts,
/// mol2ShapeOpts, mol1ConfId, mol2ConfId)`.
pub fn score_mol_with_params(
    mol1: &Mol,
    mol2: &Mol,
    mol1_shape_opts: &ShapeInputOptions,
    mol2_shape_opts: &ShapeInputOptions,
    params: &ScoreMolParams,
) -> Result<(f64, f64)> {
    Ok(pair(ffi::pubchem_score_mol(
        mol1.as_ffi(),
        mol2.as_ffi(),
        &mol1_shape_opts.to_ffi(),
        &mol2_shape_opts.to_ffi(),
        params.mol1_conf_id,
        params.mol2_conf_id,
    )?))
}

/// The scores of a shape, made without normalizing, and a molecule's first
/// conformer, where it is. Mirrors `rdShapeAlign.ScoreMol(shape, mol,
/// molShapeOpts)`.
pub fn score_mol_shape_input(
    shape: &ShapeInput,
    mol: &Mol,
    mol_shape_opts: &ShapeInputOptions,
) -> Result<(f64, f64)> {
    score_mol_shape_input_with_mol_conf_id(shape, mol, mol_shape_opts, -1)
}

/// Mirrors `rdShapeAlign.ScoreMol(shape, mol, molShapeOpts, molConfId)`.
pub fn score_mol_shape_input_with_mol_conf_id(
    shape: &ShapeInput,
    mol: &Mol,
    mol_shape_opts: &ShapeInputOptions,
    mol_conf_id: i32,
) -> Result<(f64, f64)> {
    Ok(pair(ffi::pubchem_score_mol_shape_input(
        &shape.ptr,
        mol.as_ffi(),
        &mol_shape_opts.to_ffi(),
        mol_conf_id,
    )?))
}

/// The shape score of two shapes where they are, made without
/// normalizing; the color score is 0. Mirrors `rdShapeAlign.ScoreShape`.
pub fn score_shape(shape1: &ShapeInput, shape2: &ShapeInput) -> (f64, f64) {
    score_shape_with_use_colors(shape1, shape2, false)
}

/// Mirrors `rdShapeAlign.ScoreShape(shape1, shape2, useColors)`; with
/// colors, the color score of the features of the types both shapes have.
pub fn score_shape_with_use_colors(
    shape1: &ShapeInput,
    shape2: &ShapeInput,
    use_colors: bool,
) -> (f64, f64) {
    pair(expect_ok(ffi::pubchem_score_shape(
        &shape1.ptr,
        &shape2.ptr,
        use_colors,
    )))
}
