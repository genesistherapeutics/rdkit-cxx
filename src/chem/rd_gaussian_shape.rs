// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdGaussianShape`: overlays of molecules by the
//! overlap of Gaussian shapes and pharmacophore features ("colors"), scored
//! by Tversky similarity.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_mol_file;
//! use rdkit_cxx::chem::rd_gaussian_shape::align_mol;
//!
//! let reference = mol_from_mol_file("ref.mol").unwrap();
//! let mut fit = mol_from_mol_file("fit.mol").unwrap();
//! let (combo, shape, color) = align_mol(&reference, &mut fit).unwrap();
//! println!("{combo} {shape} {color}");
//! ```
//!
//! A shape is made from all of a molecule's conformers unless a conformer is
//! given. RDKit trusts what it is given, so these refuse what it would
//! mishandle: atom indices out of range (`Error::Index`); and with
//! `Error::Value`, molecules without conformers, coordinates that are not
//! finite or exceed 1e8, radii outside 0.001 to 1000 in magnitude, shapes
//! without an atom other than hydrogen, and custom features of other counts
//! for other conformers. Overlays also need Tversky weights `sim_alpha` and
//! `sim_beta` in [0, 1] (not both 0), an `opt_param` in [0, 1], and an
//! `n_steps` of at least 10: otherwise RDKit can compare NaN scores in a
//! sort, return a score it never set, or count its steps down past 0.
//!
//! Making a shape from a molecule with `atom_radii` does not change the
//! molecule, unlike Python, where the radii stay on its atoms as
//! `BespokeRadius` properties that later shapes of it read.

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::{expect_ok, to_string};
use crate::ffi::gaussian_shape::ffi;
use crate::geometry::Point3D;

rdkit_enum! {
    /// The starting orientations of an overlay. Mirrors
    /// `rdGaussianShape.StartMode`.
    StartMode: ffi::StartMode {
        Rotate0 = ROTATE_0,
        Rotate180 = ROTATE_180,
        Rotate180Wiggle = ROTATE_180_WIGGLE,
        Rotate45 = ROTATE_45,
        Rotate0Fragment = ROTATE_0_FRAGMENT,
        Rotate180Fragment = ROTATE_180_FRAGMENT,
        Rotate45Fragment = ROTATE_45_FRAGMENT,
        ALaPubchem = A_LA_PUBCHEM,
    }
}

rdkit_enum! {
    /// What drives an overlay's optimization. Mirrors
    /// `rdGaussianShape.OptimMode`.
    OptimMode: ffi::OptimMode {
        ShapeOnly = SHAPE_ONLY,
        ShapePlusColorScore = SHAPE_PLUS_COLOR_SCORE,
        ShapePlusColor = SHAPE_PLUS_COLOR,
    }
}

/// A custom feature: its type, position, radius, and the atoms it came
/// from (which may be empty).
pub type CustomFeature = (u32, Point3D, f64, Vec<u32>);

/// What `ShapeInput::best_similarity` returns: the best scores
/// (combination, shape, color), the numbers of the two shapes that reach
/// them, and the 16 values of the fit's transformation, by rows.
pub type BestSimilarityResult = ((f64, f64, f64), u32, u32, Vec<f64>);

/// What `score_molecule_all_conformers` returns: the combination scores by
/// reference and fit conformer, the best pair's conformers, and its
/// transformation's 16 values, by rows.
pub type ConformerScores = (Vec<Vec<f64>>, i32, i32, Vec<f64>);

/// How a molecule becomes a shape. Mirrors
/// `rdGaussianShape.ShapeInputOptions` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeInputOptions {
    /// Adds color features, RDKit's pharmacophore features unless
    /// `custom_features` are given.
    pub use_colors: bool,
    /// Gives every atom a carbon's radius: faster, slightly less accurate.
    pub all_carbon_radii: bool,
    /// The atoms that make up the shape; empty for all.
    pub atom_subset: Vec<u32>,
    /// Features to use instead of the default ones, one list per conformer.
    pub custom_features: Vec<Vec<CustomFeature>>,
    /// Radii for some atoms, by atom index, instead of their elements'.
    pub atom_radii: Vec<(u32, f64)>,
    /// With several conformers, drops shapes until none is more similar to
    /// another than this; -1 for none.
    pub shape_prune_threshold: f64,
    /// Sorts the shapes by decreasing volume.
    pub sort_shapes: bool,
    /// Includes dummy atoms in the shape.
    pub include_dummies: bool,
}

impl Default for ShapeInputOptions {
    fn default() -> ShapeInputOptions {
        let data = expect_ok(ffi::gaussian_shape_input_options_default());
        ShapeInputOptions {
            use_colors: data.use_colors,
            all_carbon_radii: data.all_carbon_radii,
            atom_subset: data.atom_subset,
            custom_features: Vec::new(),
            atom_radii: Vec::new(),
            shape_prune_threshold: data.shape_prune_threshold,
            sort_shapes: data.sort_shapes,
            include_dummies: data.include_dummies,
        }
    }
}

impl ShapeInputOptions {
    fn to_ffi(&self) -> ffi::GaussianShapeInputOptionsData {
        let features = self.custom_features.iter().flatten();
        ffi::GaussianShapeInputOptionsData {
            use_colors: self.use_colors,
            custom_feature_counts: self
                .custom_features
                .iter()
                .map(|f| f.len() as u32)
                .collect(),
            custom_feature_types: features.clone().map(|f| f.0).collect(),
            custom_feature_positions: features.clone().map(|f| f.1).collect(),
            custom_feature_radii: features.clone().map(|f| f.2).collect(),
            custom_feature_atom_counts: features.clone().map(|f| f.3.len() as u32).collect(),
            custom_feature_atoms: features.flat_map(|f| f.3.iter().copied()).collect(),
            atom_subset: self.atom_subset.clone(),
            atom_radii_atoms: self.atom_radii.iter().map(|&(atom, _)| atom).collect(),
            atom_radii_values: self.atom_radii.iter().map(|&(_, radius)| radius).collect(),
            all_carbon_radii: self.all_carbon_radii,
            shape_prune_threshold: self.shape_prune_threshold,
            sort_shapes: self.sort_shapes,
            include_dummies: self.include_dummies,
        }
    }
}

/// How shapes are overlaid and scored. Mirrors
/// `rdGaussianShape.ShapeOverlayOptions` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeOverlayOptions {
    /// The starting orientations of the fit.
    pub start_mode: StartMode,
    /// What the optimization maximizes.
    pub optim_mode: OptimMode,
    /// The Tversky weight of the reference's volume outside the overlap;
    /// with `sim_beta` 1, the score is a Tanimoto similarity.
    pub sim_alpha: f64,
    /// The Tversky weight of the fit's volume outside the overlap.
    pub sim_beta: f64,
    /// The weight of the color score against the shape score.
    pub opt_param: f64,
    /// The most optimization steps.
    pub n_steps: i32,
    /// Moves the shapes to their principal axes, centered on the origin,
    /// before the overlay.
    pub normalize: bool,
    /// Ignores pairs of atoms farther apart than `dist_cutoff`.
    pub use_dist_cutoff: bool,
    /// The distance, in Å, beyond which pairs are ignored.
    pub dist_cutoff: f64,
    /// Stops when an optimization step changes the score by less.
    pub shape_convergence_criterion: f64,
}

impl Default for ShapeOverlayOptions {
    fn default() -> ShapeOverlayOptions {
        let data = expect_ok(ffi::shape_overlay_options_default());
        ShapeOverlayOptions {
            start_mode: StartMode::from_ffi(data.start_mode),
            optim_mode: OptimMode::from_ffi(data.optim_mode),
            sim_alpha: data.sim_alpha,
            sim_beta: data.sim_beta,
            opt_param: data.opt_param,
            n_steps: data.n_steps,
            normalize: data.normalize,
            use_dist_cutoff: data.use_dist_cutoff,
            dist_cutoff: data.dist_cutoff,
            shape_convergence_criterion: data.shape_convergence_criterion,
        }
    }
}

impl ShapeOverlayOptions {
    fn to_ffi(&self) -> ffi::ShapeOverlayOptionsData {
        ffi::ShapeOverlayOptionsData {
            start_mode: self.start_mode.to_ffi(),
            optim_mode: self.optim_mode.to_ffi(),
            sim_alpha: self.sim_alpha,
            sim_beta: self.sim_beta,
            opt_param: self.opt_param,
            n_steps: self.n_steps,
            normalize: self.normalize,
            use_dist_cutoff: self.use_dist_cutoff,
            dist_cutoff: self.dist_cutoff,
            shape_convergence_criterion: self.shape_convergence_criterion,
        }
    }
}

fn triple(scores: &[f64]) -> (f64, f64, f64) {
    (scores[0], scores[1], scores[2])
}

fn five(scores: &[f64]) -> (f64, f64, f64, f64, f64) {
    (scores[0], scores[1], scores[2], scores[3], scores[4])
}

/// A molecule's shapes, one per conformer, of which one, the active shape,
/// takes part in overlays. Mirrors `rdGaussianShape.ShapeInput`.
pub struct ShapeInput {
    ptr: UniquePtr<ffi::GaussianShapeInput>,
}

impl Clone for ShapeInput {
    fn clone(&self) -> ShapeInput {
        ShapeInput {
            ptr: expect_ok(ffi::gaussian_shape_input_clone(&self.ptr)),
        }
    }
}

/// Parameters for `ShapeInput::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeInputParams {
    /// The conformer; -1 for all of them.
    pub conf_id: i32,
    pub shape_opt: ShapeInputOptions,
    /// The options whose distance cutoff the shapes' volumes use.
    pub overlay_opts: ShapeOverlayOptions,
}

impl Default for ShapeInputParams {
    fn default() -> ShapeInputParams {
        ShapeInputParams {
            conf_id: -1,
            shape_opt: ShapeInputOptions::default(),
            overlay_opts: ShapeOverlayOptions::default(),
        }
    }
}

/// Parameters for `ShapeInput::shape_to_mol_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeToMolParams {
    /// Adds the color features, as xenon atoms.
    pub include_colors: bool,
    /// Builds the molecule's bonds from its SMILES; otherwise only carbon
    /// atoms at the positions.
    pub with_bonds: bool,
}

impl Default for ShapeToMolParams {
    fn default() -> ShapeToMolParams {
        ShapeToMolParams {
            include_colors: false,
            with_bonds: true,
        }
    }
}

/// Parameters for `ShapeInput::best_similarity_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct BestSimilarityParams {
    /// Stops at the first overlay that scores above this; -1 for none.
    pub threshold: f64,
    pub overlay_opts: ShapeOverlayOptions,
}

impl Default for BestSimilarityParams {
    fn default() -> BestSimilarityParams {
        BestSimilarityParams {
            threshold: -1.0,
            overlay_opts: ShapeOverlayOptions::default(),
        }
    }
}

impl ShapeInput {
    /// The shapes of all the molecule's conformers, with colors. Mirrors
    /// `ShapeInput(mol)`; the module documentation lists the errors.
    pub fn new(mol: &Mol) -> Result<ShapeInput> {
        ShapeInput::new_with_params(mol, &ShapeInputParams::default())
    }

    /// Mirrors `ShapeInput(mol, confId, shapeOpt, overlayOpts)`.
    pub fn new_with_params(mol: &Mol, params: &ShapeInputParams) -> Result<ShapeInput> {
        Ok(ShapeInput {
            ptr: ffi::gaussian_shape_input_new(
                mol.as_ffi(),
                params.conf_id,
                &params.shape_opt.to_ffi(),
                &params.overlay_opts.to_ffi(),
            )?,
        })
    }

    /// The SMILES of the molecule (or atom subset) the shape came from.
    /// Mirrors `ShapeInput.GetSmiles`.
    pub fn get_smiles(&self) -> Result<String> {
        to_string(expect_ok(ffi::gaussian_shape_input_get_smiles(&self.ptr)))
    }

    /// Makes another shape the active one; one out of range is
    /// `Error::Index`. Mirrors `ShapeInput.setActiveShape`.
    pub fn set_active_shape(&mut self, arg: u32) -> Result<()> {
        Ok(ffi::gaussian_shape_input_set_active_shape(
            self.ptr.pin_mut(),
            arg,
        )?)
    }

    /// Mirrors `ShapeInput.getActiveShape`.
    pub fn get_active_shape(&self) -> u32 {
        expect_ok(ffi::gaussian_shape_input_get_active_shape(&self.ptr))
    }

    /// The atoms in each shape. Mirrors `ShapeInput.NumAtoms`.
    pub fn num_atoms(&self) -> u32 {
        expect_ok(ffi::gaussian_shape_input_num_atoms(&self.ptr))
    }

    /// The color features in each shape. Mirrors `ShapeInput.NumFeatures`.
    pub fn num_features(&self) -> u32 {
        expect_ok(ffi::gaussian_shape_input_num_features(&self.ptr))
    }

    /// The shapes: one per conformer, unless some were pruned. Mirrors
    /// `ShapeInput.NumShapes`.
    pub fn num_shapes(&self) -> u32 {
        expect_ok(ffi::gaussian_shape_input_num_shapes(&self.ptr))
    }

    /// The volume of the active shape's atoms. Mirrors
    /// `ShapeInput.ShapeVolume`.
    pub fn shape_volume(&self) -> f64 {
        expect_ok(ffi::gaussian_shape_input_shape_volume(&self.ptr))
    }

    /// The volume of the active shape's color features. Mirrors
    /// `ShapeInput.ColorVolume`.
    pub fn color_volume(&self) -> f64 {
        expect_ok(ffi::gaussian_shape_input_color_volume(&self.ptr))
    }

    /// Moves the active shape to its principal axes, centered on the origin.
    /// Mirrors `ShapeInput.NormalizeCoords`.
    pub fn normalize_coords(&mut self) {
        expect_ok(ffi::gaussian_shape_input_normalize_coords(
            self.ptr.pin_mut(),
        ))
    }

    /// A molecule, with bonds, at the active shape's coordinates. Mirrors
    /// `ShapeInput.ShapeToMol`.
    ///
    /// The molecule is not sanitized. One whose atoms are not the shape's
    /// (as when dummies are left out) is `Error::Runtime`.
    pub fn shape_to_mol(&self) -> Result<Mol> {
        self.shape_to_mol_with_params(&ShapeToMolParams::default())
    }

    /// Mirrors `ShapeInput.ShapeToMol(includeColors, withBonds)`.
    pub fn shape_to_mol_with_params(&self, params: &ShapeToMolParams) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::gaussian_shape_input_shape_to_mol(
            &self.ptr,
            params.include_colors,
            params.with_bonds,
        )?))
    }

    /// Overlays every pair of this shape's and the fit's shapes and returns
    /// the best scores (combination, shape, color), the two shapes' numbers,
    /// and the 16 values of the fit's transformation, by rows. Mirrors
    /// `ShapeInput.BestSimilarity`, with the overlay errors of the module
    /// documentation.
    pub fn best_similarity(&mut self, fit_shape: &ShapeInput) -> Result<BestSimilarityResult> {
        self.best_similarity_with_params(fit_shape, &BestSimilarityParams::default())
    }

    /// Mirrors `ShapeInput.BestSimilarity(fitShape, threshold,
    /// overlayOpts)`. When no overlay beats the threshold, the scores are -1
    /// and the shape numbers `u32::MAX`.
    pub fn best_similarity_with_params(
        &mut self,
        fit_shape: &ShapeInput,
        params: &BestSimilarityParams,
    ) -> Result<BestSimilarityResult> {
        let result = ffi::gaussian_shape_input_best_similarity(
            self.ptr.pin_mut(),
            &fit_shape.ptr,
            params.threshold,
            &params.overlay_opts.to_ffi(),
        )?;
        Ok((
            triple(&result),
            result[3] as u32,
            result[4] as u32,
            result[5..].to_vec(),
        ))
    }

    /// The best score any overlay of a shape here with one of the fit's
    /// could reach. Mirrors `ShapeInput.MaxPossibleSimilarity`.
    pub fn max_possible_similarity(&self, fit_shape: &ShapeInput) -> f64 {
        self.max_possible_similarity_with_overlay_opts(fit_shape, &ShapeOverlayOptions::default())
    }

    /// Mirrors `ShapeInput.MaxPossibleSimilarity(fitShape, overlayOpts)`.
    pub fn max_possible_similarity_with_overlay_opts(
        &self,
        fit_shape: &ShapeInput,
        overlay_opts: &ShapeOverlayOptions,
    ) -> f64 {
        expect_ok(ffi::gaussian_shape_input_max_possible_similarity(
            &self.ptr,
            &fit_shape.ptr,
            &overlay_opts.to_ffi(),
        ))
    }
}

/// Parameters for `align_mol_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolParams {
    pub ref_opts: ShapeInputOptions,
    pub fit_opts: ShapeInputOptions,
    pub overlay_opts: ShapeOverlayOptions,
    /// The reference's conformer; -1 for all of them, of which the active
    /// shape is the first (after sorting).
    pub ref_conf_id: i32,
    /// The fit's conformer; -1 for all of them, of which the active shape
    /// is used and the first conformer moves.
    pub fit_conf_id: i32,
}

impl Default for AlignMolParams {
    fn default() -> AlignMolParams {
        AlignMolParams {
            ref_opts: ShapeInputOptions::default(),
            fit_opts: ShapeInputOptions::default(),
            overlay_opts: ShapeOverlayOptions::default(),
            ref_conf_id: -1,
            fit_conf_id: -1,
        }
    }
}

/// Overlays the fit's shape on the reference's, moves the fit's conformer
/// there, and returns the scores (combination, shape, color). Mirrors
/// `rdGaussianShape.AlignMol`, with the errors of the module documentation.
pub fn align_mol(ref_: &Mol, fit: &mut Mol) -> Result<(f64, f64, f64)> {
    align_mol_with_params(ref_, fit, &AlignMolParams::default())
}

/// Mirrors `rdGaussianShape.AlignMol(ref, fit, refOpts, fitOpts,
/// overlayOpts, refConfId, fitConfId)`.
pub fn align_mol_with_params(
    ref_: &Mol,
    fit: &mut Mol,
    params: &AlignMolParams,
) -> Result<(f64, f64, f64)> {
    Ok(triple(&ffi::gaussian_align_mol(
        ref_.as_ffi(),
        fit.pin_mut(),
        &params.ref_opts.to_ffi(),
        &params.fit_opts.to_ffi(),
        &params.overlay_opts.to_ffi(),
        params.ref_conf_id,
        params.fit_conf_id,
    )?))
}

/// Parameters for `align_mol_shape_input_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolShapeInputParams {
    pub fit_opts: ShapeInputOptions,
    pub overlay_opts: ShapeOverlayOptions,
    /// The fit's conformer; -1 for all of them.
    pub fit_conf_id: i32,
}

impl Default for AlignMolShapeInputParams {
    fn default() -> AlignMolShapeInputParams {
        AlignMolShapeInputParams {
            fit_opts: ShapeInputOptions::default(),
            overlay_opts: ShapeOverlayOptions::default(),
            fit_conf_id: -1,
        }
    }
}

/// `align_mol` onto a reference shape's active shape. Mirrors
/// `rdGaussianShape.AlignMol(refShape, fit)`.
pub fn align_mol_shape_input(ref_shape: &ShapeInput, fit: &mut Mol) -> Result<(f64, f64, f64)> {
    align_mol_shape_input_with_params(ref_shape, fit, &AlignMolShapeInputParams::default())
}

/// Mirrors `rdGaussianShape.AlignMol(refShape, fit, fitOpts, overlayOpts,
/// fitConfId)`.
pub fn align_mol_shape_input_with_params(
    ref_shape: &ShapeInput,
    fit: &mut Mol,
    params: &AlignMolShapeInputParams,
) -> Result<(f64, f64, f64)> {
    Ok(triple(&ffi::gaussian_align_mol_shape_input(
        &ref_shape.ptr,
        fit.pin_mut(),
        &params.fit_opts.to_ffi(),
        &params.overlay_opts.to_ffi(),
        params.fit_conf_id,
    )?))
}

/// Overlays the fit's active shape on the reference's, moves it there, and
/// returns the scores (combination, shape, color) and the 16 values of the
/// transformation, by rows. Mirrors `rdGaussianShape.AlignShapes`.
pub fn align_shapes(
    ref_shape: &ShapeInput,
    fit_shape: &mut ShapeInput,
) -> Result<(f64, f64, f64, Vec<f64>)> {
    align_shapes_with_overlay_opts(ref_shape, fit_shape, &ShapeOverlayOptions::default())
}

/// Mirrors `rdGaussianShape.AlignShapes(refShape, fitShape, overlayOpts)`.
pub fn align_shapes_with_overlay_opts(
    ref_shape: &ShapeInput,
    fit_shape: &mut ShapeInput,
    overlay_opts: &ShapeOverlayOptions,
) -> Result<(f64, f64, f64, Vec<f64>)> {
    let result = ffi::gaussian_align_shapes(
        &ref_shape.ptr,
        fit_shape.ptr.pin_mut(),
        &overlay_opts.to_ffi(),
    )?;
    Ok((result[0], result[1], result[2], result[3..].to_vec()))
}

/// Parameters for `score_mol_with_params`, with Python's defaults (A2):
/// `align_mol_with_params`'s. The overlay's start mode and normalization
/// are ignored, since the shapes are scored where they are.
pub type ScoreMolParams = AlignMolParams;

/// The scores (combination, shape, color) and the shape and color overlap
/// volumes of two molecules where they are. Mirrors
/// `rdGaussianShape.ScoreMol`.
pub fn score_mol(ref_: &Mol, fit: &Mol) -> Result<(f64, f64, f64, f64, f64)> {
    score_mol_with_params(ref_, fit, &ScoreMolParams::default())
}

/// Mirrors `rdGaussianShape.ScoreMol(ref, fit, refOpts, fitOpts,
/// overlayOpts, refConfId, fitConfId)`.
pub fn score_mol_with_params(
    ref_: &Mol,
    fit: &Mol,
    params: &ScoreMolParams,
) -> Result<(f64, f64, f64, f64, f64)> {
    Ok(five(&ffi::gaussian_score_mol(
        ref_.as_ffi(),
        fit.as_ffi(),
        &params.ref_opts.to_ffi(),
        &params.fit_opts.to_ffi(),
        &params.overlay_opts.to_ffi(),
        params.ref_conf_id,
        params.fit_conf_id,
    )?))
}

/// Parameters for `score_mol_shape_input_with_params`, with Python's
/// defaults (A2): `align_mol_shape_input_with_params`'s.
pub type ScoreMolShapeInputParams = AlignMolShapeInputParams;

/// `score_mol` with a reference shape's active shape. Mirrors
/// `rdGaussianShape.ScoreMol(refShape, fit)`.
pub fn score_mol_shape_input(
    ref_shape: &ShapeInput,
    fit: &Mol,
) -> Result<(f64, f64, f64, f64, f64)> {
    score_mol_shape_input_with_params(ref_shape, fit, &ScoreMolShapeInputParams::default())
}

/// Mirrors `rdGaussianShape.ScoreMol(refShape, fit, fitOpts, overlayOpts,
/// fitConfId)`.
pub fn score_mol_shape_input_with_params(
    ref_shape: &ShapeInput,
    fit: &Mol,
    params: &ScoreMolShapeInputParams,
) -> Result<(f64, f64, f64, f64, f64)> {
    Ok(five(&ffi::gaussian_score_mol_shape_input(
        &ref_shape.ptr,
        fit.as_ffi(),
        &params.fit_opts.to_ffi(),
        &params.overlay_opts.to_ffi(),
        params.fit_conf_id,
    )?))
}

/// The scores and overlap volumes of two shapes' active shapes where they
/// are. Mirrors `rdGaussianShape.ScoreShape`.
pub fn score_shape(ref_shape: &ShapeInput, fit_shape: &ShapeInput) -> (f64, f64, f64, f64, f64) {
    score_shape_with_overlay_opts(ref_shape, fit_shape, &ShapeOverlayOptions::default())
}

/// Mirrors `rdGaussianShape.ScoreShape(refShape, fitShape, overlayOpts)`.
pub fn score_shape_with_overlay_opts(
    ref_shape: &ShapeInput,
    fit_shape: &ShapeInput,
    overlay_opts: &ShapeOverlayOptions,
) -> (f64, f64, f64, f64, f64) {
    five(&expect_ok(ffi::gaussian_score_shape(
        &ref_shape.ptr,
        &fit_shape.ptr,
        &overlay_opts.to_ffi(),
    )))
}

/// Parameters for `score_molecule_all_conformers_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScoreMoleculeAllConformersParams {
    pub ref_opts: ShapeInputOptions,
    pub fit_opts: ShapeInputOptions,
    pub overlay_opts: ShapeOverlayOptions,
}

/// The combination scores of overlaying each of the fit's conformers on
/// each of the reference's (`scores[i][j]` for reference conformer `i` and
/// fit conformer `j`, counted in order), the best pair, and its
/// transformation's 16 values, by rows. The molecules do not move. Mirrors
/// `rdGaussianShape.ScoreMoleculeAllConformers`.
pub fn score_molecule_all_conformers(ref_: &Mol, fit: &Mol) -> Result<ConformerScores> {
    score_molecule_all_conformers_with_params(
        ref_,
        fit,
        &ScoreMoleculeAllConformersParams::default(),
    )
}

/// Mirrors `rdGaussianShape.ScoreMoleculeAllConformers(ref, fit, refOpts,
/// fitOpts, overlayOpts)`.
pub fn score_molecule_all_conformers_with_params(
    ref_: &Mol,
    fit: &Mol,
    params: &ScoreMoleculeAllConformersParams,
) -> Result<ConformerScores> {
    let result = ffi::gaussian_score_molecule_all_conformers(
        ref_.as_ffi(),
        fit.as_ffi(),
        &params.ref_opts.to_ffi(),
        &params.fit_opts.to_ffi(),
        &params.overlay_opts.to_ffi(),
    )?;
    let (num_ref, num_fit) = (result[0] as usize, result[1] as usize);
    let scores = result[2..2 + num_ref * num_fit]
        .chunks(num_fit.max(1))
        .map(<[f64]>::to_vec)
        .collect();
    let rest = &result[2 + num_ref * num_fit..];
    Ok((scores, rest[0] as i32, rest[1] as i32, rest[2..].to_vec()))
}
