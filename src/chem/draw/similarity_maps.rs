// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit/Chem/Draw/SimilarityMaps.py`: similarity maps, which weigh
//! each atom by how much a similarity (or a model's prediction) drops
//! without it, and draw the weights as Gaussians under the molecule
//! (Riniker and Landrum, J. Cheminf. 5, 43 (2013)).
//!
//! Python's fingerprint and prediction functions are closures here. A
//! fingerprint function takes the molecule and the atom to leave out
//! (`None` for none, Python's -1). The default metric, Dice similarity,
//! needs bit vectors; other fingerprints take a metric of their own
//! (`dice_similarity_sparse_int_vect`, say).
//!
//! Python's matplotlib colormaps are not ported: a map's colours are given
//! as three colours. The atom-pair, topological-torsion, and Morgan
//! fingerprint functions are not ported either: they are built on
//! `rdMolDescriptors`' deprecated fingerprint functions, which are not
//! bound. `get_atomic_weights_for_fingerprint_generator` weighs atoms with
//! the fingerprint generators that replace them.

use crate::chem::Mol;
use crate::chem::draw::rd_mol_draw_2d::{
    ContourAndDrawGaussiansParams, ContourParams, DrawColour, MolDraw2D,
    PrepareMolForDrawingParams, contour_and_draw_gaussians_with_params,
    prepare_mol_for_drawing_with_params,
};
use crate::chem::rd_depictor::compute_2d_coords;
use crate::chem::rd_fingerprint_generator::{
    AdditionalOutput, FingerprintGenerator64, GetFingerprintParams,
};
use crate::chem::{RDKFingerprintParams, rdk_fingerprint_with_outputs};
use crate::data_structs::{BitVect, ExplicitBitVect, UIntSparseIntVect, dice_similarity};
use crate::geometry::{Point2D, Point3D};
use crate::{Error, Result};

/// Each atom's weight: how much `metric` between the reference's and the
/// probe's fingerprints drops without the atom, with Dice similarity.
/// Mirrors `SimilarityMaps.GetAtomicWeightsForFingerprint`.
///
/// `fp_function` is called on the reference, then on the probe, then on
/// the probe without each atom in turn.
pub fn get_atomic_weights_for_fingerprint<F: BitVect>(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
) -> Result<Vec<f64>> {
    get_atomic_weights_for_fingerprint_with_metric(ref_mol, probe_mol, fp_function, dice_similarity)
}

/// As `get_atomic_weights_for_fingerprint`, with a metric. Mirrors
/// `SimilarityMaps.GetAtomicWeightsForFingerprint(refMol, probeMol,
/// fpFunction, metric)`.
pub fn get_atomic_weights_for_fingerprint_with_metric<F>(
    ref_mol: &Mol,
    probe_mol: &Mol,
    mut fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
    metric: impl Fn(&F, &F) -> Result<f64>,
) -> Result<Vec<f64>> {
    let ref_fp = fp_function(ref_mol, None)?;
    let base_similarity = metric(&ref_fp, &fp_function(probe_mol, None)?)?;
    (0..probe_mol.get_num_atoms())
        .map(|atom_id| {
            Ok(base_similarity - metric(&ref_fp, &fp_function(probe_mol, Some(atom_id))?)?)
        })
        .collect()
}

/// Each atom's weight: how much a model's prediction for the probe drops
/// without it. Mirrors `SimilarityMaps.GetAtomicWeightsForModel`.
pub fn get_atomic_weights_for_model<F>(
    probe_mol: &Mol,
    mut fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
    mut prediction_function: impl FnMut(&F) -> Result<f64>,
) -> Result<Vec<f64>> {
    let base_proba = prediction_function(&fp_function(probe_mol, None)?)?;
    (0..probe_mol.get_num_atoms())
        .map(|atom_id| {
            Ok(base_proba - prediction_function(&fp_function(probe_mol, Some(atom_id))?)?)
        })
        .collect()
}

/// The weights scaled so the largest in magnitude is 1 (unless all are 0),
/// and that magnitude. Mirrors `SimilarityMaps.GetStandardizedWeights`;
/// no weights is `Error::Value`, as Python's `max` raises.
pub fn get_standardized_weights(weights: &[f64]) -> Result<(Vec<f64>, f64)> {
    let max_abs_weight = weights
        .iter()
        .map(|w| w.abs())
        .reduce(|a, b| if b > a { b } else { a })
        .ok_or_else(|| Error::Value {
            message: "max() iterable argument is empty".to_owned(),
        })?;
    if max_abs_weight > 0.0 {
        Ok((
            weights.iter().map(|w| w / max_abs_weight).collect(),
            max_abs_weight,
        ))
    } else {
        Ok((weights.to_vec(), max_abs_weight))
    }
}

/// Parameters for `get_similarity_map_from_weights_with_params`, with
/// Python's defaults (A2). Python's `scale`, `size`, `coordScale`, `step`,
/// `colors`, and `alpha`, which it ignores, are not ported.
#[derive(Clone, Debug, PartialEq)]
pub struct GetSimilarityMapFromWeightsParams {
    /// The colours of the lowest, middle, and highest weights; `None` for
    /// `ContourParams`' map.
    pub color_map: Option<[DrawColour; 3]>,
    /// The Gaussians' width; `None` for 0.3 times the first bond's length
    /// (or the distance between the first two atoms), to two decimals.
    pub sigma: Option<f64>,
    /// The number of contour lines.
    pub contour_lines: u32,
    /// The grid's resolution.
    pub grid_resolution: f64,
    /// The grid's padding around the atoms.
    pub extra_grid_padding: f64,
    /// Fills only the grid boxes whose values reach `fill_threshold`.
    pub use_fill_threshold: bool,
    /// The magnitude a grid box's value needs to be filled.
    pub fill_threshold: f64,
    /// Takes `fill_threshold` as a fraction of the values' range.
    pub fill_threshold_is_fraction: bool,
}

impl Default for GetSimilarityMapFromWeightsParams {
    fn default() -> GetSimilarityMapFromWeightsParams {
        GetSimilarityMapFromWeightsParams {
            color_map: None,
            sigma: None,
            contour_lines: 10,
            grid_resolution: 0.1,
            extra_grid_padding: 0.5,
            use_fill_threshold: false,
            fill_threshold: 0.01,
            fill_threshold_is_fraction: true,
        }
    }
}

/// Draws the molecule over Gaussians on its atoms, of heights `weights`,
/// on `draw2d`, which it clears first and leaves unfinished. Mirrors
/// `SimilarityMaps.GetSimilarityMapFromWeights(mol, weights, draw2d)`,
/// which returns `draw2d`.
pub fn get_similarity_map_from_weights(
    mol: &Mol,
    weights: &[f64],
    draw2d: &mut MolDraw2D,
) -> Result<()> {
    get_similarity_map_from_weights_with_params(
        mol,
        weights,
        draw2d,
        &GetSimilarityMapFromWeightsParams::default(),
    )
}

/// Mirrors `SimilarityMaps.GetSimilarityMapFromWeights(mol, weights,
/// draw2d, colorMap, sigma=..., contourLines=..., gridResolution=...,
/// extraGridPadding=..., useFillThreshold=..., fillThreshold=...,
/// fillThresholdIsFraction=...)`. As in Python, it turns off the drawer's
/// `clear_background`. `Error::Value` for fewer than two atoms, as in
/// Python, and as `contour_and_draw_gaussians_with_params` and
/// `MolDraw2D::draw_molecule` (a weight for each atom, say).
pub fn get_similarity_map_from_weights_with_params(
    mol: &Mol,
    weights: &[f64],
    draw2d: &mut MolDraw2D,
    params: &GetSimilarityMapFromWeightsParams,
) -> Result<()> {
    if mol.get_num_atoms() < 2 {
        return Err(Error::Value {
            message: "too few atoms".to_owned(),
        });
    }
    let mut mol = prepare_mol_for_drawing_with_params(
        mol,
        &PrepareMolForDrawingParams {
            add_chiral_hs: false,
            ..Default::default()
        },
    )?;
    if mol.get_num_conformers() == 0 {
        compute_2d_coords(&mut mol)?;
    }
    let (locs, sigma) = {
        let conf = mol.get_conformer()?;
        let sigma = match params.sigma {
            Some(sigma) => sigma,
            None => {
                let (idx1, idx2) = if mol.get_num_bonds() > 0 {
                    let bond = mol.get_bond_with_idx(0)?;
                    (bond.get_begin_atom_idx(), bond.get_end_atom_idx())
                } else {
                    (0, 1)
                };
                let p1 = conf.get_atom_position(idx1)?;
                let p2 = conf.get_atom_position(idx2)?;
                let d = Point3D::new(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
                round2(0.3 * d.length())
            }
        };
        let locs = (0..mol.get_num_atoms())
            .map(|i| Ok(Point2D::from(conf.get_atom_position(i)?)))
            .collect::<Result<Vec<_>>>()?;
        (locs, sigma)
    };
    let sigmas = vec![sigma; mol.get_num_atoms() as usize];
    draw2d.clear_drawing()?;
    let mut ps = ContourParams {
        fill_grid: true,
        grid_resolution: params.grid_resolution,
        extra_grid_padding: params.extra_grid_padding,
        use_fill_threshold: params.use_fill_threshold,
        fill_threshold: params.fill_threshold,
        fill_threshold_is_fraction: params.fill_threshold_is_fraction,
        ..Default::default()
    };
    if let Some(colors) = &params.color_map {
        ps.set_colour_map(colors);
    }
    contour_and_draw_gaussians_with_params(
        draw2d,
        &locs,
        weights,
        &sigmas,
        &ContourAndDrawGaussiansParams {
            n_contours: params.contour_lines,
            params: ps,
            ..Default::default()
        },
    )?;
    let mut options = draw2d.draw_options();
    options.clear_background = false;
    draw2d.set_draw_options(&options)?;
    draw2d.draw_molecule(&mol)
}

/// Python's `round(x, 2)`: the double nearest the decimal that rounds `x`'s
/// exact value to two places, ties to even, as Rust's formatting rounds.
fn round2(x: f64) -> f64 {
    format!("{x:.2}").parse().expect("a formatted float parses")
}

/// Parameters for `get_similarity_map_for_fingerprint_with_params`, with
/// Python's defaults (A2): the metric, and the drawing's parameters.
#[derive(Clone, Debug)]
pub struct GetSimilarityMapForFingerprintParams<F> {
    /// The similarity metric.
    pub metric: fn(&F, &F) -> Result<f64>,
    /// The parameters of `get_similarity_map_from_weights_with_params`.
    pub map: GetSimilarityMapFromWeightsParams,
}

impl<F: BitVect> Default for GetSimilarityMapForFingerprintParams<F> {
    fn default() -> GetSimilarityMapForFingerprintParams<F> {
        GetSimilarityMapForFingerprintParams {
            metric: dice_similarity::<F>,
            map: GetSimilarityMapFromWeightsParams::default(),
        }
    }
}

/// Draws the probe's similarity map against the reference
/// (`get_atomic_weights_for_fingerprint`, standardized) on `draw2d`, and
/// returns the largest weight's magnitude. Mirrors
/// `SimilarityMaps.GetSimilarityMapForFingerprint`, which returns `draw2d`
/// with it.
pub fn get_similarity_map_for_fingerprint<F: BitVect>(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
    draw2d: &mut MolDraw2D,
) -> Result<f64> {
    get_similarity_map_for_fingerprint_with_params(
        ref_mol,
        probe_mol,
        fp_function,
        draw2d,
        &GetSimilarityMapForFingerprintParams::default(),
    )
}

/// Mirrors `SimilarityMaps.GetSimilarityMapForFingerprint(refMol,
/// probeMol, fpFunction, draw2d, metric, **kwargs)`, whose keywords are
/// `GetSimilarityMapFromWeights`'.
pub fn get_similarity_map_for_fingerprint_with_params<F>(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
    draw2d: &mut MolDraw2D,
    params: &GetSimilarityMapForFingerprintParams<F>,
) -> Result<f64> {
    let weights = get_atomic_weights_for_fingerprint_with_metric(
        ref_mol,
        probe_mol,
        fp_function,
        params.metric,
    )?;
    draw_standardized(probe_mol, &weights, draw2d, &params.map)
}

fn draw_standardized(
    probe_mol: &Mol,
    weights: &[f64],
    draw2d: &mut MolDraw2D,
    params: &GetSimilarityMapFromWeightsParams,
) -> Result<f64> {
    let (weights, max_weight) = get_standardized_weights(weights)?;
    get_similarity_map_from_weights_with_params(probe_mol, &weights, draw2d, params)?;
    Ok(max_weight)
}

/// The atoms that set each bit, as `AdditionalOutput::get_atoms_per_bit`
/// gives them.
type AtomsPerBit = std::collections::BTreeMap<u64, Vec<Vec<u32>>>;

/// The probe's fingerprint and the atoms that set each of its bits.
fn probe_fingerprint<T>(
    probe_mol: &Mol,
    get: impl FnOnce(&Mol, &GetFingerprintParams) -> Result<(T, AdditionalOutput)>,
) -> Result<(T, AtomsPerBit)> {
    let mut probe_ao = AdditionalOutput::default();
    probe_ao.allocate_atoms_per_bit();
    let (fp, outputs) = get(
        probe_mol,
        &GetFingerprintParams {
            additional_output: Some(probe_ao),
            ..Default::default()
        },
    )?;
    let atoms_per_bit = outputs.get_atoms_per_bit().cloned().unwrap_or_default();
    Ok((fp, atoms_per_bit))
}

fn bit_index(bit: u64) -> Result<u32> {
    u32::try_from(bit).map_err(|_| Error::Value {
        message: format!("bit {bit} is past the fingerprint's end"),
    })
}

/// Each atom's weight: how much Dice similarity between the reference's
/// and the probe's fingerprints from `fpg` drops without the bits the atom
/// sets. Mirrors
/// `SimilarityMaps.GetAtomicWeightsForFingerprintGenerator`.
pub fn get_atomic_weights_for_fingerprint_generator(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
) -> Result<Vec<f64>> {
    get_atomic_weights_for_fingerprint_generator_with_metric(
        ref_mol,
        probe_mol,
        fpg,
        dice_similarity::<ExplicitBitVect>,
    )
}

/// As `get_atomic_weights_for_fingerprint_generator`, with a metric.
/// Mirrors `SimilarityMaps.GetAtomicWeightsForFingerprintGenerator(refMol,
/// probeMol, fpg, metric=...)`.
pub fn get_atomic_weights_for_fingerprint_generator_with_metric(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
    metric: impl Fn(&ExplicitBitVect, &ExplicitBitVect) -> Result<f64>,
) -> Result<Vec<f64>> {
    let ref_fp = fpg.get_fingerprint(ref_mol)?;
    let (probe_fp, atoms_per_bit) =
        probe_fingerprint(probe_mol, |m, p| fpg.get_fingerprint_with_outputs(m, p))?;
    let base_similarity = metric(&ref_fp, &probe_fp)?;
    let mut atom_fps = vec![probe_fp; probe_mol.get_num_atoms() as usize];
    for (&bit, entries) in &atoms_per_bit {
        for entry in entries {
            for &aid in entry {
                atom_fps[aid as usize].un_set_bit(bit_index(bit)?)?;
            }
        }
    }
    atom_fps
        .iter()
        .map(|fp| Ok(base_similarity - metric(&ref_fp, fp)?))
        .collect()
}

/// As `get_atomic_weights_for_fingerprint_generator`, with count
/// fingerprints, each atom taking one from each count it adds to. Mirrors
/// `SimilarityMaps.GetAtomicWeightsForFingerprintGenerator(refMol,
/// probeMol, fpg, useCounts=True)`.
pub fn get_atomic_weights_for_fingerprint_generator_use_counts(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
) -> Result<Vec<f64>> {
    get_atomic_weights_for_fingerprint_generator_use_counts_with_metric(
        ref_mol,
        probe_mol,
        fpg,
        crate::data_structs::dice_similarity_sparse_int_vect::<UIntSparseIntVect>,
    )
}

/// As `get_atomic_weights_for_fingerprint_generator_use_counts`, with a
/// metric. Mirrors
/// `SimilarityMaps.GetAtomicWeightsForFingerprintGenerator(refMol,
/// probeMol, fpg, useCounts=True, metric=...)`.
pub fn get_atomic_weights_for_fingerprint_generator_use_counts_with_metric(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
    metric: impl Fn(&UIntSparseIntVect, &UIntSparseIntVect) -> Result<f64>,
) -> Result<Vec<f64>> {
    let ref_fp = fpg.get_count_fingerprint(ref_mol)?;
    let (probe_fp, atoms_per_bit) = probe_fingerprint(probe_mol, |m, p| {
        fpg.get_count_fingerprint_with_outputs(m, p)
    })?;
    let base_similarity = metric(&ref_fp, &probe_fp)?;
    let mut atom_fps = vec![probe_fp; probe_mol.get_num_atoms() as usize];
    for (&bit, entries) in &atoms_per_bit {
        let bit = bit_index(bit)?;
        for entry in entries {
            for &aid in entry {
                let fp = &mut atom_fps[aid as usize];
                let val = fp.get_val(bit)?;
                fp.set_val(bit, val - 1)?;
            }
        }
    }
    atom_fps
        .iter()
        .map(|fp| Ok(base_similarity - metric(&ref_fp, fp)?))
        .collect()
}

/// Parameters for `get_similarity_map_for_fingerprint_generator_with_params`
/// and its `use_counts` form, with Python's defaults (A2): the metric, and
/// the drawing's parameters.
#[derive(Clone, Debug)]
pub struct GetSimilarityMapForFingerprintGeneratorParams<F> {
    /// The similarity metric.
    pub metric: fn(&F, &F) -> Result<f64>,
    /// The parameters of `get_similarity_map_from_weights_with_params`.
    pub map: GetSimilarityMapFromWeightsParams,
}

impl Default for GetSimilarityMapForFingerprintGeneratorParams<ExplicitBitVect> {
    fn default() -> Self {
        GetSimilarityMapForFingerprintGeneratorParams {
            metric: dice_similarity::<ExplicitBitVect>,
            map: GetSimilarityMapFromWeightsParams::default(),
        }
    }
}

impl Default for GetSimilarityMapForFingerprintGeneratorParams<UIntSparseIntVect> {
    fn default() -> Self {
        GetSimilarityMapForFingerprintGeneratorParams {
            metric: crate::data_structs::dice_similarity_sparse_int_vect::<UIntSparseIntVect>,
            map: GetSimilarityMapFromWeightsParams::default(),
        }
    }
}

/// Draws the probe's similarity map against the reference
/// (`get_atomic_weights_for_fingerprint_generator`, standardized) on
/// `draw2d`, and returns the largest weight's magnitude. Mirrors
/// `SimilarityMaps.GetSimilarityMapForFingerprintGenerator`, which returns
/// `draw2d` with it.
pub fn get_similarity_map_for_fingerprint_generator(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
    draw2d: &mut MolDraw2D,
) -> Result<f64> {
    get_similarity_map_for_fingerprint_generator_with_params(
        ref_mol,
        probe_mol,
        fpg,
        draw2d,
        &GetSimilarityMapForFingerprintGeneratorParams::default(),
    )
}

/// Mirrors `SimilarityMaps.GetSimilarityMapForFingerprintGenerator(refMol,
/// probeMol, fpg, draw2d, metric, **kwargs)`, whose keywords are
/// `GetSimilarityMapFromWeights`'.
pub fn get_similarity_map_for_fingerprint_generator_with_params(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
    draw2d: &mut MolDraw2D,
    params: &GetSimilarityMapForFingerprintGeneratorParams<ExplicitBitVect>,
) -> Result<f64> {
    let weights = get_atomic_weights_for_fingerprint_generator_with_metric(
        ref_mol,
        probe_mol,
        fpg,
        params.metric,
    )?;
    draw_standardized(probe_mol, &weights, draw2d, &params.map)
}

/// As `get_similarity_map_for_fingerprint_generator`, with count
/// fingerprints. Mirrors
/// `SimilarityMaps.GetSimilarityMapForFingerprintGenerator(refMol,
/// probeMol, fpg, draw2d, useCounts=True)`.
pub fn get_similarity_map_for_fingerprint_generator_use_counts(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
    draw2d: &mut MolDraw2D,
) -> Result<f64> {
    get_similarity_map_for_fingerprint_generator_use_counts_with_params(
        ref_mol,
        probe_mol,
        fpg,
        draw2d,
        &GetSimilarityMapForFingerprintGeneratorParams::default(),
    )
}

/// Mirrors `SimilarityMaps.GetSimilarityMapForFingerprintGenerator(refMol,
/// probeMol, fpg, draw2d, metric, useCounts=True, **kwargs)`.
pub fn get_similarity_map_for_fingerprint_generator_use_counts_with_params(
    ref_mol: &Mol,
    probe_mol: &Mol,
    fpg: &FingerprintGenerator64,
    draw2d: &mut MolDraw2D,
    params: &GetSimilarityMapForFingerprintGeneratorParams<UIntSparseIntVect>,
) -> Result<f64> {
    let weights = get_atomic_weights_for_fingerprint_generator_use_counts_with_metric(
        ref_mol,
        probe_mol,
        fpg,
        params.metric,
    )?;
    draw_standardized(probe_mol, &weights, draw2d, &params.map)
}

/// Draws the probe's map of a model's predictions
/// (`get_atomic_weights_for_model`, standardized) on `draw2d`, and returns
/// the largest weight's magnitude. Mirrors
/// `SimilarityMaps.GetSimilarityMapForModel`, which returns `draw2d` with
/// it.
pub fn get_similarity_map_for_model<F>(
    probe_mol: &Mol,
    fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
    prediction_function: impl FnMut(&F) -> Result<f64>,
    draw2d: &mut MolDraw2D,
) -> Result<f64> {
    get_similarity_map_for_model_with_params(
        probe_mol,
        fp_function,
        prediction_function,
        draw2d,
        &GetSimilarityMapFromWeightsParams::default(),
    )
}

/// Mirrors `SimilarityMaps.GetSimilarityMapForModel(probeMol, fpFunction,
/// predictionFunction, draw2d, **kwargs)`, whose keywords are
/// `GetSimilarityMapFromWeights`'.
pub fn get_similarity_map_for_model_with_params<F>(
    probe_mol: &Mol,
    fp_function: impl FnMut(&Mol, Option<u32>) -> Result<F>,
    prediction_function: impl FnMut(&F) -> Result<f64>,
    draw2d: &mut MolDraw2D,
    params: &GetSimilarityMapFromWeightsParams,
) -> Result<f64> {
    let weights = get_atomic_weights_for_model(probe_mol, fp_function, prediction_function)?;
    draw_standardized(probe_mol, &weights, draw2d, params)
}

/// Parameters for `get_rdk_fingerprint_with_params`, with Python's
/// defaults (A2). Python's `fpType` takes only `'bv'`, and is not ported.
#[derive(Clone, Debug, PartialEq)]
pub struct GetRDKFingerprintParams {
    /// The atom whose paths to leave out; `None` for none (Python's -1).
    pub atom_id: Option<u32>,
    /// The number of bits.
    pub n_bits: u32,
    /// The fewest bonds in a path.
    pub min_path: u32,
    /// The most bonds in a path.
    pub max_path: u32,
    /// The bits each path sets.
    pub n_bits_per_hash: u32,
}

impl Default for GetRDKFingerprintParams {
    fn default() -> GetRDKFingerprintParams {
        GetRDKFingerprintParams {
            atom_id: None,
            n_bits: 2048,
            min_path: 1,
            max_path: 5,
            n_bits_per_hash: 2,
        }
    }
}

/// The RDKit fingerprint (`rdk_fingerprint`). Mirrors
/// `SimilarityMaps.GetRDKFingerprint(mol)`.
pub fn get_rdk_fingerprint(mol: &Mol) -> Result<ExplicitBitVect> {
    get_rdk_fingerprint_with_params(mol, &GetRDKFingerprintParams::default())
}

/// The RDKit fingerprint without the bits of `atom_id`'s paths. Mirrors
/// `SimilarityMaps.GetRDKFingerprint(mol, atomId, fpType, nBits, minPath,
/// maxPath, nBitsPerHash)`; `Error::Value` for an atom past the last, as
/// in Python. Python keeps the fingerprint on the molecule for later calls
/// (whatever their parameters, until `GetAtomicWeightsForFingerprint`
/// deletes it); this computes it each time.
pub fn get_rdk_fingerprint_with_params(
    mol: &Mol,
    params: &GetRDKFingerprintParams,
) -> Result<ExplicitBitVect> {
    let (mut fp, outputs) = rdk_fingerprint_with_outputs(
        mol,
        &RDKFingerprintParams {
            fp_size: params.n_bits,
            min_path: params.min_path,
            max_path: params.max_path,
            n_bits_per_hash: params.n_bits_per_hash,
            ..Default::default()
        },
    )?;
    let Some(atom_id) = params.atom_id else {
        return Ok(fp);
    };
    if atom_id >= mol.get_num_atoms() {
        return Err(Error::Value {
            message: "atom index greater than number of atoms".to_owned(),
        });
    }
    let bits = outputs.atom_bits[atom_id as usize]
        .iter()
        .map(|&bit| bit_index(bit))
        .collect::<Result<Vec<u32>>>()?;
    fp.un_set_bits_from_list(&bits)?;
    Ok(fp)
}
