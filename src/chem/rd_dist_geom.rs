// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdDistGeom`: conformers from distance geometry.
//!
//! The embedding runs on the calling thread whatever `num_threads` asks for,
//! and one embedding at a time runs in the process: RDKit's own threads race
//! on state they share, and every call swaps the process's SIGINT handler
//! through static variables. Each conformer's random seed depends only on
//! its index, so a seeded embedding gives the conformers RDKit's threads
//! would.

use std::collections::BTreeMap;

use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::{expect_ok, to_string};
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::embedder::ffi;
use crate::geometry::Point3D;
use crate::{Error, Result};

rdkit_enum! {
    /// The checks an embedding attempt can fail. Mirrors
    /// `rdDistGeom.EmbedFailureCauses`; `EmbedParameters::get_failure_counts`
    /// counts them in this order.
    EmbedFailureCauses: ffi::EmbedFailureCauses {
        InitialCoords = INITIAL_COORDS,
        FirstMinimization = FIRST_MINIMIZATION,
        CheckTetrahedralCenters = CHECK_TETRAHEDRAL_CENTERS,
        CheckChiralCenters = CHECK_CHIRAL_CENTERS,
        MinimizeFourthDimension = MINIMIZE_FOURTH_DIMENSION,
        EtkMinimization = ETK_MINIMIZATION,
        FinalChiralBounds = FINAL_CHIRAL_BOUNDS,
        FinalCenterInVolume = FINAL_CENTER_IN_VOLUME,
        LinearDoubleBond = LINEAR_DOUBLE_BOND,
        BadDoubleBondStereo = BAD_DOUBLE_BOND_STEREO,
        CheckChiralCenters2 = CHECK_CHIRAL_CENTERS2,
        ExceededTimeout = EXCEEDED_TIMEOUT,
        Minimization = MINIMIZATION,
        KtermViolation = KTERM_VIOLATION,
        Clash = CLASH,
    }
}

/// Parameters for embedding. Mirrors `rdDistGeom.EmbedParameters`;
/// `Default` reads a default-constructed C++ object, and [`etkdg_v3`] and
/// its siblings return RDKit's presets.
///
/// Embedding updates the parameters as Python's do: the failure counts, and
/// what RDKit changes while it works (a `max_iterations` of 0 becomes ten
/// times the number of atoms, random coordinates set `basin_thresh`).
#[derive(Clone, Debug, PartialEq)]
pub struct EmbedParameters {
    /// The most attempts per conformer; 0 for ten times the atoms.
    pub max_iterations: u32,
    /// Accepted for parity: the embedding always runs on the calling thread
    /// (see the module documentation).
    pub num_threads: i32,
    /// Seconds per fragment before giving up; 0 for no limit.
    pub timeout: u32,
    /// The random seed; −1 leaves RDKit's generator unseeded.
    pub random_seed: i32,
    /// Removes the molecule's conformers first.
    pub clear_confs: bool,
    /// Starts from random coordinates rather than distance geometry.
    pub use_random_coords: bool,
    /// The size of the box for random coordinates.
    pub box_size_mult: f64,
    /// Picks random coordinates for a negative eigenvalue.
    pub rand_neg_eig: bool,
    /// Fails with this many zero eigenvalues.
    pub num_zero_fail: u32,
    /// The distance-geometry minimizer's force tolerance.
    pub optimizer_force_tol: f64,
    /// The distance-geometry force field's basin threshold.
    pub basin_thresh: f64,
    /// Embeds even when triangle smoothing fails.
    pub ignore_smoothing_failures: bool,
    /// Enforces the chiral centers.
    pub enforce_chirality: bool,
    /// Uses experimental torsion-angle preferences.
    pub use_exp_torsion_angle_prefs: bool,
    /// Uses basic knowledge such as flat aromatic rings.
    pub use_basic_knowledge: bool,
    /// The torsion preferences' version: 1 or 2.
    pub et_version: u32,
    /// Prints the torsion preferences.
    pub verbose: bool,
    /// Keeps only conformers this far apart, by RMSD; −1 keeps all.
    pub prune_rms_thresh: f64,
    /// Prunes on heavy atoms only.
    pub only_heavy_atoms_for_rms: bool,
    /// Embeds each fragment separately.
    pub embed_fragments_separately: bool,
    /// Uses small-ring torsion preferences.
    pub use_small_ring_torsions: bool,
    /// Uses macrocycle torsion preferences.
    pub use_macrocycle_torsions: bool,
    /// Uses the 1-4 distance bounds heuristics for macrocycles.
    pub use_macrocycle14config: bool,
    /// Minimizes in stages rather than all at once.
    pub use_legacy_implementation: bool,
    /// Scales the bounds-matrix restraints against the others.
    pub bounds_mat_force_scaling: f64,
    /// Uses the molecule's symmetry when pruning.
    pub use_symmetry_for_pruning: bool,
    /// Constrains amide bonds to be trans.
    pub force_trans_amides: bool,
    /// Counts the checks that fail, for `get_failure_counts`.
    pub track_failures: bool,
    /// Seeds conformers so that embedding can be restarted.
    pub enable_sequential_random_seeds: bool,
    /// Treats conjugated terminal groups symmetrically when pruning.
    pub symmetrize_conjugated_terminal_groups_for_pruning: bool,
    /// Atoms to place relative to each other at these coordinates; see
    /// `set_coord_map`.
    pub coord_map: BTreeMap<u32, Point3D>,
    /// Coulomb-like interactions between atom pairs; see `set_cpci`.
    pub cpci: BTreeMap<(u32, u32), f64>,
    /// A bounds matrix to use instead of the molecule's, as rows; see
    /// `set_bounds_mat`.
    pub bounds_mat: Option<Vec<Vec<f64>>>,
    /// The failure counts of the last embedding with `track_failures`, by
    /// `EmbedFailureCauses`; see `get_failure_counts`.
    pub failures: Vec<u32>,
}

impl Default for EmbedParameters {
    fn default() -> EmbedParameters {
        EmbedParameters::preset(ffi::EmbedPreset::Default)
    }
}

impl EmbedParameters {
    fn preset(preset: ffi::EmbedPreset) -> EmbedParameters {
        let data = expect_ok(ffi::embed_parameters_preset(preset));
        let mut params = EmbedParameters {
            max_iterations: 0,
            num_threads: 1,
            timeout: 0,
            random_seed: -1,
            clear_confs: true,
            use_random_coords: false,
            box_size_mult: 2.0,
            rand_neg_eig: true,
            num_zero_fail: 1,
            optimizer_force_tol: 1e-3,
            basin_thresh: 5.0,
            ignore_smoothing_failures: false,
            enforce_chirality: true,
            use_exp_torsion_angle_prefs: false,
            use_basic_knowledge: false,
            et_version: 2,
            verbose: false,
            prune_rms_thresh: -1.0,
            only_heavy_atoms_for_rms: true,
            embed_fragments_separately: true,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_macrocycle14config: false,
            use_legacy_implementation: true,
            bounds_mat_force_scaling: 1.0,
            use_symmetry_for_pruning: true,
            force_trans_amides: true,
            track_failures: false,
            enable_sequential_random_seeds: false,
            symmetrize_conjugated_terminal_groups_for_pruning: true,
            coord_map: BTreeMap::new(),
            cpci: BTreeMap::new(),
            bounds_mat: None,
            failures: Vec::new(),
        };
        params.update(&data);
        params
    }

    /// Takes the fields RDKit reads and writes.
    fn update(&mut self, data: &ffi::EmbedParametersData) {
        self.max_iterations = data.max_iterations;
        self.num_threads = data.num_threads;
        self.timeout = data.timeout;
        self.random_seed = data.random_seed;
        self.clear_confs = data.clear_confs;
        self.use_random_coords = data.use_random_coords;
        self.box_size_mult = data.box_size_mult;
        self.rand_neg_eig = data.rand_neg_eig;
        self.num_zero_fail = data.num_zero_fail;
        self.optimizer_force_tol = data.optimizer_force_tol;
        self.basin_thresh = data.basin_thresh;
        self.ignore_smoothing_failures = data.ignore_smoothing_failures;
        self.enforce_chirality = data.enforce_chirality;
        self.use_exp_torsion_angle_prefs = data.use_exp_torsion_angle_prefs;
        self.use_basic_knowledge = data.use_basic_knowledge;
        self.et_version = data.et_version;
        self.verbose = data.verbose;
        self.prune_rms_thresh = data.prune_rms_thresh;
        self.only_heavy_atoms_for_rms = data.only_heavy_atoms_for_rms;
        self.embed_fragments_separately = data.embed_fragments_separately;
        self.use_small_ring_torsions = data.use_small_ring_torsions;
        self.use_macrocycle_torsions = data.use_macrocycle_torsions;
        self.use_macrocycle14config = data.use_macrocycle14config;
        self.use_legacy_implementation = data.use_legacy_implementation;
        self.bounds_mat_force_scaling = data.bounds_mat_force_scaling;
        self.use_symmetry_for_pruning = data.use_symmetry_for_pruning;
        self.force_trans_amides = data.force_trans_amides;
        self.track_failures = data.track_failures;
        self.enable_sequential_random_seeds = data.enable_sequential_random_seeds;
        self.symmetrize_conjugated_terminal_groups_for_pruning =
            data.symmetrize_conjugated_terminal_groups_for_pruning;
        self.failures = data.failures.clone();
    }

    /// A bounds matrix that is not square or is empty is `Error::Value`.
    fn to_ffi(&self) -> Result<ffi::EmbedParametersData> {
        let (has_bounds_mat, bounds_mat_size, bounds_mat) = match &self.bounds_mat {
            Some(rows) => (true, square_size(rows)?, rows.concat()),
            None => (false, 0, Vec::new()),
        };
        Ok(ffi::EmbedParametersData {
            max_iterations: self.max_iterations,
            num_threads: self.num_threads,
            timeout: self.timeout,
            random_seed: self.random_seed,
            clear_confs: self.clear_confs,
            use_random_coords: self.use_random_coords,
            box_size_mult: self.box_size_mult,
            rand_neg_eig: self.rand_neg_eig,
            num_zero_fail: self.num_zero_fail,
            optimizer_force_tol: self.optimizer_force_tol,
            basin_thresh: self.basin_thresh,
            ignore_smoothing_failures: self.ignore_smoothing_failures,
            enforce_chirality: self.enforce_chirality,
            use_exp_torsion_angle_prefs: self.use_exp_torsion_angle_prefs,
            use_basic_knowledge: self.use_basic_knowledge,
            et_version: self.et_version,
            verbose: self.verbose,
            prune_rms_thresh: self.prune_rms_thresh,
            only_heavy_atoms_for_rms: self.only_heavy_atoms_for_rms,
            embed_fragments_separately: self.embed_fragments_separately,
            use_small_ring_torsions: self.use_small_ring_torsions,
            use_macrocycle_torsions: self.use_macrocycle_torsions,
            use_macrocycle14config: self.use_macrocycle14config,
            use_legacy_implementation: self.use_legacy_implementation,
            bounds_mat_force_scaling: self.bounds_mat_force_scaling,
            use_symmetry_for_pruning: self.use_symmetry_for_pruning,
            force_trans_amides: self.force_trans_amides,
            track_failures: self.track_failures,
            enable_sequential_random_seeds: self.enable_sequential_random_seeds,
            symmetrize_conjugated_terminal_groups_for_pruning: self
                .symmetrize_conjugated_terminal_groups_for_pruning,
            coord_map_atoms: self.coord_map.keys().copied().collect(),
            coord_map_points: self.coord_map.values().copied().collect(),
            cpci_pairs: self
                .cpci
                .keys()
                .map(|&(first, second)| IndexPair { first, second })
                .collect(),
            cpci_values: self.cpci.values().copied().collect(),
            has_bounds_mat,
            bounds_mat_size,
            bounds_mat,
            failures: self.failures.clone(),
        })
    }

    /// Sets the bounds matrix to use, as rows, without triangle smoothing:
    /// lower bounds below the diagonal, upper bounds above it. Mirrors
    /// `EmbedParameters.SetBoundsMat`; a matrix that is not square, or is
    /// empty, is `Error::Value`.
    pub fn set_bounds_mat(&mut self, bounds_mat_arg: &[Vec<f64>]) -> Result<()> {
        square_size(bounds_mat_arg)?;
        self.bounds_mat = Some(bounds_mat_arg.to_vec());
        Ok(())
    }

    /// Sets Coulomb-like interactions between pairs of atoms, for the
    /// minimization. Mirrors `EmbedParameters.SetCPCI`.
    pub fn set_cpci(&mut self, cpci_dict: &BTreeMap<(u32, u32), f64>) {
        self.cpci = cpci_dict.clone();
    }

    /// Sets coordinates whose distances the conformers reproduce for their
    /// atoms. Mirrors `EmbedParameters.SetCoordMap`. RDKit ignores the map
    /// for molecules with more than one fragment, except that random
    /// coordinates apply its indices to each fragment.
    pub fn set_coord_map(&mut self, coord_map: &BTreeMap<u32, Point3D>) {
        self.coord_map = coord_map.clone();
    }

    /// How often each check failed in the last embedding with
    /// `track_failures`, by `EmbedFailureCauses`. Mirrors
    /// `EmbedParameters.GetFailureCounts`.
    pub fn get_failure_counts(&self) -> Vec<u32> {
        self.failures.clone()
    }
}

/// The size of a square, non-empty matrix given as rows, as
/// `setBoundsMatrix` requires.
fn square_size(rows: &[Vec<f64>]) -> Result<u32> {
    if rows.is_empty() {
        return Err(Error::Value {
            message: "The array has to have a nonzero size".to_owned(),
        });
    }
    if rows.iter().any(|row| row.len() != rows.len()) {
        return Err(Error::Value {
            message: "The array has to be square".to_owned(),
        });
    }
    u32::try_from(rows.len()).map_err(|_| Error::Value {
        message: "The array is too large".to_owned(),
    })
}

/// Splits a row-major matrix into rows of `n`.
fn rows(values: &[f64], n: usize) -> Vec<Vec<f64>> {
    if n == 0 {
        return Vec::new();
    }
    values.chunks(n).map(<[f64]>::to_vec).collect()
}

macro_rules! presets {
    ($($name:ident = $preset:ident, $python:literal, $what:literal;)*) => {
        $(
            #[doc = concat!("Parameters for ", $what, ". Mirrors `rdDistGeom.", $python, "`.")]
            pub fn $name() -> EmbedParameters {
                EmbedParameters::preset(ffi::EmbedPreset::$preset)
            }
        )*
    };
}

presets! {
    etkdg = ETKDG, "ETKDG", "ETKDG, version 1";
    etkdg_v2 = ETKDGv2, "ETKDGv2", "ETKDG, version 2";
    etkdg_v3 = ETKDGv3, "ETKDGv3", "ETKDG, version 3, with macrocycle torsions";
    sr_etkdg_v3 = SrETKDGv3, "srETKDGv3", "ETKDG, version 3, with small-ring torsions";
    etdg = ETDG, "ETDG", "ETDG";
    etdg_v2 = ETDGv2, "ETDGv2", "ETDG, version 2";
    kdg = KDG, "KDG", "KDG";
    dg = DG, "DG", "plain distance geometry";
}

/// Python's keyword defaults for `EmbedMolecule`, which differ from
/// `EmbedParameters`' own.
fn keyword_defaults() -> EmbedParameters {
    EmbedParameters {
        use_exp_torsion_angle_prefs: true,
        use_basic_knowledge: true,
        use_macrocycle_torsions: true,
        use_macrocycle14config: true,
        ..EmbedParameters::default()
    }
}

/// The parameters as RDKit's JSON. Mirrors
/// `rdDistGeom.EmbedParametersToJSON`.
pub fn embed_parameters_to_json(embed_parameters: &EmbedParameters) -> Result<String> {
    to_string(ffi::embed_parameters_to_json(&embed_parameters.to_ffi()?)?)
}

/// Embeds a conformer with Python's keyword defaults, returning its id or −1
/// on failure. Mirrors `rdDistGeom.EmbedMolecule(mol)`, whose defaults are
/// `EmbedParameters::default()` with `use_exp_torsion_angle_prefs`,
/// `use_basic_knowledge`, `use_macrocycle_torsions`, and
/// `use_macrocycle14config` set. A molecule without atoms is `Error::Value`.
pub fn embed_molecule(mol: &mut Mol) -> Result<i32> {
    embed_molecule_with_params(mol, &mut keyword_defaults())
}

/// Mirrors `rdDistGeom.EmbedMolecule(mol, params)`, and its keyword form,
/// whose arguments are the fields of `params` (`maxAttempts` is
/// `max_iterations`, `forceTol` `optimizer_force_tol`, and
/// `printExpTorsionAngles` `verbose`).
///
/// `Error::Value` for a molecule without atoms, above 16383 atoms in a
/// fragment (where the minimizer overflows), for distances whose squares
/// RDKit would sum past `i32::MAX` (with distance geometry, about 4295 atoms
/// in a fragment), and for coordinates or bounds that are not finite;
/// `Error::Index` for an atom in the coordinate map or CPCI past the end.
pub fn embed_molecule_with_params(mol: &mut Mol, params: &mut EmbedParameters) -> Result<i32> {
    let ids = embed_multiple_confs_with_params(mol, 1, params)?;
    Ok(ids.first().copied().unwrap_or(-1))
}

/// Embeds ten conformers with Python's keyword defaults, returning their
/// ids. Mirrors `rdDistGeom.EmbedMultipleConfs(mol)`, whose defaults are
/// those of `embed_molecule` with `only_heavy_atoms_for_rms` unset.
pub fn embed_multiple_confs(mol: &mut Mol) -> Result<Vec<i32>> {
    let mut params = EmbedParameters {
        only_heavy_atoms_for_rms: false,
        ..keyword_defaults()
    };
    embed_multiple_confs_with_params(mol, 10, &mut params)
}

/// Mirrors `rdDistGeom.EmbedMultipleConfs(mol, numConfs, params)`, and its
/// keyword form, as `embed_molecule_with_params`; a timeout is a list
/// ending in −1.
pub fn embed_multiple_confs_with_params(
    mol: &mut Mol,
    num_confs: u32,
    params: &mut EmbedParameters,
) -> Result<Vec<i32>> {
    let mut data = params.to_ffi()?;
    let ids = ffi::embed_multiple_confs(mol.pin_mut(), num_confs, &mut data)?;
    params.update(&data);
    Ok(ids)
}

/// Parameters for `get_molecule_bounds_matrix_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetMoleculeBoundsMatrixParams {
    /// Sets 1-5 bounds from the topology.
    pub set15bounds: bool,
    /// Scales the van der Waals lower bounds down.
    pub scale_vdw: bool,
    /// Triangle-smooths the matrix.
    pub do_triangle_smoothing: bool,
    /// Uses macrocycle heuristics for 1-4 bounds.
    pub use_macrocycle14config: bool,
    /// Keeps amides trans.
    pub force_trans_amides: bool,
    /// Sets 1-4 bounds from the topology.
    pub set14bounds: bool,
    /// Sets 1-3 bounds from the topology.
    pub set13bounds: bool,
}

impl Default for GetMoleculeBoundsMatrixParams {
    fn default() -> GetMoleculeBoundsMatrixParams {
        GetMoleculeBoundsMatrixParams {
            set15bounds: true,
            scale_vdw: false,
            do_triangle_smoothing: true,
            use_macrocycle14config: false,
            force_trans_amides: true,
            set14bounds: true,
            set13bounds: true,
        }
    }
}

/// Parameters for `get_molecule_bounds_matrix_embed_parameters_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetMoleculeBoundsMatrixEmbedParametersParams {
    /// Triangle-smooths the matrix.
    pub do_triangle_smoothing: bool,
    /// Scales the van der Waals lower bounds down.
    pub scale_vdw: bool,
    /// Sets 1-5 bounds from the topology.
    pub set15bounds: bool,
    /// Sets 1-4 bounds from the topology.
    pub set14bounds: bool,
    /// Sets 1-3 bounds from the topology.
    pub set13bounds: bool,
}

impl Default for GetMoleculeBoundsMatrixEmbedParametersParams {
    fn default() -> GetMoleculeBoundsMatrixEmbedParametersParams {
        GetMoleculeBoundsMatrixEmbedParametersParams {
            do_triangle_smoothing: true,
            scale_vdw: false,
            set15bounds: true,
            set14bounds: true,
            set13bounds: true,
        }
    }
}

/// The molecule's distance bounds, as rows: lower bounds below the diagonal,
/// upper bounds above it. Mirrors `rdDistGeom.GetMoleculeBoundsMatrix`.
/// `Error::Value` for a molecule without atoms and above 46340 atoms (U3).
pub fn get_molecule_bounds_matrix(mol: &Mol) -> Result<Vec<Vec<f64>>> {
    get_molecule_bounds_matrix_with_params(mol, &GetMoleculeBoundsMatrixParams::default())
}

/// Mirrors `rdDistGeom.GetMoleculeBoundsMatrix(mol, set15bounds, scaleVDW,
/// doTriangleSmoothing, useMacrocycle14config, forceTransAmides,
/// set14bounds, set13bounds)`.
pub fn get_molecule_bounds_matrix_with_params(
    mol: &Mol,
    params: &GetMoleculeBoundsMatrixParams,
) -> Result<Vec<Vec<f64>>> {
    let embed_parameters = EmbedParameters {
        use_macrocycle14config: params.use_macrocycle14config,
        force_trans_amides: params.force_trans_amides,
        ..EmbedParameters::default()
    };
    get_molecule_bounds_matrix_embed_parameters_with_params(
        mol,
        &embed_parameters,
        &GetMoleculeBoundsMatrixEmbedParametersParams {
            do_triangle_smoothing: params.do_triangle_smoothing,
            scale_vdw: params.scale_vdw,
            set15bounds: params.set15bounds,
            set14bounds: params.set14bounds,
            set13bounds: params.set13bounds,
        },
    )
}

/// The molecule's distance bounds under `embed_parameters`. Mirrors
/// `rdDistGeom.GetMoleculeBoundsMatrix(mol, embedParameters)` (A5).
pub fn get_molecule_bounds_matrix_embed_parameters(
    mol: &Mol,
    embed_parameters: &EmbedParameters,
) -> Result<Vec<Vec<f64>>> {
    get_molecule_bounds_matrix_embed_parameters_with_params(
        mol,
        embed_parameters,
        &GetMoleculeBoundsMatrixEmbedParametersParams::default(),
    )
}

/// Mirrors `rdDistGeom.GetMoleculeBoundsMatrix(mol, embedParameters,
/// doTriangleSmoothing, scaleVDW, set15bounds, set14bounds, set13bounds)`.
pub fn get_molecule_bounds_matrix_embed_parameters_with_params(
    mol: &Mol,
    embed_parameters: &EmbedParameters,
    params: &GetMoleculeBoundsMatrixEmbedParametersParams,
) -> Result<Vec<Vec<f64>>> {
    let values = ffi::get_molecule_bounds_matrix(
        mol.as_ffi(),
        &embed_parameters.to_ffi()?,
        params.do_triangle_smoothing,
        params.scale_vdw,
        params.set15bounds,
        params.set14bounds,
        params.set13bounds,
    )?;
    Ok(rows(&values, mol.get_num_atoms() as usize))
}

/// An experimental torsion's terms.
/// A Rust addition, for the keys of Python's dictionary.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum TorsionTerms {
    /// A cosine series: Python's `V` and `signs`.
    Cosine {
        /// The force constants, Python's `V`.
        v: Vec<f64>,
        /// The signs of the terms.
        signs: Vec<i32>,
    },
}

/// A bond matched by an experimental torsion pattern: one of the
/// dictionaries `rdDistGeom.GetExperimentalTorsions` returns.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct ExperimentalTorsion {
    /// Python's `bondIndex`.
    pub bond_index: u32,
    /// Python's `torsionIndex`, the pattern's index.
    pub torsion_index: u64,
    /// Python's `smarts`, the pattern.
    pub smarts: String,
    /// Python's `atomIndices`, the four atoms of the torsion.
    pub atom_indices: Vec<u32>,
    /// The terms.
    pub terms: TorsionTerms,
}

/// The bonds that match experimental torsion patterns, with Python's keyword
/// defaults. Mirrors `rdDistGeom.GetExperimentalTorsions(mol)`, whose
/// defaults are those of `embed_molecule`.
pub fn get_experimental_torsions(mol: &Mol) -> Result<Vec<ExperimentalTorsion>> {
    get_experimental_torsions_with_params(mol, &keyword_defaults())
}

/// Mirrors `rdDistGeom.GetExperimentalTorsions(mol, embedParams)`, and its
/// keyword form, whose arguments are fields of `embed_params`
/// (`printExpTorsionAngles` is `verbose`). RDKit gives `Error::Value` for an
/// `et_version` other than 1 or 2 when it reads the torsion preferences.
pub fn get_experimental_torsions_with_params(
    mol: &Mol,
    embed_params: &EmbedParameters,
) -> Result<Vec<ExperimentalTorsion>> {
    let torsions = ffi::get_experimental_torsions(
        mol.as_ffi(),
        embed_params.use_exp_torsion_angle_prefs,
        embed_params.use_small_ring_torsions,
        embed_params.use_macrocycle_torsions,
        embed_params.use_basic_knowledge,
        embed_params.et_version,
        embed_params.verbose,
    )?;
    torsions
        .into_iter()
        .map(|t| {
            Ok(ExperimentalTorsion {
                bond_index: t.bond_index,
                torsion_index: t.torsion_index,
                smarts: to_string(t.smarts)?,
                atom_indices: t.atom_indices,
                terms: TorsionTerms::Cosine {
                    v: t.v,
                    signs: t.signs,
                },
            })
        })
        .collect()
}

/// The failure causes in the order RDKit checks them, with each one's index
/// in the failure counts. Mirrors `rdDistGeom.OrderedEmbedFailureCauses()`.
pub fn ordered_embed_failure_causes() -> Vec<(EmbedFailureCauses, usize)> {
    ordered_embed_failure_causes_with_legacy_implementation(true)
}

/// Mirrors
/// `rdDistGeom.OrderedEmbedFailureCauses(legacyImplementation)`; the
/// all-in-one list, like Python's, gives `LinearDoubleBond` the index 9.
pub fn ordered_embed_failure_causes_with_legacy_implementation(
    legacy_implementation: bool,
) -> Vec<(EmbedFailureCauses, usize)> {
    use EmbedFailureCauses as E;
    if legacy_implementation {
        vec![
            (E::InitialCoords, 0),
            (E::FirstMinimization, 1),
            (E::CheckTetrahedralCenters, 2),
            (E::CheckChiralCenters, 3),
            (E::MinimizeFourthDimension, 4),
            (E::EtkMinimization, 5),
            (E::LinearDoubleBond, 8),
            (E::CheckChiralCenters2, 10),
            (E::FinalChiralBounds, 6),
            (E::FinalCenterInVolume, 7),
            (E::BadDoubleBondStereo, 9),
            (E::ExceededTimeout, 11),
        ]
    } else {
        vec![
            (E::InitialCoords, 0),
            (E::Minimization, 12),
            (E::KtermViolation, 13),
            (E::CheckTetrahedralCenters, 2),
            (E::CheckChiralCenters, 3),
            (E::FinalChiralBounds, 6),
            (E::FinalCenterInVolume, 7),
            (E::BadDoubleBondStereo, 9),
            (E::Clash, 14),
            (E::LinearDoubleBond, 9),
            (E::ExceededTimeout, 11),
        ]
    }
}
