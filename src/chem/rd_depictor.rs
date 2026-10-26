// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdDepictor`: 2D coordinates for molecules.
//!
//! RDKit keeps its depiction settings in process-wide globals that every
//! coordinate generation reads: whether CoordGen is preferred, the ring
//! system templates, and the bond length, which Python's `Compute2DCoords`
//! changes for one call. The wrappers serialize 2D depiction across threads
//! (CONVENTIONS.md U4), so a call sees one set of settings throughout.
//!
//! Depiction functions refuse, with `Error::Value`, coordinates they read
//! (coordinate maps, a reference's or the molecule's conformer, ring system
//! templates) and bond lengths that are not finite or are beyond 1e5 in
//! magnitude: RDKit sorts angles it derives from them and casts some to
//! `int`. Through RDKit's own depictor they also refuse molecules above
//! 46340 atoms (its distance matrix is computed in `int`) or, with an
//! explicit hydrogen, above 21474 (so are its atom ranks), a square-planar
//! atom without neighbors, and a trigonal-bipyramidal atom with more than
//! three equatorial neighbors (RDKit reads past its tables for both).
//! Through CoordGen, when it is preferred, they refuse molecules whose
//! Morgan scores overflow CoordGen's `int`, such as a chain of 40 carbons,
//! or whose formal charges do (U9); `force_rdkit` depicts those with
//! RDKit's depictor.

use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::path::Path;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::depictor::ffi;
use crate::geometry::Point2D;
use crate::util::path_bytes;

/// Whether RDKit was built with CoordGen. Mirrors
/// `rdDepictor.IsCoordGenSupportAvailable`.
pub fn is_coord_gen_support_available() -> bool {
    expect_ok(ffi::is_coord_gen_support_available())
}

/// Sets whether 2D depiction uses CoordGen rather than RDKit's own
/// depictor, for every thread. Mirrors `rdDepictor.SetPreferCoordGen`; like
/// it, does nothing without CoordGen support.
pub fn set_prefer_coord_gen(val: bool) {
    expect_ok(ffi::set_prefer_coord_gen(val));
}

/// Whether 2D depiction uses CoordGen. Mirrors
/// `rdDepictor.GetPreferCoordGen`; always false without CoordGen support.
pub fn get_prefer_coord_gen() -> bool {
    expect_ok(ffi::get_prefer_coord_gen())
}

/// Sets whether 2D depiction prefers CoordGen until dropped, then restores
/// the preference it replaced. Mirrors `rdDepictor.UsingCoordGen`, a context
/// manager in Python.
///
/// The preference is process-wide, so while the guard lives, depictions on
/// every thread use it. Setting and restoring it each take the depiction
/// lock, so a depiction on another thread sees one preference for its whole
/// call; the guard does not hold the lock in between, so the thread that
/// made it can depict. Guards restore in the order they are dropped, so
/// guards that overlap on several threads can leave another guard's value
/// behind, as in Python. The guard is neither `Send` nor `Sync`, so it is
/// dropped on the thread that made it.
///
/// ```no_run
/// use rdkit_cxx::chem::mol_from_smiles;
/// use rdkit_cxx::chem::rd_depictor::{UsingCoordGen, compute_2d_coords};
///
/// let mut mol = mol_from_smiles("c1ccccc1O").unwrap();
/// {
///     let _coord_gen = UsingCoordGen::new(true);
///     compute_2d_coords(&mut mol).unwrap();
/// }
/// ```
///
/// ```compile_fail,E0277
/// use rdkit_cxx::chem::rd_depictor::UsingCoordGen;
///
/// fn send<T: Send>(_: T) {}
/// send(UsingCoordGen::new(true));
/// ```
pub struct UsingCoordGen {
    initial: bool,
    _not_send: PhantomData<*const ()>,
}

impl std::fmt::Debug for UsingCoordGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UsingCoordGen").finish_non_exhaustive()
    }
}

impl UsingCoordGen {
    /// Prefers CoordGen, or not, until the guard is dropped. Mirrors
    /// `UsingCoordGen(temp_state)` and its `__enter__`.
    pub fn new(temp_state: bool) -> UsingCoordGen {
        UsingCoordGen {
            initial: expect_ok(ffi::replace_prefer_coord_gen(temp_state)),
            _not_send: PhantomData,
        }
    }
}

impl Drop for UsingCoordGen {
    fn drop(&mut self) {
        expect_ok(ffi::set_prefer_coord_gen(self.initial));
    }
}

/// Replaces the ring system templates with those in a file: one CXSMILES per
/// line, each a single ring system with 2D coordinates. Mirrors
/// `rdDepictor.SetRingSystemTemplates`. `Error::Value` (Python's
/// `ValueError`) for a file that cannot be read or holds an invalid
/// template, which leaves the templates as they were, and for template
/// coordinates that are not finite or beyond 1e5.
pub fn set_ring_system_templates(template_path: impl AsRef<Path>) -> Result<()> {
    Ok(ffi::set_ring_system_templates(path_bytes(
        template_path.as_ref(),
    ))?)
}

/// Adds the ring system templates in a file, which take precedence over
/// those already loaded. Mirrors `rdDepictor.AddRingSystemTemplates`; errors
/// as `set_ring_system_templates`.
pub fn add_ring_system_templates(template_path: impl AsRef<Path>) -> Result<()> {
    Ok(ffi::add_ring_system_templates(path_bytes(
        template_path.as_ref(),
    ))?)
}

/// Replaces the ring system templates with RDKit's defaults. Mirrors
/// `rdDepictor.LoadDefaultRingSystemTemplates`.
pub fn load_default_ring_system_templates() {
    expect_ok(ffi::load_default_ring_system_templates());
}

/// How a depiction is constrained to a reference. Mirrors
/// `rdDepictor.ConstrainedDepictionParams`, plus the `confId` and `refPatt`
/// arguments of `GenerateDepictionMatching2DStructure` (A3).
#[derive(Clone, Debug)]
pub struct ConstrainedDepictionParams {
    /// Without a match to the reference, depicts the molecule
    /// unconstrained instead of failing.
    pub accept_failure: bool,
    /// Uses RDKit's depictor even when CoordGen is preferred.
    pub force_rdkit: bool,
    /// Ignores terminal dummy atoms of the reference that match an
    /// implicit hydrogen or are attached to a query atom.
    pub allow_r_groups: bool,
    /// Keeps the coordinates of conformer `existing_conf_id` (or computes
    /// unconstrained ones) and aligns them rigidly to the reference, rather
    /// than building the rest of the molecule around the reference's
    /// coordinates.
    pub align_only: bool,
    /// Updates or clears existing mol block wedging as needed, rather than
    /// always keeping it.
    pub adjust_mol_block_wedging: bool,
    /// The molecule's conformer to align (`align_only`) or to compare
    /// wedging against, or -1 for the first.
    pub existing_conf_id: i32,
    /// Uses ring system templates for complex ring systems.
    pub use_ring_templates: bool,
    /// The reference's conformer, or -1 for the first (Python's `confId`).
    pub conf_id: i32,
    /// A query that maps the molecule's atoms to the reference's (Python's
    /// `refPatt`); the reference itself when `None`. The atom-map forms
    /// ignore it.
    pub ref_patt: Option<Mol>,
}

impl Default for ConstrainedDepictionParams {
    fn default() -> ConstrainedDepictionParams {
        let data = expect_ok(ffi::constrained_depiction_params_default());
        ConstrainedDepictionParams {
            accept_failure: data.accept_failure,
            force_rdkit: data.force_rdkit,
            allow_r_groups: data.allow_r_groups,
            align_only: data.align_only,
            adjust_mol_block_wedging: data.adjust_mol_block_wedging,
            existing_conf_id: data.existing_conf_id,
            use_ring_templates: data.use_ring_templates,
            conf_id: data.conf_id,
            ref_patt: None,
        }
    }
}

impl ConstrainedDepictionParams {
    fn to_ffi(&self) -> ffi::ConstrainedDepictionParamsData {
        ffi::ConstrainedDepictionParamsData {
            accept_failure: self.accept_failure,
            force_rdkit: self.force_rdkit,
            allow_r_groups: self.allow_r_groups,
            align_only: self.align_only,
            adjust_mol_block_wedging: self.adjust_mol_block_wedging,
            existing_conf_id: self.existing_conf_id,
            use_ring_templates: self.use_ring_templates,
            conf_id: self.conf_id,
        }
    }
}

/// Parameters for `compute_2d_coords_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct Compute2DCoordsParams {
    /// Orients the depiction canonically (ignored with a coordinate map).
    pub canon_orient: bool,
    /// Removes the molecule's conformers first, rather than adding one.
    pub clear_confs: bool,
    /// Fixed 2D positions for some atoms.
    pub coord_map: BTreeMap<u32, Point2D>,
    /// The rotatable bonds flipped at random in each sample.
    pub n_flips_per_sample: u32,
    /// The number of random samples; with `n_flips_per_sample`, replaces
    /// RDKit's bond flips for collisions.
    pub n_sample: u32,
    /// The seed for sampling; nonpositive leaves RDKit's generator as it is.
    pub sample_seed: i32,
    /// Nonzero also permutes the bonds around degree-4 atoms while sampling
    /// (an `int` in Python).
    pub permute_deg4_nodes: i32,
    /// The bond length for this call; nonpositive keeps RDKit's 1.5.
    pub bond_length: f64,
    /// Uses RDKit's depictor even when CoordGen is preferred.
    pub force_rdkit: bool,
    /// Uses ring system templates for complex ring systems.
    pub use_ring_templates: bool,
    /// Opens angles along paths between colliding atoms.
    pub use_path_angle_expansion: bool,
}

impl Default for Compute2DCoordsParams {
    fn default() -> Compute2DCoordsParams {
        Compute2DCoordsParams {
            canon_orient: true,
            clear_confs: true,
            coord_map: BTreeMap::new(),
            n_flips_per_sample: 0,
            n_sample: 0,
            sample_seed: 0,
            permute_deg4_nodes: 0,
            bond_length: -1.0,
            force_rdkit: false,
            use_ring_templates: false,
            use_path_angle_expansion: true,
        }
    }
}

/// Gives the molecule a conformer with 2D coordinates and returns its id.
/// Mirrors `rdDepictor.Compute2DCoords`. `Error::Value` for a molecule the
/// depictor would mishandle (see the module documentation).
pub fn compute_2d_coords(mol: &mut Mol) -> Result<u32> {
    compute_2d_coords_with_params(mol, &Compute2DCoordsParams::default())
}

/// Mirrors `rdDepictor.Compute2DCoords(mol, canonOrient, clearConfs,
/// coordMap, nFlipsPerSample, nSample, sampleSeed, permuteDeg4Nodes,
/// bondLength, forceRDKit, useRingTemplates, usePathAngleExpansion)`.
/// `Error::Index` for a coordinate map atom past the end, and `Error::Value`
/// for coordinates or a bond length that are not finite or beyond 1e5.
pub fn compute_2d_coords_with_params(mol: &mut Mol, params: &Compute2DCoordsParams) -> Result<u32> {
    let data = ffi::Compute2DCoordsData {
        canon_orient: params.canon_orient,
        clear_confs: params.clear_confs,
        coord_map_keys: params.coord_map.keys().copied().collect(),
        coord_map_values: params.coord_map.values().copied().collect(),
        n_flips_per_sample: params.n_flips_per_sample,
        n_sample: params.n_sample,
        sample_seed: params.sample_seed,
        permute_deg4_nodes: params.permute_deg4_nodes,
        bond_length: params.bond_length,
        force_rdkit: params.force_rdkit,
        use_ring_templates: params.use_ring_templates,
        use_path_angle_expansion: params.use_path_angle_expansion,
    };
    Ok(ffi::compute_2d_coords(mol.pin_mut(), &data)?)
}

/// Parameters for `compute_2d_coords_mimic_distmat_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct Compute2DCoordsMimicDistmatParams {
    /// Orients the depiction canonically.
    pub canon_orient: bool,
    /// Removes the molecule's conformers first, rather than adding one.
    pub clear_confs: bool,
    /// The weight, from 0 to 1, of matching the distances rather than
    /// spreading the atoms out.
    pub weight_dist_mat: f64,
    /// The rotatable bonds flipped at random in each sample.
    pub n_flips_per_sample: u32,
    /// The number of random samples.
    pub n_sample: u32,
    /// The seed for sampling; nonpositive leaves RDKit's generator as it is.
    pub sample_seed: i32,
    /// Also permutes the bonds around degree-4 atoms while sampling.
    pub permute_deg4_nodes: bool,
    /// The bond length for this call; nonpositive keeps RDKit's 1.5.
    pub bond_length: f64,
    /// Accepted for Python's signature; RDKit always uses its own depictor
    /// here.
    pub force_rdkit: bool,
}

impl Default for Compute2DCoordsMimicDistmatParams {
    fn default() -> Compute2DCoordsMimicDistmatParams {
        Compute2DCoordsMimicDistmatParams {
            canon_orient: false,
            clear_confs: true,
            weight_dist_mat: 0.5,
            n_flips_per_sample: 3,
            n_sample: 100,
            sample_seed: 100,
            permute_deg4_nodes: true,
            bond_length: -1.0,
            force_rdkit: false,
        }
    }
}

/// Gives the molecule a conformer whose interatomic distances mimic
/// `dist_mat` and returns its id. `dist_mat` is the condensed lower
/// triangle, entry `j * (j - 1) / 2 + i` for atoms `i < j`; negative entries
/// are ignored. Mirrors `rdDepictor.Compute2DCoordsMimicDistmat`.
/// `Error::Value` for a matrix of the wrong size.
pub fn compute_2d_coords_mimic_distmat(mol: &mut Mol, dist_mat: &[f64]) -> Result<u32> {
    compute_2d_coords_mimic_distmat_with_params(
        mol,
        dist_mat,
        &Compute2DCoordsMimicDistmatParams::default(),
    )
}

/// Mirrors `rdDepictor.Compute2DCoordsMimicDistmat(mol, distMat,
/// canonOrient, clearConfs, weightDistMat, nFlipsPerSample, nSample,
/// sampleSeed, permuteDeg4Nodes, bondLength, forceRDKit)`.
pub fn compute_2d_coords_mimic_distmat_with_params(
    mol: &mut Mol,
    dist_mat: &[f64],
    params: &Compute2DCoordsMimicDistmatParams,
) -> Result<u32> {
    let data = ffi::Compute2DCoordsMimicDistmatData {
        canon_orient: params.canon_orient,
        clear_confs: params.clear_confs,
        weight_dist_mat: params.weight_dist_mat,
        n_flips_per_sample: params.n_flips_per_sample,
        n_sample: params.n_sample,
        sample_seed: params.sample_seed,
        permute_deg4_nodes: params.permute_deg4_nodes,
        bond_length: params.bond_length,
        force_rdkit: params.force_rdkit,
    };
    Ok(ffi::compute_2d_coords_mimic_distmat(
        mol.pin_mut(),
        dist_mat,
        &data,
    )?)
}

/// Depicts the molecule with the atoms that match the reference at the
/// reference's 2D coordinates, and returns the (reference atom, molecule
/// atom) pairs it matched. The molecule comes back with one conformer.
/// Mirrors `rdDepictor.GenerateDepictionMatching2DStructure(mol,
/// reference)`. `Error::Value` (Python's `ValueError`) when the
/// molecule does not match the reference, or for a missing reference
/// conformer.
pub fn generate_depiction_matching_2d_structure(
    mol: &mut Mol,
    reference: &Mol,
) -> Result<Vec<(u32, u32)>> {
    generate_depiction_matching_2d_structure_with_params(
        mol,
        reference,
        &ConstrainedDepictionParams::default(),
    )
}

/// Mirrors `rdDepictor.GenerateDepictionMatching2DStructure(mol, reference,
/// confId, refPatt, params)` and its keyword form (`acceptFailure`,
/// `forceRDKit`, `allowRGroups`).
pub fn generate_depiction_matching_2d_structure_with_params(
    mol: &mut Mol,
    reference: &Mol,
    params: &ConstrainedDepictionParams,
) -> Result<Vec<(u32, u32)>> {
    let pattern = params.ref_patt.as_ref();
    let pairs = ffi::generate_depiction_matching_2d_structure(
        mol.pin_mut(),
        reference.as_ffi(),
        pattern.unwrap_or(reference).as_ffi(),
        pattern.is_some(),
        &params.to_ffi(),
    )?;
    Ok(pairs.iter().map(|p| (p.first, p.second)).collect())
}

/// Depicts the molecule with the atoms in `atom_map`'s (reference atom,
/// molecule atom) pairs at the reference's coordinates. Mirrors
/// `rdDepictor.GenerateDepictionMatching2DStructure(mol, reference,
/// atomMap)` (A5). `Error::Index` for an atom past the end of either
/// molecule, and `Error::Value` (Python's `ValueError`) for a map
/// longer than the reference.
pub fn generate_depiction_matching_2d_structure_atom_map(
    mol: &mut Mol,
    reference: &Mol,
    atom_map: &[(u32, u32)],
) -> Result<()> {
    generate_depiction_matching_2d_structure_atom_map_with_params(
        mol,
        reference,
        atom_map,
        &ConstrainedDepictionParams::default(),
    )
}

/// Mirrors `rdDepictor.GenerateDepictionMatching2DStructure(mol, reference,
/// atomMap, confId, params)` and `(mol, reference, atomMap, confId,
/// forceRDKit)`. Ignores `params.ref_patt`.
pub fn generate_depiction_matching_2d_structure_atom_map_with_params(
    mol: &mut Mol,
    reference: &Mol,
    atom_map: &[(u32, u32)],
    params: &ConstrainedDepictionParams,
) -> Result<()> {
    let pairs: Vec<IndexPair> = atom_map
        .iter()
        .map(|&(first, second)| IndexPair { first, second })
        .collect();
    Ok(ffi::generate_depiction_matching_2d_structure_atom_map(
        mol.pin_mut(),
        reference.as_ffi(),
        &pairs,
        &params.to_ffi(),
    )?)
}

/// Parameters for `generate_depiction_matching_3d_structure_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug)]
pub struct GenerateDepictionMatching3DStructureParams {
    /// The reference's conformer, or -1 for the first.
    pub conf_id: i32,
    /// A query that maps a subset of the molecule's atoms to the
    /// reference's; all atoms, in order, when `None`.
    pub ref_patt: Option<Mol>,
    /// Without a match, depicts the molecule unconstrained instead of
    /// failing.
    pub accept_failure: bool,
    /// Uses RDKit's depictor even when CoordGen is preferred.
    pub force_rdkit: bool,
}

impl Default for GenerateDepictionMatching3DStructureParams {
    fn default() -> GenerateDepictionMatching3DStructureParams {
        GenerateDepictionMatching3DStructureParams {
            conf_id: -1,
            ref_patt: None,
            accept_failure: false,
            force_rdkit: false,
        }
    }
}

/// Depicts the molecule so that its 2D distances mimic the 3D distances of
/// the reference, by default the same molecule with 3D coordinates. Mirrors
/// `rdDepictor.GenerateDepictionMatching3DStructure`. `Error::Value`
/// (Python's `ValueError`) when the reference is smaller than the
/// molecule or the pattern matches neither, and for a missing reference
/// conformer; `Error::Runtime` when a pattern maps to a reference atom past
/// the reference's conformer.
pub fn generate_depiction_matching_3d_structure(mol: &mut Mol, reference: &Mol) -> Result<()> {
    generate_depiction_matching_3d_structure_with_params(
        mol,
        reference,
        &GenerateDepictionMatching3DStructureParams::default(),
    )
}

/// Mirrors `rdDepictor.GenerateDepictionMatching3DStructure(mol, reference,
/// confId, refPatt, acceptFailure, forceRDKit)`.
pub fn generate_depiction_matching_3d_structure_with_params(
    mol: &mut Mol,
    reference: &Mol,
    params: &GenerateDepictionMatching3DStructureParams,
) -> Result<()> {
    let pattern = params.ref_patt.as_ref();
    Ok(ffi::generate_depiction_matching_3d_structure(
        mol.pin_mut(),
        reference.as_ffi(),
        params.conf_id,
        pattern.unwrap_or(reference).as_ffi(),
        pattern.is_some(),
        params.accept_failure,
        params.force_rdkit,
    )?)
}

/// Parameters for `straighten_depiction_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct StraightenDepictionParams {
    /// The conformer to rotate, or -1 for the first.
    pub conf_id: i32,
    /// Rotates as little as possible, to bond angles of 0, 30, 60, or 90
    /// degrees with the x axis, rather than to 30 or 90.
    pub minimize_rotation: bool,
}

impl Default for StraightenDepictionParams {
    fn default() -> StraightenDepictionParams {
        StraightenDepictionParams {
            conf_id: -1,
            minimize_rotation: false,
        }
    }
}

/// Rotates the depiction so that most bonds make 30 or 90 degrees with the
/// x axis. Mirrors `rdDepictor.StraightenDepiction`. Does nothing to a
/// molecule without bonds; otherwise `Error::Value` for a missing conformer
/// or coordinates that are not finite or beyond 1e5.
pub fn straighten_depiction(mol: &mut Mol) -> Result<()> {
    straighten_depiction_with_params(mol, &StraightenDepictionParams::default())
}

/// Mirrors `rdDepictor.StraightenDepiction(mol, confId, minimizeRotation)`.
pub fn straighten_depiction_with_params(
    mol: &mut Mol,
    params: &StraightenDepictionParams,
) -> Result<()> {
    Ok(ffi::straighten_depiction(
        mol.pin_mut(),
        params.conf_id,
        params.minimize_rotation,
    )?)
}

/// Parameters for `normalize_depiction_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct NormalizeDepictionParams {
    /// The conformer to normalize, or -1 for the first.
    pub conf_id: i32,
    /// Positive aligns the main axis with x, negative with y, and 0 leaves
    /// the orientation alone.
    pub canonicalize: i32,
    /// The scale to apply; negative scales bonds to RDKit's standard length.
    pub scale_factor: f64,
}

impl Default for NormalizeDepictionParams {
    fn default() -> NormalizeDepictionParams {
        NormalizeDepictionParams {
            conf_id: -1,
            canonicalize: 1,
            scale_factor: -1.0,
        }
    }
}

/// Centers, orients, and scales the depiction, and returns the scale it
/// applied (-1 for none). Mirrors `rdDepictor.NormalizeDepiction`. Returns
/// -1 for a molecule without bonds; otherwise `Error::Value` for a missing
/// conformer or coordinates that are not finite or beyond 1e5.
pub fn normalize_depiction(mol: &mut Mol) -> Result<f64> {
    normalize_depiction_with_params(mol, &NormalizeDepictionParams::default())
}

/// Mirrors `rdDepictor.NormalizeDepiction(mol, confId, canonicalize,
/// scaleFactor)`.
pub fn normalize_depiction_with_params(
    mol: &mut Mol,
    params: &NormalizeDepictionParams,
) -> Result<f64> {
    Ok(ffi::normalize_depiction(
        mol.pin_mut(),
        params.conf_id,
        params.canonicalize,
        params.scale_factor,
    )?)
}
