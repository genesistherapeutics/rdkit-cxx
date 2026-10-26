// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.AllChem`, a Python module that gathers RDKit's
//! chemistry modules in one namespace, as this one re-exports them, and
//! adds helpers built on them: transforming, comparing, and constraining
//! conformers, grid shapes and volumes, library enumeration, and bond
//! orders from a template.
//!
//! Python also patches `Compute2DCoords` and `ComputeGasteigerCharges` onto
//! `Mol` as methods; here they are `rd_depictor::compute_2d_coords` and
//! `rd_partial_charges::compute_gasteiger_charges`.

pub use crate::chem::chemical_features::*;
pub use crate::chem::enumerate_stereoisomers::{StereoEnumerationOptions, enumerate_stereoisomers};
pub use crate::chem::rd_chem_reactions::*;
pub use crate::chem::rd_depictor::*;
pub use crate::chem::rd_dist_geom::*;
pub use crate::chem::rd_fingerprint_generator::*;
pub use crate::chem::rd_force_field_helpers::*;
pub use crate::chem::rd_mol_align::*;
pub use crate::chem::rd_mol_descriptors::*;
pub use crate::chem::rd_mol_enumerator::*;
pub use crate::chem::rd_mol_transforms::*;
pub use crate::chem::rd_partial_charges::*;
pub use crate::chem::rd_reduced_graphs::*;
pub use crate::chem::rd_shape_helpers::*;
pub use crate::chem::rd_sln_parse::*;
pub use crate::chem::rdqueries::*;
pub use crate::chem::*;
// Names more than one of those modules define are the last import's, as in
// Python, where each `import *` replaces the names before it.
pub use crate::chem::rd_chem_reactions::SanitizeFlags;

use crate::force_field::rd_force_field::{ForceField, MinimizeParams};
use crate::geometry::rd_geometry::{UniformGrid3D, UniformGrid3DParams};
use crate::py::float_repr;
use crate::{Error, Result};

/// Parameters for `transform_mol_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct TransformMolParams {
    /// The conformer to transform; -1 for the first.
    pub conf_id: i32,
    /// Keeps the other conformers; otherwise only the one with id
    /// `conf_id` (0 for -1) stays, with id 0.
    pub keep_confs: bool,
}

impl Default for TransformMolParams {
    fn default() -> TransformMolParams {
        TransformMolParams {
            conf_id: -1,
            keep_confs: false,
        }
    }
}

/// Applies the transform, a 4x4 matrix as rows, to the molecule's first
/// conformer and removes its other conformers. Mirrors
/// `AllChem.TransformMol(mol, tform)`.
pub fn transform_mol(mol: &mut Mol, tform: &[Vec<f64>]) -> Result<()> {
    transform_mol_with_params(mol, tform, &TransformMolParams::default())
}

/// As [`transform_mol`], for the conformer `params.conf_id`. Without
/// `keep_confs`, Python keeps the conformer with id 0 when `conf_id` is -1,
/// which need not be the first conformer it transformed, and so does this
/// port; a molecule without that conformer is `Error::Value`. Mirrors
/// `AllChem.TransformMol(mol, tform, confId, keepConfs)`.
pub fn transform_mol_with_params(
    mol: &mut Mol,
    tform: &[Vec<f64>],
    params: &TransformMolParams,
) -> Result<()> {
    transform_conformer(&mut mol.get_conformer_with_id_mut(params.conf_id)?, tform)?;
    if !params.keep_confs {
        let conf_id = if params.conf_id == -1 {
            0
        } else {
            params.conf_id
        };
        let all_conf_ids: Vec<u32> = mol.get_conformers().map(|c| c.get_id()).collect();
        for cid in all_conf_ids {
            if i64::from(cid) != i64::from(conf_id) {
                mol.remove_conformer(cid);
            }
        }
        mol.get_conformer_with_id_mut(conf_id)?.set_id(0);
    }
    Ok(())
}

/// Parameters for `compute_mol_shape_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputeMolShapeParams {
    /// The conformer to encode; -1 for the first.
    pub conf_id: i32,
    /// The grid's dimensions.
    pub box_dim: (f64, f64, f64),
    /// The grid's spacing.
    pub spacing: f64,
    /// `EncodeShape`'s other arguments, which Python passes as keywords;
    /// their `conf_id` is replaced by `conf_id`.
    pub encode_shape: EncodeShapeParams,
}

impl Default for ComputeMolShapeParams {
    fn default() -> ComputeMolShapeParams {
        ComputeMolShapeParams {
            conf_id: -1,
            box_dim: (20.0, 20.0, 20.0),
            spacing: 0.5,
            encode_shape: EncodeShapeParams::default(),
        }
    }
}

/// The molecule's shape, encoded on a 20 Å grid with 0.5 Å spacing.
/// Mirrors `AllChem.ComputeMolShape(mol)`.
pub fn compute_mol_shape(mol: &Mol) -> Result<UniformGrid3D> {
    compute_mol_shape_with_params(mol, &ComputeMolShapeParams::default())
}

/// As [`compute_mol_shape`]. Mirrors `AllChem.ComputeMolShape(mol, confId,
/// boxDim, spacing, **kwargs)`.
pub fn compute_mol_shape_with_params(
    mol: &Mol,
    params: &ComputeMolShapeParams,
) -> Result<UniformGrid3D> {
    let (x, y, z) = params.box_dim;
    let mut res = UniformGrid3D::new_with_params(
        x,
        y,
        z,
        &UniformGrid3DParams {
            spacing: params.spacing,
            ..Default::default()
        },
    )?;
    encode_shape_with_params(
        mol,
        &mut res,
        &EncodeShapeParams {
            conf_id: params.conf_id,
            ..params.encode_shape.clone()
        },
    )?;
    Ok(res)
}

/// Parameters for `compute_mol_volume_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputeMolVolumeParams {
    /// The conformer to measure; -1 for the first.
    pub conf_id: i32,
    /// The grid's spacing.
    pub grid_spacing: f64,
    /// The margin around the molecule.
    pub box_margin: f64,
}

impl Default for ComputeMolVolumeParams {
    fn default() -> ComputeMolVolumeParams {
        ComputeMolVolumeParams {
            conf_id: -1,
            grid_spacing: 0.2,
            box_margin: 2.0,
        }
    }
}

/// The volume of the molecule's first conformer, from a grid encoding of
/// its shape: the volume of the grid points inside it. Mirrors
/// `AllChem.ComputeMolVolume(mol)`.
pub fn compute_mol_volume(mol: &Mol) -> Result<f64> {
    compute_mol_volume_with_params(mol, &ComputeMolVolumeParams::default())
}

/// As [`compute_mol_volume`]. Mirrors `AllChem.ComputeMolVolume(mol,
/// confId, gridSpacing, boxMargin)`.
pub fn compute_mol_volume_with_params(mol: &Mol, params: &ComputeMolVolumeParams) -> Result<f64> {
    let mut mol = mol.clone();
    let (lo, hi) = {
        let mut conf = mol.get_conformer_with_id_mut(params.conf_id)?;
        canonicalize_conformer_with_params(
            &mut conf,
            &CanonicalizeConformerParams {
                ignore_hs: false,
                ..Default::default()
            },
        )?;
        compute_conf_box(&conf)?
    };
    let margin = 2.0 * params.box_margin;
    let mut shape = UniformGrid3D::new_with_params(
        hi.x - lo.x + margin,
        hi.y - lo.y + margin,
        hi.z - lo.z + margin,
        &UniformGrid3DParams {
            spacing: params.grid_spacing,
            ..Default::default()
        },
    )?;
    encode_shape_with_params(
        &mol,
        &mut shape,
        &EncodeShapeParams {
            conf_id: params.conf_id,
            ignore_hs: false,
            vdw_scale: 1.0,
            ..Default::default()
        },
    )?;
    let occ_vect = shape.get_occupancy_vect();
    let mut voxels = 0u32;
    for i in 0..occ_vect.len() {
        if occ_vect.get_val(i)? == 3 {
            voxels += 1;
        }
    }
    Ok(params.grid_spacing.powi(3) * f64::from(voxels))
}

/// Parameters for `get_conformer_rms_with_params` and
/// `get_conformer_rms_matrix_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetConformerRMSParams {
    /// The atoms to align on; `None` (or empty, as in Python) for all.
    pub atom_ids: Option<Vec<u32>>,
    /// Takes the conformers as aligned already.
    pub prealigned: bool,
}

impl GetConformerRMSParams {
    fn atom_ids(&self) -> Vec<u32> {
        self.atom_ids.clone().unwrap_or_default()
    }
}

/// The RMS distance between two conformers, over all atoms, after aligning
/// the second to the first (which leaves it aligned). Mirrors
/// `AllChem.GetConformerRMS(mol, confId1, confId2)`.
pub fn get_conformer_rms(mol: &mut Mol, conf_id1: i32, conf_id2: i32) -> Result<f64> {
    get_conformer_rms_with_params(mol, conf_id1, conf_id2, &GetConformerRMSParams::default())
}

/// As [`get_conformer_rms`], aligning on `params.atom_ids`, or not at all
/// with `prealigned`; the distance is over all atoms either way. Mirrors
/// `AllChem.GetConformerRMS(mol, confId1, confId2, atomIds, prealigned)`.
pub fn get_conformer_rms_with_params(
    mol: &mut Mol,
    conf_id1: i32,
    conf_id2: i32,
    params: &GetConformerRMSParams,
) -> Result<f64> {
    if !params.prealigned {
        align_mol_conformers_with_params(
            mol,
            &AlignMolConformersParams {
                atom_ids: params.atom_ids(),
                conf_ids: vec![conf_id1, conf_id2],
                ..Default::default()
            },
        )?;
    }
    let conf1 = mol.get_conformer_with_id(conf_id1)?;
    let conf2 = mol.get_conformer_with_id(conf_id2)?;
    let mut ssr = 0.0;
    for i in 0..mol.get_num_atoms() {
        let d = conf1
            .get_atom_position(i)?
            .distance(&conf2.get_atom_position(i)?);
        ssr += d * d;
    }
    ssr /= f64::from(mol.get_num_atoms());
    Ok(ssr.sqrt())
}

/// The RMS distances between every pair of the molecule's conformers, as the
/// lower triangle of the matrix, row by row (`[d(1,0), d(2,0), d(2,1), ...]`
/// for conformers in the molecule's order), after aligning every conformer
/// to the first (which leaves them aligned). Mirrors
/// `AllChem.GetConformerRMSMatrix(mol)`.
pub fn get_conformer_rms_matrix(mol: &mut Mol) -> Result<Vec<f64>> {
    get_conformer_rms_matrix_with_params(mol, &GetConformerRMSParams::default())
}

/// As [`get_conformer_rms_matrix`]. Mirrors
/// `AllChem.GetConformerRMSMatrix(mol, atomIds, prealigned)`.
pub fn get_conformer_rms_matrix_with_params(
    mol: &mut Mol,
    params: &GetConformerRMSParams,
) -> Result<Vec<f64>> {
    let conf_ids: Vec<i32> = mol.get_conformers().map(|c| c.get_id() as i32).collect();
    let rmsvals = if !params.prealigned {
        align_mol_conformers_with_params(
            mol,
            &AlignMolConformersParams {
                atom_ids: params.atom_ids(),
                ..Default::default()
            },
        )?
    } else {
        let mut rmsvals = Vec::new();
        for &conf_id in conf_ids.iter().skip(1) {
            rmsvals.push(get_conformer_rms_with_params(
                mol,
                conf_ids[0],
                conf_id,
                params,
            )?);
        }
        rmsvals
    };
    let prealigned = GetConformerRMSParams {
        prealigned: true,
        ..params.clone()
    };
    let mut cmat = Vec::new();
    for i in 1..conf_ids.len() {
        cmat.push(rmsvals[i - 1]);
        for j in 1..i {
            cmat.push(get_conformer_rms_with_params(
                mol,
                conf_ids[i],
                conf_ids[j],
                &prealigned,
            )?);
        }
    }
    Ok(cmat)
}

/// The products of one combination of reactants and the reactants. Mirrors
/// the `ProductReactants` named tuple that
/// `AllChem.EnumerateLibraryFromReaction(..., returnReactants=True)` yields.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct ProductReactants {
    /// One set of the reaction's products.
    pub products: Vec<Mol>,
    /// The reactants they came from, one from each sidechain set.
    pub reactants: Vec<Mol>,
}

/// The virtual library of a reaction and its sidechain sets, made as it is
/// asked for: each set of products, for each combination of reactants, one
/// from each set, the last set varying fastest. Built by
/// [`enumerate_library_from_reaction`].
pub struct EnumerateLibraryFromReaction<'a> {
    reaction: &'a mut ChemicalReaction,
    sidechain_sets: &'a [Vec<Mol>],
    /// The next combination, as an index into each set; `None` when done.
    next: Option<Vec<usize>>,
    pending: std::collections::VecDeque<ProductReactants>,
}

impl std::fmt::Debug for EnumerateLibraryFromReaction<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EnumerateLibraryFromReaction")
            .finish_non_exhaustive()
    }
}

impl EnumerateLibraryFromReaction<'_> {
    fn advance(&mut self) -> Option<Vec<usize>> {
        let current = self.next.take()?;
        let mut next = current.clone();
        for depth in (0..next.len()).rev() {
            next[depth] += 1;
            if next[depth] < self.sidechain_sets[depth].len() {
                self.next = Some(next);
                break;
            }
            next[depth] = 0;
        }
        Some(current)
    }

    /// The next set of products with its reactants. Mirrors a step of
    /// `AllChem.EnumerateLibraryFromReaction(..., returnReactants=True)`.
    pub fn next_with_reactants(&mut self) -> Option<Result<ProductReactants>> {
        while self.pending.is_empty() {
            let combination = self.advance()?;
            let chains: Vec<Mol> = combination
                .iter()
                .zip(self.sidechain_sets)
                .map(|(&i, set)| set[i].clone())
                .collect();
            match self.reaction.run_reactants(&chains) {
                Ok(prod_sets) => {
                    self.pending
                        .extend(prod_sets.into_iter().map(|products| ProductReactants {
                            products,
                            reactants: chains.clone(),
                        }))
                }
                Err(e) => {
                    self.next = None;
                    return Some(Err(e));
                }
            }
        }
        self.pending.pop_front().map(Ok)
    }
}

/// The next set of products.
impl Iterator for EnumerateLibraryFromReaction<'_> {
    type Item = Result<Vec<Mol>>;

    fn next(&mut self) -> Option<Result<Vec<Mol>>> {
        self.next_with_reactants()
            .map(|res| res.map(|pr| pr.products))
    }
}

/// The reaction's virtual library from the sidechain sets, one set per
/// reactant template, made as it is iterated. A number of sets other than
/// the reaction's reactant templates is `Error::Value`, as in Python, which
/// raises it at the first step; a reaction without reactant templates is
/// `Error::Index`, as Python's `IndexError`. Python's sets may be lazy
/// iterators; here they are slices. Mirrors
/// `AllChem.EnumerateLibraryFromReaction(reaction, sidechainSets)`.
pub fn enumerate_library_from_reaction<'a>(
    reaction: &'a mut ChemicalReaction,
    sidechain_sets: &'a [Vec<Mol>],
) -> Result<EnumerateLibraryFromReaction<'a>> {
    let n = reaction.get_num_reactant_templates() as usize;
    if sidechain_sets.len() != n {
        return Err(Error::Value {
            message: format!("{} sidechains provided, {n} required", sidechain_sets.len()),
        });
    }
    if n == 0 {
        return Err(Error::index(0, 0));
    }
    let next = sidechain_sets
        .iter()
        .all(|set| !set.is_empty())
        .then(|| vec![0; n]);
    Ok(EnumerateLibraryFromReaction {
        reaction,
        sidechain_sets,
        next,
        pending: Default::default(),
    })
}

/// As [`enumerate_library_from_reaction`], yielding each set of products
/// with its reactants. Mirrors
/// `AllChem.EnumerateLibraryFromReaction(reaction, sidechainSets,
/// returnReactants=True)` (A4).
pub fn enumerate_library_from_reaction_return_reactants<'a>(
    reaction: &'a mut ChemicalReaction,
    sidechain_sets: &'a [Vec<Mol>],
) -> Result<impl Iterator<Item = Result<ProductReactants>> + 'a> {
    let mut library = enumerate_library_from_reaction(reaction, sidechain_sets)?;
    Ok(std::iter::from_fn(move || library.next_with_reactants()))
}

/// A function that makes the force field `constrained_embed` cleans up with,
/// for the molecule's conformer 0. Mirrors `ConstrainedEmbed`'s
/// `getForceField` argument, which Python calls as `getForceField(mol,
/// confId=0)`; any function or closure will do.
pub type GetForceField = dyn Fn(&mut Mol) -> Result<ForceField<'_>>;

/// UFF's force field for the molecule's conformer 0, `constrained_embed`'s
/// default. Mirrors `rdForceFieldHelpers.UFFGetMoleculeForceField(mol,
/// confId=0)`.
fn uff_force_field(mol: &mut Mol) -> Result<ForceField<'_>> {
    uff_get_molecule_force_field_with_params(
        mol,
        &UFFGetMoleculeForceFieldParams {
            conf_id: 0,
            ..Default::default()
        },
    )
}

/// Parameters for `constrained_embed_with_params`, with Python's defaults
/// (A2).
#[derive(Clone)]
pub struct ConstrainedEmbedParams<'a> {
    /// Pulls the matched atoms toward the core's positions with extra
    /// points while cleaning up; otherwise constrains their distances.
    pub use_tethers: bool,
    /// The core's conformer that gives the constraints; -1 for the first.
    pub core_conf_id: i32,
    /// The embedder's random seed.
    pub randomseed: i32,
    /// Makes the force field for the cleanup.
    pub get_force_field: &'a GetForceField,
    /// The embedder's other arguments, which Python passes as keywords;
    /// their `coord_map` and `random_seed` are replaced. The default is
    /// `EmbedMolecule`'s keyword defaults.
    pub embed_params: EmbedParameters,
}

impl std::fmt::Debug for ConstrainedEmbedParams<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ConstrainedEmbedParams")
            .field("use_tethers", &self.use_tethers)
            .field("core_conf_id", &self.core_conf_id)
            .field("randomseed", &self.randomseed)
            .field("embed_params", &self.embed_params)
            .finish_non_exhaustive()
    }
}

impl Default for ConstrainedEmbedParams<'_> {
    fn default() -> Self {
        ConstrainedEmbedParams {
            use_tethers: true,
            core_conf_id: -1,
            randomseed: 2342,
            get_force_field: &uff_force_field,
            embed_params: EmbedParameters {
                use_exp_torsion_angle_prefs: true,
                use_basic_knowledge: true,
                use_macrocycle_torsions: true,
                use_macrocycle14config: true,
                ..EmbedParameters::default()
            },
        }
    }
}

/// Embeds the molecule with the atoms that match the core at the core's
/// coordinates, cleans it up with UFF, and aligns it to the core, setting
/// the alignment's RMS as its `EmbedRMS` property, written as Python writes
/// a `float`. A molecule that does not match the core, or that does not
/// embed, is `Error::Value`, as in Python. Python returns the molecule it
/// changed. Mirrors `AllChem.ConstrainedEmbed(mol, core)`.
pub fn constrained_embed(mol: &mut Mol, core: &Mol) -> Result<()> {
    constrained_embed_with_params(mol, core, &ConstrainedEmbedParams::default())
}

/// As [`constrained_embed`]. Mirrors `AllChem.ConstrainedEmbed(mol, core,
/// useTethers, coreConfId, randomseed, getForceField, **kwargs)`.
pub fn constrained_embed_with_params(
    mol: &mut Mol,
    core: &Mol,
    params: &ConstrainedEmbedParams,
) -> Result<()> {
    let match_ = mol.get_substruct_match(core)?;
    if match_.is_empty() {
        return Err(Error::Value {
            message: "molecule doesn't match the core".to_owned(),
        });
    }
    let core_conf = core.get_conformer_with_id(params.core_conf_id)?;
    let mut coord_map = std::collections::BTreeMap::new();
    for (i, &idx_i) in match_.iter().enumerate() {
        coord_map.insert(idx_i, core_conf.get_atom_position(i as u32)?);
    }
    let mut embed_params = EmbedParameters {
        coord_map: coord_map.clone(),
        random_seed: params.randomseed,
        ..params.embed_params.clone()
    };
    let ci = embed_molecule_with_params(mol, &mut embed_params)?;
    if ci < 0 {
        return Err(Error::Value {
            message: "Could not embed molecule.".to_owned(),
        });
    }
    let alg_map = GetAlignmentTransformParams {
        atom_map: match_
            .iter()
            .enumerate()
            .map(|(i, &j)| (j, i as u32))
            .collect(),
        ..Default::default()
    };
    let rms = if !params.use_tethers {
        {
            let mut ff = (params.get_force_field)(mol)?;
            for (i, &idx_i) in match_.iter().enumerate() {
                for &idx_j in &match_[i + 1..] {
                    let d = coord_map[&idx_i].distance(&coord_map[&idx_j]);
                    ff.add_distance_constraint(idx_i, idx_j, d, d, 100.0)?;
                }
            }
            ff.initialize()?;
            let mut n = 4;
            let mut more = ff.minimize()?;
            while more != 0 && n != 0 {
                more = ff.minimize()?;
                n -= 1;
            }
        }
        align_mol_with_params(mol, core, &alg_map)?
    } else {
        align_mol_with_params(mol, core, &alg_map)?;
        {
            let mut ff = (params.get_force_field)(mol)?;
            // Python tethers to the core's default conformer whatever
            // `coreConfId` is.
            let conf = core.get_conformer()?;
            for i in 0..core.get_num_atoms() {
                let p = conf.get_atom_position(i)?;
                let p_idx = ff.add_extra_point_with_fixed(p.x, p.y, p.z, true) - 1;
                let target = *match_
                    .get(i as usize)
                    .ok_or(Error::index(i, match_.len()))?;
                ff.add_distance_constraint(p_idx, target, 0.0, 0.0, 100.0)?;
            }
            ff.initialize()?;
            let minimize = MinimizeParams {
                energy_tol: 1e-4,
                force_tol: 1e-3,
                ..Default::default()
            };
            let mut n = 4;
            let mut more = ff.minimize_with_params(&minimize)?;
            while more != 0 && n != 0 {
                more = ff.minimize_with_params(&minimize)?;
                n -= 1;
            }
        }
        align_mol_with_params(mol, core, &alg_map)?
    };
    mol.set_prop("EmbedRMS", &float_repr(rms));
    Ok(())
}

/// A copy of the molecule with the bond orders, aromaticity, charges, and
/// explicit hydrogens of a template that matches it with every bond single
/// and every charge zero, sanitized; the molecule unchanged when the
/// template matches it as it is. The template should have no explicit
/// hydrogens. A template that matches neither way is `Error::Value`, as in
/// Python; when it matches several ways, the first is used, with a warning
/// to RDKit's log, as in Python. Mirrors
/// `AllChem.AssignBondOrdersFromTemplate`.
pub fn assign_bond_orders_from_template(refmol: &Mol, mol: &Mol) -> Result<Mol> {
    let mut refmol2 = refmol.clone();
    let mut mol2 = mol.clone();
    if !mol2.get_substruct_match(&refmol2)?.is_empty() {
        return Ok(mol2);
    }
    for idx in 0..mol2.get_num_bonds() {
        let mut b = mol2.get_bond_with_idx_mut(idx)?;
        if b.get_bond_type() != BondType::Single {
            b.set_bond_type(BondType::Single);
            b.set_is_aromatic(false);
        }
    }
    for idx in 0..refmol2.get_num_bonds() {
        let mut b = refmol2.get_bond_with_idx_mut(idx)?;
        b.set_bond_type(BondType::Single);
        b.set_is_aromatic(false);
    }
    for idx in 0..refmol2.get_num_atoms() {
        refmol2.get_atom_with_idx_mut(idx)?.set_formal_charge(0);
    }
    for idx in 0..mol2.get_num_atoms() {
        mol2.get_atom_with_idx_mut(idx)?.set_formal_charge(0);
    }
    let matches = mol2.get_substruct_matches_with_params(
        &refmol2,
        &SubstructMatchParameters {
            uniquify: false,
            ..Default::default()
        },
    )?;
    let Some(matching) = matches.first() else {
        return Err(Error::Value {
            message: "No matching found".to_owned(),
        });
    };
    if matches.len() > 1 {
        // Python logs through RDLogger, which writes to RDKit's log.
        crate::rd_base::log_message(
            "rdApp.warning",
            "WARNING: More than one matching pattern found - picking one\n",
        );
    }
    let at = |i: u32| -> Result<u32> {
        matching
            .get(i as usize)
            .copied()
            .ok_or(Error::index(i, matching.len()))
    };
    for b in refmol.get_bonds() {
        let (atom1, atom2) = (at(b.get_begin_atom_idx())?, at(b.get_end_atom_idx())?);
        let idx = mol2
            .get_bond_between_atoms(atom1, atom2)?
            .map(|b2| b2.get_idx())
            .ok_or_else(|| Error::Value {
                message: format!("no bond between atoms {atom1} and {atom2}"),
            })?;
        let mut b2 = mol2.get_bond_with_idx_mut(idx)?;
        b2.set_bond_type(b.get_bond_type());
        b2.set_is_aromatic(b.get_is_aromatic());
    }
    for a in refmol.get_atoms() {
        let mut a2 = mol2.get_atom_with_idx_mut(at(a.get_idx())?)?;
        a2.set_hybridization(a.get_hybridization());
        a2.set_is_aromatic(a.get_is_aromatic());
        a2.set_num_explicit_hs(a.get_num_explicit_hs());
        a2.set_formal_charge(a.get_formal_charge());
    }
    sanitize_mol(&mut mol2)?;
    Ok(mol2)
}
