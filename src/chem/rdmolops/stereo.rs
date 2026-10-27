// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Stereochemistry perception, wedging, and enhanced stereo.

use crate::chem::rdmolops::SanitizeFlags;
use crate::chem::{BondMut, Conformer, Mol, StereoGroupType};
use crate::error::{decode, expect_ok};
use crate::ffi::chirality::ffi;
use crate::{Error, Result};

/// Parameters for `assign_stereochemistry_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AssignStereochemistryParams {
    /// Removes chirality and bond stereo that can't hold.
    pub clean_it: bool,
    /// Assigns again even if the molecule says it has been.
    pub force: bool,
    /// Marks atoms that could be stereocenters with `_ChiralityPossible`.
    pub flag_possible_stereo_centers: bool,
}

/// Assigns CIP labels and bond stereo from the chiral tags and bond
/// directions. Mirrors `rdmolops.AssignStereochemistry`.
pub fn assign_stereochemistry(mol: &mut Mol) -> Result<()> {
    assign_stereochemistry_with_params(mol, &AssignStereochemistryParams::default())
}

/// Mirrors `rdmolops.AssignStereochemistry(mol, cleanIt, force,
/// flagPossibleStereoCenters)`.
pub fn assign_stereochemistry_with_params(
    mol: &mut Mol,
    params: &AssignStereochemistryParams,
) -> Result<()> {
    Ok(ffi::assign_stereochemistry(
        mol.pin_mut(),
        params.clean_it,
        params.force,
        params.flag_possible_stereo_centers,
    )?)
}

/// Parameters for the functions that read stereo from a conformer, with
/// Python's defaults (A2).
macro_rules! conformer_stereo_params {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq)]
        pub struct $name {
            /// The conformer; -1 for the first.
            pub conf_id: i32,
            /// Replaces the stereo the molecule already has.
            pub replace_existing_tags: bool,
        }

        impl Default for $name {
            fn default() -> $name {
                $name {
                    conf_id: -1,
                    replace_existing_tags: true,
                }
            }
        }
    };
}

conformer_stereo_params!(
    /// Parameters for `assign_stereochemistry_from_3d_with_params`.
    AssignStereochemistryFrom3DParams
);
conformer_stereo_params!(
    /// Parameters for `assign_chiral_types_from_bond_dirs_with_params`.
    AssignChiralTypesFromBondDirsParams
);
conformer_stereo_params!(
    /// Parameters for `assign_atom_chiral_tags_from_structure_with_params`.
    AssignAtomChiralTagsFromStructureParams
);

/// Assigns chiral tags, CIP labels, and bond stereo from a 3D conformer; a
/// 2D one changes nothing. Mirrors `rdmolops.AssignStereochemistryFrom3D`.
pub fn assign_stereochemistry_from_3d(mol: &mut Mol) -> Result<()> {
    assign_stereochemistry_from_3d_with_params(mol, &Default::default())
}

/// Mirrors `rdmolops.AssignStereochemistryFrom3D(mol, confId,
/// replaceExistingTags)`.
pub fn assign_stereochemistry_from_3d_with_params(
    mol: &mut Mol,
    params: &AssignStereochemistryFrom3DParams,
) -> Result<()> {
    Ok(ffi::assign_stereochemistry_from_3d(
        mol.pin_mut(),
        params.conf_id,
        params.replace_existing_tags,
    )?)
}

/// Assigns chiral tags from the wedged bonds and a conformer's 2D
/// coordinates. Mirrors `rdmolops.AssignChiralTypesFromBondDirs`.
pub fn assign_chiral_types_from_bond_dirs(mol: &mut Mol) -> Result<()> {
    assign_chiral_types_from_bond_dirs_with_params(mol, &Default::default())
}

/// Mirrors `rdmolops.AssignChiralTypesFromBondDirs(mol, confId,
/// replaceExistingTags)`.
pub fn assign_chiral_types_from_bond_dirs_with_params(
    mol: &mut Mol,
    params: &AssignChiralTypesFromBondDirsParams,
) -> Result<()> {
    Ok(ffi::assign_chiral_types_from_bond_dirs(
        mol.pin_mut(),
        params.conf_id,
        params.replace_existing_tags,
    )?)
}

/// Assigns chiral tags from a 3D conformer, to every atom with four
/// neighbors whether or not it is a stereocenter. Mirrors
/// `rdmolops.AssignAtomChiralTagsFromStructure`.
pub fn assign_atom_chiral_tags_from_structure(mol: &mut Mol) -> Result<()> {
    assign_atom_chiral_tags_from_structure_with_params(mol, &Default::default())
}

/// Mirrors `rdmolops.AssignAtomChiralTagsFromStructure(mol, confId,
/// replaceExistingTags)`.
pub fn assign_atom_chiral_tags_from_structure_with_params(
    mol: &mut Mol,
    params: &AssignAtomChiralTagsFromStructureParams,
) -> Result<()> {
    Ok(ffi::assign_atom_chiral_tags_from_structure(
        mol.pin_mut(),
        params.conf_id,
        params.replace_existing_tags,
    )?)
}

/// Assigns chiral tags from the mol-file parities in `molParity` atom
/// properties. Mirrors `rdmolops.AssignAtomChiralTagsFromMolParity`.
pub fn assign_atom_chiral_tags_from_mol_parity(mol: &mut Mol) -> Result<()> {
    assign_atom_chiral_tags_from_mol_parity_with_replace_existing_tags(mol, true)
}

/// Mirrors `rdmolops.AssignAtomChiralTagsFromMolParity(mol,
/// replaceExistingTags)`.
pub fn assign_atom_chiral_tags_from_mol_parity_with_replace_existing_tags(
    mol: &mut Mol,
    replace_existing_tags: bool,
) -> Result<()> {
    Ok(ffi::assign_atom_chiral_tags_from_mol_parity(
        mol.pin_mut(),
        replace_existing_tags,
    )?)
}

/// Marks the double bonds that could be cis or trans `STEREOANY`. Mirrors
/// `rdmolops.FindPotentialStereoBonds`.
pub fn find_potential_stereo_bonds(mol: &mut Mol) -> Result<()> {
    find_potential_stereo_bonds_with_clean_it(mol, false)
}

/// With `clean_it`, also removes stereo from double bonds that can't have
/// it. Mirrors `rdmolops.FindPotentialStereoBonds(mol, cleanIt)`.
pub fn find_potential_stereo_bonds_with_clean_it(mol: &mut Mol, clean_it: bool) -> Result<()> {
    Ok(ffi::find_potential_stereo_bonds(mol.pin_mut(), clean_it)?)
}

/// Removes all stereochemistry: chiral tags, bond stereo, and stereo
/// groups. Mirrors `rdmolops.RemoveStereochemistry`.
pub fn remove_stereochemistry(mol: &mut Mol) -> Result<()> {
    Ok(ffi::remove_stereochemistry(mol.pin_mut())?)
}

/// Removes chiral tags that came from 3D coordinates rather than from the
/// input. Mirrors `rdmolops.RemoveNonExplicit3DChirality`.
pub fn remove_non_explicit_3d_chirality(mol: &mut Mol) -> Result<()> {
    Ok(ffi::remove_non_explicit_3d_chirality(mol.pin_mut())?)
}

/// Sets the directions of single bonds next to double bonds from a
/// conformer's coordinates, which `assign_stereochemistry` then reads.
/// Mirrors `rdmolops.DetectBondStereochemistry`.
pub fn detect_bond_stereochemistry(mol: &mut Mol) -> Result<()> {
    detect_bond_stereochemistry_with_conf_id(mol, -1)
}

/// Mirrors `rdmolops.DetectBondStereochemistry(mol, confId)`.
pub fn detect_bond_stereochemistry_with_conf_id(mol: &mut Mol, conf_id: i32) -> Result<()> {
    Ok(ffi::detect_bond_stereochemistry(mol.pin_mut(), conf_id)?)
}

/// Sets the directions of single bonds next to double bonds from the bond
/// stereo. Mirrors `rdmolops.SetDoubleBondNeighborDirections`.
pub fn set_double_bond_neighbor_directions(mol: &mut Mol) -> Result<()> {
    Ok(ffi::set_double_bond_neighbor_directions(mol.pin_mut())?)
}

/// Uses `conf`'s coordinates, which need a position for every atom; a copy
/// of one of the molecule's (`Conformer::from`) will do. Mirrors
/// `rdmolops.SetDoubleBondNeighborDirections(mol, conf)`.
pub fn set_double_bond_neighbor_directions_with_conf(
    mol: &mut Mol,
    conf: &Conformer,
) -> Result<()> {
    Ok(ffi::set_double_bond_neighbor_directions_with_conf(
        mol.pin_mut(),
        conf.raw(),
    )?)
}

/// Sets bond stereo from the directions of the neighboring single bonds.
/// Mirrors `rdmolops.SetBondStereoFromDirections`.
pub fn set_bond_stereo_from_directions(mol: &mut Mol) -> Result<()> {
    Ok(ffi::set_bond_stereo_from_directions(mol.pin_mut())?)
}

/// Parameters for `find_meso_centers_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindMesoCentersParams {
    /// Distinguishes atoms by isotope.
    pub include_isotopes: bool,
    /// Distinguishes atoms by atom map number.
    pub include_atom_maps: bool,
}

impl Default for FindMesoCentersParams {
    fn default() -> FindMesoCentersParams {
        FindMesoCentersParams {
            include_isotopes: true,
            include_atom_maps: false,
        }
    }
}

/// The pairs of stereocenters that make the molecule meso. Mirrors
/// `rdmolops.FindMesoCenters`.
pub fn find_meso_centers(mol: &Mol) -> Result<Vec<(u32, u32)>> {
    find_meso_centers_with_params(mol, &FindMesoCentersParams::default())
}

/// Mirrors `rdmolops.FindMesoCenters(mol, includeIsotopes, includeAtomMaps)`.
pub fn find_meso_centers_with_params(
    mol: &Mol,
    params: &FindMesoCentersParams,
) -> Result<Vec<(u32, u32)>> {
    Ok(ffi::find_meso_centers(
        mol.as_ffi(),
        params.include_isotopes,
        params.include_atom_maps,
    )?
    .into_iter()
    .map(|pair| (pair.first, pair.second))
    .collect())
}

/// When every specified stereocenter is in one AND or OR group, notes the
/// molecule as "AND enantiomer" or "OR enantiomer" and removes the group.
/// Mirrors `rdmolops.SimplifyEnhancedStereo`.
pub fn simplify_enhanced_stereo(mol: &mut Mol) -> Result<()> {
    simplify_enhanced_stereo_with_remove_affected_stereo_groups(mol, true)
}

/// Mirrors `rdmolops.SimplifyEnhancedStereo(mol, removeAffectedStereoGroups)`.
pub fn simplify_enhanced_stereo_with_remove_affected_stereo_groups(
    mol: &mut Mol,
    remove_affected_stereo_groups: bool,
) -> Result<()> {
    Ok(ffi::simplify_enhanced_stereo(
        mol.pin_mut(),
        remove_affected_stereo_groups,
    )?)
}

/// Adds stereo groups from the mol file's chiral flag: an ABS group when it
/// is 1, and groups of `zero_flag_group_type` when it is 0. Mirrors
/// `rdmolops.TranslateChiralFlagToStereoGroups`.
pub fn translate_chiral_flag_to_stereo_groups(mol: &mut Mol) -> Result<()> {
    translate_chiral_flag_to_stereo_groups_with_zero_flag_group_type(
        mol,
        StereoGroupType::StereoAnd,
    )
}

/// Mirrors `rdmolops.TranslateChiralFlagToStereoGroups(mol,
/// zeroFlagGroupType)`.
pub fn translate_chiral_flag_to_stereo_groups_with_zero_flag_group_type(
    mol: &mut Mol,
    zero_flag_group_type: StereoGroupType,
) -> Result<()> {
    Ok(ffi::translate_chiral_flag_to_stereo_groups(
        mol.pin_mut(),
        zero_flag_group_type.to_ffi(),
    )?)
}

/// Parameters for `add_stereo_annotations_with_params`, with Python's
/// defaults (A2). `{cip}` and `{id}` stand for the CIP label and the group
/// id.
#[derive(Clone, Debug, PartialEq)]
pub struct AddStereoAnnotationsParams {
    /// The label for atoms in an ABS group.
    pub abs_label: String,
    /// The label for atoms in an OR group.
    pub or_label: String,
    /// The label for atoms in an AND group.
    pub and_label: String,
    /// The label for other stereocenters.
    pub cip_label: String,
    /// The label for stereo bonds.
    pub bond_label: String,
}

impl Default for AddStereoAnnotationsParams {
    fn default() -> AddStereoAnnotationsParams {
        AddStereoAnnotationsParams {
            abs_label: "abs ({cip})".to_owned(),
            or_label: "or{id}".to_owned(),
            and_label: "and{id}".to_owned(),
            cip_label: "({cip})".to_owned(),
            bond_label: "({cip})".to_owned(),
        }
    }
}

/// Labels atoms and bonds with their R/S, E/Z, and stereo groups, in the
/// `atomNote` and `bondNote` properties drawings show. Mirrors
/// `rdmolops.AddStereoAnnotations`.
pub fn add_stereo_annotations(mol: &mut Mol) -> Result<()> {
    add_stereo_annotations_with_params(mol, &AddStereoAnnotationsParams::default())
}

/// Mirrors `rdmolops.AddStereoAnnotations(mol, absLabel, orLabel, andLabel,
/// cipLabel, bondLabel)`.
pub fn add_stereo_annotations_with_params(
    mol: &mut Mol,
    params: &AddStereoAnnotationsParams,
) -> Result<()> {
    Ok(ffi::add_stereo_annotations(
        mol.pin_mut(),
        &params.abs_label,
        &params.or_label,
        &params.and_label,
        &params.cip_label,
        &params.bond_label,
    )?)
}

/// Thresholds for `AddWavyBondsForStereoAny`'s `add_when_impossible`: how
/// constrained a single bond may be and still be made wavy. Mirrors
/// `rdmolops.StereoBondThresholds`.
#[derive(Clone, Copy, Debug)]
pub struct StereoBondThresholds;

impl StereoBondThresholds {
    /// Next to a double bond without stereo.
    pub const DBL_BOND_NO_STEREO: u32 = 1000;
    /// Next to a double bond with stereo.
    pub const DBL_BOND_SPECIFIED_STEREO: u32 = 10000;
    /// To an atom with specified chirality.
    pub const CHIRAL_ATOM: u32 = 100000;
    /// With its direction already set.
    pub const DIRECTION_SET: u32 = 1000000;
}

/// Parameters for `add_wavy_bonds_for_stereo_any_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AddWavyBondsForStereoAnyParams {
    /// Clears the double bonds' stereo once their neighbors are wavy.
    pub clear_double_bond_flags: bool,
    /// The `StereoBondThresholds` up to which a bond is made wavy when no
    /// unconstrained one is available.
    pub add_when_impossible: u32,
}

impl Default for AddWavyBondsForStereoAnyParams {
    fn default() -> AddWavyBondsForStereoAnyParams {
        AddWavyBondsForStereoAnyParams {
            clear_double_bond_flags: true,
            add_when_impossible: StereoBondThresholds::DBL_BOND_NO_STEREO,
        }
    }
}

/// Makes a single bond next to each `STEREOANY` double bond wavy. Mirrors
/// `rdmolops.AddWavyBondsForStereoAny`.
pub fn add_wavy_bonds_for_stereo_any(mol: &mut Mol) -> Result<()> {
    add_wavy_bonds_for_stereo_any_with_params(mol, &AddWavyBondsForStereoAnyParams::default())
}

/// Mirrors `rdmolops.AddWavyBondsForStereoAny(mol, clearDoubleBondFlags,
/// addWhenImpossible)`.
pub fn add_wavy_bonds_for_stereo_any_with_params(
    mol: &mut Mol,
    params: &AddWavyBondsForStereoAnyParams,
) -> Result<()> {
    Ok(ffi::add_wavy_bonds_for_stereo_any(
        mol.pin_mut(),
        params.clear_double_bond_flags,
        params.add_when_impossible,
    )?)
}

/// Wedges or dashes a single bond from the chiral atom `from_atom_idx`, by
/// the bond's molecule's conformer `conf_id`. Python passes the conformer
/// itself, which must belong to the bond's molecule; a `BondMut` borrows
/// that molecule, so the conformer goes by id. Mirrors `rdmolops.WedgeBond`.
pub fn wedge_bond(bond: &mut BondMut<'_>, from_atom_idx: u32, conf_id: i32) -> Result<()> {
    Ok(ffi::wedge_bond(bond.raw_mut(), from_atom_idx, conf_id)?)
}

/// How `wedge_mol_bonds` wedges. Mirrors `rdmolops.BondWedgingParameters`.
#[derive(Clone, Debug, PartialEq)]
pub struct BondWedgingParameters {
    /// Wedges two bonds at a chiral center when possible: not ring bonds,
    /// not bonds to other chiral centers, and not bonds more than 120
    /// degrees apart.
    pub wedge_two_bonds_if_possible: bool,
}

impl Default for BondWedgingParameters {
    fn default() -> BondWedgingParameters {
        let data = expect_ok(ffi::bond_wedging_parameters_default());
        BondWedgingParameters {
            wedge_two_bonds_if_possible: data.wedge_two_bonds_if_possible,
        }
    }
}

/// Wedges the molecule's bonds for its chiral centers, as in mol files, by
/// `conformer`'s coordinates, one per atom; a copy of one of the molecule's
/// (`Conformer::from`) will do. Mirrors `rdmolops.WedgeMolBonds`.
pub fn wedge_mol_bonds(mol: &mut Mol, conformer: &Conformer) -> Result<()> {
    wedge_mol_bonds_with_params(mol, conformer, &BondWedgingParameters::default())
}

/// Mirrors `rdmolops.WedgeMolBonds(mol, conformer, params)`.
pub fn wedge_mol_bonds_with_params(
    mol: &mut Mol,
    conformer: &Conformer,
    params: &BondWedgingParameters,
) -> Result<()> {
    Ok(ffi::wedge_mol_bonds(
        mol.pin_mut(),
        conformer.raw(),
        &ffi::BondWedgingParametersData {
            wedge_two_bonds_if_possible: params.wedge_two_bonds_if_possible,
        },
    )?)
}

/// Replaces the bond wedging with the wedges the mol block had. Mirrors
/// `rdmolops.ReapplyMolBlockWedging`.
pub fn reapply_mol_block_wedging(mol: &mut Mol) -> Result<()> {
    reapply_mol_block_wedging_with_all_bond_types(mol, true)
}

/// Mirrors `rdmolops.ReapplyMolBlockWedging(mol, allBondTypes)`:
/// `all_bond_types` (Python's default, true) reapplies wedges to bonds other
/// than single and aromatic ones too.
pub fn reapply_mol_block_wedging_with_all_bond_types(
    mol: &mut Mol,
    all_bond_types: bool,
) -> Result<()> {
    Ok(ffi::reapply_mol_block_wedging(
        mol.pin_mut(),
        all_bond_types,
    )?)
}

/// Whether stereo perception from 3D recognizes non-tetrahedral chirality.
/// Mirrors `rdmolops.GetAllowNontetrahedralChirality`.
pub fn get_allow_nontetrahedral_chirality() -> bool {
    expect_ok(ffi::get_allow_nontetrahedral_chirality())
}

/// Mirrors `rdmolops.SetAllowNontetrahedralChirality`.
///
/// # Safety
///
/// RDKit keeps this setting in an environment variable, which it sets with
/// `setenv`, so this has [`std::env::set_var`]'s contract: no other thread
/// may read or write the environment while it runs. RDKit reads it during
/// stereo perception (CONVENTIONS.md U5).
pub unsafe fn set_allow_nontetrahedral_chirality(val: bool) {
    expect_ok(ffi::set_allow_nontetrahedral_chirality(val))
}

/// Whether stereo perception uses RDKit's legacy code. Mirrors
/// `rdmolops.GetUseLegacyStereoPerception`.
pub fn get_use_legacy_stereo_perception() -> bool {
    expect_ok(ffi::get_use_legacy_stereo_perception())
}

/// Mirrors `rdmolops.SetUseLegacyStereoPerception`.
///
/// # Safety
///
/// As for [`set_allow_nontetrahedral_chirality`]: RDKit reads this
/// environment variable during sanitization and stereo perception.
pub unsafe fn set_use_legacy_stereo_perception(val: bool) {
    expect_ok(ffi::set_use_legacy_stereo_perception(val))
}

/// Parameters for `expand_attachment_points_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandAttachmentPointsParams {
    /// Adds the dummy atoms as queries that match any atom.
    pub add_as_queries: bool,
    /// Gives the dummy atoms coordinates.
    pub add_coords: bool,
}

impl Default for ExpandAttachmentPointsParams {
    fn default() -> ExpandAttachmentPointsParams {
        ExpandAttachmentPointsParams {
            add_as_queries: true,
            add_coords: true,
        }
    }
}

/// Turns `molAttchpt` attachment points into dummy atoms. Mirrors
/// `rdmolops.ExpandAttachmentPoints`.
pub fn expand_attachment_points(mol: &mut Mol) -> Result<()> {
    expand_attachment_points_with_params(mol, &ExpandAttachmentPointsParams::default())
}

/// Mirrors `rdmolops.ExpandAttachmentPoints(mol, addAsQueries, addCoords)`.
pub fn expand_attachment_points_with_params(
    mol: &mut Mol,
    params: &ExpandAttachmentPointsParams,
) -> Result<()> {
    Ok(ffi::expand_attachment_points(
        mol.pin_mut(),
        params.add_as_queries,
        params.add_coords,
    )?)
}

/// Turns dummy atoms made by `expand_attachment_points` back into
/// attachment points. Mirrors `rdmolops.CollapseAttachmentPoints`.
pub fn collapse_attachment_points(mol: &mut Mol) -> Result<()> {
    collapse_attachment_points_with_marked_only(mol, true)
}

/// Without `marked_only`, collapses any suitable dummy atom. Mirrors
/// `rdmolops.CollapseAttachmentPoints(mol, markedOnly)`.
pub fn collapse_attachment_points_with_marked_only(mol: &mut Mol, marked_only: bool) -> Result<()> {
    Ok(ffi::collapse_attachment_points(mol.pin_mut(), marked_only)?)
}

/// The problems sanitization would find, as the errors it would return,
/// without changing the molecule. Mirrors `rdmolops.DetectChemistryProblems`.
pub fn detect_chemistry_problems(mol: &Mol) -> Result<Vec<Error>> {
    detect_chemistry_problems_with_sanitize_ops(mol, SanitizeFlags::SANITIZE_ALL)
}

/// Mirrors `rdmolops.DetectChemistryProblems(mol, sanitizeOps)`.
pub fn detect_chemistry_problems_with_sanitize_ops(
    mol: &Mol,
    sanitize_ops: SanitizeFlags,
) -> Result<Vec<Error>> {
    let problems = ffi::detect_chemistry_problems(mol.as_ffi(), sanitize_ops.bits())?;
    Ok(problems
        .iter()
        .map(|problem| decode(&problem.to_string_lossy()))
        .collect())
}
