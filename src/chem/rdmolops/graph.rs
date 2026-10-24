// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdmolops's graph matrices, paths, subgraphs, and subsets.

use std::collections::BTreeMap;

use super::mol_list_into_map;
use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::expect_ok;
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::graph::ffi;

/// Splits a row-major matrix into rows of `n`.
fn rows<T: Copy>(values: &[T], n: usize) -> Vec<Vec<T>> {
    if n == 0 {
        return Vec::new();
    }
    values.chunks(n).map(<[T]>::to_vec).collect()
}

fn to_map(pairs: Vec<IndexPair>) -> BTreeMap<u32, u32> {
    pairs.into_iter().map(|p| (p.first, p.second)).collect()
}

/// Parameters for `get_distance_matrix_with_params`, with Python's defaults
/// (A2). Python's `force` and `prefix` control RDKit's matrix cache, which
/// the crate never keeps (CONVENTIONS.md U6).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetDistanceMatrixParams {
    /// Weights each bond by the inverse of its order, 2/3 for aromatic bonds.
    pub use_bo: bool,
    /// Puts 6 / atomic number on the diagonal.
    pub use_atom_wts: bool,
}

/// The topological distance between each pair of atoms, as rows. Mirrors
/// `rdmolops.GetDistanceMatrix`. `Error::Value` above 46340 atoms, where
/// RDKit's `int` arithmetic overflows (U3).
pub fn get_distance_matrix(mol: &Mol) -> Result<Vec<Vec<f64>>> {
    get_distance_matrix_with_params(mol, &GetDistanceMatrixParams::default())
}

/// Mirrors `rdmolops.GetDistanceMatrix(mol, useBO, useAtomWts)`.
pub fn get_distance_matrix_with_params(
    mol: &Mol,
    params: &GetDistanceMatrixParams,
) -> Result<Vec<Vec<f64>>> {
    let values = ffi::get_distance_matrix(mol.as_ffi(), params.use_bo, params.use_atom_wts)?;
    Ok(rows(&values, mol.get_num_atoms() as usize))
}

/// Parameters for `get_3d_distance_matrix_with_params`, with Python's
/// defaults (A2). Python's `force` and `prefix` control RDKit's matrix
/// cache, which the crate never keeps (U6).
#[derive(Clone, Debug, PartialEq)]
pub struct Get3DDistanceMatrixParams {
    /// The conformer; −1 is the default conformer.
    pub conf_id: i32,
    /// Puts 6 / atomic number on the diagonal.
    pub use_atom_wts: bool,
}

impl Default for Get3DDistanceMatrixParams {
    fn default() -> Get3DDistanceMatrixParams {
        Get3DDistanceMatrixParams {
            conf_id: -1,
            use_atom_wts: false,
        }
    }
}

/// The Euclidean distance between each pair of atoms in a conformer, as
/// rows. Mirrors `rdmolops.Get3DDistanceMatrix`. `Error::Value` for a
/// missing conformer and above 46340 atoms (U3).
pub fn get_3d_distance_matrix(mol: &Mol) -> Result<Vec<Vec<f64>>> {
    get_3d_distance_matrix_with_params(mol, &Get3DDistanceMatrixParams::default())
}

/// Mirrors `rdmolops.Get3DDistanceMatrix(mol, confId, useAtomWts)`.
pub fn get_3d_distance_matrix_with_params(
    mol: &Mol,
    params: &Get3DDistanceMatrixParams,
) -> Result<Vec<Vec<f64>>> {
    let values = ffi::get_3d_distance_matrix(mol.as_ffi(), params.conf_id, params.use_atom_wts)?;
    Ok(rows(&values, mol.get_num_atoms() as usize))
}

/// 1 for each bonded pair of atoms and 0 elsewhere, as rows. Mirrors
/// `rdmolops.GetAdjacencyMatrix`. `Error::Value` above 46340 atoms (U3).
pub fn get_adjacency_matrix(mol: &Mol) -> Result<Vec<Vec<i32>>> {
    get_adjacency_matrix_with_empty_val(mol, 0)
}

/// As `get_adjacency_matrix`, with `empty_val` for each unbonded pair, the
/// diagonal included. Mirrors `rdmolops.GetAdjacencyMatrix(mol,
/// emptyVal=...)`, whose `memset` keeps only the low byte of `emptyVal`;
/// here every unbonded entry is `empty_val`.
pub fn get_adjacency_matrix_with_empty_val(mol: &Mol, empty_val: i32) -> Result<Vec<Vec<i32>>> {
    let values = ffi::get_adjacency_matrix(mol.as_ffi(), false, empty_val)?;
    // Each entry is 1.0 or empty_val, so the conversion is exact.
    let values: Vec<i32> = values.into_iter().map(|v| v as i32).collect();
    Ok(rows(&values, mol.get_num_atoms() as usize))
}

/// Each bonded pair's bond order (its valence contribution) and 0
/// elsewhere, as rows. Mirrors `rdmolops.GetAdjacencyMatrix(mol,
/// useBO=True)` (A4).
pub fn get_adjacency_matrix_use_bo(mol: &Mol) -> Result<Vec<Vec<f64>>> {
    get_adjacency_matrix_use_bo_with_empty_val(mol, 0)
}

/// As `get_adjacency_matrix_use_bo`, with `empty_val` for each unbonded
/// pair, the diagonal included. Mirrors `rdmolops.GetAdjacencyMatrix(mol,
/// useBO=True, emptyVal=...)`.
pub fn get_adjacency_matrix_use_bo_with_empty_val(
    mol: &Mol,
    empty_val: i32,
) -> Result<Vec<Vec<f64>>> {
    let values = ffi::get_adjacency_matrix(mol.as_ffi(), true, empty_val)?;
    Ok(rows(&values, mol.get_num_atoms() as usize))
}

/// Parameters for `find_all_paths_of_length_n_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindAllPathsOfLengthNParams {
    /// Counts and returns bonds; otherwise atoms.
    pub use_bonds: bool,
    /// Includes hydrogens.
    pub use_hs: bool,
    /// Only paths that start at this atom; −1 for all.
    pub rooted_at_atom: i32,
    /// Only paths that are shortest paths between their ends.
    pub only_shortest_paths: bool,
}

impl Default for FindAllPathsOfLengthNParams {
    fn default() -> FindAllPathsOfLengthNParams {
        FindAllPathsOfLengthNParams {
            use_bonds: true,
            use_hs: false,
            rooted_at_atom: -1,
            only_shortest_paths: false,
        }
    }
}

/// Every linear path of `length` bonds, as bond indices. Mirrors
/// `rdmolops.FindAllPathsOfLengthN`. `Error::Value` above 46340 atoms (U3).
pub fn find_all_paths_of_length_n(mol: &Mol, length: u32) -> Result<Vec<Vec<u32>>> {
    find_all_paths_of_length_n_with_params(mol, length, &FindAllPathsOfLengthNParams::default())
}

/// Mirrors `rdmolops.FindAllPathsOfLengthN(mol, length, useBonds, useHs,
/// rootedAtAtom, onlyShortestPaths)`.
pub fn find_all_paths_of_length_n_with_params(
    mol: &Mol,
    length: u32,
    params: &FindAllPathsOfLengthNParams,
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(ffi::find_all_paths_of_length_n(
        mol.as_ffi(),
        length,
        params.use_bonds,
        params.use_hs,
        params.rooted_at_atom,
        params.only_shortest_paths,
    )?))
}

/// Parameters for `find_all_subgraphs_of_length_n_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindAllSubgraphsOfLengthNParams {
    /// Includes hydrogens.
    pub use_hs: bool,
    /// Only subgraphs that contain this atom; −1 for all.
    pub rooted_at_atom: i32,
}

impl Default for FindAllSubgraphsOfLengthNParams {
    fn default() -> FindAllSubgraphsOfLengthNParams {
        FindAllSubgraphsOfLengthNParams {
            use_hs: false,
            rooted_at_atom: -1,
        }
    }
}

/// Every connected subgraph of `length` bonds, as bond indices. Mirrors
/// `rdmolops.FindAllSubgraphsOfLengthN`.
pub fn find_all_subgraphs_of_length_n(mol: &Mol, length: u32) -> Result<Vec<Vec<u32>>> {
    find_all_subgraphs_of_length_n_with_params(
        mol,
        length,
        &FindAllSubgraphsOfLengthNParams::default(),
    )
}

/// Mirrors `rdmolops.FindAllSubgraphsOfLengthN(mol, length, useHs,
/// rootedAtAtom)`.
pub fn find_all_subgraphs_of_length_n_with_params(
    mol: &Mol,
    length: u32,
    params: &FindAllSubgraphsOfLengthNParams,
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(ffi::find_all_subgraphs_of_length_n(
        mol.as_ffi(),
        length,
        params.use_hs,
        params.rooted_at_atom,
    )?))
}

/// Parameters for `find_all_subgraphs_of_length_m_to_n_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindAllSubgraphsOfLengthMToNParams {
    /// Includes hydrogens.
    pub use_hs: bool,
    /// Only subgraphs that contain this atom; −1 for all.
    pub rooted_at_atom: i32,
}

impl Default for FindAllSubgraphsOfLengthMToNParams {
    fn default() -> FindAllSubgraphsOfLengthMToNParams {
        FindAllSubgraphsOfLengthMToNParams {
            use_hs: false,
            rooted_at_atom: -1,
        }
    }
}

/// The connected subgraphs of each length from `min` to `max` bonds: entry
/// `i` holds those of `min + i` bonds. Mirrors
/// `rdmolops.FindAllSubgraphsOfLengthMToN`. `Error::Value` when `min >
/// max`.
pub fn find_all_subgraphs_of_length_m_to_n(
    mol: &Mol,
    min: u32,
    max: u32,
) -> Result<Vec<Vec<Vec<u32>>>> {
    find_all_subgraphs_of_length_m_to_n_with_params(
        mol,
        min,
        max,
        &FindAllSubgraphsOfLengthMToNParams::default(),
    )
}

/// Mirrors `rdmolops.FindAllSubgraphsOfLengthMToN(mol, min, max, useHs,
/// rootedAtAtom)`.
pub fn find_all_subgraphs_of_length_m_to_n_with_params(
    mol: &Mol,
    min: u32,
    max: u32,
    params: &FindAllSubgraphsOfLengthMToNParams,
) -> Result<Vec<Vec<Vec<u32>>>> {
    let mut counts = Vec::new();
    let lists = ffi::find_all_subgraphs_of_length_m_to_n(
        mol.as_ffi(),
        min,
        max,
        params.use_hs,
        params.rooted_at_atom,
        &mut counts,
    )?;
    let mut lists = IndexList::into_nested(lists).into_iter();
    // Lengths past the number of bonds have no subgraphs, and C++ reports
    // none for them.
    let mut groups: Vec<Vec<Vec<u32>>> = counts
        .into_iter()
        .map(|n| lists.by_ref().take(n as usize).collect())
        .collect();
    groups.resize((max - min) as usize + 1, Vec::new());
    Ok(groups)
}

/// Parameters for `find_unique_subgraphs_of_length_n_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindUniqueSubgraphsOfLengthNParams {
    /// Includes hydrogens.
    pub use_hs: bool,
    /// Uses bond orders to tell subgraphs apart.
    pub use_bo: bool,
    /// Only subgraphs that contain this atom; −1 for all.
    pub rooted_at_atom: i32,
}

impl Default for FindUniqueSubgraphsOfLengthNParams {
    fn default() -> FindUniqueSubgraphsOfLengthNParams {
        FindUniqueSubgraphsOfLengthNParams {
            use_hs: false,
            use_bo: true,
            rooted_at_atom: -1,
        }
    }
}

/// One connected subgraph of `length` bonds for each distinct
/// substructure, as bond indices. Mirrors
/// `rdmolops.FindUniqueSubgraphsOfLengthN`.
pub fn find_unique_subgraphs_of_length_n(mol: &Mol, length: u32) -> Result<Vec<Vec<u32>>> {
    find_unique_subgraphs_of_length_n_with_params(
        mol,
        length,
        &FindUniqueSubgraphsOfLengthNParams::default(),
    )
}

/// Mirrors `rdmolops.FindUniqueSubgraphsOfLengthN(mol, length, useHs, useBO,
/// rootedAtAtom)`.
pub fn find_unique_subgraphs_of_length_n_with_params(
    mol: &Mol,
    length: u32,
    params: &FindUniqueSubgraphsOfLengthNParams,
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(
        ffi::find_unique_subgraphs_of_length_n(
            mol.as_ffi(),
            length,
            params.use_hs,
            params.use_bo,
            params.rooted_at_atom,
        )?,
    ))
}

/// Parameters for `find_atom_environment_of_radius_n_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindAtomEnvironmentOfRadiusNParams {
    /// Includes hydrogens.
    pub use_hs: bool,
    /// Returns nothing unless the environment reaches the full radius.
    pub enforce_size: bool,
}

impl Default for FindAtomEnvironmentOfRadiusNParams {
    fn default() -> FindAtomEnvironmentOfRadiusNParams {
        FindAtomEnvironmentOfRadiusNParams {
            use_hs: false,
            enforce_size: true,
        }
    }
}

/// What Python's `atomMap` output argument receives (A9).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FindAtomEnvironmentOfRadiusNOutputs {
    /// Each atom in the environment and its distance from the root.
    pub atom_map: BTreeMap<u32, u32>,
}

/// The bonds within `radius` bonds of `rooted_at_atom`. Mirrors
/// `rdmolops.FindAtomEnvironmentOfRadiusN`. `Error::Value` for an atom
/// index past the end.
pub fn find_atom_environment_of_radius_n(
    mol: &Mol,
    radius: u32,
    rooted_at_atom: u32,
) -> Result<Vec<u32>> {
    find_atom_environment_of_radius_n_with_params(
        mol,
        radius,
        rooted_at_atom,
        &FindAtomEnvironmentOfRadiusNParams::default(),
    )
}

/// Mirrors `rdmolops.FindAtomEnvironmentOfRadiusN(mol, radius, rootedAtAtom,
/// useHs, enforceSize)`.
pub fn find_atom_environment_of_radius_n_with_params(
    mol: &Mol,
    radius: u32,
    rooted_at_atom: u32,
    params: &FindAtomEnvironmentOfRadiusNParams,
) -> Result<Vec<u32>> {
    Ok(find_atom_environment_of_radius_n_with_outputs(mol, radius, rooted_at_atom, params)?.0)
}

/// The environment and Python's `atomMap` (A9). Mirrors
/// `rdmolops.FindAtomEnvironmentOfRadiusN(..., atomMap)`.
pub fn find_atom_environment_of_radius_n_with_outputs(
    mol: &Mol,
    radius: u32,
    rooted_at_atom: u32,
    params: &FindAtomEnvironmentOfRadiusNParams,
) -> Result<(Vec<u32>, FindAtomEnvironmentOfRadiusNOutputs)> {
    let mut atom_map = Vec::new();
    let path = ffi::find_atom_environment_of_radius_n(
        mol.as_ffi(),
        radius,
        rooted_at_atom,
        params.use_hs,
        params.enforce_size,
        &mut atom_map,
    )?;
    Ok((
        path,
        FindAtomEnvironmentOfRadiusNOutputs {
            atom_map: to_map(atom_map),
        },
    ))
}

/// What Python's `atomMap` output argument receives (A9).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PathToSubmolOutputs {
    /// Each atom's index in the molecule and its index in the submolecule.
    pub atom_map: BTreeMap<u32, u32>,
}

/// The submolecule the bonds and their atoms make. Mirrors
/// `rdmolops.PathToSubmol`; bond indices past the end are skipped.
pub fn path_to_submol(mol: &Mol, path: &[u32]) -> Result<Mol> {
    path_to_submol_with_use_query(mol, path, false)
}

/// Mirrors `rdmolops.PathToSubmol(mol, path, useQuery)`; with `use_query`,
/// the submolecule copies the atoms and bonds as queries.
pub fn path_to_submol_with_use_query(mol: &Mol, path: &[u32], use_query: bool) -> Result<Mol> {
    Ok(path_to_submol_with_outputs(mol, path, use_query)?.0)
}

/// The submolecule and Python's `atomMap` (A9). Mirrors
/// `rdmolops.PathToSubmol(mol, path, useQuery, atomMap)`.
pub fn path_to_submol_with_outputs(
    mol: &Mol,
    path: &[u32],
    use_query: bool,
) -> Result<(Mol, PathToSubmolOutputs)> {
    let mut atom_map = Vec::new();
    let submol = ffi::path_to_submol(mol.as_ffi(), path, use_query, &mut atom_map)?;
    Ok((
        Mol::from_ffi(submol),
        PathToSubmolOutputs {
            atom_map: to_map(atom_map),
        },
    ))
}

/// Parameters for `split_mol_by_pdb_residues_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SplitMolByPDBResiduesParams {
    /// Only these residue names, or with `negate_list` all but these; `None`
    /// for all.
    pub white_list: Option<Vec<String>>,
    /// Inverts `white_list`.
    pub negate_list: bool,
}

/// A molecule for each PDB residue name, with the atoms that have it; atoms
/// without PDB residue information go under `""`. Mirrors
/// `rdmolops.SplitMolByPDBResidues`.
pub fn split_mol_by_pdb_residues(mol: &Mol) -> Result<BTreeMap<String, Mol>> {
    split_mol_by_pdb_residues_with_params(mol, &SplitMolByPDBResiduesParams::default())
}

/// Mirrors `rdmolops.SplitMolByPDBResidues(mol, whiteList, negateList)`.
pub fn split_mol_by_pdb_residues_with_params(
    mol: &Mol,
    params: &SplitMolByPDBResiduesParams,
) -> Result<BTreeMap<String, Mol>> {
    let white_list = params.white_list.as_deref();
    mol_list_into_map(ffi::split_mol_by_pdb_residues(
        mol.as_ffi(),
        white_list.is_some(),
        white_list.unwrap_or_default(),
        params.negate_list,
    )?)
}

/// Parameters for `split_mol_by_pdb_chain_id_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SplitMolByPDBChainIdParams {
    /// Only these chain ids, or with `negate_list` all but these; `None` for
    /// all.
    pub white_list: Option<Vec<String>>,
    /// Inverts `white_list`.
    pub negate_list: bool,
}

/// A molecule for each PDB chain id, with the atoms that have it; atoms
/// without PDB residue information go under `""`. Mirrors
/// `rdmolops.SplitMolByPDBChainId`.
pub fn split_mol_by_pdb_chain_id(mol: &Mol) -> Result<BTreeMap<String, Mol>> {
    split_mol_by_pdb_chain_id_with_params(mol, &SplitMolByPDBChainIdParams::default())
}

/// Mirrors `rdmolops.SplitMolByPDBChainId(mol, whiteList, negateList)`.
pub fn split_mol_by_pdb_chain_id_with_params(
    mol: &Mol,
    params: &SplitMolByPDBChainIdParams,
) -> Result<BTreeMap<String, Mol>> {
    let white_list = params.white_list.as_deref();
    mol_list_into_map(ffi::split_mol_by_pdb_chain_id(
        mol.as_ffi(),
        white_list.is_some(),
        white_list.unwrap_or_default(),
        params.negate_list,
    )?)
}

/// A copy with its atoms reordered: atom `i` of the copy is atom
/// `new_order[i]` of the molecule. Mirrors `rdmolops.RenumberAtoms`.
/// `Error::Value` unless `new_order` is a permutation of the atom indices.
pub fn renumber_atoms(mol: &Mol, new_order: &[u32]) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::renumber_atoms(mol.as_ffi(), new_order)?))
}

rdkit_enum! {
    /// How `copy_mol_subset` reads its indices. Mirrors
    /// `rdmolops.SubsetMethod`.
    SubsetMethod: ffi::SubsetMethod {
        BondsBetweenAtoms = BONDS_BETWEEN_ATOMS,
        Bonds = BONDS,
    }
}

/// What `copy_mol_subset` copies. Mirrors `rdmolops.SubsetOptions` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct SubsetOptions {
    /// Sanitizes the subset.
    pub sanitize: bool,
    /// Clears the subset's computed properties.
    pub clear_computed_props: bool,
    /// Copies the atoms and bonds as queries.
    pub copy_as_query: bool,
    /// Copies the conformers' coordinates.
    pub copy_coordinates: bool,
    /// The conformer to copy; −1 for all. Python's value is unsigned, so its
    /// default reads 4294967295.
    pub conformer_idx: i32,
    /// How the indices are read.
    pub method: SubsetMethod,
}

impl SubsetOptions {
    fn to_ffi(&self) -> ffi::SubsetOptionsData {
        ffi::SubsetOptionsData {
            sanitize: self.sanitize,
            clear_computed_props: self.clear_computed_props,
            copy_as_query: self.copy_as_query,
            copy_coordinates: self.copy_coordinates,
            conformer_idx: self.conformer_idx,
            method: self.method.to_ffi(),
        }
    }
}

impl Default for SubsetOptions {
    fn default() -> SubsetOptions {
        let data = expect_ok(ffi::subset_options_default());
        SubsetOptions {
            sanitize: data.sanitize,
            clear_computed_props: data.clear_computed_props,
            copy_as_query: data.copy_as_query,
            copy_coordinates: data.copy_coordinates,
            conformer_idx: data.conformer_idx,
            method: SubsetMethod::from_ffi(data.method),
        }
    }
}

/// Where `copy_mol_subset` put each atom and bond. Mirrors
/// `rdmolops.SubsetInfo`, which Python passes as an output argument (A9).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SubsetInfo {
    /// Each copied atom's index in the molecule and in the subset.
    pub atom_mapping: BTreeMap<u32, u32>,
    /// Each copied bond's index in the molecule and in the subset.
    pub bond_mapping: BTreeMap<u32, u32>,
}

impl SubsetInfo {
    fn ffi() -> ffi::SubsetInfoData {
        ffi::SubsetInfoData {
            atom_mapping: Vec::new(),
            bond_mapping: Vec::new(),
        }
    }

    fn from_ffi(data: ffi::SubsetInfoData) -> SubsetInfo {
        SubsetInfo {
            atom_mapping: to_map(data.atom_mapping),
            bond_mapping: to_map(data.bond_mapping),
        }
    }
}

/// A copy of the atoms (or with `SubsetMethod::Bonds`, the bonds) in `path`
/// and what connects them. Mirrors `rdmolops.CopyMolSubset(mol, path)`;
/// indices past the end are skipped.
pub fn copy_mol_subset(mol: &Mol, path: &[u32]) -> Result<Mol> {
    copy_mol_subset_with_params(mol, path, &SubsetOptions::default())
}

/// Mirrors `rdmolops.CopyMolSubset(mol, path, options)` (A3).
pub fn copy_mol_subset_with_params(
    mol: &Mol,
    path: &[u32],
    options: &SubsetOptions,
) -> Result<Mol> {
    Ok(copy_mol_subset_with_outputs(mol, path, options)?.0)
}

/// The subset and its `SubsetInfo` (A9). Mirrors
/// `rdmolops.CopyMolSubset(mol, path, subsetInfo, options)`.
pub fn copy_mol_subset_with_outputs(
    mol: &Mol,
    path: &[u32],
    options: &SubsetOptions,
) -> Result<(Mol, SubsetInfo)> {
    let mut info = SubsetInfo::ffi();
    let subset = ffi::copy_mol_subset(mol.as_ffi(), path, &options.to_ffi(), &mut info)?;
    Ok((Mol::from_ffi(subset), SubsetInfo::from_ffi(info)))
}

/// A copy of exactly these atoms and bonds; a bond is kept only with both
/// of its atoms. Mirrors `rdmolops.CopyMolSubset(mol, atomIndices,
/// bondIndices)` (A5). `Error::Index` for an index past the end.
pub fn copy_mol_subset_atoms_bonds(
    mol: &Mol,
    atom_indices: &[u32],
    bond_indices: &[u32],
) -> Result<Mol> {
    copy_mol_subset_atoms_bonds_with_params(
        mol,
        atom_indices,
        bond_indices,
        &SubsetOptions::default(),
    )
}

/// Mirrors `rdmolops.CopyMolSubset(mol, atomIndices, bondIndices, options)`
/// (A3).
pub fn copy_mol_subset_atoms_bonds_with_params(
    mol: &Mol,
    atom_indices: &[u32],
    bond_indices: &[u32],
    options: &SubsetOptions,
) -> Result<Mol> {
    Ok(copy_mol_subset_atoms_bonds_with_outputs(mol, atom_indices, bond_indices, options)?.0)
}

/// The subset and its `SubsetInfo` (A9). Mirrors
/// `rdmolops.CopyMolSubset(mol, atomIndices, bondIndices, subsetInfo,
/// options)`.
pub fn copy_mol_subset_atoms_bonds_with_outputs(
    mol: &Mol,
    atom_indices: &[u32],
    bond_indices: &[u32],
    options: &SubsetOptions,
) -> Result<(Mol, SubsetInfo)> {
    let mut info = SubsetInfo::ffi();
    let subset = ffi::copy_mol_subset_atoms_bonds(
        mol.as_ffi(),
        atom_indices,
        bond_indices,
        &options.to_ffi(),
        &mut info,
    )?;
    Ok((Mol::from_ffi(subset), SubsetInfo::from_ffi(info)))
}
