// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdmolops's chemical transformations: substructure replacement, core and
//! sidechain replacement, fragmentation, molzip, and recursive queries.

use std::collections::BTreeMap;
use std::path::Path;

use cxx::UniquePtr;

use super::{mol_list_into_map, mol_list_into_vec, with_borrowed_mols, with_borrowed_mols_shared};
use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::BondType;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::expect_ok;
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::mol::ffi::RWMol;
use crate::ffi::transforms::ffi;
use crate::geometry::Point3D;

fn optional_mol(ptr: UniquePtr<RWMol>) -> Option<Mol> {
    (!ptr.is_null()).then(|| Mol::from_ffi(ptr))
}

fn to_pairs(pairs: &[(u32, u32)]) -> Vec<IndexPair> {
    pairs
        .iter()
        .map(|&(first, second)| IndexPair { first, second })
        .collect()
}

/// Parameters for `delete_substructs_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DeleteSubstructsParams {
    /// Deletes only matches that are whole fragments.
    pub only_frags: bool,
    /// Uses chirality in the match.
    pub use_chirality: bool,
}

/// A copy without the atoms that match the query. Mirrors
/// `rdmolops.DeleteSubstructs`.
pub fn delete_substructs(mol: &Mol, query: &Mol) -> Result<Mol> {
    delete_substructs_with_params(mol, query, &DeleteSubstructsParams::default())
}

/// Mirrors `rdmolops.DeleteSubstructs(mol, query, onlyFrags, useChirality)`.
pub fn delete_substructs_with_params(
    mol: &Mol,
    query: &Mol,
    params: &DeleteSubstructsParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::delete_substructs(
        mol.as_ffi(),
        query.as_ffi(),
        params.only_frags,
        params.use_chirality,
    )?))
}

/// Parameters for `replace_substructs_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReplaceSubstructsParams {
    /// Replaces every match in one molecule; otherwise one molecule per
    /// match.
    pub replace_all: bool,
    /// The replacement's atom that bonds to the molecule.
    pub replacement_connection_point: u32,
    /// Uses chirality in the match.
    pub use_chirality: bool,
}

/// Copies with the query's matches replaced, one per match. Mirrors
/// `rdmolops.ReplaceSubstructs`. `Error::Runtime` for a
/// `replacement_connection_point` past the replacement's end.
pub fn replace_substructs(mol: &Mol, query: &Mol, replacement: &Mol) -> Result<Vec<Mol>> {
    replace_substructs_with_params(mol, query, replacement, &ReplaceSubstructsParams::default())
}

/// Mirrors `rdmolops.ReplaceSubstructs(mol, query, replacement, replaceAll,
/// replacementConnectionPoint, useChirality)`.
pub fn replace_substructs_with_params(
    mol: &Mol,
    query: &Mol,
    replacement: &Mol,
    params: &ReplaceSubstructsParams,
) -> Result<Vec<Mol>> {
    Ok(mol_list_into_vec(ffi::replace_substructs(
        mol.as_ffi(),
        query.as_ffi(),
        replacement.as_ffi(),
        params.replace_all,
        params.replacement_connection_point,
        params.use_chirality,
    )?))
}

/// The core that matches `core_query`, with a dummy atom for each
/// sidechain, or `None` without a match. Mirrors
/// `rdmolops.ReplaceSidechains`.
pub fn replace_sidechains(mol: &Mol, core_query: &Mol) -> Result<Option<Mol>> {
    replace_sidechains_with_use_chirality(mol, core_query, false)
}

/// Mirrors `rdmolops.ReplaceSidechains(mol, coreQuery, useChirality)`.
pub fn replace_sidechains_with_use_chirality(
    mol: &Mol,
    core_query: &Mol,
    use_chirality: bool,
) -> Result<Option<Mol>> {
    Ok(optional_mol(ffi::replace_sidechains(
        mol.as_ffi(),
        core_query.as_ffi(),
        use_chirality,
    )?))
}

/// Parameters for `replace_core_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ReplaceCoreParams {
    /// Replaces the dummy atoms the core matched too.
    pub replace_dummies: bool,
    /// Labels each sidechain's dummy with the core atom's index; otherwise
    /// they are numbered in order.
    pub label_by_index: bool,
    /// Returns `None` unless each sidechain attaches at a core dummy.
    pub require_dummy_match: bool,
}

impl Default for ReplaceCoreParams {
    fn default() -> ReplaceCoreParams {
        ReplaceCoreParams {
            replace_dummies: true,
            label_by_index: false,
            require_dummy_match: false,
        }
    }
}

/// The sidechains left when the matched core is removed, each attached to
/// a dummy atom, or `None` when they cannot be. `matches` pairs core atom
/// indices with molecule atom indices. Mirrors `rdmolops.ReplaceCore(mol,
/// core, matches)`. `Error::Value` for an index past either end.
pub fn replace_core(mol: &Mol, core: &Mol, matches: &[(u32, u32)]) -> Result<Option<Mol>> {
    replace_core_with_params(mol, core, matches, &ReplaceCoreParams::default())
}

/// Mirrors `rdmolops.ReplaceCore(mol, core, matches, replaceDummies,
/// labelByIndex, requireDummyMatch)`.
pub fn replace_core_with_params(
    mol: &Mol,
    core: &Mol,
    matches: &[(u32, u32)],
    params: &ReplaceCoreParams,
) -> Result<Option<Mol>> {
    Ok(optional_mol(ffi::replace_core_with_matches(
        mol.as_ffi(),
        core.as_ffi(),
        &to_pairs(matches),
        params.replace_dummies,
        params.label_by_index,
        params.require_dummy_match,
    )?))
}

/// Parameters for `replace_core_query_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ReplaceCoreQueryParams {
    /// Replaces the dummy atoms the core matched too.
    pub replace_dummies: bool,
    /// Labels each sidechain's dummy with the core atom's index; otherwise
    /// they are numbered in order.
    pub label_by_index: bool,
    /// Returns `None` unless each sidechain attaches at a core dummy.
    pub require_dummy_match: bool,
    /// Uses chirality in the match.
    pub use_chirality: bool,
}

impl Default for ReplaceCoreQueryParams {
    fn default() -> ReplaceCoreQueryParams {
        ReplaceCoreQueryParams {
            replace_dummies: true,
            label_by_index: false,
            require_dummy_match: false,
            use_chirality: false,
        }
    }
}

/// As `replace_core`, matching `core_query` itself. Mirrors
/// `rdmolops.ReplaceCore(mol, coreQuery)` (A5).
pub fn replace_core_query(mol: &Mol, core_query: &Mol) -> Result<Option<Mol>> {
    replace_core_query_with_params(mol, core_query, &ReplaceCoreQueryParams::default())
}

/// Mirrors `rdmolops.ReplaceCore(mol, coreQuery, replaceDummies,
/// labelByIndex, requireDummyMatch, useChirality)`.
pub fn replace_core_query_with_params(
    mol: &Mol,
    core_query: &Mol,
    params: &ReplaceCoreQueryParams,
) -> Result<Option<Mol>> {
    Ok(optional_mol(ffi::replace_core(
        mol.as_ffi(),
        core_query.as_ffi(),
        params.replace_dummies,
        params.label_by_index,
        params.require_dummy_match,
        params.use_chirality,
    )?))
}

/// The Murcko scaffold: the rings and the linkers between them. Mirrors
/// `rdmolops.MurckoDecompose`.
pub fn murcko_decompose(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::murcko_decompose(mol.as_ffi())?))
}

/// One molecule holding both, unbonded. Mirrors `rdmolops.CombineMols`.
pub fn combine_mols(mol1: &Mol, mol2: &Mol) -> Result<Mol> {
    combine_mols_with_offset(mol1, mol2, Point3D::default())
}

/// As `combine_mols`, with `mol2`'s coordinates moved by `offset`. Mirrors
/// `rdmolops.CombineMols(mol1, mol2, offset)`.
pub fn combine_mols_with_offset(mol1: &Mol, mol2: &Mol, offset: Point3D) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::combine_mols(
        mol1.as_ffi(),
        mol2.as_ffi(),
        &offset,
    )?))
}

/// The match whose core atoms have the most heavy-atom substituents; each
/// match gives a molecule atom index for each core atom. Mirrors
/// `rdmolops.GetMostSubstitutedCoreMatch`. `Error::Value` for no matches or
/// an empty one.
pub fn get_most_substituted_core_match(
    mol: &Mol,
    core: &Mol,
    matches: &[Vec<u32>],
) -> Result<Vec<u32>> {
    Ok(ffi::get_most_substituted_core_match(
        mol.as_ffi(),
        core.as_ffi(),
        &IndexList::from_nested(matches),
    )?)
}

/// The matches, most substituted first. Mirrors
/// `rdmolops.SortMatchesByDegreeOfCoreSubstitution`. `Error::Value` for no
/// matches or an empty one.
pub fn sort_matches_by_degree_of_core_substitution(
    mol: &Mol,
    core: &Mol,
    matches: &[Vec<u32>],
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(
        ffi::sort_matches_by_degree_of_core_substitution(
            mol.as_ffi(),
            core.as_ffi(),
            &IndexList::from_nested(matches),
        )?,
    ))
}

/// Adds `query` as a recursive query on an atom, combined with its existing
/// query. Mirrors `rdmolops.AddRecursiveQuery`. `Error::Value` for an atom
/// index past the end.
pub fn add_recursive_query(mol: &mut Mol, query: &Mol, atom_idx: u32) -> Result<()> {
    add_recursive_query_with_preserve_existing_query(mol, query, atom_idx, true)
}

/// Mirrors `rdmolops.AddRecursiveQuery(mol, query, atomIdx,
/// preserveExistingQuery)`; without `preserve_existing_query`, the recursive
/// query replaces the atom's query.
pub fn add_recursive_query_with_preserve_existing_query(
    mol: &mut Mol,
    query: &Mol,
    atom_idx: u32,
    preserve_existing_query: bool,
) -> Result<()> {
    Ok(ffi::add_recursive_query(
        mol.pin_mut(),
        query.as_ffi(),
        atom_idx,
        preserve_existing_query,
    )?)
}

/// For each atom with the property `prop_name`, adds the queries it names
/// (comma-separated, lowercased) as recursive queries. Mirrors
/// `rdmolops.MolAddRecursiveQueries`. `Error::Key` for a name with no query
/// that is not SMARTS either.
pub fn mol_add_recursive_queries(
    mol: &mut Mol,
    queries: &BTreeMap<String, Mol>,
    prop_name: &str,
) -> Result<()> {
    with_borrowed_mols_shared(
        queries.iter().map(|(name, query)| (name.as_str(), query)),
        |list| ffi::add_recursive_queries(mol.pin_mut(), list, prop_name),
    )?;
    Ok(())
}

/// Parameters for `parse_mol_query_def_file_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ParseMolQueryDefFileParams {
    /// Lowercases the names.
    pub standardize: bool,
    /// The characters that separate columns.
    pub delimiter: String,
    /// Lines that start with this are skipped.
    pub comment: String,
    /// The column with each query's name.
    pub name_column: u32,
    /// The column with each query's SMARTS.
    pub smarts_column: u32,
}

impl Default for ParseMolQueryDefFileParams {
    fn default() -> ParseMolQueryDefFileParams {
        ParseMolQueryDefFileParams {
            standardize: true,
            delimiter: "\t".to_owned(),
            comment: "//".to_owned(),
            name_column: 0,
            smarts_column: 1,
        }
    }
}

/// The named SMARTS queries in a file; lines whose SMARTS do not parse are
/// skipped. Mirrors `rdmolops.ParseMolQueryDefFile(fileobj)` with a file
/// name. `Error::Os` if the file cannot be opened.
pub fn parse_mol_query_def_file(path: impl AsRef<Path>) -> Result<BTreeMap<String, Mol>> {
    parse_mol_query_def_file_with_params(path, &ParseMolQueryDefFileParams::default())
}

/// Mirrors `rdmolops.ParseMolQueryDefFile(fileobj, standardize, delimiter,
/// comment, nameColumn, smartsColumn)` with a file name.
pub fn parse_mol_query_def_file_with_params(
    path: impl AsRef<Path>,
    params: &ParseMolQueryDefFileParams,
) -> Result<BTreeMap<String, Mol>> {
    mol_list_into_map(ffi::parse_mol_query_def_file(
        path.as_ref().as_os_str().as_encoded_bytes(),
        params.standardize,
        &params.delimiter,
        &params.comment,
        params.name_column,
        params.smarts_column,
    )?)
}

/// As `parse_mol_query_def_file`, from the file's text. Mirrors
/// `rdmolops.ParseMolQueryDefFile(fileobj)` with a file object (A5).
pub fn parse_mol_query_def_file_text(text: &str) -> Result<BTreeMap<String, Mol>> {
    parse_mol_query_def_file_text_with_params(text, &ParseMolQueryDefFileParams::default())
}

/// Mirrors `rdmolops.ParseMolQueryDefFile(fileobj, standardize, delimiter,
/// comment, nameColumn, smartsColumn)` with a file object.
pub fn parse_mol_query_def_file_text_with_params(
    text: &str,
    params: &ParseMolQueryDefFileParams,
) -> Result<BTreeMap<String, Mol>> {
    mol_list_into_map(ffi::parse_mol_query_def_text(
        text,
        params.standardize,
        &params.delimiter,
        &params.comment,
        params.name_column,
        params.smarts_column,
    )?)
}

/// Parameters for `fragment_on_bonds_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentOnBondsParams {
    /// Caps each broken bond's ends with dummy atoms.
    pub add_dummies: bool,
    /// The isotope labels of the two dummies for each bond; by default each
    /// dummy is labeled with the index of the atom it replaces.
    pub dummy_labels: Option<Vec<(u32, u32)>>,
    /// The type of each new bond to a dummy; by default the broken bond's.
    /// Python treats an empty list as `None`.
    pub bond_types: Option<Vec<BondType>>,
}

impl Default for FragmentOnBondsParams {
    fn default() -> FragmentOnBondsParams {
        FragmentOnBondsParams {
            add_dummies: true,
            dummy_labels: None,
            bond_types: None,
        }
    }
}

/// What Python's `cutsPerAtom` output argument receives (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct FragmentOnBondsOutputs {
    /// The number of broken bonds at each atom.
    pub cuts_per_atom: Vec<u32>,
}

struct FragmentArgs {
    dummy_labels: Vec<IndexPair>,
    bond_types: Vec<crate::ffi::bond::ffi::BondType>,
}

fn fragment_args(
    dummy_labels: &Option<Vec<(u32, u32)>>,
    bond_types: &Option<Vec<BondType>>,
) -> FragmentArgs {
    FragmentArgs {
        dummy_labels: dummy_labels.as_deref().map(to_pairs).unwrap_or_default(),
        bond_types: bond_types.iter().flatten().map(|t| t.to_ffi()).collect(),
    }
}

/// A copy with the bonds broken. Mirrors `rdmolops.FragmentOnBonds`.
/// `Error::Value` for no bonds, a repeated bond, or option lists shorter
/// than the bonds; `Error::Index` for a bond index past the end.
pub fn fragment_on_bonds(mol: &Mol, bond_indices: &[u32]) -> Result<Mol> {
    fragment_on_bonds_with_params(mol, bond_indices, &FragmentOnBondsParams::default())
}

/// Mirrors `rdmolops.FragmentOnBonds(mol, bondIndices, addDummies,
/// dummyLabels, bondTypes)`.
pub fn fragment_on_bonds_with_params(
    mol: &Mol,
    bond_indices: &[u32],
    params: &FragmentOnBondsParams,
) -> Result<Mol> {
    Ok(fragment_on_bonds_with_outputs(mol, bond_indices, params)?.0)
}

/// The fragmented copy and Python's `cutsPerAtom` (A9). Mirrors
/// `rdmolops.FragmentOnBonds(..., cutsPerAtom)`.
pub fn fragment_on_bonds_with_outputs(
    mol: &Mol,
    bond_indices: &[u32],
    params: &FragmentOnBondsParams,
) -> Result<(Mol, FragmentOnBondsOutputs)> {
    let args = fragment_args(&params.dummy_labels, &params.bond_types);
    let mut cuts_per_atom = Vec::new();
    let result = ffi::fragment_on_bonds(
        mol.as_ffi(),
        bond_indices,
        params.add_dummies,
        params.dummy_labels.is_some(),
        &args.dummy_labels,
        params.bond_types.is_some(),
        &args.bond_types,
        &mut cuts_per_atom,
    )?;
    Ok((
        Mol::from_ffi(result),
        FragmentOnBondsOutputs { cuts_per_atom },
    ))
}

/// Parameters for `fragment_on_some_bonds_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentOnSomeBondsParams {
    /// How many of the bonds each fragmentation breaks.
    pub num_to_break: u32,
    /// Caps each broken bond's ends with dummy atoms.
    pub add_dummies: bool,
    /// The isotope labels of the two dummies for each bond; by default each
    /// dummy is labeled with the index of the atom it replaces.
    pub dummy_labels: Option<Vec<(u32, u32)>>,
    /// The type of each new bond to a dummy; by default the broken bond's.
    /// Python treats an empty list as `None`.
    pub bond_types: Option<Vec<BondType>>,
}

impl Default for FragmentOnSomeBondsParams {
    fn default() -> FragmentOnSomeBondsParams {
        FragmentOnSomeBondsParams {
            num_to_break: 1,
            add_dummies: true,
            dummy_labels: None,
            bond_types: None,
        }
    }
}

fn fragment_on_some_bonds_impl(
    mol: &Mol,
    bond_indices: &[u32],
    params: &FragmentOnSomeBondsParams,
    return_cuts_per_atom: bool,
) -> Result<(Vec<Mol>, Vec<Vec<u32>>)> {
    let args = fragment_args(&params.dummy_labels, &params.bond_types);
    let mut cuts_per_atom = Vec::new();
    let frags = ffi::fragment_on_some_bonds(
        mol.as_ffi(),
        bond_indices,
        params.num_to_break,
        params.add_dummies,
        params.dummy_labels.is_some(),
        &args.dummy_labels,
        params.bond_types.is_some(),
        &args.bond_types,
        return_cuts_per_atom,
        &mut cuts_per_atom,
    )?;
    Ok((
        mol_list_into_vec(frags),
        IndexList::into_nested(cuts_per_atom),
    ))
}

/// A copy for each way to break `num_to_break` of the bonds. Mirrors
/// `rdmolops.FragmentOnSomeBonds`. `Error::Value` for no bonds, more than
/// 63, a repeated bond, or option lists shorter than the bonds;
/// `Error::Index` for a bond index past the end.
pub fn fragment_on_some_bonds(mol: &Mol, bond_indices: &[u32]) -> Result<Vec<Mol>> {
    fragment_on_some_bonds_with_params(mol, bond_indices, &FragmentOnSomeBondsParams::default())
}

/// Mirrors `rdmolops.FragmentOnSomeBonds(mol, bondIndices, numToBreak,
/// addDummies, dummyLabels, bondTypes)`.
pub fn fragment_on_some_bonds_with_params(
    mol: &Mol,
    bond_indices: &[u32],
    params: &FragmentOnSomeBondsParams,
) -> Result<Vec<Mol>> {
    Ok(fragment_on_some_bonds_impl(mol, bond_indices, params, false)?.0)
}

/// The copies and, for each, the number of broken bonds at each atom.
/// Mirrors `rdmolops.FragmentOnSomeBonds(mol, bondIndices,
/// returnCutsPerAtom=True)` (A4).
pub fn fragment_on_some_bonds_return_cuts_per_atom(
    mol: &Mol,
    bond_indices: &[u32],
) -> Result<(Vec<Mol>, Vec<Vec<u32>>)> {
    fragment_on_some_bonds_return_cuts_per_atom_with_params(
        mol,
        bond_indices,
        &FragmentOnSomeBondsParams::default(),
    )
}

/// Mirrors `rdmolops.FragmentOnSomeBonds(mol, bondIndices, numToBreak,
/// addDummies, dummyLabels, bondTypes, returnCutsPerAtom=True)`.
pub fn fragment_on_some_bonds_return_cuts_per_atom_with_params(
    mol: &Mol,
    bond_indices: &[u32],
    params: &FragmentOnSomeBondsParams,
) -> Result<(Vec<Mol>, Vec<Vec<u32>>)> {
    fragment_on_some_bonds_impl(mol, bond_indices, params, true)
}

/// A copy broken at its BRICS bonds, each end capped with a dummy labeled
/// by its BRICS environment. Mirrors `rdmolops.FragmentOnBRICSBonds`.
pub fn fragment_on_brics_bonds(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::fragment_on_brics_bonds(mol.as_ffi())?))
}

rdkit_enum! {
    /// How molzip finds the atoms to join. Mirrors `rdmolops.MolzipLabel`.
    MolzipLabel: ffi::MolzipLabel {
        AtomMapNumber = AtomMapNumber,
        Isotope = Isotope,
        FragmentOnBonds = FragmentOnBonds,
        AtomType = AtomType,
    }
}

/// How molzip joins molecules. Mirrors `rdmolops.MolzipParams` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct MolzipParams {
    /// How the atoms to join are labeled.
    pub label: MolzipLabel,
    /// The atom symbols that mark attachment points for
    /// `MolzipLabel::AtomType`.
    pub atom_symbols: Vec<String>,
    /// Checks valences after zipping.
    pub enforce_valence_rules: bool,
    /// Gives the result, and for `molzip_fragments` the inputs, 2D
    /// coordinates; those calls return `Error::Value` for a result the
    /// depictor would mishandle (see `chem::rd_depictor`).
    pub generate_coordinates: bool,
    /// Aligns the inputs' coordinates along the new bonds.
    pub align_coordinates: bool,
}

impl MolzipParams {
    /// Sets the symbols for `MolzipLabel::AtomType`. Mirrors
    /// `MolzipParams.setAtomSymbols`.
    pub fn set_atom_symbols<S: AsRef<str>>(&mut self, symbols: &[S]) {
        self.atom_symbols = symbols.iter().map(|s| s.as_ref().to_owned()).collect();
    }

    fn to_ffi(&self) -> ffi::MolzipParamsData {
        ffi::MolzipParamsData {
            label: self.label.to_ffi(),
            atom_symbols: self.atom_symbols.clone(),
            enforce_valence_rules: self.enforce_valence_rules,
            generate_coordinates: self.generate_coordinates,
            align_coordinates: self.align_coordinates,
        }
    }
}

impl Default for MolzipParams {
    fn default() -> MolzipParams {
        let data = expect_ok(ffi::molzip_params_default());
        MolzipParams {
            label: MolzipLabel::from_ffi(data.label),
            atom_symbols: data.atom_symbols,
            enforce_valence_rules: data.enforce_valence_rules,
            generate_coordinates: data.generate_coordinates,
            align_coordinates: data.align_coordinates,
        }
    }
}

/// `a` and `b` joined at their matching labeled atoms, or `None` when
/// RDKit cannot join them. Mirrors `rdmolops.molzip(a, b)`.
/// `Error::Value` with `MolzipLabel::FragmentOnBonds` for a dummy without
/// exactly one neighbor, which RDKit would dereference as null.
pub fn molzip(a: &Mol, b: &Mol) -> Result<Option<Mol>> {
    molzip_with_params(a, b, &MolzipParams::default())
}

/// Mirrors `rdmolops.molzip(a, b, params)`.
pub fn molzip_with_params(a: &Mol, b: &Mol, params: &MolzipParams) -> Result<Option<Mol>> {
    Ok(optional_mol(ffi::molzip(
        a.as_ffi(),
        b.as_ffi(),
        &params.to_ffi(),
    )?))
}

/// The fragments within one molecule joined at their matching labeled
/// atoms. Mirrors `rdmolops.molzip(a)` (A5).
pub fn molzip_single(a: &Mol) -> Result<Option<Mol>> {
    molzip_single_with_params(a, &MolzipParams::default())
}

/// Mirrors `rdmolops.molzip(a, params)`.
pub fn molzip_single_with_params(a: &Mol, params: &MolzipParams) -> Result<Option<Mol>> {
    Ok(optional_mol(ffi::molzip_single(
        a.as_ffi(),
        &params.to_ffi(),
    )?))
}

/// The molecules of an R-group decomposition, core first, joined at their
/// matching labeled atoms; `None` for no molecules. With
/// `generate_coordinates`, each input gets the coordinates it has in the
/// result. Mirrors `rdmolops.molzipFragments`. `Error::Value` for a labeled
/// atom without exactly one neighbor, or a linker (two labeled atoms bonded
/// to each other): RDKit dereferences a null pointer for both when zipping
/// a list.
pub fn molzip_fragments(mols: &mut [Mol]) -> Result<Option<Mol>> {
    molzip_fragments_with_params(mols, &MolzipParams::default())
}

/// Mirrors `rdmolops.molzipFragments(mols, params)`.
pub fn molzip_fragments_with_params(
    mols: &mut [Mol],
    params: &MolzipParams,
) -> Result<Option<Mol>> {
    let data = params.to_ffi();
    Ok(optional_mol(with_borrowed_mols(
        mols.iter_mut().map(|mol| ("", mol)),
        |list| ffi::molzip_fragments(list, &data),
    )?))
}

/// An R-group decomposition row joined at its matching labeled atoms; the
/// row needs a `"Core"`. As `molzip_fragments`, the inputs may gain
/// coordinates. Mirrors `rdmolops.molzip(row)` (A5). `Error::Value`
/// without a core, and as `molzip_fragments`.
pub fn molzip_row(row: &mut BTreeMap<String, Mol>) -> Result<Option<Mol>> {
    molzip_row_with_params(row, &MolzipParams::default())
}

/// Mirrors `rdmolops.molzip(row, params)`.
pub fn molzip_row_with_params(
    row: &mut BTreeMap<String, Mol>,
    params: &MolzipParams,
) -> Result<Option<Mol>> {
    let data = params.to_ffi();
    Ok(optional_mol(with_borrowed_mols(
        row.iter_mut().map(|(name, mol)| (name.as_str(), mol)),
        |list| ffi::molzip_row(list, &data),
    )?))
}
