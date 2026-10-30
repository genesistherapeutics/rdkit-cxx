// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdRGroupDecomposition`: decomposing molecules into
//! a core and the R groups attached to it.
//!
//! RDKit's GA, and its parallel runs, share a random number generator and
//! a scoring static across the process, so decompositions process one at a
//! time, and parallel runs run one after another on the calling thread
//! (CONVENTIONS.md U4); `ga_parallel_runs` is accepted for parity.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smiles;
//! use rdkit_cxx::chem::rd_r_group_decomposition::r_group_decompose_as_smiles;
//!
//! let cores = vec![mol_from_smiles("c1ccccc1")?];
//! let mols = vec![mol_from_smiles("c1ccccc1Cl")?, mol_from_smiles("CCCC")?];
//! let (rows, unmatched) = r_group_decompose_as_smiles(&cores, &mols)?;
//! assert_eq!(rows[0]["R1"], "Cl[*:1]");
//! assert_eq!(unmatched, [1]);
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

use std::collections::BTreeMap;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::{Mol, SubstructMatchParameters};
use crate::chem::rdmolfiles::mol_to_smiles;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::error::expect_ok;
use crate::ffi::molops::ffi::{MolList, mol_list_name, mol_list_size, mol_list_take};
use crate::ffi::rgroup::ffi;

macro_rules! flags {
    ($(#[$meta:meta])* $type:ident { $($name:ident),* $(,)? }) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $type(u32);

        #[allow(non_upper_case_globals)]
        impl $type {
            $(
                #[doc = concat!("Mirrors `", stringify!($type), ".", stringify!($name), "`.")]
                pub const $name: $type = $type(ffi::$type::$name.repr);
            )*

            const NAMES: &[(&str, $type)] = &[$((stringify!($name), $type::$name)),*];
        }

        crate::flags::flag_ops!($type);
    };
}

flags!(
    /// Where a core's R-group labels are: bits, combined. Mirrors
    /// `rdRGroupDecomposition.RGroupLabels` (CONVENTIONS.md F1); the
    /// constants keep Python's spelling.
    RGroupLabels {
        IsotopeLabels,
        AtomMapLabels,
        AtomIndexLabels,
        RelabelDuplicateLabels,
        MDLRGroupLabels,
        DummyAtomLabels,
        AutoDetect,
    }
);

flags!(
    /// How molecules' matches are chosen. Mirrors
    /// `rdRGroupDecomposition.RGroupMatching` (CONVENTIONS.md F1); the
    /// constants keep Python's spelling.
    RGroupMatching {
        Greedy,
        GreedyChunks,
        Exhaustive,
        NoSymmetrization,
        GA,
    }
);

flags!(
    /// Where the results' R-group labels are written: bits, combined.
    /// Mirrors `rdRGroupDecomposition.RGroupLabelling` (CONVENTIONS.md F1);
    /// the constants keep Python's spelling.
    RGroupLabelling {
        AtomMap,
        Isotope,
        MDLRGroup,
    }
);

flags!(
    /// How cores are aligned to each other. Mirrors
    /// `rdRGroupDecomposition.RGroupCoreAlignment` (CONVENTIONS.md F1),
    /// whose deprecated `None` is `NoAlignment`; the constants keep Python's
    /// spelling.
    RGroupCoreAlignment {
        NoAlignment,
        MCS,
    }
);

flags!(
    /// How a decomposition is scored. Mirrors
    /// `rdRGroupDecomposition.RGroupScore` (CONVENTIONS.md F1); the
    /// constants keep Python's spelling.
    RGroupScore {
        Match,
        FingerprintVariance,
    }
);

/// How a decomposition labels and matches. Mirrors
/// `rdRGroupDecomposition.RGroupDecompositionParameters`.
#[derive(Clone, Debug, PartialEq)]
pub struct RGroupDecompositionParameters {
    /// Where the cores' labels are. Mirrors
    /// `RGroupDecompositionParameters.labels`.
    pub labels: RGroupLabels,
    /// Mirrors `RGroupDecompositionParameters.matchingStrategy`.
    pub matching_strategy: RGroupMatching,
    /// Mirrors `RGroupDecompositionParameters.scoreMethod`.
    pub score_method: RGroupScore,
    /// Where the results' labels go. Mirrors
    /// `RGroupDecompositionParameters.rgroupLabelling`.
    pub rgroup_labelling: RGroupLabelling,
    /// Mirrors `RGroupDecompositionParameters.alignment`.
    pub alignment: RGroupCoreAlignment,
    /// The molecules `GreedyChunks` processes at a time, above 0 for it.
    /// Mirrors `RGroupDecompositionParameters.chunkSize`.
    pub chunk_size: u32,
    /// Decomposes only at the cores' labeled R groups. Mirrors
    /// `RGroupDecompositionParameters.onlyMatchAtRGroups`.
    pub only_match_at_r_groups: bool,
    /// Drops labeled R groups that are only hydrogens. Mirrors
    /// `RGroupDecompositionParameters.removeAllHydrogenRGroups`.
    pub remove_all_hydrogen_r_groups: bool,
    /// Drops them and their labels from the core too. Mirrors
    /// `RGroupDecompositionParameters.removeAllHydrogenRGroupsAndLabels`.
    pub remove_all_hydrogen_r_groups_and_labels: bool,
    /// Removes hydrogens from the results. Mirrors
    /// `RGroupDecompositionParameters.removeHydrogensPostMatch`.
    pub remove_hydrogens_post_match: bool,
    /// Allows labeled R groups of degree 2 or more. Mirrors
    /// `RGroupDecompositionParameters.allowNonTerminalRGroups`.
    pub allow_non_terminal_r_groups: bool,
    /// Allows several R groups on an unlabeled core atom. Mirrors
    /// `RGroupDecompositionParameters.allowMultipleRGroupsOnUnlabelled`.
    pub allow_multiple_r_groups_on_unlabelled: bool,
    /// Matches the cores' tautomers. RDKit's tautomer matcher reads the
    /// matches of a recursive SMARTS (`$(...)`) as the core's, past their
    /// end, so a core with one is `Error::Value`. Mirrors
    /// `RGroupDecompositionParameters.doTautomers`.
    pub do_tautomers: bool,
    /// Enumerates the cores' link nodes, position variation, and repeat
    /// units. Link nodes RDKit's enumerator would read past an empty list
    /// on (as `rd_generalized_substruct` checks them) are `Error::Value`.
    /// Mirrors `RGroupDecompositionParameters.doEnumeration`.
    pub do_enumeration: bool,
    /// Seconds, none when not positive. Mirrors
    /// `RGroupDecompositionParameters.timeout`.
    pub timeout: f64,
    /// The GA's population, or RDKit's guess for -1. Mirrors
    /// `RGroupDecompositionParameters.gaPopulationSize`.
    pub ga_population_size: i32,
    /// Mirrors `RGroupDecompositionParameters.gaMaximumOperations`.
    pub ga_maximum_operations: i32,
    /// Mirrors
    /// `RGroupDecompositionParameters.gaNumberOperationsWithoutImprovement`.
    pub ga_number_operations_without_improvement: i32,
    /// The GA's seed: -1 for its default, -2 for a random one. Mirrors
    /// `RGroupDecompositionParameters.gaRandomSeed`.
    pub ga_random_seed: i32,
    /// Mirrors `RGroupDecompositionParameters.gaNumberRuns`.
    pub ga_number_runs: i32,
    /// Accepted for parity: the runs go one after another on the calling
    /// thread (see the module documentation). Mirrors
    /// `RGroupDecompositionParameters.gaParallelRuns`.
    pub ga_parallel_runs: bool,
    /// How cores are matched. Mirrors
    /// `RGroupDecompositionParameters.substructMatchParams`, which Python
    /// reads only, as a reference that edits the parameters'.
    pub substruct_match_params: SubstructMatchParameters,
}

/// RDKit's defaults. Mirrors `RGroupDecompositionParameters()`.
impl Default for RGroupDecompositionParameters {
    fn default() -> RGroupDecompositionParameters {
        let d = expect_ok(ffi::r_group_decomposition_parameters_default());
        RGroupDecompositionParameters {
            labels: RGroupLabels(d.labels),
            matching_strategy: RGroupMatching(d.matching_strategy),
            score_method: RGroupScore(d.score_method),
            rgroup_labelling: RGroupLabelling(d.rgroup_labelling),
            alignment: RGroupCoreAlignment(d.alignment),
            chunk_size: d.chunk_size,
            only_match_at_r_groups: d.only_match_at_r_groups,
            remove_all_hydrogen_r_groups: d.remove_all_hydrogen_r_groups,
            remove_all_hydrogen_r_groups_and_labels: d.remove_all_hydrogen_r_groups_and_labels,
            remove_hydrogens_post_match: d.remove_hydrogens_post_match,
            allow_non_terminal_r_groups: d.allow_non_terminal_r_groups,
            allow_multiple_r_groups_on_unlabelled: d.allow_multiple_r_groups_on_unlabelled,
            do_tautomers: d.do_tautomers,
            do_enumeration: d.do_enumeration,
            timeout: d.timeout,
            ga_population_size: d.ga_population_size,
            ga_maximum_operations: d.ga_maximum_operations,
            ga_number_operations_without_improvement: d.ga_number_operations_without_improvement,
            ga_random_seed: d.ga_random_seed,
            ga_number_runs: d.ga_number_runs,
            ga_parallel_runs: d.ga_parallel_runs,
            substruct_match_params: SubstructMatchParameters::from_ffi(d.substruct_match_params),
        }
    }
}

impl RGroupDecompositionParameters {
    fn to_ffi(&self) -> ffi::RGroupDecompositionParametersData {
        ffi::RGroupDecompositionParametersData {
            labels: self.labels.bits(),
            matching_strategy: self.matching_strategy.bits(),
            score_method: self.score_method.bits(),
            rgroup_labelling: self.rgroup_labelling.bits(),
            alignment: self.alignment.bits(),
            chunk_size: self.chunk_size,
            only_match_at_r_groups: self.only_match_at_r_groups,
            remove_all_hydrogen_r_groups: self.remove_all_hydrogen_r_groups,
            remove_all_hydrogen_r_groups_and_labels: self.remove_all_hydrogen_r_groups_and_labels,
            remove_hydrogens_post_match: self.remove_hydrogens_post_match,
            allow_non_terminal_r_groups: self.allow_non_terminal_r_groups,
            allow_multiple_r_groups_on_unlabelled: self.allow_multiple_r_groups_on_unlabelled,
            do_tautomers: self.do_tautomers,
            do_enumeration: self.do_enumeration,
            timeout: self.timeout,
            ga_population_size: self.ga_population_size,
            ga_maximum_operations: self.ga_maximum_operations,
            ga_number_operations_without_improvement: self.ga_number_operations_without_improvement,
            ga_random_seed: self.ga_random_seed,
            ga_number_runs: self.ga_number_runs,
            ga_parallel_runs: self.ga_parallel_runs,
            substruct_match_params: self.substruct_match_params.to_ffi(),
        }
    }
}

/// One molecule's decomposition: its core and R groups by label (`Core`,
/// `R1`, ...).
/// A Rust addition, naming Python's dictionary.
pub type RGroupRow<T> = BTreeMap<String, T>;

/// Each label's molecules, one per decomposed molecule (an empty molecule
/// where one has no such group). A Rust addition, naming Python's dictionary.
pub type RGroupColumns<T> = BTreeMap<String, Vec<T>>;

/// Moves the named molecules out of a table, in order, with their labels.
fn take_table(table: ffi::RGroupTable) -> Result<Vec<Vec<(String, Mol)>>> {
    let mut list: UniquePtr<MolList> = table.mols;
    let size = expect_ok(mol_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    let mut entries = Vec::with_capacity(size);
    for i in 0..size {
        let name = expect_ok(mol_list_name(&list, i))
            .to_str()
            .map_err(|err| crate::Error::Value {
                message: format!("invalid UTF-8 in a string from RDKit: {err}"),
            })?
            .to_owned();
        entries.push((
            name,
            Mol::from_ffi(expect_ok(mol_list_take(list.pin_mut(), i))),
        ));
    }
    let mut entries = entries.into_iter();
    Ok(table
        .sizes
        .iter()
        .map(|&n| entries.by_ref().take(n as usize).collect())
        .collect())
}

fn rows(table: ffi::RGroupTable) -> Result<Vec<RGroupRow<Mol>>> {
    Ok(take_table(table)?
        .into_iter()
        .map(|row| row.into_iter().collect())
        .collect())
}

fn columns(table: ffi::RGroupTable) -> Result<RGroupColumns<Mol>> {
    Ok(take_table(table)?
        .into_iter()
        .filter_map(|column| {
            let label = column.first()?.0.clone();
            Some((label, column.into_iter().map(|(_, mol)| mol).collect()))
        })
        .collect())
}

/// RDKit's `MolToSmiles(mol, isomericSmiles=true)`, as Python's `asSmiles`.
fn smiles_rows(rows: Vec<RGroupRow<Mol>>) -> Result<Vec<RGroupRow<String>>> {
    rows.into_iter()
        .map(|row| {
            row.into_iter()
                .map(|(label, mol)| Ok((label, mol_to_smiles(&mol)?)))
                .collect()
        })
        .collect()
}

fn smiles_columns(columns: RGroupColumns<Mol>) -> Result<RGroupColumns<String>> {
    columns
        .into_iter()
        .map(|(label, mols)| {
            Ok((
                label,
                mols.iter().map(mol_to_smiles).collect::<Result<_>>()?,
            ))
        })
        .collect()
}

/// A decomposition of molecules into the cores' R groups, added one at a
/// time and then processed. Mirrors
/// `rdRGroupDecomposition.RGroupDecomposition`.
///
/// Every method takes `&mut self`, since the getters rewrite the
/// decomposition's matched cores (CONVENTIONS.md S4); results are copies.
/// `RGroupDecomposition` is `Send` but not `Sync`.
pub struct RGroupDecomposition {
    /// Never null.
    ptr: UniquePtr<ffi::RGroupDecomposition>,
}

impl RGroupDecomposition {
    /// A decomposition for the cores, from most to least specific, with
    /// RDKit's parameters. Mirrors `RGroupDecomposition(cores)`, for one
    /// core or a list of them.
    pub fn new(cores: &[Mol]) -> Result<RGroupDecomposition> {
        RGroupDecomposition::new_with_params(cores, &RGroupDecompositionParameters::default())
    }

    /// A decomposition with parameters. `GreedyChunks` with a `chunk_size`
    /// of 0, a core label (an atom map number or MDL R label on a dummy
    /// atom, or a `tempRlabel` property) of 2^30 or more, and, for MCS
    /// alignment, a core bond of
    /// a type past `BondType::Zero` (from a pickle) are `Error::Value`: RDKit
    /// divides by the chunk size, adds one to labels in `int`, and indexes
    /// a table by bond types. A core RDKit cannot prepare is
    /// `Error::Runtime`. Mirrors `RGroupDecomposition(cores, params)`.
    pub fn new_with_params(
        cores: &[Mol],
        params: &RGroupDecompositionParameters,
    ) -> Result<RGroupDecomposition> {
        let data = params.to_ffi();
        let ptr = with_borrowed_mols_shared(cores.iter().map(|core| ("", core)), |list| {
            ffi::r_group_decomposition_new(list, &data)
        })?;
        assert!(!ptr.is_null(), "RDKit returned a null decomposition");
        Ok(RGroupDecomposition { ptr })
    }

    /// Adds a molecule (a copy), returning its index, -1 when no core
    /// matches it, or -2 when it has no R groups; with `Greedy` and
    /// `GreedyChunks`, it may process the molecules so far. Mirrors
    /// `RGroupDecomposition.Add`.
    pub fn add(&mut self, mol: &Mol) -> Result<i32> {
        Ok(ffi::r_group_decomposition_add(
            self.ptr.pin_mut(),
            mol.as_ffi(),
        )?)
    }

    /// The index of the core that matches the molecule, or -1. Mirrors
    /// `RGroupDecomposition.GetMatchingCoreIdx`.
    pub fn get_matching_core_idx(&mut self, mol: &Mol) -> Result<i32> {
        Ok(self.get_matching_core_idx_with_outputs(mol)?.0)
    }

    /// As `get_matching_core_idx`, with the core's matches, which Python's
    /// `matches` output argument receives (A9). Mirrors
    /// `RGroupDecomposition.GetMatchingCoreIdx(mol, matches)`.
    pub fn get_matching_core_idx_with_outputs(
        &mut self,
        mol: &Mol,
    ) -> Result<(i32, GetMatchingCoreIdxOutputs)> {
        let result =
            ffi::r_group_decomposition_get_matching_core_idx(self.ptr.pin_mut(), mol.as_ffi())?;
        let matches = result
            .matches
            .iter()
            .map(|pairs| {
                pairs
                    .indices
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|&[core, atom]| (core, atom))
                    .collect()
            })
            .collect();
        Ok((result.core_idx, GetMatchingCoreIdxOutputs { matches }))
    }

    /// Processes the molecules, returning whether it could; it must come
    /// before the getters. Mirrors `RGroupDecomposition.Process`.
    pub fn process(&mut self) -> bool {
        expect_ok(ffi::r_group_decomposition_process(self.ptr.pin_mut()))
    }

    /// As `process`, with the best score. Mirrors
    /// `RGroupDecomposition.ProcessAndScore`.
    pub fn process_and_score(&mut self) -> (bool, f64) {
        let result = expect_ok(ffi::r_group_decomposition_process_and_score(
            self.ptr.pin_mut(),
        ));
        (result.success, result.score)
    }

    /// The labels the results use (`Core`, `R1`, ...), sorted; molecules
    /// added since processing are `Error::Runtime`. Mirrors
    /// `RGroupDecomposition.GetRGroupLabels`.
    pub fn get_r_group_labels(&mut self) -> Result<Vec<String>> {
        ffi::r_group_decomposition_get_r_group_labels(self.ptr.pin_mut())?
            .iter()
            .map(|label| {
                label
                    .to_str()
                    .map(str::to_owned)
                    .map_err(|err| crate::Error::Value {
                        message: format!("invalid UTF-8 in a string from RDKit: {err}"),
                    })
            })
            .collect()
    }

    /// Each decomposed molecule's core and R groups; molecules added since
    /// processing are `Error::Runtime`. Mirrors
    /// `RGroupDecomposition.GetRGroupsAsRows`.
    pub fn get_r_groups_as_rows(&mut self) -> Result<Vec<RGroupRow<Mol>>> {
        rows(ffi::r_group_decomposition_get_r_groups_as_rows(
            self.ptr.pin_mut(),
        )?)
    }

    /// As `get_r_groups_as_rows`, with isomeric SMILES. Mirrors
    /// `RGroupDecomposition.GetRGroupsAsRows(asSmiles=True)`.
    pub fn get_r_groups_as_rows_as_smiles(&mut self) -> Result<Vec<RGroupRow<String>>> {
        smiles_rows(self.get_r_groups_as_rows()?)
    }

    /// The cores and R groups by label, one per decomposed molecule.
    /// Mirrors `RGroupDecomposition.GetRGroupsAsColumns`.
    pub fn get_r_groups_as_columns(&mut self) -> Result<RGroupColumns<Mol>> {
        columns(ffi::r_group_decomposition_get_r_groups_as_columns(
            self.ptr.pin_mut(),
        )?)
    }

    /// As `get_r_groups_as_columns`, with isomeric SMILES. Mirrors
    /// `RGroupDecomposition.GetRGroupsAsColumns(asSmiles=True)`.
    pub fn get_r_groups_as_columns_as_smiles(&mut self) -> Result<RGroupColumns<String>> {
        smiles_columns(self.get_r_groups_as_columns()?)
    }
}

impl std::fmt::Debug for RGroupDecomposition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RGroupDecomposition")
            .finish_non_exhaustive()
    }
}

/// What Python's `matches` output argument receives (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct GetMatchingCoreIdxOutputs {
    /// Each of the core's matches, as (core atom, molecule atom) pairs
    /// for the molecule's heavy atoms.
    pub matches: Vec<Vec<(u32, u32)>>,
}

fn decompose(
    cores: &[Mol],
    mols: &[Mol],
    options: &RGroupDecompositionParameters,
    as_rows: bool,
) -> Result<(ffi::RGroupTable, Vec<u32>)> {
    let data = options.to_ffi();
    let result = with_borrowed_mols_shared(cores.iter().map(|core| ("", core)), |cores| {
        with_borrowed_mols_shared(mols.iter().map(|mol| ("", mol)), |mols| {
            ffi::r_group_decompose(cores, mols, &data, as_rows)
        })
    })?;
    Ok((result.table, result.unmatched))
}

/// Decomposes the molecules over the cores, from most to least specific,
/// returning each decomposed molecule's core and R groups and the indices
/// of the molecules no core matched. Past the timeout, after a molecule,
/// is `Error::Runtime`. Mirrors `rdRGroupDecomposition.RGroupDecompose`.
pub fn r_group_decompose(cores: &[Mol], mols: &[Mol]) -> Result<(Vec<RGroupRow<Mol>>, Vec<u32>)> {
    r_group_decompose_with_options(cores, mols, &RGroupDecompositionParameters::default())
}

/// Mirrors `rdRGroupDecomposition.RGroupDecompose(cores, mols,
/// options=options)`; the options are checked as
/// `RGroupDecomposition::new_with_params` checks them.
pub fn r_group_decompose_with_options(
    cores: &[Mol],
    mols: &[Mol],
    options: &RGroupDecompositionParameters,
) -> Result<(Vec<RGroupRow<Mol>>, Vec<u32>)> {
    let (table, unmatched) = decompose(cores, mols, options, true)?;
    Ok((rows(table)?, unmatched))
}

/// As `r_group_decompose`, with isomeric SMILES. Mirrors
/// `rdRGroupDecomposition.RGroupDecompose(asSmiles=True)`.
pub fn r_group_decompose_as_smiles(
    cores: &[Mol],
    mols: &[Mol],
) -> Result<(Vec<RGroupRow<String>>, Vec<u32>)> {
    r_group_decompose_as_smiles_with_options(cores, mols, &RGroupDecompositionParameters::default())
}

/// Mirrors `rdRGroupDecomposition.RGroupDecompose(asSmiles=True,
/// options=options)`.
pub fn r_group_decompose_as_smiles_with_options(
    cores: &[Mol],
    mols: &[Mol],
    options: &RGroupDecompositionParameters,
) -> Result<(Vec<RGroupRow<String>>, Vec<u32>)> {
    let (rows, unmatched) = r_group_decompose_with_options(cores, mols, options)?;
    Ok((smiles_rows(rows)?, unmatched))
}

/// As `r_group_decompose`, by label. Mirrors
/// `rdRGroupDecomposition.RGroupDecompose(asRows=False)`.
pub fn r_group_decompose_as_columns(
    cores: &[Mol],
    mols: &[Mol],
) -> Result<(RGroupColumns<Mol>, Vec<u32>)> {
    r_group_decompose_as_columns_with_options(
        cores,
        mols,
        &RGroupDecompositionParameters::default(),
    )
}

/// Mirrors `rdRGroupDecomposition.RGroupDecompose(asRows=False,
/// options=options)`.
pub fn r_group_decompose_as_columns_with_options(
    cores: &[Mol],
    mols: &[Mol],
    options: &RGroupDecompositionParameters,
) -> Result<(RGroupColumns<Mol>, Vec<u32>)> {
    let (table, unmatched) = decompose(cores, mols, options, false)?;
    Ok((columns(table)?, unmatched))
}

/// As `r_group_decompose_as_columns`, with isomeric SMILES. Mirrors
/// `rdRGroupDecomposition.RGroupDecompose(asSmiles=True, asRows=False)`.
pub fn r_group_decompose_as_columns_as_smiles(
    cores: &[Mol],
    mols: &[Mol],
) -> Result<(RGroupColumns<String>, Vec<u32>)> {
    r_group_decompose_as_columns_as_smiles_with_options(
        cores,
        mols,
        &RGroupDecompositionParameters::default(),
    )
}

/// Mirrors `rdRGroupDecomposition.RGroupDecompose(asSmiles=True,
/// asRows=False, options=options)`.
pub fn r_group_decompose_as_columns_as_smiles_with_options(
    cores: &[Mol],
    mols: &[Mol],
    options: &RGroupDecompositionParameters,
) -> Result<(RGroupColumns<String>, Vec<u32>)> {
    let (columns, unmatched) = r_group_decompose_as_columns_with_options(cores, mols, options)?;
    Ok((smiles_columns(columns)?, unmatched))
}

/// Parameters for `relabel_mapped_dummies_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelabelMappedDummiesParams {
    /// Where to read the labels from, in priority order: atom map numbers,
    /// isotopes, MDL R labels.
    pub input_labels: RGroupLabelling,
    /// Where to write them.
    pub output_labels: RGroupLabelling,
}

impl Default for RelabelMappedDummiesParams {
    fn default() -> RelabelMappedDummiesParams {
        RelabelMappedDummiesParams {
            input_labels: RGroupLabelling::AtomMap
                | RGroupLabelling::Isotope
                | RGroupLabelling::MDLRGroup,
            output_labels: RGroupLabelling::MDLRGroup,
        }
    }
}

/// Moves dummy atoms' R-group labels to MDL R labels, so that drawings show
/// them as R groups. Mirrors `rdRGroupDecomposition.RelabelMappedDummies`.
pub fn relabel_mapped_dummies(mol: &mut Mol) -> Result<()> {
    relabel_mapped_dummies_with_params(mol, &RelabelMappedDummiesParams::default())
}

/// A dummy atom whose atom map number is `i32::MIN` is `Error::Value`
/// when atom map numbers are read: RDKit takes its absolute value. Writing
/// an atom map number of 1000 or more is `Error::Runtime`. Mirrors
/// `rdRGroupDecomposition.RelabelMappedDummies(mol, inputLabels,
/// outputLabels)`.
pub fn relabel_mapped_dummies_with_params(
    mol: &mut Mol,
    params: &RelabelMappedDummiesParams,
) -> Result<()> {
    Ok(ffi::relabel_mapped_dummies(
        mol.pin_mut(),
        params.input_labels.bits(),
        params.output_labels.bits(),
    )?)
}
