// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Query adjustment.

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::adjust_query::ffi;

/// Which atoms or bonds an `AdjustQueryParameters` flag skips, as bits.
/// Mirrors `rdmolops.AdjustQueryWhichFlags` (CONVENTIONS.md F1).
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct AdjustQueryWhichFlags(u32);

macro_rules! adjust_query_which_flags {
    ($($name:ident),* $(,)?) => {
        impl AdjustQueryWhichFlags {
            $(
                #[doc = concat!("Mirrors `AdjustQueryWhichFlags.", stringify!($name), "`.")]
                pub const $name: AdjustQueryWhichFlags =
                    AdjustQueryWhichFlags(ffi::AdjustQueryWhichFlags::$name.repr);
            )*

            const NAMES: &[(&str, AdjustQueryWhichFlags)] =
                &[$((stringify!($name), AdjustQueryWhichFlags::$name)),*];
        }
    };
}

adjust_query_which_flags!(
    ADJUST_IGNORENONE,
    ADJUST_IGNORECHAINS,
    ADJUST_IGNORERINGS,
    ADJUST_IGNOREDUMMIES,
    ADJUST_IGNORENONDUMMIES,
    ADJUST_IGNOREMAPPED,
    ADJUST_IGNOREALL,
);

crate::flags::flag_ops!(AdjustQueryWhichFlags);

/// How `adjust_query_properties` changes a query. Mirrors
/// `rdmolops.AdjustQueryParameters`.
#[derive(Clone, Debug, PartialEq)]
pub struct AdjustQueryParameters {
    /// Adds degree queries.
    pub adjust_degree: bool,
    /// Which atoms get degree queries.
    pub adjust_degree_flags: AdjustQueryWhichFlags,
    /// Adds heavy-atom degree queries.
    pub adjust_heavy_degree: bool,
    /// Which atoms get heavy-atom degree queries.
    pub adjust_heavy_degree_flags: AdjustQueryWhichFlags,
    /// Adds ring-count queries.
    pub adjust_ring_count: bool,
    /// Which atoms get ring-count queries.
    pub adjust_ring_count_flags: AdjustQueryWhichFlags,
    /// Makes dummy atoms without isotopes match any atom.
    pub make_dummies_queries: bool,
    /// Perceives aromaticity.
    pub aromatize_if_possible: bool,
    /// Makes bonds match any bond.
    pub make_bonds_generic: bool,
    /// Which bonds match any bond.
    pub make_bonds_generic_flags: AdjustQueryWhichFlags,
    /// Makes atoms match any atom.
    pub make_atoms_generic: bool,
    /// Which atoms match any atom.
    pub make_atoms_generic_flags: AdjustQueryWhichFlags,
    /// Adds ring-or-chain queries.
    pub adjust_ring_chain: bool,
    /// Which atoms get ring-or-chain queries.
    pub adjust_ring_chain_flags: AdjustQueryWhichFlags,
    /// Removes stereo from double bonds without the `molStereoCare` property.
    pub use_stereo_care_for_bonds: bool,
    /// Makes conjugated five-membered rings' bonds match single, double, or aromatic.
    pub adjust_conjugated_five_rings: bool,
    /// Uses the MDL aromaticity model for five-membered rings.
    pub set_mdl_five_ring_aromaticity: bool,
    /// Makes single bonds between aromatic atoms and degree-one neighbors match single or aromatic.
    pub adjust_single_bonds_to_degree_one_neighbors: bool,
    /// Makes single bonds between aromatic or conjugated atoms match single or aromatic.
    pub adjust_single_bonds_between_aromatic_atoms: bool,
}

impl AdjustQueryParameters {
    pub(crate) fn from_ffi(data: ffi::AdjustQueryParametersData) -> AdjustQueryParameters {
        AdjustQueryParameters {
            adjust_degree: data.adjust_degree,
            adjust_degree_flags: AdjustQueryWhichFlags(data.adjust_degree_flags),
            adjust_heavy_degree: data.adjust_heavy_degree,
            adjust_heavy_degree_flags: AdjustQueryWhichFlags(data.adjust_heavy_degree_flags),
            adjust_ring_count: data.adjust_ring_count,
            adjust_ring_count_flags: AdjustQueryWhichFlags(data.adjust_ring_count_flags),
            make_dummies_queries: data.make_dummies_queries,
            aromatize_if_possible: data.aromatize_if_possible,
            make_bonds_generic: data.make_bonds_generic,
            make_bonds_generic_flags: AdjustQueryWhichFlags(data.make_bonds_generic_flags),
            make_atoms_generic: data.make_atoms_generic,
            make_atoms_generic_flags: AdjustQueryWhichFlags(data.make_atoms_generic_flags),
            adjust_ring_chain: data.adjust_ring_chain,
            adjust_ring_chain_flags: AdjustQueryWhichFlags(data.adjust_ring_chain_flags),
            use_stereo_care_for_bonds: data.use_stereo_care_for_bonds,
            adjust_conjugated_five_rings: data.adjust_conjugated_five_rings,
            set_mdl_five_ring_aromaticity: data.set_mdl_five_ring_aromaticity,
            adjust_single_bonds_to_degree_one_neighbors: data
                .adjust_single_bonds_to_degree_one_neighbors,
            adjust_single_bonds_between_aromatic_atoms: data
                .adjust_single_bonds_between_aromatic_atoms,
        }
    }

    pub(crate) fn to_ffi(&self) -> ffi::AdjustQueryParametersData {
        ffi::AdjustQueryParametersData {
            adjust_degree: self.adjust_degree,
            adjust_degree_flags: self.adjust_degree_flags.bits(),
            adjust_heavy_degree: self.adjust_heavy_degree,
            adjust_heavy_degree_flags: self.adjust_heavy_degree_flags.bits(),
            adjust_ring_count: self.adjust_ring_count,
            adjust_ring_count_flags: self.adjust_ring_count_flags.bits(),
            make_dummies_queries: self.make_dummies_queries,
            aromatize_if_possible: self.aromatize_if_possible,
            make_bonds_generic: self.make_bonds_generic,
            make_bonds_generic_flags: self.make_bonds_generic_flags.bits(),
            make_atoms_generic: self.make_atoms_generic,
            make_atoms_generic_flags: self.make_atoms_generic_flags.bits(),
            adjust_ring_chain: self.adjust_ring_chain,
            adjust_ring_chain_flags: self.adjust_ring_chain_flags.bits(),
            use_stereo_care_for_bonds: self.use_stereo_care_for_bonds,
            adjust_conjugated_five_rings: self.adjust_conjugated_five_rings,
            set_mdl_five_ring_aromaticity: self.set_mdl_five_ring_aromaticity,
            adjust_single_bonds_to_degree_one_neighbors: self
                .adjust_single_bonds_to_degree_one_neighbors,
            adjust_single_bonds_between_aromatic_atoms: self
                .adjust_single_bonds_between_aromatic_atoms,
        }
    }

    /// Parameters that change nothing, to turn on one adjustment at a time.
    /// Mirrors `AdjustQueryParameters.NoAdjustments`.
    pub fn no_adjustments() -> AdjustQueryParameters {
        AdjustQueryParameters::from_ffi(expect_ok(ffi::adjust_query_parameters_no_adjustments()))
    }
}

impl Default for AdjustQueryParameters {
    fn default() -> AdjustQueryParameters {
        AdjustQueryParameters::from_ffi(expect_ok(ffi::adjust_query_parameters_default()))
    }
}

/// A copy of a query with RDKit's default adjustments: degree queries on
/// ring atoms that aren't dummies, dummies that match any atom, and
/// aromaticity. Mirrors `rdmolops.AdjustQueryProperties`.
pub fn adjust_query_properties(mol: &Mol) -> Result<Mol> {
    adjust_query_properties_with_params(mol, &AdjustQueryParameters::default())
}

/// Mirrors `rdmolops.AdjustQueryProperties(mol, params)`.
pub fn adjust_query_properties_with_params(
    mol: &Mol,
    params: &AdjustQueryParameters,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::adjust_query_properties(
        mol.as_ffi(),
        &params.to_ffi(),
    )?))
}

/// As `adjust_query_properties`, and turns the generic groups in the
/// molecule's atom labels and substance groups into queries. Mirrors
/// `rdmolops.AdjustQueryPropertiesWithGenericGroups`.
pub fn adjust_query_properties_with_generic_groups(mol: &Mol) -> Result<Mol> {
    adjust_query_properties_with_generic_groups_with_params(mol, &AdjustQueryParameters::default())
}

/// Mirrors `rdmolops.AdjustQueryPropertiesWithGenericGroups(mol, params)`.
pub fn adjust_query_properties_with_generic_groups_with_params(
    mol: &Mol,
    params: &AdjustQueryParameters,
) -> Result<Mol> {
    Ok(Mol::from_ffi(
        ffi::adjust_query_properties_with_generic_groups(mol.as_ffi(), &params.to_ffi())?,
    ))
}

/// Parameters for `set_generic_queries_from_properties_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SetGenericQueriesFromPropertiesParams {
    /// Reads generic groups from atom labels.
    pub use_atom_labels: bool,
    /// Reads generic groups from substance groups.
    pub use_s_groups: bool,
}

impl Default for SetGenericQueriesFromPropertiesParams {
    fn default() -> SetGenericQueriesFromPropertiesParams {
        SetGenericQueriesFromPropertiesParams {
            use_atom_labels: true,
            use_s_groups: true,
        }
    }
}

/// Records the generic groups (`ARY`, `ALK`, and the like) the molecule's
/// atom labels and substance groups name, for generic substructure matches.
/// Mirrors `rdmolops.SetGenericQueriesFromProperties`.
pub fn set_generic_queries_from_properties(mol: &mut Mol) -> Result<()> {
    set_generic_queries_from_properties_with_params(mol, &Default::default())
}

/// Mirrors `rdmolops.SetGenericQueriesFromProperties(mol, useAtomLabels,
/// useSGroups)`.
pub fn set_generic_queries_from_properties_with_params(
    mol: &mut Mol,
    params: &SetGenericQueriesFromPropertiesParams,
) -> Result<()> {
    Ok(ffi::set_generic_queries_from_properties(
        mol.pin_mut(),
        params.use_atom_labels,
        params.use_s_groups,
    )?)
}

/// Turns the generic-group queries back into substance groups. Mirrors
/// `rdmolops.ConvertGenericQueriesToSubstanceGroups`.
pub fn convert_generic_queries_to_substance_groups(mol: &mut Mol) -> Result<()> {
    Ok(ffi::convert_generic_queries_to_substance_groups(
        mol.pin_mut(),
    )?)
}
