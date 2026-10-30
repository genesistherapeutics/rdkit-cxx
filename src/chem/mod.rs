// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem`, which re-exports `rdchem`, `rdmolfiles`, and
//! `rdmolops` (and `rdCIPLabeler` and `rdMolInterchange`), and defines
//! Python functions of its own (`FindMolChiralCenters`), ported here; its
//! other submodules are reached by path, as in Python.

pub mod all_chem;
pub mod atom_pairs;
pub mod brics;
pub mod build_fragment_catalog;
pub mod chem_utils;
pub mod chemical_features;
pub mod crippen;
pub mod descriptors;
pub mod descriptors_3d;
pub mod draw;
pub mod e_state;
pub mod enumerate_heterocycles;
pub mod enumerate_stereoisomers;
pub mod feat_finder_cli;
pub mod feat_maps;
pub mod features;
pub mod filter_catalog;
pub mod fraggle;
pub mod fragment_catalog;
pub mod fragment_matcher;
pub mod fragments;
pub mod functional_groups;
pub mod graph_descriptors;
pub mod graphs;
pub mod inchi;
mod init;
pub mod lipinski;
pub mod maccs_keys;
pub mod mol_key;
pub mod mol_standardize;
pub mod mol_surf;
pub mod pharm_2d;
pub mod pharm_3d;
pub mod qed;
pub mod randomize;
pub mod rd_abbreviations;
pub mod rd_chem_reactions;
pub mod rd_chemical_features;
pub mod rd_cip_labeler;
pub mod rd_conformer_parser;
pub mod rd_coord_gen;
pub mod rd_depictor;
pub mod rd_deprotect;
pub mod rd_determine_bonds;
pub mod rd_dist_geom;
pub mod rd_fingerprint_generator;
pub mod rd_fmcs;
pub mod rd_force_field_helpers;
pub mod rd_free_sasa;
pub mod rd_generalized_substruct;
pub mod rd_mhfp_fingerprint;
pub mod rd_mif;
pub mod rd_mmpa;
pub mod rd_mol_align;
pub mod rd_mol_catalog;
pub mod rd_mol_chemical_features;
pub mod rd_mol_descriptors;
pub mod rd_mol_enumerator;
pub mod rd_mol_hash;
pub mod rd_mol_interchange;
pub mod rd_mol_processing;
pub mod rd_mol_transforms;
pub mod rd_partial_charges;
pub mod rd_r_group_decomposition;
pub mod rd_rascal_mces;
pub mod rd_reduced_graphs;
pub mod rd_shape_align;
pub mod rd_shape_helpers;
pub mod rd_sln_parse;
pub mod rd_substruct_library;
pub mod rd_tautomer_query;
pub mod rdchem;
pub mod rdfiltercatalog;
pub mod rdfragcatalog;
pub mod rdinchi;
pub mod rdmolfiles;
pub mod rdmolops;
pub mod rdqueries;
pub mod rdtrajectory;
pub mod recap;
pub mod reduced_graphs;
pub mod registration_hash;
pub mod salt_remover;
pub mod satis;
pub mod scaffolds;
pub mod simple_enum;
pub mod spacial_score;
pub mod subshape;
pub mod template_align;
pub mod torsion_fingerprints;

pub use inchi::{
    InchiLogLevel, MolFromInchiParams, MolToInchiParams, inchi_available, inchi_to_inchi_key,
    mol_from_inchi, mol_from_inchi_with_params, mol_to_inchi, mol_to_inchi_and_aux_info,
    mol_to_inchi_and_aux_info_with_params, mol_to_inchi_key, mol_to_inchi_key_with_options,
    mol_to_inchi_with_params,
};
pub use init::*;
pub use rd_cip_labeler::*;
pub use rd_mol_interchange::*;
pub use rdchem::*;
pub use rdmolfiles::*;
pub use rdmolops::*;
