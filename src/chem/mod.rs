// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem`, which re-exports `rdchem`, `rdmolfiles`, and
//! `rdmolops` (and `rdCIPLabeler` and `rdMolInterchange`), and defines
//! Python functions of its own (`FindMolChiralCenters`), ported here; its
//! other submodules are reached by path, as in Python.

pub mod atom_pairs;
pub mod draw;
pub mod feat_finder_cli;
pub mod feat_maps;
pub mod features;
mod init;
pub mod lipinski;
pub mod maccs_keys;
pub mod mol_standardize;
pub mod pharm_2d;
pub mod pharm_3d;
pub(crate) mod py_math;
pub mod rd_abbreviations;
pub mod rd_chem_draw;
pub mod rd_chem_reactions;
pub mod rd_chemical_features;
pub mod rd_cip_labeler;
pub mod rd_conformer_parser;
pub mod rd_coord_gen;
pub mod rd_depictor;
pub mod rd_deprotect;
pub mod rd_determine_bonds;
pub mod rd_dist_geom;
pub mod rd_enumerate_stereoisomers;
pub mod rd_fingerprint_generator;
pub mod rd_fmcs;
pub mod rd_force_field_helpers;
pub mod rd_free_sasa;
pub mod rd_gaussian_shape;
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
pub mod rd_synthon_space_search;
pub mod rd_tautomer_query;
pub mod rdchem;
pub mod rdfiltercatalog;
pub mod rdfragcatalog;
pub mod rdinchi;
pub mod rdmolfiles;
pub mod rdmolops;
pub mod rdqueries;
pub mod registration_hash;
pub mod scaffolds;
pub mod subshape;
pub mod torsion_fingerprints;

pub use init::*;
pub use rd_cip_labeler::*;
pub use rd_mol_interchange::*;
pub use rdchem::*;
pub use rdmolfiles::*;
pub use rdmolops::*;
