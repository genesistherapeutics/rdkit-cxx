// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Private cxx bridges, one per wrapper header (CONVENTIONS.md B3).

pub(crate) mod adjust_query;
pub(crate) mod alignment;
pub(crate) mod atom;
pub(crate) mod bond;
pub(crate) mod chem_reactions;
pub(crate) mod chirality;
pub(crate) mod cip_labeler;
pub(crate) mod conformer;
pub(crate) mod conformer_parser;
pub(crate) mod coord_gen;
pub(crate) mod data_structs;
pub(crate) mod dclv;
pub(crate) mod depictor;
pub(crate) mod deprotect;
pub(crate) mod descriptors;
pub(crate) mod descriptors_3d;
pub(crate) mod determine_bonds;
pub(crate) mod distance_geometry;
pub(crate) mod embedder;
pub(crate) mod enumerate;
pub(crate) mod errors;
pub(crate) mod fingerprint_generator;
pub(crate) mod fingerprints;
pub(crate) mod force_field;
pub(crate) mod force_field_helpers;
pub(crate) mod fpb;
pub(crate) mod geometry;
pub(crate) mod graph;
pub(crate) mod grid;
pub(crate) mod mhfp;
pub(crate) mod mmpa;
pub(crate) mod mol;
pub(crate) mod mol_align;
pub(crate) mod mol_bundle;
pub(crate) mod mol_hash;
pub(crate) mod mol_transforms;
pub(crate) mod molfiles;
pub(crate) mod molops;
pub(crate) mod monomer_info;
pub(crate) mod partial_charges;
pub(crate) mod perception;
pub(crate) mod periodic_table;
pub(crate) mod png;
pub(crate) mod properties;
pub(crate) mod props;
pub(crate) mod query;
pub(crate) mod rd_base;
pub(crate) mod rdqueries;
pub(crate) mod reduced_graphs;
pub(crate) mod resonance;
pub(crate) mod ring_info;
pub(crate) mod shape_helpers;
pub(crate) mod sln_parse;
pub(crate) mod standardize;
pub(crate) mod standardize_charge;
pub(crate) mod standardize_fragment;
pub(crate) mod standardize_normalize;
pub(crate) mod standardize_pipeline;
pub(crate) mod standardize_tautomer;
pub(crate) mod standardize_validate;
pub(crate) mod stereo;
pub(crate) mod substance_group;
pub(crate) mod substruct;
pub(crate) mod suppliers;
pub(crate) mod transforms;
pub(crate) mod writers;
