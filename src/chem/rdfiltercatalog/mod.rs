// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdfiltercatalog`: catalogs of substructure filters,
//! such as RDKit's built-in PAINS, Brenk, NIH, ZINC, and ChEMBL filters, and
//! the matchers they are made of.
//!
//! A `FilterCatalog` holds `FilterCatalogEntry`s, each a named matcher with
//! properties. Matchers are a `SmartsMatcher`, the combinations in
//! `filter_match_ops`, an `ExclusionList`, a tree of them in a
//! `FilterHierarchyMatcher`, or a `CustomFilterMatcher` written in Rust;
//! each dereferences to `FilterMatcher`.
//!
//! RDKit shares matchers and entries between catalogs, entries, and match
//! results. Here a matcher or entry that is shared is never changed: each
//! concrete matcher owns its own, every matcher another holds is a copy,
//! and changing a shared entry changes a copy of it. Catalog pickles are
//! Boost.Serialization archives, so `FilterCatalog::from_pickle` is an
//! `unsafe fn` (CONVENTIONS.md U2).

mod catalog;
pub mod filter_match_ops;
mod matchers;

pub use catalog::*;
pub use matchers::{
    CustomFilterMatcher, ExclusionList, FilterHierarchyMatcher, FilterHierarchyMatcherMut,
    FilterMatch, FilterMatcher, SmartsMatcher, SmartsMatcherParams,
};
