/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdfiltercatalog: catalogs of substructure filters (PAINS, Brenk, and the
// others RDKit ships), the matchers they are built from, and matchers whose
// work is a Rust callback.
//
// Matchers and catalog entries are boost::shared_ptrs that RDKit shares
// between catalogs, entries, and match results. A matcher reachable from
// more than one place is never changed: Rust's concrete matchers own theirs
// alone (RDKit copies a matcher into whatever holds it), and an entry is
// copied before it is changed unless its handle is its only owner.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
#include <GraphMol/FilterCatalog/FilterMatcherBase.h>
#include <GraphMol/FilterCatalog/FilterMatchers.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Defined by the bridges (src/ffi/filter_catalog.rs, chirality.rs).
struct FilterMatcherCallback;
struct IndexPair;

using FilterCatalog = RDKit::FilterCatalog;
using FilterCatalogParams = RDKit::FilterCatalogParams;

// Any matcher.
struct FilterMatcherHandle {
  boost::shared_ptr<RDKit::FilterMatcherBase> ptr;
};

// An entry, maybe shared with catalogs and other handles; changed only
// through entry_mut, which copies it first unless this handle is its only
// owner.
struct FilterCatalogEntryHandle {
  boost::shared_ptr<RDKit::FilterCatalogEntry> ptr;
};

// Matches crossing the bridge, either way.
struct FilterMatchList {
  std::vector<RDKit::FilterMatch> matches;
};

// Copies of matchers, for ExclusionList.SetExclusionPatterns.
struct FilterMatcherList {
  std::vector<boost::shared_ptr<RDKit::FilterMatcherBase>> matchers;
};

// Entries a call returns, in rows (RunFilterCatalog) or one row.
struct FilterCatalogEntryList {
  std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry>> entries;
  std::vector<std::size_t> row_ends;
};

// Matchers.
std::unique_ptr<FilterMatcherHandle> filter_matcher_new(
    rust::Str name, rust::Box<FilterMatcherCallback> callback);
std::unique_ptr<FilterMatcherHandle> filter_matcher_copy(
    const FilterMatcherHandle &matcher);
std::unique_ptr<std::string> filter_matcher_get_name(
    const FilterMatcherHandle &matcher);
bool filter_matcher_is_valid(const FilterMatcherHandle &matcher);
bool filter_matcher_has_match(const FilterMatcherHandle &matcher,
                              const RWMol &mol);
std::unique_ptr<FilterMatchList> filter_matcher_get_matches(
    const FilterMatcherHandle &matcher, const RWMol &mol);

std::unique_ptr<FilterMatcherHandle> smarts_matcher_new(rust::Str name);
std::unique_ptr<FilterMatcherHandle> smarts_matcher_from_mol(const RWMol &mol);
std::unique_ptr<FilterMatcherHandle> smarts_matcher_new_mol(
    rust::Str name, const RWMol &mol, std::uint32_t min_count,
    std::uint32_t max_count);
std::unique_ptr<FilterMatcherHandle> smarts_matcher_new_smarts(
    rust::Str name, rust::Str smarts, std::uint32_t min_count,
    std::uint32_t max_count);
// A copy of the pattern, or null when there is none.
std::unique_ptr<RWMol> smarts_matcher_get_pattern(
    const FilterMatcherHandle &matcher);
void smarts_matcher_set_pattern_mol(FilterMatcherHandle &matcher,
                                    const RWMol &mol);
void smarts_matcher_set_pattern_smarts(FilterMatcherHandle &matcher,
                                       rust::Str smarts);
std::uint32_t smarts_matcher_get_min_count(const FilterMatcherHandle &matcher);
void smarts_matcher_set_min_count(FilterMatcherHandle &matcher,
                                  std::uint32_t count);
std::uint32_t smarts_matcher_get_max_count(const FilterMatcherHandle &matcher);
void smarts_matcher_set_max_count(FilterMatcherHandle &matcher,
                                  std::uint32_t count);

std::unique_ptr<FilterMatcherHandle> exclusion_list_new();
void exclusion_list_add_pattern(FilterMatcherHandle &list,
                                const FilterMatcherHandle &base);
std::unique_ptr<FilterMatcherList> filter_matcher_list_new();
void filter_matcher_list_push(FilterMatcherList &list,
                              const FilterMatcherHandle &matcher);
void exclusion_list_set_exclusion_patterns(FilterMatcherHandle &list,
                                           const FilterMatcherList &patterns);

std::unique_ptr<FilterMatcherHandle> filter_hierarchy_matcher_new();
std::unique_ptr<FilterMatcherHandle> filter_hierarchy_matcher_from_matcher(
    const FilterMatcherHandle &matcher);
void filter_hierarchy_matcher_set_pattern(FilterMatcherHandle &node,
                                          const FilterMatcherHandle &matcher);
// The child, as it is in the node's tree.
std::unique_ptr<FilterMatcherHandle> filter_hierarchy_matcher_add_child(
    FilterMatcherHandle &node, const FilterMatcherHandle &hierarchy);

std::unique_ptr<FilterMatcherHandle> filter_match_ops_and(
    const FilterMatcherHandle &arg1, const FilterMatcherHandle &arg2);
std::unique_ptr<FilterMatcherHandle> filter_match_ops_or(
    const FilterMatcherHandle &arg1, const FilterMatcherHandle &arg2);
std::unique_ptr<FilterMatcherHandle> filter_match_ops_not(
    const FilterMatcherHandle &arg1);

std::size_t filter_match_list_size(const FilterMatchList &list);
std::unique_ptr<FilterMatcherHandle> filter_match_list_matcher(
    const FilterMatchList &list, std::size_t i);
rust::Vec<IndexPair> filter_match_list_atom_pairs(const FilterMatchList &list,
                                                  std::size_t i);
void filter_match_list_push(FilterMatchList &list,
                            const FilterMatcherHandle &matcher,
                            rust::Slice<const IndexPair> atom_pairs);

// Entries.
std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_new();
std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_new_matcher(
    rust::Str name, const FilterMatcherHandle &matcher);
std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_share(
    const FilterCatalogEntryHandle &entry);
bool filter_catalog_entry_same(const FilterCatalogEntryHandle &a,
                               const FilterCatalogEntryHandle &b);
bool filter_catalog_entry_is_valid(const FilterCatalogEntryHandle &entry);
std::unique_ptr<std::string> filter_catalog_entry_get_description(
    const FilterCatalogEntryHandle &entry);
void filter_catalog_entry_set_description(FilterCatalogEntryHandle &entry,
                                          rust::Str description);
std::unique_ptr<FilterMatchList> filter_catalog_entry_get_filter_matches(
    const FilterCatalogEntryHandle &entry, const RWMol &mol);
bool filter_catalog_entry_has_filter_match(
    const FilterCatalogEntryHandle &entry, const RWMol &mol);
rust::Vec<std::uint8_t> filter_catalog_entry_serialize(
    const FilterCatalogEntryHandle &entry);
std::unique_ptr<std::vector<std::string>> filter_catalog_entry_get_prop_list(
    const FilterCatalogEntryHandle &entry);
void filter_catalog_entry_set_prop(FilterCatalogEntryHandle &entry,
                                   rust::Str key, rust::Str val);
std::unique_ptr<std::string> filter_catalog_entry_get_prop(
    const FilterCatalogEntryHandle &entry, rust::Str key);
void filter_catalog_entry_clear_prop(FilterCatalogEntryHandle &entry,
                                     rust::Str key);

std::size_t filter_catalog_entry_list_size(const FilterCatalogEntryList &list);
std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_list_get(
    const FilterCatalogEntryList &list, std::size_t i);
rust::Vec<std::size_t> filter_catalog_entry_list_row_ends(
    const FilterCatalogEntryList &list);

// Parameters and catalogs. Catalog sets are FilterCatalogs bits.
std::unique_ptr<FilterCatalogParams> filter_catalog_params_new();
std::unique_ptr<FilterCatalogParams> filter_catalog_params_from_catalogs(
    std::uint32_t catalogs);
bool filter_catalog_params_add_catalog(FilterCatalogParams &params,
                                       std::uint32_t catalogs);
// FilterCatalogs' values in a fixed order, for a test that pins them.
rust::Vec<std::uint32_t> filter_catalogs_values();

std::unique_ptr<FilterCatalog> filter_catalog_new();
std::unique_ptr<FilterCatalog> filter_catalog_from_params(
    const FilterCatalogParams &params);
std::unique_ptr<FilterCatalog> filter_catalog_from_catalogs(
    std::uint32_t catalogs);
// The caller promises bytes RDKit wrote (CONVENTIONS.md U2).
std::unique_ptr<FilterCatalog> filter_catalog_from_pickle(
    rust::Slice<const std::uint8_t> pickle);
rust::Vec<std::uint8_t> filter_catalog_serialize(const FilterCatalog &catalog);
void filter_catalog_add_entry(FilterCatalog &catalog,
                              const FilterCatalogEntryHandle &entry);
bool filter_catalog_remove_entry(FilterCatalog &catalog, std::uint32_t idx);
bool filter_catalog_remove_entry_entry(FilterCatalog &catalog,
                                       const FilterCatalogEntryHandle &entry);
std::uint32_t filter_catalog_get_num_entries(const FilterCatalog &catalog);
std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_get_entry(
    const FilterCatalog &catalog, std::uint32_t idx);
bool filter_catalog_has_match(const FilterCatalog &catalog, const RWMol &mol);
// Null when nothing matches.
std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_get_first_match(
    const FilterCatalog &catalog, const RWMol &mol);
std::unique_ptr<FilterCatalogEntryList> filter_catalog_get_matches(
    const FilterCatalog &catalog, const RWMol &mol);
std::unique_ptr<FilterMatchList> filter_catalog_get_filter_matches(
    const FilterCatalog &catalog, const RWMol &mol);

bool filter_catalog_can_serialize();
std::unique_ptr<FilterCatalogEntryList> run_filter_catalog(
    const FilterCatalog &catalog, rust::Slice<const rust::Str> smiles,
    std::int32_t num_threads);
std::unique_ptr<FilterCatalog> get_functional_group_hierarchy();
std::unique_ptr<MolList> get_flattened_functional_group_hierarchy(
    bool normalized);

}  // namespace rdkit_cxx
