// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;
use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdfiltercatalog::FilterMatch;
use crate::chem::rdfiltercatalog::matchers::{FilterMatcher, matches_from_ffi};
use crate::chem::rdmolops::mol_list_into_map;
use crate::error::{expect_ok, to_string};
use crate::ffi::filter_catalog::ffi;

/// A named matcher with string properties, such as a PAINS filter with its
/// reference and scope. Mirrors `rdfiltercatalog.FilterCatalogEntry`.
///
/// Catalogs share their entries with the entries they return, and a clone
/// shares its entry too, as Python's objects do; changing an entry that is
/// shared changes a copy, so a catalog's entries never change. It is
/// `Send + Sync`.
pub struct FilterCatalogEntry {
    /// Never null.
    ptr: UniquePtr<ffi::FilterCatalogEntryHandle>,
}

impl FilterCatalogEntry {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::FilterCatalogEntryHandle>) -> FilterCatalogEntry {
        assert!(!ptr.is_null(), "RDKit returned a null entry");
        FilterCatalogEntry { ptr }
    }

    fn raw(&self) -> &ffi::FilterCatalogEntryHandle {
        self.ptr.as_ref().unwrap()
    }

    /// An entry described by `name`, with a copy of `matcher`. Mirrors
    /// `FilterCatalogEntry(name, matcher)`. `Error::Value` for a
    /// `FilterHierarchyMatcher` without a pattern, or a matcher holding one,
    /// whose `IsValid` dereferences null in Python.
    pub fn new(name: &str, matcher: &FilterMatcher) -> Result<FilterCatalogEntry> {
        Ok(FilterCatalogEntry::from_ffi(
            ffi::filter_catalog_entry_new_matcher(name, matcher.raw())?,
        ))
    }

    /// Whether the entry has a valid matcher. Mirrors
    /// `FilterCatalogEntry.IsValid`.
    ///
    /// # Panics
    ///
    /// When a `CustomFilterMatcher`'s `is_valid` panics.
    pub fn is_valid(&self) -> bool {
        expect_ok(ffi::filter_catalog_entry_is_valid(self.raw()))
    }

    /// Mirrors `FilterCatalogEntry.GetDescription`. `Error::Value` for one,
    /// from a pickle, that is not UTF-8 (and likewise for the other
    /// strings).
    pub fn get_description(&self) -> Result<String> {
        to_string(expect_ok(ffi::filter_catalog_entry_get_description(
            self.raw(),
        )))
    }

    /// Mirrors `FilterCatalogEntry.SetDescription`.
    pub fn set_description(&mut self, description: &str) {
        expect_ok(ffi::filter_catalog_entry_set_description(
            self.ptr.pin_mut(),
            description,
        ))
    }

    /// The matches of the entry's matcher, empty when it does not match or
    /// is not valid. Mirrors `FilterCatalogEntry.GetFilterMatches`.
    /// `Error::Runtime` when a `CustomFilterMatcher` panics.
    pub fn get_filter_matches(&self, mol: &Mol) -> Result<Vec<FilterMatch>> {
        Ok(matches_from_ffi(
            &*ffi::filter_catalog_entry_get_filter_matches(self.raw(), mol.as_ffi())?,
        ))
    }

    /// Whether the entry's matcher is valid and matches. Mirrors
    /// `FilterCatalogEntry.HasFilterMatch`. `Error::Runtime` when a
    /// `CustomFilterMatcher` panics.
    pub fn has_filter_match(&self, mol: &Mol) -> Result<bool> {
        Ok(ffi::filter_catalog_entry_has_filter_match(
            self.raw(),
            mol.as_ffi(),
        )?)
    }

    /// The entry as a Boost.Serialization archive. Mirrors
    /// `FilterCatalogEntry.Serialize`. `Error::Value` when a `SmartsMatcher`
    /// in it has no pattern, which Python's dereferences, and for a
    /// `CustomFilterMatcher`, whose class Boost does not know.
    /// `Error::Value` for a pattern a pickle cannot hold (`Mol::to_binary`).
    pub fn serialize(&self) -> Result<Vec<u8>> {
        Ok(ffi::filter_catalog_entry_serialize(self.raw())?)
    }

    /// The properties' names, the description's among them. Mirrors
    /// `FilterCatalogEntry.GetPropList`.
    pub fn get_prop_list(&self) -> Result<Vec<String>> {
        expect_ok(ffi::filter_catalog_entry_get_prop_list(self.raw()))
            .iter()
            .map(|key| {
                key.to_str()
                    .map(str::to_owned)
                    .map_err(|err| crate::Error::Value {
                        message: format!("invalid UTF-8 in a string from RDKit: {err}"),
                    })
            })
            .collect()
    }

    /// Mirrors `FilterCatalogEntry.SetProp`.
    pub fn set_prop(&mut self, key: &str, val: &str) {
        expect_ok(ffi::filter_catalog_entry_set_prop(
            self.ptr.pin_mut(),
            key,
            val,
        ))
    }

    /// Mirrors `FilterCatalogEntry.GetProp`. `Error::Key` for a missing
    /// property.
    pub fn get_prop(&self, key: &str) -> Result<String> {
        to_string(ffi::filter_catalog_entry_get_prop(self.raw(), key)?)
    }

    /// Removes a property, if the entry has it. Mirrors
    /// `FilterCatalogEntry.ClearProp`.
    pub fn clear_prop(&mut self, key: &str) {
        expect_ok(ffi::filter_catalog_entry_clear_prop(
            self.ptr.pin_mut(),
            key,
        ))
    }

    /// Whether two handles share one entry, as `remove_entry_filter_catalog_entry`
    /// compares them. Mirrors Python's `is` on entries.
    pub fn is_same_entry(&self, other: &FilterCatalogEntry) -> bool {
        expect_ok(ffi::filter_catalog_entry_same(self.raw(), other.raw()))
    }
}

/// An entry without a matcher, which is not valid. Mirrors
/// `FilterCatalogEntry()`.
impl Default for FilterCatalogEntry {
    fn default() -> FilterCatalogEntry {
        FilterCatalogEntry::from_ffi(expect_ok(ffi::filter_catalog_entry_new()))
    }
}

/// Shares the entry, which a change then copies.
impl Clone for FilterCatalogEntry {
    fn clone(&self) -> FilterCatalogEntry {
        FilterCatalogEntry::from_ffi(expect_ok(ffi::filter_catalog_entry_share(self.raw())))
    }
}

impl fmt::Debug for FilterCatalogEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("FilterCatalogEntry")
            .field(&self.get_description())
            .finish()
    }
}

/// Sets of RDKit's built-in filter catalogs, as bits. Mirrors
/// `rdfiltercatalog.FilterCatalogParams.FilterCatalogs` (CONVENTIONS.md F1);
/// the constants keep Python's spelling.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct FilterCatalogs(u32);

macro_rules! filter_catalogs {
    ($($name:ident = $value:expr),* $(,)?) => {
        #[allow(non_upper_case_globals)]
        impl FilterCatalogs {
            $(
                #[doc = concat!("Mirrors `FilterCatalogs.", stringify!($name), "`.")]
                pub const $name: FilterCatalogs = FilterCatalogs($value);
            )*

            const NAMES: &[(&str, FilterCatalogs)] =
                &[$((stringify!($name), FilterCatalogs::$name)),*];
        }
    };
}

// In the order filter_catalogs_values returns them; a test checks the values
// against FilterCatalog.h.
filter_catalogs! {
    PAINS_A = 1 << 1,
    PAINS_B = 1 << 2,
    PAINS_C = 1 << 3,
    PAINS = 0b1110,
    BRENK = 1 << 4,
    NIH = 1 << 5,
    ZINC = 1 << 6,
    CHEMBL_Glaxo = 1 << 7,
    CHEMBL_Dundee = 1 << 8,
    CHEMBL_BMS = 1 << 9,
    CHEMBL_SureChEMBL = 1 << 10,
    CHEMBL_MLSMR = 1 << 11,
    CHEMBL_Inpharmatica = 1 << 12,
    CHEMBL_LINT = 1 << 13,
    CHEMBL = 0b11_1111_1000_0000,
    ALL = 0b11_1111_1111_1110,
}

crate::flags::flag_ops!(FilterCatalogs);

#[cfg(test)]
pub(crate) fn filter_catalogs_values_for_tests() -> (Vec<u32>, Vec<u32>) {
    let rust = FilterCatalogs::NAMES
        .iter()
        .map(|(_, catalogs)| catalogs.bits())
        .collect();
    (rust, expect_ok(ffi::filter_catalogs_values()))
}

/// The built-in catalogs a `FilterCatalog` is filled from. Mirrors
/// `rdfiltercatalog.FilterCatalogParams`.
pub struct FilterCatalogParams {
    /// Never null.
    ptr: UniquePtr<ffi::FilterCatalogParams>,
}

impl FilterCatalogParams {
    /// Parameters with the catalogs `catalogs` has. Mirrors
    /// `FilterCatalogParams(catalogs)`.
    pub fn new(catalogs: FilterCatalogs) -> FilterCatalogParams {
        FilterCatalogParams {
            ptr: expect_ok(ffi::filter_catalog_params_from_catalogs(catalogs.bits())),
        }
    }

    /// Adds the catalogs `catalogs` has, returning whether it has any; bits
    /// no constant names are ignored, as RDKit ignores them. Mirrors
    /// `FilterCatalogParams.AddCatalog`.
    pub fn add_catalog(&mut self, catalogs: FilterCatalogs) -> bool {
        expect_ok(ffi::filter_catalog_params_add_catalog(
            self.ptr.pin_mut(),
            catalogs.bits(),
        ))
    }
}

/// Parameters without catalogs. Mirrors `FilterCatalogParams()`.
impl Default for FilterCatalogParams {
    fn default() -> FilterCatalogParams {
        FilterCatalogParams {
            ptr: expect_ok(ffi::filter_catalog_params_new()),
        }
    }
}

impl fmt::Debug for FilterCatalogParams {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FilterCatalogParams")
            .finish_non_exhaustive()
    }
}

fn entries_from_ffi(list: &ffi::FilterCatalogEntryList) -> Vec<FilterCatalogEntry> {
    (0..expect_ok(ffi::filter_catalog_entry_list_size(list)))
        .map(|i| {
            FilterCatalogEntry::from_ffi(expect_ok(ffi::filter_catalog_entry_list_get(list, i)))
        })
        .collect()
}

/// A list of filters to run on molecules. Mirrors
/// `rdfiltercatalog.FilterCatalog`.
///
/// It is `Send + Sync`: its entries are not changed in place, and RDKit runs
/// one catalog's matchers from several threads (CONVENTIONS.md S5). Its
/// matching, like every substructure match, holds a process-wide lock
/// (U4), so threads match one at a time.
pub struct FilterCatalog {
    /// Never null.
    ptr: UniquePtr<ffi::FilterCatalog>,
}

impl FilterCatalog {
    fn raw(&self) -> &ffi::FilterCatalog {
        self.ptr.as_ref().unwrap()
    }

    /// A catalog of the entries of the parameters' built-in catalogs.
    /// Mirrors `FilterCatalog(params)`.
    pub fn new(params: &FilterCatalogParams) -> FilterCatalog {
        FilterCatalog {
            ptr: expect_ok(ffi::filter_catalog_from_params(
                params.ptr.as_ref().unwrap(),
            )),
        }
    }

    /// A catalog of the entries of built-in catalogs. Mirrors
    /// `FilterCatalog(catalogs)`.
    pub fn from_catalogs(catalogs: FilterCatalogs) -> FilterCatalog {
        FilterCatalog {
            ptr: expect_ok(ffi::filter_catalog_from_catalogs(catalogs.bits())),
        }
    }

    /// A catalog from `serialize`'s archive. Mirrors `FilterCatalog(pickle)`.
    /// `Error::Value` for an archive whose matchers RDKit would dereference
    /// null in or recurse through forever (an empty `FilterHierarchyMatcher`,
    /// a cycle), or whose patterns' atoms have monomer info of a type their
    /// class lacks (U10), which Python's reads.
    ///
    /// # Safety
    ///
    /// `pickle` must be bytes `serialize` (this crate's or RDKit's)
    /// returned. RDKit reads it as a Boost.Serialization archive of
    /// polymorphic matchers, which trusts its class and object ids; other
    /// bytes may make it read out of bounds (CONVENTIONS.md U2).
    pub unsafe fn from_pickle(pickle: &[u8]) -> Result<FilterCatalog> {
        Ok(FilterCatalog {
            ptr: ffi::filter_catalog_from_pickle(pickle)?,
        })
    }

    /// The catalog's entries as a Boost.Serialization archive. Mirrors
    /// `FilterCatalog.Serialize`. `Error::Value` when a `SmartsMatcher` in
    /// it has no pattern, which Python's dereferences, and for a
    /// `CustomFilterMatcher`, whose class Boost does not know;
    /// `Error::Runtime` without Boost serialization
    /// (`filter_catalog_can_serialize`).
    /// `Error::Value` for a pattern a pickle cannot hold (`Mol::to_binary`).
    pub fn serialize(&self) -> Result<Vec<u8>> {
        Ok(ffi::filter_catalog_serialize(self.raw())?)
    }

    /// Adds a copy of an entry. Mirrors `FilterCatalog.AddEntry`.
    pub fn add_entry(&mut self, entry: &FilterCatalogEntry) {
        expect_ok(ffi::filter_catalog_add_entry(
            self.ptr.pin_mut(),
            entry.raw(),
        ))
    }

    /// Removes the entry at an index, returning whether there was one; later
    /// entries move down. Mirrors `FilterCatalog.RemoveEntry(idx)`.
    pub fn remove_entry(&mut self, idx: u32) -> bool {
        expect_ok(ffi::filter_catalog_remove_entry(self.ptr.pin_mut(), idx))
    }

    /// Removes an entry the catalog returned, returning whether it found
    /// it: entries are found by identity, so not a copy `add_entry` made from
    /// it, nor one changed since (which copies it). Mirrors
    /// `FilterCatalog.RemoveEntry(entry)`.
    pub fn remove_entry_filter_catalog_entry(&mut self, entry: &FilterCatalogEntry) -> bool {
        expect_ok(ffi::filter_catalog_remove_entry_entry(
            self.ptr.pin_mut(),
            entry.raw(),
        ))
    }

    /// Mirrors `FilterCatalog.GetNumEntries`.
    pub fn get_num_entries(&self) -> u32 {
        expect_ok(ffi::filter_catalog_get_num_entries(self.raw()))
    }

    /// The entry at an index, shared with the catalog. Mirrors
    /// `FilterCatalog.GetEntry`. `Error::Index` for an index past the
    /// entries.
    pub fn get_entry(&self, idx: u32) -> Result<FilterCatalogEntry> {
        Ok(FilterCatalogEntry::from_ffi(ffi::filter_catalog_get_entry(
            self.raw(),
            idx,
        )?))
    }

    /// Mirrors `FilterCatalog.GetEntryWithIdx`, as `get_entry`.
    pub fn get_entry_with_idx(&self, idx: u32) -> Result<FilterCatalogEntry> {
        self.get_entry(idx)
    }

    /// Whether any valid entry matches. Mirrors `FilterCatalog.HasMatch`.
    /// `Error::Runtime` when a `CustomFilterMatcher` panics (and likewise for
    /// the other matching methods).
    pub fn has_match(&self, mol: &Mol) -> Result<bool> {
        Ok(ffi::filter_catalog_has_match(self.raw(), mol.as_ffi())?)
    }

    /// The first entry that matches, or `None`. Mirrors
    /// `FilterCatalog.GetFirstMatch`.
    pub fn get_first_match(&self, mol: &Mol) -> Result<Option<FilterCatalogEntry>> {
        let entry = ffi::filter_catalog_get_first_match(self.raw(), mol.as_ffi())?;
        Ok((!entry.is_null()).then(|| FilterCatalogEntry::from_ffi(entry)))
    }

    /// The entries that match. Mirrors `FilterCatalog.GetMatches`.
    pub fn get_matches(&self, mol: &Mol) -> Result<Vec<FilterCatalogEntry>> {
        Ok(entries_from_ffi(&*ffi::filter_catalog_get_matches(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// The matches of every entry that matches. Mirrors
    /// `FilterCatalog.GetFilterMatches`.
    pub fn get_filter_matches(&self, mol: &Mol) -> Result<Vec<FilterMatch>> {
        Ok(matches_from_ffi(&*ffi::filter_catalog_get_filter_matches(
            self.raw(),
            mol.as_ffi(),
        )?))
    }
}

/// A catalog without entries. Mirrors `FilterCatalog()`.
impl Default for FilterCatalog {
    fn default() -> FilterCatalog {
        FilterCatalog {
            ptr: expect_ok(ffi::filter_catalog_new()),
        }
    }
}

impl fmt::Debug for FilterCatalog {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FilterCatalog")
            .field("num_entries", &self.get_num_entries())
            .finish()
    }
}

/// Whether catalogs can be pickled, which needs an RDKit built with Boost
/// serialization. Mirrors `rdfiltercatalog.FilterCatalogCanSerialize`.
pub fn filter_catalog_can_serialize() -> bool {
    expect_ok(ffi::filter_catalog_can_serialize())
}

/// The entries of a catalog that each SMILES matches, on one thread.
/// Mirrors `rdfiltercatalog.RunFilterCatalog(filterCatalog, smiles)`.
pub fn run_filter_catalog<S: AsRef<str>>(
    filter_catalog: &FilterCatalog,
    smiles: &[S],
) -> Result<Vec<Vec<FilterCatalogEntry>>> {
    run_filter_catalog_with_num_threads(filter_catalog, smiles, 1)
}

/// `num_threads` threads, or with zero or less, the hardware's concurrency
/// plus `num_threads` (and at least one). A SMILES that does not parse gets
/// one entry, "no valid RDKit molecule", without a matcher; one that fails
/// sanitization is the call's error, as in Python. The threads are the
/// crate's (RDKit's own count threads in 32 bits), and only parse in
/// parallel: every substructure match holds the crate's matching lock, so
/// matches take turns (CONVENTIONS.md U4). Mirrors
/// `rdfiltercatalog.RunFilterCatalog(filterCatalog, smiles, numThreads)`.
pub fn run_filter_catalog_with_num_threads<S: AsRef<str>>(
    filter_catalog: &FilterCatalog,
    smiles: &[S],
    num_threads: i32,
) -> Result<Vec<Vec<FilterCatalogEntry>>> {
    let smiles: Vec<&str> = smiles.iter().map(AsRef::as_ref).collect();
    let list = ffi::run_filter_catalog(filter_catalog.raw(), &smiles, num_threads)?;
    let mut entries = entries_from_ffi(&list).into_iter();
    let mut start = 0;
    Ok(expect_ok(ffi::filter_catalog_entry_list_row_ends(&list))
        .into_iter()
        .map(|end| {
            let row = entries.by_ref().take(end - start).collect();
            start = end;
            row
        })
        .collect())
}

/// RDKit's functional groups as a catalog of hierarchies, one per family of
/// groups. Mirrors `rdfiltercatalog.GetFunctionalGroupHierarchy`, which
/// returns RDKit's global catalog; this one shares its entries.
pub fn get_functional_group_hierarchy() -> FilterCatalog {
    FilterCatalog {
        ptr: expect_ok(ffi::get_functional_group_hierarchy()),
    }
}

/// RDKit's functional groups by name ("Halogen.Bromine.Aliphatic"), as
/// copies of their query molecules. Mirrors
/// `rdfiltercatalog.GetFlattenedFunctionalGroupHierarchy`.
pub fn get_flattened_functional_group_hierarchy() -> BTreeMap<String, Mol> {
    get_flattened_functional_group_hierarchy_with_normalized(false)
}

/// With `normalized`, the names are lower case. Mirrors
/// `rdfiltercatalog.GetFlattenedFunctionalGroupHierarchy(normalized)`.
pub fn get_flattened_functional_group_hierarchy_with_normalized(
    normalized: bool,
) -> BTreeMap<String, Mol> {
    mol_list_into_map(expect_ok(ffi::get_flattened_functional_group_hierarchy(
        normalized,
    )))
    .expect("RDKit's functional group names are UTF-8")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_values_match_rdkit() {
        let (rust, cxx) = filter_catalogs_values_for_tests();
        assert_eq!(rust, cxx);
        assert_eq!(
            format!("{:?}", FilterCatalogs::PAINS_A | FilterCatalogs::PAINS_B),
            "FilterCatalogs(0x6)"
        );
        assert_eq!(
            format!("{:?}", FilterCatalogs::PAINS),
            "FilterCatalogs::PAINS"
        );
    }
}
