// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::marker::PhantomData;
use std::ops::Deref;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::filter_catalog::{FilterMatcherCallback, ffi};

/// What a matcher written in Rust does. Mirrors subclassing
/// `rdfiltercatalog.FilterMatcher` (`PythonFilterMatcher`) in Python and
/// overriding `IsValid`, `HasMatch`, and `GetMatches`.
///
/// RDKit's catalogs share their matchers between threads, `run_filter_catalog`
/// among them, so a matcher is `Send + Sync`. Each method gets its own copy
/// of the molecule, and runs without the lock matching holds (CONVENTIONS.md
/// U4), so it may itself call functions that match. A panic is caught where
/// RDKit calls the method: the call that made it returns `Error::Runtime`,
/// or, for `is_valid` called by `FilterMatcher::is_valid`, panics.
pub trait CustomFilterMatcher: Send + Sync + 'static {
    /// Whether the matcher can match. Mirrors `FilterMatcher.IsValid`.
    fn is_valid(&self) -> bool {
        true
    }

    /// Mirrors `FilterMatcher.HasMatch`.
    fn has_match(&self, mol: &Mol) -> bool;

    /// Whether the molecule matches, appending the matches to `matches`
    /// (which count only when it does). Mirrors `FilterMatcher.GetMatches`,
    /// which Python's subclasses define as `GetMatches(mol, vect)`.
    fn get_matches(&self, mol: &Mol, matches: &mut Vec<FilterMatch>) -> bool;
}

impl<M: CustomFilterMatcher + ?Sized> CustomFilterMatcher for Box<M> {
    fn is_valid(&self) -> bool {
        (**self).is_valid()
    }

    fn has_match(&self, mol: &Mol) -> bool {
        (**self).has_match(mol)
    }

    fn get_matches(&self, mol: &Mol, matches: &mut Vec<FilterMatch>) -> bool {
        (**self).get_matches(mol, matches)
    }
}

/// Any matcher: one of RDKit's, or one made from a `CustomFilterMatcher`.
/// Mirrors `rdfiltercatalog.FilterMatcher` (`FilterMatcherBase`).
///
/// The concrete matchers (`SmartsMatcher`, `ExclusionList`,
/// `FilterHierarchyMatcher`, and `filter_match_ops`) dereference to it, as
/// they subclass it in Python. It is `Send + Sync`. A clone is RDKit's copy
/// of the matcher, sharing what the matcher holds, which no one changes.
pub struct FilterMatcher {
    /// Never null.
    ptr: UniquePtr<ffi::FilterMatcherHandle>,
}

impl FilterMatcher {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::FilterMatcherHandle>) -> FilterMatcher {
        assert!(!ptr.is_null(), "RDKit returned a null matcher");
        FilterMatcher { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::FilterMatcherHandle {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> std::pin::Pin<&mut ffi::FilterMatcherHandle> {
        self.ptr.pin_mut()
    }

    /// A matcher named `name` that runs `matcher`. Python subclasses the
    /// class and overrides its methods; here they are `matcher`'s. Mirrors
    /// `FilterMatcher(name)`.
    pub fn new<M: CustomFilterMatcher>(name: &str, matcher: M) -> FilterMatcher {
        let callback = Box::new(FilterMatcherCallback(Box::new(matcher)));
        FilterMatcher::from_ffi(expect_ok(ffi::filter_matcher_new(name, callback)))
    }

    /// The matcher's name; composite matchers name their parts. Mirrors
    /// `FilterMatcher.GetName` and `__str__`. `Error::Value` for a name, from
    /// a pickle, that is not UTF-8.
    pub fn get_name(&self) -> Result<String> {
        to_string(expect_ok(ffi::filter_matcher_get_name(self.raw())))
    }

    /// Whether the matcher can match: a `SmartsMatcher` needs a pattern,
    /// and a composite valid parts. Mirrors `FilterMatcher.IsValid`. False
    /// for a `FilterHierarchyMatcher` without a pattern, whose `IsValid`
    /// dereferences null in Python.
    ///
    /// # Panics
    ///
    /// When a `CustomFilterMatcher`'s `is_valid` panics.
    pub fn is_valid(&self) -> bool {
        expect_ok(ffi::filter_matcher_is_valid(self.raw()))
    }

    /// Mirrors `FilterMatcher.HasMatch`. `Error::Runtime` for an invalid
    /// matcher, as in Python, and when a `CustomFilterMatcher` panics;
    /// `Error::Value` for a `FilterHierarchyMatcher` without a pattern,
    /// where Python's dereferences null.
    pub fn has_match(&self, mol: &Mol) -> Result<bool> {
        Ok(ffi::filter_matcher_has_match(self.raw(), mol.as_ffi())?)
    }

    /// The matches of the matcher's parts, empty when it does not match.
    /// Mirrors `FilterMatcher.GetMatches`. Errors as for `has_match`.
    pub fn get_matches(&self, mol: &Mol) -> Result<Vec<FilterMatch>> {
        Ok(matches_from_ffi(&*ffi::filter_matcher_get_matches(
            self.raw(),
            mol.as_ffi(),
        )?))
    }
}

impl Clone for FilterMatcher {
    fn clone(&self) -> FilterMatcher {
        FilterMatcher::from_ffi(expect_ok(ffi::filter_matcher_copy(self.raw())))
    }
}

impl fmt::Debug for FilterMatcher {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("FilterMatcher")
            .field(&self.get_name())
            .finish()
    }
}

/// The matcher's name, as Python's `str`; invalid UTF-8 is replaced.
impl fmt::Display for FilterMatcher {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get_name() {
            Ok(name) => f.write_str(&name),
            Err(_) => f.write_str("<a name that is not UTF-8>"),
        }
    }
}

pub(crate) fn matches_from_ffi(list: &ffi::FilterMatchList) -> Vec<FilterMatch> {
    (0..expect_ok(ffi::filter_match_list_size(list)))
        .map(|i| FilterMatch {
            filter_match: FilterMatcher::from_ffi(expect_ok(ffi::filter_match_list_matcher(
                list, i,
            ))),
            atom_pairs: expect_ok(ffi::filter_match_list_atom_pairs(list, i))
                .into_iter()
                .map(|pair| (pair.first, pair.second))
                .collect(),
        })
        .collect()
}

/// A matcher that matched, and the (query atom, molecule atom) pairs of its
/// match. Mirrors `rdfiltercatalog.FilterMatch`.
///
/// Matches from `filter_match_ops::Not` and `ExclusionList` have no atom
/// pairs, since what they match is an absence.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct FilterMatch {
    /// Mirrors `FilterMatch.filterMatch`.
    pub filter_match: FilterMatcher,
    /// Mirrors `FilterMatch.atomPairs`.
    pub atom_pairs: Vec<(u32, u32)>,
}

impl FilterMatch {
    /// A match of a copy of `filter`. Mirrors `FilterMatch(filter,
    /// atomPairs)`, which holds the matcher itself.
    pub fn new(filter: &FilterMatcher, atom_pairs: &[(u32, u32)]) -> FilterMatch {
        FilterMatch {
            filter_match: filter.clone(),
            atom_pairs: atom_pairs.to_vec(),
        }
    }
}

/// Implements what the concrete matchers share: `Deref` to `FilterMatcher`,
/// `Clone`, `Debug`, and conversion into a `FilterMatcher`.
macro_rules! matcher_type {
    ($name:ident) => {
        impl Deref for $name {
            type Target = FilterMatcher;

            fn deref(&self) -> &FilterMatcher {
                &self.matcher
            }
        }

        impl Clone for $name {
            fn clone(&self) -> $name {
                $name {
                    matcher: self.matcher.clone(),
                }
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_tuple(stringify!($name))
                    .field(&self.get_name())
                    .finish()
            }
        }

        impl From<$name> for FilterMatcher {
            fn from(matcher: $name) -> FilterMatcher {
                matcher.matcher
            }
        }
    };
}

pub(crate) use matcher_type;

/// The `minCount` and `maxCount` of `SmartsMatcher(name, mol, minCount,
/// maxCount)` and `SmartsMatcher(name, smarts, minCount, maxCount)`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmartsMatcherParams {
    /// The fewest matches that make the matcher match.
    pub min_count: u32,
    /// The most matches that make the matcher match; `u32::MAX` for any
    /// number.
    pub max_count: u32,
}

impl Default for SmartsMatcherParams {
    fn default() -> SmartsMatcherParams {
        SmartsMatcherParams {
            min_count: 1,
            max_count: u32::MAX,
        }
    }
}

/// A substructure pattern that matches a molecule with between a minimum
/// and a maximum number of matches. Mirrors `rdfiltercatalog.SmartsMatcher`.
pub struct SmartsMatcher {
    matcher: FilterMatcher,
}

matcher_type!(SmartsMatcher);

impl SmartsMatcher {
    fn wrap(ptr: UniquePtr<ffi::FilterMatcherHandle>) -> SmartsMatcher {
        SmartsMatcher {
            matcher: FilterMatcher::from_ffi(ptr),
        }
    }

    /// A matcher without a pattern, which is not valid until it has one, with
    /// counts 0 and `u32::MAX`. Mirrors `SmartsMatcher(name)`.
    pub fn new(name: &str) -> SmartsMatcher {
        SmartsMatcher::wrap(expect_ok(ffi::smarts_matcher_new(name)))
    }

    /// A matcher of a copy of a query molecule, named "Unnamed
    /// SmartsMatcher", with the default counts. Mirrors
    /// `SmartsMatcher(rhs)`.
    pub fn from_rhs(rhs: &Mol) -> SmartsMatcher {
        SmartsMatcher::wrap(expect_ok(ffi::smarts_matcher_from_mol(rhs.as_ffi())))
    }

    /// A matcher of a copy of a query molecule. Mirrors
    /// `SmartsMatcher(name, mol)`.
    pub fn from_mol(name: &str, mol: &Mol) -> SmartsMatcher {
        SmartsMatcher::from_mol_with_params(name, mol, &SmartsMatcherParams::default())
    }

    /// Mirrors `SmartsMatcher(name, mol, minCount, maxCount)`.
    pub fn from_mol_with_params(
        name: &str,
        mol: &Mol,
        params: &SmartsMatcherParams,
    ) -> SmartsMatcher {
        SmartsMatcher::wrap(expect_ok(ffi::smarts_matcher_new_mol(
            name,
            mol.as_ffi(),
            params.min_count,
            params.max_count,
        )))
    }

    /// A matcher of a SMARTS, with hydrogens merged; one that does not parse
    /// gives a matcher that is not valid, as in Python. Mirrors
    /// `SmartsMatcher(name, smarts)`.
    pub fn from_smarts(name: &str, smarts: &str) -> SmartsMatcher {
        SmartsMatcher::from_smarts_with_params(name, smarts, &SmartsMatcherParams::default())
    }

    /// Mirrors `SmartsMatcher(name, smarts, minCount, maxCount)`.
    pub fn from_smarts_with_params(
        name: &str,
        smarts: &str,
        params: &SmartsMatcherParams,
    ) -> SmartsMatcher {
        SmartsMatcher::wrap(expect_ok(ffi::smarts_matcher_new_smarts(
            name,
            smarts,
            params.min_count,
            params.max_count,
        )))
    }

    /// A copy of the pattern, or `None` without one. Mirrors
    /// `SmartsMatcher.GetPattern`, which returns the pattern itself.
    pub fn get_pattern(&self) -> Option<Mol> {
        let mol = expect_ok(ffi::smarts_matcher_get_pattern(self.raw()));
        (!mol.is_null()).then(|| Mol::from_ffi(mol))
    }

    /// Matches a copy of a query molecule. Mirrors
    /// `SmartsMatcher.SetPattern(pat: Mol)`.
    pub fn set_pattern(&mut self, pat: &Mol) {
        expect_ok(ffi::smarts_matcher_set_pattern_mol(
            self.matcher.raw_mut(),
            pat.as_ffi(),
        ))
    }

    /// Matches a SMARTS; one that does not parse leaves no pattern. The
    /// counts stay. Mirrors `SmartsMatcher.SetPattern(pat: str)`.
    pub fn set_pattern_str(&mut self, pat: &str) {
        expect_ok(ffi::smarts_matcher_set_pattern_smarts(
            self.matcher.raw_mut(),
            pat,
        ))
    }

    /// Mirrors `SmartsMatcher.GetMinCount`.
    pub fn get_min_count(&self) -> u32 {
        expect_ok(ffi::smarts_matcher_get_min_count(self.raw()))
    }

    /// Mirrors `SmartsMatcher.SetMinCount`.
    pub fn set_min_count(&mut self, count: u32) {
        expect_ok(ffi::smarts_matcher_set_min_count(
            self.matcher.raw_mut(),
            count,
        ))
    }

    /// `u32::MAX` for any number. Mirrors `SmartsMatcher.GetMaxCount`.
    pub fn get_max_count(&self) -> u32 {
        expect_ok(ffi::smarts_matcher_get_max_count(self.raw()))
    }

    /// Mirrors `SmartsMatcher.SetMaxCount`.
    pub fn set_max_count(&mut self, count: u32) {
        expect_ok(ffi::smarts_matcher_set_max_count(
            self.matcher.raw_mut(),
            count,
        ))
    }
}

/// Matches molecules that none of its patterns match. Mirrors
/// `rdfiltercatalog.ExclusionList`.
pub struct ExclusionList {
    matcher: FilterMatcher,
}

matcher_type!(ExclusionList);

/// A list without patterns, which matches every molecule. Mirrors
/// `ExclusionList()`.
impl Default for ExclusionList {
    fn default() -> ExclusionList {
        ExclusionList {
            matcher: FilterMatcher::from_ffi(expect_ok(ffi::exclusion_list_new())),
        }
    }
}

impl ExclusionList {
    /// Adds a copy of a pattern. Mirrors `ExclusionList.AddPattern`.
    /// `Error::Runtime` for an invalid one, as in Python; `Error::Value` for
    /// a `FilterHierarchyMatcher` without a pattern, or one that holds one,
    /// whose `IsValid` dereferences null in Python (and likewise wherever a
    /// matcher takes a copy of another).
    pub fn add_pattern(&mut self, base: &FilterMatcher) -> Result<()> {
        Ok(ffi::exclusion_list_add_pattern(
            self.matcher.raw_mut(),
            base.raw(),
        )?)
    }

    /// Replaces the patterns with copies of these. Mirrors
    /// `ExclusionList.SetExclusionPatterns`, which, unlike `add_pattern`,
    /// takes invalid ones.
    pub fn set_exclusion_patterns(&mut self, list: &[&FilterMatcher]) -> Result<()> {
        let mut patterns = expect_ok(ffi::filter_matcher_list_new());
        for pattern in list {
            ffi::filter_matcher_list_push(patterns.pin_mut(), pattern.raw())?;
        }
        Ok(ffi::exclusion_list_set_exclusion_patterns(
            self.matcher.raw_mut(),
            &patterns,
        )?)
    }
}

/// A tree of matchers: a node matches when its pattern does, with the
/// matches of its matching children if any match, else its own. Mirrors
/// `rdfiltercatalog.FilterHierarchyMatcher`.
///
/// A node made by `default` has no pattern. RDKit's methods dereference it,
/// so it is not valid, its matching is `Error::Value`, and no other matcher
/// or entry takes it; `set_pattern` gives it one.
pub struct FilterHierarchyMatcher {
    matcher: FilterMatcher,
}

matcher_type!(FilterHierarchyMatcher);

/// A node without a pattern. Mirrors `FilterHierarchyMatcher()`.
impl Default for FilterHierarchyMatcher {
    fn default() -> FilterHierarchyMatcher {
        FilterHierarchyMatcher {
            matcher: FilterMatcher::from_ffi(expect_ok(ffi::filter_hierarchy_matcher_new())),
        }
    }
}

impl FilterHierarchyMatcher {
    /// A node matching a copy of `matcher`, without children. Mirrors
    /// `FilterHierarchyMatcher(matcher)`. `Error::Value` as for
    /// `ExclusionList::add_pattern`.
    pub fn new(matcher: &FilterMatcher) -> Result<FilterHierarchyMatcher> {
        Ok(FilterHierarchyMatcher {
            matcher: FilterMatcher::from_ffi(ffi::filter_hierarchy_matcher_from_matcher(
                matcher.raw(),
            )?),
        })
    }

    /// Matches a copy of `matcher`. Mirrors
    /// `FilterHierarchyMatcher.SetPattern`. `Error::Runtime` for an invalid
    /// matcher, as in Python; `Error::Value` as for
    /// `ExclusionList::add_pattern`.
    pub fn set_pattern(&mut self, matcher: &FilterMatcher) -> Result<()> {
        set_pattern(&mut self.matcher, matcher)
    }

    /// Adds a copy of `hierarchy` as a child, returning the child in this
    /// tree, which borrows the tree to add children of its own. Mirrors
    /// `FilterHierarchyMatcher.AddChild`. `Error::Runtime` for a hierarchy
    /// without a valid pattern, as in Python.
    pub fn add_child(
        &mut self,
        hierarchy: &FilterHierarchyMatcher,
    ) -> Result<FilterHierarchyMatcherMut<'_>> {
        add_child(&mut self.matcher, hierarchy)
    }
}

fn set_pattern(node: &mut FilterMatcher, matcher: &FilterMatcher) -> Result<()> {
    Ok(ffi::filter_hierarchy_matcher_set_pattern(
        node.raw_mut(),
        matcher.raw(),
    )?)
}

fn add_child<'m>(
    node: &'m mut FilterMatcher,
    hierarchy: &FilterHierarchyMatcher,
) -> Result<FilterHierarchyMatcherMut<'m>> {
    let child = ffi::filter_hierarchy_matcher_add_child(node.raw_mut(), hierarchy.raw())?;
    Ok(FilterHierarchyMatcherMut {
        matcher: FilterMatcher::from_ffi(child),
        _tree: PhantomData,
    })
}

/// A child in a `FilterHierarchyMatcher`'s tree, which borrows the tree
/// (CONVENTIONS.md W2). Mirrors the `FilterHierarchyMatcher` that
/// `FilterHierarchyMatcher.AddChild` returns.
pub struct FilterHierarchyMatcherMut<'m> {
    /// Shares the child with the tree, which only this view changes while it
    /// lives.
    matcher: FilterMatcher,
    _tree: PhantomData<&'m mut FilterHierarchyMatcher>,
}

impl Deref for FilterHierarchyMatcherMut<'_> {
    type Target = FilterMatcher;

    fn deref(&self) -> &FilterMatcher {
        &self.matcher
    }
}

impl fmt::Debug for FilterHierarchyMatcherMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("FilterHierarchyMatcherMut")
            .field(&self.get_name())
            .finish()
    }
}

impl FilterHierarchyMatcherMut<'_> {
    /// Mirrors `FilterHierarchyMatcher.SetPattern`, as
    /// `FilterHierarchyMatcher::set_pattern`.
    pub fn set_pattern(&mut self, matcher: &FilterMatcher) -> Result<()> {
        set_pattern(&mut self.matcher, matcher)
    }

    /// Mirrors `FilterHierarchyMatcher.AddChild`, as
    /// `FilterHierarchyMatcher::add_child`.
    pub fn add_child(
        &mut self,
        hierarchy: &FilterHierarchyMatcher,
    ) -> Result<FilterHierarchyMatcherMut<'_>> {
        add_child(&mut self.matcher, hierarchy)
    }
}
