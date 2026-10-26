// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Query atoms and bonds.
//!
//! Python's `QueryAtom` and `QueryBond` are subclasses; here every atom and
//! bond type has the query methods, and those Python defines only on the
//! subclasses (`expand_query`, `set_query`) return `Error::Value` for a plain
//! atom or bond.

use crate::Result;
use crate::chem::rdchem::atom::expect_ok_index;
use crate::chem::rdchem::{AsAtom, AsBond, Atom, AtomMut, AtomRef, Bond, BondMut, BondRef, Mol};
use crate::error::{expect_ok, to_string};
use crate::ffi::query::ffi;

/// How `expand_query` combines two queries. Mirrors
/// `rdchem.CompositeQueryType`.
///
/// RDKit's enum has no fixed underlying type, so cxx cannot check it; it
/// crosses as `i32`, and tests check each variant's effect.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum CompositeQueryType {
    /// Mirrors `COMPOSITE_AND`.
    CompositeAnd,
    /// Mirrors `COMPOSITE_OR`.
    CompositeOr,
    /// Mirrors `COMPOSITE_XOR`.
    CompositeXor,
}

impl CompositeQueryType {
    fn to_ffi(self) -> i32 {
        match self {
            CompositeQueryType::CompositeAnd => 0,
            CompositeQueryType::CompositeOr => 1,
            CompositeQueryType::CompositeXor => 2,
        }
    }
}

/// Parameters for `expand_query_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandQueryParams {
    /// How the queries combine.
    pub how: CompositeQueryType,
    /// Keeps this query ahead of the other one.
    pub maintain_order: bool,
}

impl Default for ExpandQueryParams {
    fn default() -> ExpandQueryParams {
        ExpandQueryParams {
            how: CompositeQueryType::CompositeAnd,
            maintain_order: true,
        }
    }
}

/// Stamps the atom query getters onto a type with
/// `fn raw(&self) -> &ffi::Atom` (CONVENTIONS.md W3).
macro_rules! atom_query_getters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// A description of the query, one line per node, or empty.
            /// Mirrors `Atom.DescribeQuery`.
            pub fn describe_query(&self) -> Result<String> {
                to_string(expect_ok(ffi::atom_describe_query(self.raw())))
            }

            /// The query's type label, or empty. Mirrors `Atom.GetQueryType`.
            pub fn get_query_type(&self) -> Result<String> {
                to_string(expect_ok(ffi::atom_get_query_type(self.raw())))
            }

            /// Whether `other` matches this atom, or this atom's query. A
            /// query that needs `other`'s molecule or ring information fails
            /// without them. Mirrors `Atom.Match` (N8).
            pub fn match_(&self, other: &dyn AsAtom) -> Result<bool> {
                Ok(ffi::atom_match(self.raw(), other.atom())?)
            }
        }
    };
}

/// Stamps the atom query setters onto a type with
/// `fn raw_mut(&mut self) -> Pin<&mut ffi::Atom>` (W3).
macro_rules! atom_query_setters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Combines `other`'s query with this one, by AND. Both must be
            /// query atoms. Mirrors `QueryAtom.ExpandQuery`.
            pub fn expand_query(&mut self, other: &dyn AsAtom) -> Result<()> {
                self.expand_query_with_params(other, &ExpandQueryParams::default())
            }

            /// Mirrors `QueryAtom.ExpandQuery(other, how, maintainOrder)`.
            pub fn expand_query_with_params(
                &mut self,
                other: &dyn AsAtom,
                params: &ExpandQueryParams,
            ) -> Result<()> {
                Ok(ffi::atom_expand_query(
                    self.raw_mut(),
                    other.atom(),
                    params.how.to_ffi(),
                    params.maintain_order,
                )?)
            }

            /// Replaces this query with a copy of `other`'s. Both must be
            /// query atoms. Mirrors `QueryAtom.SetQuery`.
            pub fn set_query(&mut self, other: &dyn AsAtom) -> Result<()> {
                Ok(ffi::atom_set_query(self.raw_mut(), other.atom())?)
            }
        }
    };
}

atom_query_getters!(Atom);
atom_query_getters!(AtomRef<'m>);
atom_query_getters!(AtomMut<'m>);
atom_query_setters!(Atom);
atom_query_setters!(AtomMut<'m>);

/// Stamps the bond query getters onto a type with
/// `fn raw(&self) -> &ffi::Bond` (W3).
macro_rules! bond_query_getters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// A description of the query, one line per node, or empty.
            /// Mirrors `Bond.DescribeQuery`.
            pub fn describe_query(&self) -> Result<String> {
                to_string(expect_ok(ffi::bond_describe_query(self.raw())))
            }

            /// Whether `other` matches this bond, or this bond's query.
            /// Mirrors `Bond.Match` (N8).
            pub fn match_(&self, other: &dyn AsBond) -> Result<bool> {
                Ok(ffi::bond_match(self.raw(), other.bond())?)
            }
        }
    };
}

/// Stamps the bond query setters onto a type with
/// `fn raw_mut(&mut self) -> Pin<&mut ffi::Bond>` (W3).
macro_rules! bond_query_setters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Combines `other`'s query with this one, by AND. Both must be
            /// query bonds. Mirrors `QueryBond.ExpandQuery`.
            pub fn expand_query(&mut self, other: &dyn AsBond) -> Result<()> {
                self.expand_query_with_params(other, &ExpandQueryParams::default())
            }

            /// Mirrors `QueryBond.ExpandQuery(other, how, maintainOrder)`.
            pub fn expand_query_with_params(
                &mut self,
                other: &dyn AsBond,
                params: &ExpandQueryParams,
            ) -> Result<()> {
                Ok(ffi::bond_expand_query(
                    self.raw_mut(),
                    other.bond(),
                    params.how.to_ffi(),
                    params.maintain_order,
                )?)
            }

            /// Replaces this query with a copy of `other`'s. Both must be
            /// query bonds. Mirrors `QueryBond.SetQuery`.
            pub fn set_query(&mut self, other: &dyn AsBond) -> Result<()> {
                Ok(ffi::bond_set_query(self.raw_mut(), other.bond())?)
            }
        }
    };
}

bond_query_getters!(Bond);
bond_query_getters!(BondRef<'m>);
bond_query_getters!(BondMut<'m>);
bond_query_setters!(Bond);
bond_query_setters!(BondMut<'m>);

impl Mol {
    /// The atoms `qa`'s query matches; `qa` must be a query atom. Mirrors
    /// `Mol.GetAtomsMatchingQuery`.
    pub fn get_atoms_matching_query(&self, qa: &dyn AsAtom) -> Result<Vec<AtomRef<'_>>> {
        Ok(ffi::mol_get_atoms_matching_query(self.as_ffi(), qa.atom())?
            .into_iter()
            .map(|idx| expect_ok_index(self.get_atom_with_idx(idx)))
            .collect())
    }
}
