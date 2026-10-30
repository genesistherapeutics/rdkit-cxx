// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdfiltercatalog.FilterMatchOps`: matchers that
//! combine others.
//!
//! Each holds copies of its arguments. `Error::Value` for an argument that is
//! a `FilterHierarchyMatcher` without a pattern, or holds one, whose `IsValid`
//! dereferences null in Python.

use std::fmt;
use std::ops::Deref;

use crate::Result;
use crate::chem::rdfiltercatalog::FilterMatcher;
use crate::chem::rdfiltercatalog::matchers::matcher_type;
use crate::ffi::filter_catalog::ffi;

/// Matches when both arguments match, with both's matches. Mirrors
/// `FilterMatchOps.And`.
pub struct And {
    matcher: FilterMatcher,
}

matcher_type!(And);

impl And {
    /// Mirrors `And(arg1, arg2)`.
    pub fn new(arg1: &FilterMatcher, arg2: &FilterMatcher) -> Result<And> {
        Ok(And {
            matcher: FilterMatcher::from_ffi(ffi::filter_match_ops_and(arg1.raw(), arg2.raw())?),
        })
    }
}

/// Matches when either argument matches, with both's matches. Mirrors
/// `FilterMatchOps.Or`.
pub struct Or {
    matcher: FilterMatcher,
}

matcher_type!(Or);

impl Or {
    /// Mirrors `Or(arg1, arg2)`.
    pub fn new(arg1: &FilterMatcher, arg2: &FilterMatcher) -> Result<Or> {
        Ok(Or {
            matcher: FilterMatcher::from_ffi(ffi::filter_match_ops_or(arg1.raw(), arg2.raw())?),
        })
    }
}

/// Matches when its argument does not, with no matches. Mirrors
/// `FilterMatchOps.Not`.
pub struct Not {
    matcher: FilterMatcher,
}

matcher_type!(Not);

impl Not {
    /// Mirrors `Not(arg1)`.
    pub fn new(arg1: &FilterMatcher) -> Result<Not> {
        Ok(Not {
            matcher: FilterMatcher::from_ffi(ffi::filter_match_ops_not(arg1.raw())?),
        })
    }
}
