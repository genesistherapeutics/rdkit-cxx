// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FilterCatalog`, a Python module that re-exports
//! `rdfiltercatalog`, as this one does, and defines `FilterMatcher`, a base
//! class for matchers written in Python.
//!
//! Python subclasses its `FilterMatcher` and overrides `IsValid`,
//! `HasMatch`, and `GetMatches`. In Rust a matcher implements
//! [`CustomFilterMatcher`], and [`FilterMatcher::new`] wraps it, as
//! `rdfiltercatalog` binds it; Python's `FilterMatcher` is that type.

pub use crate::chem::rdfiltercatalog::*;
