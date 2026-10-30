// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.EState`, a Python package: Kier and Hall's
//! EState descriptors (L. B. Kier and L. H. Hall, _Molecular Structure
//! Description: The Electrotopological State_, Academic Press (1999)). It
//! re-exports `atom_types` and `e_state`, as Python's package does.

pub mod atom_types;
#[allow(clippy::module_inception)]
pub mod e_state;
pub mod e_state_vsa;
pub mod fingerprinter;

pub use atom_types::*;
pub use e_state::*;
