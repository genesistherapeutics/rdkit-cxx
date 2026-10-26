// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.MolKey`, a package of pure-Python modules for
//! molecule keys from InChIs, ported here.

pub mod inchi_info;
#[allow(clippy::module_inception)]
pub mod mol_key;
