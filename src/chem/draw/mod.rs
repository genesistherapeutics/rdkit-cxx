// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.Draw`'s C++ module, `rdMolDraw2D`. The rest of
//! `rdkit.Chem.Draw` is Python, of which the functions that draw molecules
//! to images are ported here.

mod init;
pub mod rd_mol_draw_2d;

pub use init::*;
