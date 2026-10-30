// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.Draw`: its C++ module, `rdMolDraw2D`, and ports of
//! its Python code, the package's functions (re-exported here) and
//! `SimilarityMaps`.
//!
//! The package's functions draw a molecule, a grid of them, or a
//! reaction to PNG bytes or SVG text in one call ([`mol_to_image`],
//! [`mols_to_grid_image`]). For more control, a drawer from
//! [`rd_mol_draw_2d`] draws one or more molecules with
//! [`MolDrawOptions`](rd_mol_draw_2d::MolDrawOptions) and highlights.
//! Python's PIL images and matplotlib colormaps are not ported; images are
//! PNG bytes.
//!
//! ```no_run
//! use rdkit_cxx::chem::draw::mol_to_image;
//! use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDraw2DSVG;
//! use rdkit_cxx::chem::mol_from_smiles;
//!
//! let mol = mol_from_smiles("c1ccccc1O")?;
//! std::fs::write("phenol.png", mol_to_image(&mol)?).unwrap();
//!
//! let mut drawer = MolDraw2DSVG::new(300, 300)?;
//! let mut options = drawer.draw_options();
//! options.add_atom_indices = true;
//! drawer.set_draw_options(&options)?;
//! drawer.draw_molecule(&mol)?;
//! drawer.finish_drawing();
//! let svg = drawer.get_drawing_text()?;
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

mod init;
pub mod rd_mol_draw_2d;
pub mod similarity_maps;

pub use init::*;
