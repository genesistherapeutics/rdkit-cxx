// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Shared by the drawing examples: where they write their images, how
//! Jupyter shows a molecule, and what they check in the images.
//!
//! Jupyter shows a molecule as a PNG or an SVG, by a display setting
//! (`IPythonConsole.ipython_useSVG`). The examples draw SVGs, whose content
//! they can check: each atom and bond index or atom note is a
//! `class='note'` element for each character, each CIP label character a
//! `class='CIP_Code'` element, and each atom label character a
//! `class='atom-N'` element.

#![allow(dead_code)]

use std::collections::BTreeSet;
use std::path::PathBuf;

use rdkit_cxx::Result;
use rdkit_cxx::chem::Mol;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{DrawMoleculeParams, MolDraw2DSVG, MolDrawOptions};
use rdkit_cxx::chem::draw::{MolToFileParams, mol_to_file_with_params};

/// The size Jupyter draws molecules at, unless a recipe sets it. Mirrors
/// `IPythonConsole.molSize`.
pub const MOL_SIZE: (i32, i32) = (450, 150);

/// How Jupyter shows a molecule: Python's private `Draw._moltoSVG` at
/// `size`, with the console's `drawOptions`. Mirrors
/// `IPythonConsole._toSVG`. `Draw.MolToFile` writes an SVG with the same
/// `_moltoSVG`, so this draws through it.
pub fn show(mol: &Mol, size: (i32, i32), options: &MolDrawOptions) -> Result<String> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/cookbook");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("show-{}.svg", std::process::id()));
    mol_to_file_with_params(
        mol,
        &path,
        &MolToFileParams {
            size,
            options: Some(options.clone()),
            ..Default::default()
        },
    )?;
    let svg = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    Ok(svg)
}

/// As [`show`], for a molecule holding the atoms its last substructure
/// searches matched (Python's `__sssAtoms`): `_moltoSVG` highlights them,
/// and the bonds between them (it passes no bonds to highlight, which
/// `MolToFile` does), with the drawer preparing the molecule as `_moltoSVG`
/// does before drawing.
pub fn show_with_highlights(
    mol: &Mol,
    size: (i32, i32),
    options: &MolDrawOptions,
    highlight_atoms: &[u32],
) -> Result<String> {
    let mut drawer = MolDraw2DSVG::new(size.0, size.1)?;
    drawer.set_draw_options(options)?;
    drawer.draw_molecule_with_params(
        mol,
        &DrawMoleculeParams {
            highlight_atoms: highlight_atoms.to_vec(),
            highlight_bonds: None,
            ..Default::default()
        },
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Writes an image under `target/cookbook`, returning its path.
pub fn write_image(name: &str, data: impl AsRef<[u8]>) -> std::io::Result<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/cookbook");
    std::fs::create_dir_all(&dir)?;
    let path = dir.join(name);
    std::fs::write(&path, data)?;
    println!("wrote {}", path.display());
    Ok(path)
}

/// The number of the SVG's elements of a class.
pub fn count_class(svg: &str, class: &str) -> usize {
    svg.matches(&format!("class='{class}'")).count()
}

/// The atoms the SVG draws a label for.
pub fn labelled_atoms(svg: &str) -> BTreeSet<u32> {
    svg.split("class='atom-")
        .skip(1)
        .filter_map(|rest| {
            let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
            digits.parse().ok()
        })
        .collect()
}

/// The colors the SVG fills or strokes with.
pub fn colours(svg: &str) -> BTreeSet<String> {
    ["fill:#", "stroke:#"]
        .iter()
        .flat_map(|prefix| {
            svg.split(prefix)
                .skip(1)
                .map(|rest| format!("#{}", &rest[..6]))
        })
        .collect()
}

/// A PNG's width and height, from its header.
pub fn png_size(png: &[u8]) -> (u32, u32) {
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n", "not a PNG");
    let word = |at: usize| u32::from_be_bytes(png[at..at + 4].try_into().unwrap());
    (word(16), word(20))
}
