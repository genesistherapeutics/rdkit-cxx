// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Shared by the drawing examples: where they write their images, how
//! Jupyter shows a molecule, and what they check in the images.
//!
//! Jupyter shows a molecule as a PNG or an SVG, by a display setting
//! (`IPythonConsole.ipython_useSVG`). The examples draw SVGs, whose content
//! they can check: each atom and bond index, atom note, or CIP label is a
//! `class='note'` element for each character (RDKit 2024.03 draws CIP labels
//! as atom notes), and each atom label character a `class='atom-N'`
//! element.

#![allow(dead_code)]

use std::collections::BTreeSet;
use std::path::PathBuf;

use rdkit_cxx::Result;
use rdkit_cxx::chem::Mol;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::draw::{MolToSVGParams, mol_to_svg_with_params};

/// The size Jupyter draws molecules at, unless a recipe sets it. Mirrors
/// `IPythonConsole.molSize`.
pub const MOL_SIZE: (i32, i32) = (450, 150);

/// How Jupyter shows a molecule: `Draw.MolToSVG` at `size`, with the
/// console's `drawOptions`. Mirrors `IPythonConsole._toSVG`.
pub fn show(mol: &Mol, size: (i32, i32), options: &MolDrawOptions) -> Result<String> {
    mol_to_svg_with_params(
        mol,
        &MolToSVGParams {
            size,
            draw_options: Some(options.clone()),
            ..Default::default()
        },
    )
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
