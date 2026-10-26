// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The ported rdkit.Chem.Draw functions against Python's in RDKit 2026.09.1,
// whose SVGs they reproduce byte for byte (by SHA-1). PNGs come from Cairo
// and zlib, whose bytes vary by build, so only their sizes are checked.

use rdkit_cxx::chem::draw::rd_mol_draw_2d::{DrawColour, MolDrawOptions};
use rdkit_cxx::chem::draw::{
    DrawMolWithMatchesParams, MolToImageParams, MolToSVGParams, MolsToGridImageParams,
    draw_mol_with_matches_as_svg, draw_mol_with_matches_as_svg_with_params,
    mol_to_image_with_params, mol_to_svg, mol_to_svg_with_params, mols_to_grid_image,
    mols_to_grid_image_use_svg, mols_to_grid_image_use_svg_with_params,
};
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles};

fn sha1(text: &str) -> String {
    let mut h = sha1_smol::Sha1::new();
    h.update(text.as_bytes());
    h.digest().to_string()
}

/// A PNG's width and height, from its IHDR chunk.
fn png_size(png: &[u8]) -> (u32, u32) {
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    let word = |at: usize| u32::from_be_bytes(png[at..at + 4].try_into().unwrap());
    (word(16), word(20))
}

fn aspirin() -> Mol {
    mol_from_smiles("c1cc(C(=O)O)c(OC(=O)C)cc1").unwrap()
}

fn atom_indices() -> MolDrawOptions {
    MolDrawOptions {
        add_atom_indices: true,
        ..Default::default()
    }
}

fn grid_mols() -> Vec<Mol> {
    ["CCO", "c1ccccc1O", "C[C@H](N)C(=O)O", "CC(=O)Nc1ccc(O)cc1"]
        .iter()
        .map(|smiles| mol_from_smiles(smiles).unwrap())
        .collect()
}

#[test]
fn mol_to_svg_matches_python() {
    let m = aspirin();
    assert_eq!(
        sha1(&mol_to_svg(&m).unwrap()),
        "9afc250ec287966afa6f8f9dabc8d34d4d54852d"
    );
    let highlighted = MolToSVGParams {
        size: (250, 200),
        highlight_atoms: vec![3, 4, 5],
        legend: "aspirin".to_owned(),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mol_to_svg_with_params(&m, &highlighted).unwrap()),
        "4bb4f2b5f55657c1ca32363c1e233f5186108092"
    );
    let with_options = MolToSVGParams {
        draw_options: Some(atom_indices()),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mol_to_svg_with_params(&m, &with_options).unwrap()),
        "df4bd1f8d3333e23d8c37e45a7925c175f8eeb96"
    );
}

#[test]
fn mols_to_grid_image_matches_python() {
    let ms = grid_mols();
    assert_eq!(
        sha1(&mols_to_grid_image_use_svg(&ms).unwrap()),
        "bba40b87101112bc4c668bca0bad71ae7eaa4386"
    );
    // A legend too many is cut, as in Python.
    let params = MolsToGridImageParams {
        mols_per_row: 2,
        sub_img_size: (150, 150),
        legends: Some(["a", "b", "c", "d", "e"].map(str::to_owned).to_vec()),
        highlight_atom_lists: Some(vec![vec![0], vec![6], vec![1, 2], vec![]]),
        draw_options: Some(atom_indices()),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mols_to_grid_image_use_svg_with_params(&ms, &params).unwrap()),
        "0be1e27b5ed1c03659ceabdc3fb2bffbeb0dca27"
    );
    assert_eq!(png_size(&mols_to_grid_image(&ms).unwrap()), (600, 400));
    let no_rows = MolsToGridImageParams {
        mols_per_row: 0,
        ..Default::default()
    };
    assert!(mols_to_grid_image_use_svg_with_params(&ms, &no_rows).is_err());
}

#[test]
fn draw_mol_with_matches_matches_python() {
    let m = aspirin();
    let q = mol_from_smarts("C(=O)O").unwrap();
    let matches = m.get_substruct_matches(&q).unwrap();
    assert_eq!(matches, [[3, 4, 5], [8, 9, 7]]);
    assert_eq!(
        sha1(&draw_mol_with_matches_as_svg(&m, &matches).unwrap()),
        "181dbb9601892df9de9c0c6f21a2635d05445ffd"
    );
    let params = DrawMolWithMatchesParams {
        qry: Some(q),
        colors: Some(vec![
            DrawColour::new(1.0, 0.0, 0.0),
            DrawColour::new(0.0, 0.0, 1.0),
        ]),
        ..Default::default()
    };
    assert_eq!(
        sha1(&draw_mol_with_matches_as_svg_with_params(&m, &matches, &params).unwrap()),
        "31c5bfbeb48c183f055caab54701aefc318d22f0"
    );
    let no_colors = DrawMolWithMatchesParams {
        colors: Some(Vec::new()),
        ..Default::default()
    };
    assert!(draw_mol_with_matches_as_svg_with_params(&m, &matches, &no_colors).is_err());
    assert!(draw_mol_with_matches_as_svg_with_params(&m, &[], &no_colors).is_ok());
}

#[test]
fn mol_to_image_is_a_png_of_its_size() {
    let params = MolToImageParams {
        size: (250, 200),
        ..Default::default()
    };
    assert_eq!(
        png_size(&mol_to_image_with_params(&aspirin(), &params).unwrap()),
        (250, 200)
    );
}
