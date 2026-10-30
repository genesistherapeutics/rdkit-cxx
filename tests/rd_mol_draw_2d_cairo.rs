// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// MolDraw2DCairo. Expected values are from Python's rdMolDraw2D; the
// pixels depend on the fonts, so the tests check the PNG's size and
// metadata, not its bytes. Each refused input has a test.

mod common;

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::*;
use rdkit_cxx::chem::rd_chem_reactions::{
    ChemicalReaction, reaction_from_png_string, reaction_from_smarts, reaction_to_smarts,
};
use rdkit_cxx::chem::{
    Mol, mol_add_recursive_queries, mol_from_png_string, mol_from_smiles, mol_to_smiles,
    mols_from_png_string,
};

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn assert_value_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    let err = result.unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}

fn assert_runtime_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    let err = result.unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
}

// The PNG signature, and the image's size from its IHDR chunk.
fn png_size(png: &[u8]) -> (i32, i32) {
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    assert_eq!(&png[12..16], b"IHDR");
    let be = |b: &[u8]| i32::from_be_bytes(b.try_into().unwrap());
    (be(&png[16..20]), be(&png[20..24]))
}

fn assert_near(actual: (i32, i32), expected: (i32, i32)) {
    assert!(
        (actual.0 - expected.0).abs() <= 1 && (actual.1 - expected.1).abs() <= 1,
        "{actual:?} {expected:?}"
    );
}

#[test]
fn draws_png_with_metadata() {
    let mut d = MolDraw2DCairo::new(300, 200).unwrap();
    assert_eq!((d.width(), d.height()), (300, 200));
    d.draw_molecule(&mol("c1ccccc1O")).unwrap();
    d.finish_drawing();
    let png = d.get_drawing_text().unwrap();
    assert_eq!(png_size(&png), (300, 200));
    // Python's MolFromPNGString reads the molecule back: Oc1ccccc1.
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        mol_to_smiles(&unsafe { mol_from_png_string(&png) }.unwrap()).unwrap(),
        "Oc1ccccc1"
    );

    // Each molecule drawn adds its own (MOL1, SMILES1, rdkitPKL1 for the
    // second).
    d.draw_molecule(&mol("CCN")).unwrap();
    // SAFETY: a pickle RDKit wrote.
    let smiles: Vec<String> = unsafe { mols_from_png_string(&d.get_drawing_text().unwrap()) }
        .unwrap()
        .into_iter()
        .map(|m| mol_to_smiles(&m.unwrap()).unwrap())
        .collect();
    assert_eq!(smiles, ["Oc1ccccc1", "CCN"]);

    let mut plain = MolDraw2DCairo::new(300, 200).unwrap();
    let opts = MolDrawOptions {
        include_metadata: false,
        ..Default::default()
    };
    plain.set_draw_options(&opts).unwrap();
    plain.draw_molecule(&mol("c1ccccc1O")).unwrap();
    let png = plain.get_drawing_text().unwrap();
    assert_eq!(png_size(&png), (300, 200));
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { mol_from_png_string(&png) }.is_err());

    // Nothing drawn yet: the background alone, with no metadata.
    let blank = MolDraw2DCairo::new(20, 10).unwrap();
    assert_eq!(png_size(&blank.get_drawing_text().unwrap()), (20, 10));
}

// The metadata holds each molecule's and reaction's pickle, which RDKit
// would write with a count past 255 truncated and could not read back.
#[test]
fn metadata_needs_picklable_molecules() {
    let mut query = mol("CC");
    query
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("query", &vec!["a"; 256].join(","));
    let queries = BTreeMap::from([("a".to_owned(), mol("C"))]);
    mol_add_recursive_queries(&mut query, &queries, "query").unwrap();
    let mut d = MolDraw2DCairo::new(300, 200).unwrap();
    assert_value_error(d.draw_molecule(&query));
    let params = MolDraw2DCairoParams {
        panel_width: 150,
        panel_height: 200,
        ..Default::default()
    };
    let mut d = MolDraw2DCairo::new_with_params(300, 200, &params).unwrap();
    assert_value_error(d.draw_molecules(&[mol("CCO"), query.clone()]));
    let mut rxn = ChemicalReaction::default();
    rxn.add_reactant_template(&query);
    rxn.add_product_template(&mol("CC"));
    let mut d = MolDraw2DCairo::new(300, 200).unwrap();
    assert_value_error(d.draw_reaction(&rxn));

    let mut plain = MolDraw2DCairo::new(300, 200).unwrap();
    let opts = MolDrawOptions {
        include_metadata: false,
        ..Default::default()
    };
    plain.set_draw_options(&opts).unwrap();
    plain.draw_molecule(&query).unwrap();
    // SVG metadata holds no pickles.
    MolDraw2DSVG::new(300, 200)
        .unwrap()
        .draw_molecule(&query)
        .unwrap();
}

#[test]
fn flexible_canvases() {
    let mut d = MolDraw2DCairo::new(-1, -1).unwrap();
    // Before its first drawing a flexible canvas has no surface: Python
    // raises RuntimeError (no draw context) for these, and its ClearDrawing
    // would set the canvas up at -1 by -1 pixels, for good.
    assert_runtime_error(d.get_drawing_text());
    assert_runtime_error(d.set_colour(DrawColour::new(0.0, 0.0, 0.0)));
    assert_value_error(d.clear_drawing());
    d.draw_molecule(&mol("c1ccccc1O")).unwrap();
    // Python: 112 by 57.
    assert_near((d.width(), d.height()), (112, 57));
    assert_eq!(
        png_size(&d.get_drawing_text().unwrap()),
        (d.width(), d.height())
    );
    d.clear_drawing().unwrap();
    d.set_colour(DrawColour::new(0.0, 0.0, 0.0)).unwrap();

    let mut svg = MolDraw2DSVG::new(-1, -1).unwrap();
    assert_value_error(svg.clear_drawing());

    let no_freetype = MolDraw2DCairoParams {
        no_freetype: true,
        ..Default::default()
    };
    let mut d = MolDraw2DCairo::new_with_params(-1, -1, &no_freetype).unwrap();
    d.draw_molecule(&mol("c1ccccc1O")).unwrap();
    // Python: 118 by 57, from the system's fonts.
    assert!((d.width() - 118).abs() <= 8 && (d.height() - 57).abs() <= 8);
    assert_eq!(
        png_size(&d.get_drawing_text().unwrap()),
        (d.width(), d.height())
    );
    // Without FreeType the text must be ASCII: Cairo measures it a byte at a
    // time, and a byte of UTF-8 alone stops its drawing.
    let params = DrawMoleculeParams {
        legend: "caf\u{e9}".into(),
        ..Default::default()
    };
    let mut d = MolDraw2DCairo::new_with_params(300, 200, &no_freetype).unwrap();
    assert_value_error(d.draw_molecule_with_params(&mol("CO"), &params));
    MolDraw2DCairo::new(300, 200)
        .unwrap()
        .draw_molecule_with_params(&mol("CO"), &params)
        .unwrap();

    let mut d = MolDraw2DCairo::new(-1, -1).unwrap();
    let atoms = BTreeMap::from([(0, vec![red(), blue()])]);
    let none = BTreeMap::new();
    d.draw_molecule_with_highlights(
        &mol("c1ccccc1O"),
        "",
        &atoms,
        &none,
        &BTreeMap::new(),
        &BTreeMap::new(),
    )
    .unwrap();
    // Python: 112 by 63.
    assert_near((d.width(), d.height()), (112, 63));
}

fn red() -> DrawColour {
    DrawColour::new(1.0, 0.0, 0.0)
}

fn blue() -> DrawColour {
    DrawColour::new(0.0, 0.0, 1.0)
}

#[test]
fn canvas_sizes() {
    for (w, h) in [(0, 10), (40000, 10), (-1, 10), (-2, -2)] {
        assert_value_error(MolDraw2DCairo::new(w, h));
    }
    // Cairo's image surfaces hold at most 2^31 bytes, 4 a pixel; Python's
    // drawer makes one in an error state, from which GetDrawingText gets no
    // PNG.
    assert_value_error(MolDraw2DCairo::new(32767, 32767));
    assert_value_error(MolDraw2DCairo::new(23171, 23171));
    let d = MolDraw2DCairo::new(1, 1).unwrap();
    assert_eq!(png_size(&d.get_drawing_text().unwrap()), (1, 1));
}

#[test]
fn flexible_canvases_cairo_cannot_make() {
    // At 1000 pixels an angstrom, Python's flexible canvas for C40 is 55729
    // by 882 pixels, more than Cairo makes, and its GetDrawingText raises
    // after reading a PNG signature that is not there.
    let mut chain = mol(&"C".repeat(40));
    rdkit_cxx::chem::rd_depictor::compute_2d_coords(&mut chain).unwrap();
    let huge = MolDrawOptions {
        scaling_factor: 1000.0,
        ..Default::default()
    };
    let mut d = MolDraw2DCairo::new(-1, -1).unwrap();
    d.set_draw_options(&huge).unwrap();
    assert_value_error(d.draw_molecule(&chain));
    // A C50 ring, 26638 by 26704 pixels in Python: each side fits, but not
    // the bytes.
    let ring = mol(&format!("C1{}1", "C".repeat(49)));
    let (w, h) = d.get_mol_size(&ring).unwrap();
    assert_near((w, h), (26638, 26704));
    assert_value_error(d.draw_molecule(&ring));
    // The drawer is as it was.
    assert_eq!((d.width(), d.height()), (-1, -1));
    d.set_draw_options(&MolDrawOptions::default()).unwrap();
    d.draw_molecule(&chain).unwrap();
    png_size(&d.get_drawing_text().unwrap());

    // An SVG canvas has no such limit.
    let mut svg = MolDraw2DSVG::new(-1, -1).unwrap();
    svg.set_draw_options(&huge).unwrap();
    svg.draw_molecule(&chain).unwrap();
    assert_near((svg.width(), svg.height()), (55729, 882));

    // Python's flexible canvas for C30>>C30 at this scale is 86154 by 1322.
    let c30 = "C".repeat(30);
    let rxn = reaction_from_smarts(&format!("{c30}>>{c30}")).unwrap();
    let mut d = MolDraw2DCairo::new(-1, -1).unwrap();
    d.set_draw_options(&huge).unwrap();
    assert_value_error(d.draw_reaction(&rxn));
    let mut svg = MolDraw2DSVG::new(-1, -1).unwrap();
    svg.set_draw_options(&huge).unwrap();
    svg.draw_reaction(&rxn).unwrap();
    assert_near((svg.width(), svg.height()), (86154, 1322));
}

#[test]
fn reactions_and_grids() {
    let rxn = reaction_from_smarts("[C:1]=[O:2].[N:3]>>[C:1](-[N:3])-[O:2]").unwrap();
    let mut d = MolDraw2DCairo::new(-1, -1).unwrap();
    d.draw_reaction(&rxn).unwrap();
    // Python: 238 by 26, with the reaction in the metadata.
    assert_near((d.width(), d.height()), (238, 26));
    let png = d.get_drawing_text().unwrap();
    assert_eq!(png_size(&png), (d.width(), d.height()));
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        reaction_to_smarts(&unsafe { reaction_from_png_string(&png) }.unwrap()).unwrap(),
        "[C:1]=[O:2].[N:3]>>[C:1](-[N:3])-[O:2]"
    );

    let mut d = MolDraw2DCairo::new_with_params(
        400,
        200,
        &MolDraw2DCairoParams {
            panel_width: 200,
            panel_height: 200,
            ..Default::default()
        },
    )
    .unwrap();
    d.draw_molecules(&[mol("c1ccccc1O"), mol("CCN")]).unwrap();
    // SAFETY: a pickle RDKit wrote.
    let smiles: Vec<String> = unsafe { mols_from_png_string(&d.get_drawing_text().unwrap()) }
        .unwrap()
        .into_iter()
        .map(|m| mol_to_smiles(&m.unwrap()).unwrap())
        .collect();
    assert_eq!(smiles, ["Oc1ccccc1", "CCN"]);
}

#[test]
fn writes_png_files() {
    let mut d = MolDraw2DCairo::new(100, 100).unwrap();
    d.draw_molecule(&mol("CO")).unwrap();
    let path = common::temp_path("drawing.png");
    d.write_drawing_text(&path).unwrap();
    assert_eq!(std::fs::read(&path).unwrap(), d.get_drawing_text().unwrap());
    std::fs::remove_file(&path).unwrap();
    // RDKit logs an error and Python returns None.
    let err = d.write_drawing_text("/nonexistent/dir/x.png").unwrap_err();
    assert!(matches!(err, Error::Os { .. }), "{err:?}");
    assert_value_error(d.write_drawing_text("a\0b.png"));
}

#[test]
fn cairo_drawers_are_send() {
    fn is_send<T: Send>() {}
    is_send::<MolDraw2DCairo>();
    // A drawer drawn on one thread finishes on another.
    let mut d = MolDraw2DCairo::new(100, 100).unwrap();
    d.draw_molecule(&mol("CO")).unwrap();
    let png = std::thread::spawn(move || {
        d.draw_molecule(&mol("CCO")).unwrap();
        d.get_drawing_text().unwrap()
    })
    .join()
    .unwrap();
    assert_eq!(png_size(&png), (100, 100));
}
