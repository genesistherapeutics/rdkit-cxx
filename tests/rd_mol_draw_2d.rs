// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// MolDraw2DSVG and the functions that draw molecules. Expected values are
// from Python's rdMolDraw2D; SVG text depends on the fonts, so the tests
// check its structure, not its bytes. Each refused input has a test.

mod common;

use std::collections::BTreeMap;

use common::assert_close_to;
use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::*;
use rdkit_cxx::chem::rd_depictor::compute_2d_coords;
use rdkit_cxx::chem::{
    Mol, create_mol_data_substance_group, create_mol_substance_group,
    get_mol_substance_group_with_idx_mut, mol_from_smiles,
};
use rdkit_cxx::geometry::Point3D;

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn depicted(smiles: &str) -> Mol {
    let mut m = mol(smiles);
    compute_2d_coords(&mut m).unwrap();
    m
}

fn svg() -> MolDraw2DSVG {
    MolDraw2DSVG::new(300, 300).unwrap()
}

fn assert_value_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    let err = result.unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}

fn assert_index_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    let err = result.unwrap_err();
    assert!(matches!(err, Error::Index { .. }), "{err:?}");
}

fn with_legend(legend: &str) -> DrawMoleculeParams {
    DrawMoleculeParams {
        legend: legend.to_string(),
        ..Default::default()
    }
}

#[test]
fn drawer_defaults_match_python() {
    let d = svg();
    assert_eq!((d.width(), d.height()), (300, 300));
    assert_close_to(d.font_size(), 0.6, 1e-12);
    assert_eq!(d.line_width(), 2.0);
    assert!(d.fill_polys());
    assert!(!d.flexi_mode());
    assert_eq!((d.offset().x, d.offset().y), (0.0, 0.0));
    assert_eq!(d.draw_options(), MolDrawOptions::default());
    assert_eq!(
        DrawMoleculeParams::default(),
        DrawMoleculeParams {
            highlight_atoms: vec![],
            highlight_bonds: None,
            highlight_atom_colors: BTreeMap::new(),
            highlight_bond_colors: BTreeMap::new(),
            highlight_atom_radii: BTreeMap::new(),
            conf_id: -1,
            legend: String::new(),
        }
    );
    let p = MolDraw2DSVGParams::default();
    assert_eq!(
        (p.panel_width, p.panel_height, p.no_freetype),
        (-1, -1, false)
    );
    let p = TagAtomsParams::default();
    assert_eq!((p.radius, p.events.len()), (0.2, 0));
    let p = PrepareMolForDrawingParams::default();
    assert_eq!(
        (
            p.kekulize,
            p.add_chiral_hs,
            p.wedge_bonds,
            p.force_coords,
            p.wavy_bonds
        ),
        (true, true, true, false, false)
    );
    let p = MolToSVGParams::default();
    assert_eq!(
        (p.width, p.height, p.line_width_mult, p.include_atom_circles),
        (300, 300, 1, true)
    );
    assert_eq!((p.kekulize, p.conf_id), (true, -1));
    let p = PrepareAndDrawMoleculeParams::default();
    assert_eq!((p.conf_id, p.kekulize, p.legend.as_str()), (-1, true, ""));
    assert_eq!(GetMolSizeParams::default().conf_id, -1);
}

#[test]
fn draws_svg() {
    let mut d = svg();
    let phenol = mol("c1ccccc1O");
    d.draw_molecule_with_params(&phenol, &with_legend("phenol"))
        .unwrap();
    d.finish_drawing();
    let text = d.get_drawing_text().unwrap();
    assert!(text.starts_with("<?xml version='1.0' encoding='iso-8859-1'?>"));
    assert!(text.ends_with("</svg>\n"));
    assert!(text.contains("width='300px'"));
    assert!(text.contains("atom-6") && text.contains("bond-6"));
    // RDKit's depiction drew a copy: the molecule keeps no conformer.
    assert_eq!(phenol.get_num_conformers(), 0);

    // Without FreeType, text is SVG text: the O, the H, and the legend.
    let params = MolDraw2DSVGParams {
        no_freetype: true,
        ..Default::default()
    };
    let mut d = MolDraw2DSVG::new_with_params(300, 300, &params).unwrap();
    d.draw_molecule_with_params(&phenol, &with_legend("x"))
        .unwrap();
    d.finish_drawing();
    assert_eq!(d.get_drawing_text().unwrap().matches("<text").count(), 3);
}

#[test]
fn mol_size_and_flexible_canvas() {
    let ethanol = mol("CCO");
    let mut flexi = MolDraw2DSVG::new(-1, -1).unwrap();
    let (w, h) = flexi.get_mol_size(&ethanol).unwrap();
    assert!((w - 70).abs() <= 1 && (h - 21).abs() <= 1, "{w} {h}");
    let mut fixed = MolDraw2DSVG::new(250, 200).unwrap();
    assert_eq!(fixed.get_mol_size(&ethanol).unwrap(), (250, 200));
    flexi.draw_molecule(&ethanol).unwrap();
    assert_eq!((flexi.width(), flexi.height()), (w, h));
}

#[test]
fn setters() {
    let mut d = svg();
    d.set_font_size(20.0).unwrap();
    assert_eq!(d.font_size(), 20.0);
    d.set_line_width(3.0).unwrap();
    assert_eq!(d.line_width(), 3.0);
    assert_eq!(d.draw_options().bond_line_width, 3.0);
    d.set_fill_polys(false);
    assert!(!d.fill_polys());
    d.set_flexi_mode(true);
    assert!(d.flexi_mode());
    d.set_offset(5, 7).unwrap();
    assert_eq!((d.offset().x, d.offset().y), (5.0, 7.0));
    d.set_colour(DrawColour::new(0.1, 0.2, 0.3)).unwrap();
    let mut opts = d.draw_options();
    opts.add_atom_indices = true;
    d.set_draw_options(&opts).unwrap();
    assert!(d.draw_options().add_atom_indices);
    d.clear_drawing().unwrap();

    update_drawer_params_from_json(&mut d, r#"{"addBondIndices": true, "bondLineWidth": 4}"#)
        .unwrap();
    assert!(d.draw_options().add_bond_indices);
    assert_eq!(d.line_width(), 4.0);
    set_dark_mode_mol_draw_2d(&mut d);
    assert_eq!(
        d.draw_options().get_background_colour(),
        DrawColour::new(0.0, 0.0, 0.0)
    );
    set_monochrome_mode_mol_draw_2d(
        &mut d,
        DrawColour::new(0.1, 0.2, 0.3),
        DrawColour::new(0.9, 0.9, 0.9),
    )
    .unwrap();
    let opts = d.draw_options();
    assert_eq!(opts.get_background_colour(), DrawColour::new(0.9, 0.9, 0.9));
    assert_eq!(
        opts.get_atom_palette(),
        BTreeMap::from([(-1, DrawColour::new(0.1, 0.2, 0.3))])
    );
}

#[test]
fn setters_check_their_range() {
    let mut d = svg();
    assert_value_error(d.set_font_size(0.0));
    assert_value_error(d.set_font_size(f64::NAN));
    assert_value_error(d.set_line_width(100.0));
    assert_value_error(d.set_offset(100_000, 0));
    assert_value_error(d.set_colour(DrawColour::new(1.5, 0.0, 0.0)));
    let mut opts = MolDrawOptions {
        padding: 0.6,
        ..Default::default()
    };
    assert_value_error(d.set_draw_options(&opts));
    opts.padding = 0.05;
    opts.font_file = "a\0b".into();
    assert_value_error(d.set_draw_options(&opts));
    assert_value_error(update_drawer_params_from_json(
        &mut d,
        r#"{"scalingFactor": 0}"#,
    ));
    assert_value_error(set_monochrome_mode_mol_draw_2d(
        &mut d,
        DrawColour::new(2.0, 0.0, 0.0),
        DrawColour::new(1.0, 1.0, 1.0),
    ));
    // Nothing changed.
    assert_eq!(d.draw_options(), MolDrawOptions::default());
    assert_eq!(d.line_width(), 2.0);
}

#[test]
fn canvas_sizes() {
    for (w, h) in [(0, 0), (40_000, 10), (-1, 300), (300, -1), (-5, -5)] {
        assert_value_error(MolDraw2DSVG::new(w, h));
    }
    let panels = MolDraw2DSVGParams {
        panel_width: 100,
        ..Default::default()
    };
    assert_value_error(MolDraw2DSVG::new_with_params(-1, -1, &panels));
    assert_value_error(MolDraw2DSVG::new_with_params(
        300,
        300,
        &MolDraw2DSVGParams {
            panel_height: 0,
            ..Default::default()
        },
    ));
    // Too small to draw in: DrawMol's int sizes come out zero.
    let mut d = MolDraw2DSVG::new(1, 1).unwrap();
    assert_value_error(d.draw_molecule(&mol("CCO")));
    let mut d = MolDraw2DSVG::new(10, 10).unwrap();
    d.draw_molecule(&mol("CCO")).unwrap();
    // A drawer that refused an input can still draw.
    let mut d = svg();
    assert_value_error(d.draw_molecule_with_params(&mol("CCO"), &with_legend("<sub>2</sub>")));
    d.draw_molecule(&mol("CCO")).unwrap();
}

#[test]
fn highlights() {
    let mut d = svg();
    let m = mol("c1ccccc1O");
    let params = DrawMoleculeParams {
        highlight_atoms: vec![0, 1, 6],
        highlight_atom_colors: BTreeMap::from([(0, DrawColour::new(1.0, 0.0, 0.0))]),
        highlight_atom_radii: BTreeMap::from([(1, 0.5)]),
        ..Default::default()
    };
    d.draw_molecule_with_params(&m, &params).unwrap();
    let params = DrawMoleculeParams {
        highlight_atoms: vec![0, 1],
        highlight_bonds: Some(vec![3]),
        highlight_bond_colors: BTreeMap::from([(3, DrawColour::new(0.0, 0.0, 1.0))]),
        ..Default::default()
    };
    d.draw_molecule_with_params(&m, &params).unwrap();

    let bad_atom = DrawMoleculeParams {
        highlight_atoms: vec![7],
        ..Default::default()
    };
    assert_index_error(d.draw_molecule_with_params(&m, &bad_atom));
    let bad_bond = DrawMoleculeParams {
        highlight_bonds: Some(vec![7]),
        ..Default::default()
    };
    assert_index_error(d.draw_molecule_with_params(&m, &bad_bond));
    let bad_colour = DrawMoleculeParams {
        highlight_atom_colors: BTreeMap::from([(0, DrawColour::new(f64::NAN, 0.0, 0.0))]),
        ..Default::default()
    };
    assert_value_error(d.draw_molecule_with_params(&m, &bad_colour));
    let bad_radius = DrawMoleculeParams {
        highlight_atom_radii: BTreeMap::from([(0, 100.0)]),
        ..Default::default()
    };
    assert_value_error(d.draw_molecule_with_params(&m, &bad_radius));
    let bad_size = GetMolSizeParams {
        highlight_bonds: vec![99],
        ..Default::default()
    };
    assert_index_error(d.get_mol_size_with_params(&m, &bad_size));
}

#[test]
fn conformers() {
    let mut d = svg();
    let m = depicted("CCO");
    let params = DrawMoleculeParams {
        conf_id: 5,
        ..Default::default()
    };
    // Checked before RDKit draws, so the drawer stays usable.
    assert_value_error(d.draw_molecule_with_params(&m, &params));
    d.draw_molecule(&m).unwrap();
    for v in [1e6, f64::NAN, f64::INFINITY] {
        let mut far = m.clone();
        far.get_conformer_mut()
            .unwrap()
            .set_atom_position(0, Point3D::new(v, 0.0, 0.0))
            .unwrap();
        assert_value_error(d.draw_molecule(&far));
    }
}

#[test]
fn degenerate_layouts_keep_the_drawer_usable() {
    // Two coincident ring-fusion atoms: RDKit throws while laying out the
    // double bonds, which the wrappers' layout of it finds before RDKit
    // draws. (Coincident atoms of two bonds are refused before that.)
    let mut m = depicted("c1ccc2ccccc2c1");
    let p = m.get_conformer().unwrap().get_atom_position(3).unwrap();
    m.get_conformer_mut()
        .unwrap()
        .set_atom_position(8, p)
        .unwrap();
    let mut d = svg();
    assert!(matches!(
        d.draw_molecule(&m).unwrap_err(),
        Error::Runtime { .. }
    ));
    d.draw_molecule(&mol("CCO")).unwrap();
    d.finish_drawing();
    assert!(d.get_drawing_text().unwrap().ends_with("</svg>\n"));
}

// Text with markup RDKit reads at index -1 (adjustStringRectsForSuperSubScript)
// or 0 of nothing (alignString).
const BAD_MARKUP: [&str; 4] = [
    "<sub>2</sub>",
    "<sup>2</sup>",
    "<sub></sub>",
    "<sub>x</sub>y",
];

#[test]
fn legends() {
    let m = mol("CCO");
    let mut d = svg();
    for legend in [
        "",
        "ethanol",
        "a\nb",
        "x<sub>2</sub>",
        "<sup>1</sup>x",
        "\n\n",
        "é",
        "x<sub>",
    ] {
        d.draw_molecule_with_params(&m, &with_legend(legend))
            .unwrap();
    }
    for legend in BAD_MARKUP {
        assert_value_error(d.draw_molecule_with_params(&m, &with_legend(legend)));
    }
    assert_value_error(d.draw_molecule_with_params(&m, &with_legend("a\n<sub>b")));
    assert_value_error(d.draw_molecule_with_params(&m, &with_legend(&"x".repeat(10_001))));

    // Without FreeType, RDKit indexes its width table by char.
    let params = MolDraw2DSVGParams {
        no_freetype: true,
        ..Default::default()
    };
    let mut d = MolDraw2DSVG::new_with_params(300, 300, &params).unwrap();
    d.draw_molecule_with_params(&m, &with_legend("plain"))
        .unwrap();
    assert_value_error(d.draw_molecule_with_params(&m, &with_legend("é")));
    assert_value_error(d.draw_molecule_with_params(&m, &with_legend("\t")));
}

#[test]
fn atom_labels_and_notes() {
    let mut d = svg();
    let m = mol("CCO");
    let mut opts = MolDrawOptions {
        atom_labels: BTreeMap::from([(0, "R<sub>1</sub>".to_string())]),
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    d.draw_molecule(&m).unwrap();
    for label in BAD_MARKUP.into_iter().chain(["a:b>c<d"]) {
        opts.atom_labels = BTreeMap::from([(0, label.to_string())]);
        d.set_draw_options(&opts).unwrap();
        assert_value_error(d.draw_molecule(&m));
    }
    // Labels RDKit does not draw are not checked.
    opts.no_atom_labels = true;
    d.set_draw_options(&opts).unwrap();
    d.draw_molecule(&m).unwrap();
    d.set_draw_options(&MolDrawOptions::default()).unwrap();

    for key in ["atomLabel", "_displayLabel", "_displayLabelW", "atomNote"] {
        let mut labelled = m.clone();
        labelled
            .get_atom_with_idx_mut(1)
            .unwrap()
            .set_prop(key, "<sub>x</sub>");
        assert_value_error(d.draw_molecule(&labelled));
        labelled
            .get_atom_with_idx_mut(1)
            .unwrap()
            .set_prop(key, "ok");
        d.draw_molecule(&labelled).unwrap();
    }
    let mut noted = m.clone();
    noted
        .get_bond_with_idx_mut(0)
        .unwrap()
        .set_prop("bondNote", "<sup>x</sup>");
    assert_value_error(d.draw_molecule(&noted));
    let mut noted = m.clone();
    noted.set_prop("molNote", "<sub>x");
    assert_value_error(d.draw_molecule(&noted));

    // A dummy's label goes through isupper, so it must be ASCII.
    let mut dummy = mol("*C");
    dummy
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("dummyLabel", "Ω");
    assert_value_error(d.draw_molecule(&dummy));
    dummy
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("dummyLabel", "R");
    d.draw_molecule(&dummy).unwrap();
    // A charged dummy with an empty label is a superscript alone.
    let mut charged = mol("[*+]");
    charged
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("dummyLabel", "");
    assert_value_error(d.draw_molecule(&charged));
}

#[test]
fn stereo_annotations() {
    let mut d = svg();
    let opts = MolDrawOptions {
        add_stereo_annotation: true,
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    let mut m = mol("C[C@H](F)Cl");
    m.get_atom_with_idx_mut(1)
        .unwrap()
        .set_prop("_CIPCode", "S");
    d.draw_molecule(&m).unwrap();

    // RDKit 2025.09 labels stereo groups "abs", "or1", "and1", always.
    let grouped = mol("C[C@H](F)Cl |o1:1|");
    d.draw_molecule(&grouped).unwrap();
    let grouped = mol("C[C@H](F)Cl |&1:1|");
    d.draw_molecule(&grouped).unwrap();
}

#[test]
fn atom_regions() {
    let mut d = svg();
    // A region of one atom is not drawn.
    let mut opts = MolDrawOptions {
        atom_regions: vec![vec![0, 1], vec![99]],
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    d.draw_molecule(&mol("CCO")).unwrap();
    opts.atom_regions = vec![vec![0, 99]];
    d.set_draw_options(&opts).unwrap();
    assert_index_error(d.draw_molecule(&mol("CCO")));
}

#[test]
fn ring_information_must_be_cycles() {
    // add_ring refuses a ring that is not a cycle, so no drawing meets one.
    let mut m = depicted("CC=CC");
    assert_value_error(m.get_ring_info_mut().add_ring(&[1, 2], &[0, 1]));
    svg().draw_molecule(&m).unwrap();
}

#[test]
fn substance_groups() {
    let mut d = svg();
    let mut m = depicted("CCO");
    {
        let mut sg = create_mol_data_substance_group(&mut m, "pH", "7");
        sg.add_atom_with_idx(0).unwrap();
        sg.set_prop("FIELDDISP", "    1.0000    2.0000    DA    ALL  0       0");
    }
    d.draw_molecule(&m).unwrap();
    // RDKit reads FIELDDISP[25].
    let mut short = depicted("CCO");
    {
        let mut sg = create_mol_data_substance_group(&mut short, "pH", "7");
        sg.add_atom_with_idx(0).unwrap();
        sg.set_prop("FIELDDISP", "    1.0000    2.0000");
    }
    assert_value_error(d.draw_molecule(&short));
    let mut far = depicted("CCO");
    {
        let mut sg = create_mol_data_substance_group(&mut far, "pH", "7");
        sg.add_atom_with_idx(0).unwrap();
        sg.set_prop("FIELDDISP", "     1e300    2.0000    DA    ALL  0       0");
    }
    assert_value_error(d.draw_molecule(&far));
    let mut text = depicted("CCO");
    {
        let mut sg = create_mol_data_substance_group(&mut text, "pH", "<sub>7");
        sg.add_atom_with_idx(0).unwrap();
    }
    assert_value_error(d.draw_molecule(&text));

    // Brackets: coordinates and labels.
    let bracket = [
        Point3D::new(-1.0, 1.0, 0.0),
        Point3D::new(-1.0, -1.0, 0.0),
        Point3D::new(0.0, 0.0, 0.0),
    ];
    let mut sru = depicted("CCO");
    {
        let mut sg = create_mol_substance_group(&mut sru, "SRU");
        sg.add_atom_with_idx(1).unwrap();
        sg.add_bracket(&bracket).unwrap();
        sg.set_prop("CONNECT", "ht");
        sg.set_prop("LABEL", "n");
    }
    d.draw_molecule(&sru).unwrap();
    let mut empty = sru.clone();
    get_mol_substance_group_with_idx_mut(&mut empty, 0)
        .unwrap()
        .set_prop("CONNECT", "");
    assert_value_error(d.draw_molecule(&empty));
    let mut typed = depicted("CCO");
    {
        let mut sg = create_mol_substance_group(&mut typed, "SRU");
        sg.add_bracket(&bracket).unwrap();
        sg.set_prop("TYPE", "SRÜ");
    }
    assert_value_error(d.draw_molecule(&typed));
    let mut far = depicted("CCO");
    {
        let mut sg = create_mol_substance_group(&mut far, "SRU");
        let mut b = bracket;
        b[0].x = f64::NAN;
        sg.add_bracket(&b).unwrap();
    }
    assert_value_error(d.draw_molecule(&far));
}

#[test]
fn link_nodes() {
    let mut d = svg();
    let mut m = depicted("CCC");
    // Atom 2 links to atoms 1 and 3 (one-based); atoms 1 and 2 coincide.
    m.set_prop("_molLinkNodes", "1 3 2 2 1 2 3");
    d.draw_molecule(&m).unwrap();
    let p = m.get_conformer().unwrap().get_atom_position(1).unwrap();
    m.get_conformer_mut()
        .unwrap()
        .set_atom_position(0, p)
        .unwrap();
    m.get_conformer_mut()
        .unwrap()
        .set_atom_position(2, Point3D::new(p.x - 1.5, p.y, 0.0))
        .unwrap();
    assert_value_error(d.draw_molecule(&m));
}

#[test]
fn bitmap_fonts_are_refused() {
    // FreeType gives a bitmap font no units per em, which RDKit divides by.
    let path = common::temp_path("tiny.bdf");
    std::fs::write(
        &path,
        "STARTFONT 2.1\nFONT tiny\nSIZE 8 75 75\nFONTBOUNDINGBOX 8 8 0 0\n\
         STARTPROPERTIES 2\nFONT_ASCENT 8\nFONT_DESCENT 0\nENDPROPERTIES\nCHARS 1\n\
         STARTCHAR A\nENCODING 65\nSWIDTH 1000 0\nDWIDTH 8 0\nBBX 8 8 0 0\nBITMAP\n\
         FF\n81\n81\nFF\n81\n81\n81\n81\nENDCHAR\nENDFONT\n",
    )
    .unwrap();
    let mut d = svg();
    let mut opts = MolDrawOptions {
        font_file: path.to_str().unwrap().into(),
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    assert_value_error(d.draw_molecule(&mol("CO")));
    // A font that does not load falls back to the built-in one, as in RDKit,
    // whose warning names the default font from $RDBASE, and fails without
    // it; the drawer can still draw.
    opts.font_file = "/nonexistent.ttf".into();
    d.set_draw_options(&opts).unwrap();
    match std::env::var_os("RDBASE") {
        Some(_) => d.draw_molecule(&mol("CO")).unwrap(),
        None => assert!(matches!(
            d.draw_molecule(&mol("CO")).unwrap_err(),
            Error::Runtime { .. }
        )),
    }
    opts.font_file = "BuiltinRobotoRegular".into();
    d.set_draw_options(&opts).unwrap();
    d.draw_molecule(&mol("CO")).unwrap();
}

#[test]
fn comic_mode() {
    let mut d = svg();
    let opts = MolDrawOptions {
        comic_mode: true,
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    d.draw_molecule(&mol("c1ccccc1O")).unwrap();
}

#[test]
fn tags_and_metadata() {
    let m = depicted("c1ccccc1O");
    let mut d = svg();
    // Nothing drawn yet: RDKit's own check, and the wrappers'.
    assert!(matches!(
        d.tag_atoms(&m).unwrap_err(),
        Error::Runtime { .. }
    ));
    assert_value_error(d.add_molecule_metadata(&m));
    d.draw_molecule(&m).unwrap();
    let params = TagAtomsParams {
        events: BTreeMap::from([("onclick".to_string(), "select".to_string())]),
        ..Default::default()
    };
    d.tag_atoms_with_params(&m, &params).unwrap();
    d.add_molecule_metadata(&m).unwrap();
    assert_value_error(d.add_molecule_metadata_with_conf_id(&m, 3));
    d.finish_drawing();
    let text = d.get_drawing_text().unwrap();
    assert_eq!(text.matches("atom-selector").count(), 7);
    assert_eq!(text.matches("bond-selector").count(), 7);
    assert_eq!(text.matches("<rdkit:atom").count(), 7);
    assert_eq!(text.matches("<rdkit:bond").count(), 7);
    assert!(text.contains("onclick='select(0);'"));
    // A molecule with more atoms than the drawing: RDKit's own check.
    let mut d = svg();
    d.draw_molecule(&m).unwrap();
    assert!(matches!(
        d.tag_atoms(&depicted("c1ccccc1OCCC")).unwrap_err(),
        Error::Runtime { .. }
    ));
}

#[test]
fn module_functions() {
    let ethanol = depicted("CCO");
    assert_close_to(mean_bond_length(&ethanol).unwrap(), 1.5, 1e-9);
    assert_eq!(mean_bond_length(&mol("C")).unwrap(), 0.0);
    assert_value_error(mean_bond_length(&mol("CC")));
    assert_value_error(mean_bond_length_with_conf_id(&ethanol, 4));

    let bicycle = mol("C1C[C@H]2CC[C@@H]1C2");
    let prepared = prepare_mol_for_drawing(&bicycle).unwrap();
    assert_eq!(
        (prepared.get_num_atoms(), prepared.get_num_conformers()),
        (9, 1)
    );
    assert_eq!(bicycle.get_num_atoms(), 7);
    let benzene = prepare_mol_for_drawing(&mol("c1ccccc1")).unwrap();
    let orders: Vec<f64> = benzene
        .get_bonds()
        .map(|b| b.get_bond_type_as_double().unwrap())
        .collect();
    assert_eq!(orders, [1.0, 2.0, 1.0, 2.0, 1.0, 2.0]);
    let params = PrepareMolForDrawingParams {
        kekulize: false,
        add_chiral_hs: false,
        ..Default::default()
    };
    let plain = prepare_mol_for_drawing_with_params(&bicycle, &params).unwrap();
    assert_eq!(plain.get_num_atoms(), 7);

    let text = mol_to_svg(&mol("CCO")).unwrap();
    assert!(text.contains("width='300px'") && text.ends_with("</svg>\n"));
    let params = MolToSVGParams {
        width: 200,
        height: 100,
        highlight_atoms: vec![0, 1],
        ..Default::default()
    };
    let text = mol_to_svg_with_params(&mol("CCO"), &params).unwrap();
    assert!(text.contains("width='200px'") && text.contains("height='100px'"));
    assert_eq!(text.matches("ellipse").count(), 2);
    for params in [
        MolToSVGParams {
            width: 0,
            ..Default::default()
        },
        MolToSVGParams {
            width: u32::MAX,
            ..Default::default()
        },
        MolToSVGParams {
            line_width_mult: 100,
            ..Default::default()
        },
    ] {
        assert_value_error(mol_to_svg_with_params(&mol("CCO"), &params));
    }

    let mut d = svg();
    prepare_and_draw_molecule(&mut d, &bicycle).unwrap();
    let params = PrepareAndDrawMoleculeParams {
        legend: "bicycle".into(),
        highlight_atoms: vec![0],
        ..Default::default()
    };
    prepare_and_draw_molecule_with_params(&mut d, &bicycle, &params).unwrap();
    // The option is restored.
    assert!(d.draw_options().prepare_mols_before_drawing);
}

#[test]
fn drawers_are_send() {
    fn is_send<T: Send>() {}
    is_send::<MolDraw2DSVG>();
    is_send::<MolDraw2D>();
}

#[test]
fn flexible_canvases_too_small_to_scale() {
    // RDKit's scale for these comes out zero, and with a fixed font size it
    // then divides by it and converts the NaN size to int.
    let mut d = MolDraw2DSVG::new(-1, -1).unwrap();
    let mut opts = MolDrawOptions {
        fixed_font_size: 20,
        scaling_factor: 1.0,
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    assert_value_error(d.draw_molecule(&mol("[Fe]")));
    assert_value_error(d.get_mol_size(&mol("[Fe]")));
    opts.scaling_factor = 20.0;
    d.set_draw_options(&opts).unwrap();
    d.draw_molecule(&mol("[Fe]")).unwrap();
    assert!(d.width() > 0 && d.height() > 0);

    // A molecule a thousandth of the usual size.
    let mut tiny = depicted("CC");
    let n = tiny.get_num_atoms();
    let mut conf = tiny.get_conformer_mut().unwrap();
    for i in 0..n {
        let p = conf.get_atom_position(i).unwrap();
        conf.set_atom_position(i, Point3D::new(p.x * 1e-3, p.y * 1e-3, 0.0))
            .unwrap();
    }
    // RDKit leaves out bonds shorter than 0.01, so nothing is left to size
    // the drawing by, on any canvas.
    let mut d = MolDraw2DSVG::new(-1, -1).unwrap();
    assert_value_error(d.draw_molecule(&tiny));
    assert_value_error(svg().draw_molecule(&tiny));
}

#[test]
fn extents_must_include_something_drawn() {
    let mut d = svg();
    let opts = MolDrawOptions {
        drawing_extents_include: DrawElement::HIGHLIGHTS,
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    assert_value_error(d.draw_molecule(&mol("CCO")));
    let params = DrawMoleculeParams {
        highlight_atoms: vec![0],
        ..Default::default()
    };
    d.draw_molecule_with_params(&mol("CCO"), &params).unwrap();
}

#[test]
fn attachment_points_need_bond_length() {
    let mut d = svg();
    let mut opts = MolDrawOptions {
        dummies_are_attachments: true,
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    let mut m = depicted("*CC");
    d.draw_molecule(&m).unwrap();
    let n = m.get_num_atoms();
    let mut conf = m.get_conformer_mut().unwrap();
    for i in 0..n {
        let p = conf.get_atom_position(i).unwrap();
        conf.set_atom_position(i, Point3D::new(p.x * 1e-5, p.y * 1e-5, 0.0))
            .unwrap();
    }
    assert_value_error(d.draw_molecule(&m));
    // A wavy line one unit long, with waves too wide to fit two: RDKit's
    // drawWavyLine refuses fewer than two segments mid-drawing.
    let mut big = depicted("*CC");
    let n = big.get_num_atoms();
    let mut conf = big.get_conformer_mut().unwrap();
    for i in 0..n {
        let p = conf.get_atom_position(i).unwrap();
        conf.set_atom_position(i, Point3D::new(p.x * 4.0, p.y * 4.0, 0.0))
            .unwrap();
    }
    opts.multiple_bond_offset = 1.0;
    d.set_draw_options(&opts).unwrap();
    assert_value_error(d.draw_molecule(&big));
    d.draw_molecule(&mol("CCO")).unwrap();
}

#[test]
fn comic_mode_needs_bonds_with_length() {
    let mut d = svg();
    let opts = MolDrawOptions {
        comic_mode: true,
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    let mut m = depicted("CCO");
    let p = m.get_conformer().unwrap().get_atom_position(0).unwrap();
    m.get_conformer_mut()
        .unwrap()
        .set_atom_position(1, p)
        .unwrap();
    assert_value_error(d.draw_molecule(&m));
}

/// RDKit 2026.03 reads an uninitialized bond type for an atom of two bonds
/// with a neighbor at its position (isLinearAtom); an atom of three draws.
#[test]
fn atoms_of_two_bonds_need_their_neighbors_apart() {
    let at_atom_0 = |smiles: &str| {
        let mut m = depicted(smiles);
        let p = m.get_conformer().unwrap().get_atom_position(0).unwrap();
        m.get_conformer_mut()
            .unwrap()
            .set_atom_position(1, p)
            .unwrap();
        m
    };
    assert_value_error(svg().draw_molecule(&at_atom_0("CCO")));
    svg().draw_molecule(&at_atom_0("CC(C)C")).unwrap();
}

/// RDKit's hand-drawn lines never end for a line of no length, and a small
/// canvas gives a reaction's plus signs and arrow no length.
#[test]
fn comic_mode_needs_reaction_lines_with_length() {
    let draw = |smarts: &'static str, w: i32, h: i32| {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let mut d = MolDraw2DSVG::new(w, h).unwrap();
            let opts = MolDrawOptions {
                comic_mode: true,
                ..Default::default()
            };
            d.set_draw_options(&opts).unwrap();
            let rxn = rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts(smarts).unwrap();
            tx.send(d.draw_reaction(&rxn).is_ok()).unwrap();
        });
        rx.recv_timeout(std::time::Duration::from_secs(60))
            .expect("the drawing hung")
    };
    // A plus of no length hangs only when its ends' random shifts match,
    // one time in four, so 40x40 can finish in RDKit too.
    let plus = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]";
    for (w, h) in [(40, 10), (40, 8), (10, 10), (3, 3), (40, 40)] {
        assert!(!draw(plus, w, h), "{w}x{h}");
    }
    for (w, h) in [(10, 10), (3, 3)] {
        assert!(!draw("[C:1]=O>>[C:1]O", w, h), "{w}x{h}");
    }
    for (w, h) in [(100, 100), (300, 300)] {
        assert!(draw(plus, w, h), "{w}x{h}");
    }
}

fn red() -> DrawColour {
    DrawColour::new(1.0, 0.0, 0.0)
}

fn blue() -> DrawColour {
    DrawColour::new(0.0, 0.0, 1.0)
}

fn green() -> DrawColour {
    DrawColour::new(0.0, 1.0, 0.0)
}

#[test]
fn multi_colour_highlights() {
    let m = mol("c1ccccc1O");
    let atoms = BTreeMap::from([
        (0, vec![red(), blue()]),
        (1, vec![green()]),
        (6, vec![red()]),
    ]);
    let bonds = BTreeMap::from([(0, vec![red()]), (5, vec![blue(), green()])]);
    let radii = BTreeMap::from([(0, 0.4)]);
    let multipliers = BTreeMap::from([(0, 2)]);
    // Python's counts of each color in the SVG.
    for (style, counts) in [
        (MultiColourHighlightStyle::CircleAndLine, (9, 4, 4)),
        (MultiColourHighlightStyle::Lasso, (5, 1, 1)),
    ] {
        let mut d = svg();
        let opts = MolDrawOptions {
            multi_colour_highlight_style: style,
            ..Default::default()
        };
        d.set_draw_options(&opts).unwrap();
        d.draw_molecule_with_highlights(&m, "mch", &atoms, &bonds, &radii, &multipliers)
            .unwrap();
        d.finish_drawing();
        let text = d.get_drawing_text().unwrap();
        let count = |c: &str| text.matches(c).count();
        assert_eq!(
            (count("#FF0000"), count("#0000FF"), count("#00FF00")),
            counts,
            "{style:?}"
        );
    }
    // RDKit draws the first conformer whatever confId is.
    let mut d = svg();
    d.draw_molecule_with_highlights_with_conf_id(&m, "", &atoms, &bonds, &radii, &multipliers, 7)
        .unwrap();
}

#[test]
fn multi_colour_highlights_are_checked() {
    let m = mol("c1ccccc1O");
    let none = BTreeMap::new();
    let no_radii = BTreeMap::new();
    let no_multipliers = BTreeMap::new();
    for style in [
        MultiColourHighlightStyle::CircleAndLine,
        MultiColourHighlightStyle::Lasso,
    ] {
        let mut d = svg();
        let opts = MolDrawOptions {
            multi_colour_highlight_style: style,
            ..Default::default()
        };
        d.set_draw_options(&opts).unwrap();
        // DrawMolMCHLasso sets bits of an unchecked atom index.
        let far = BTreeMap::from([(99, vec![red()])]);
        assert_index_error(d.draw_molecule_with_highlights(
            &m,
            "",
            &far,
            &none,
            &no_radii,
            &no_multipliers,
        ));
        assert_index_error(d.draw_molecule_with_highlights(
            &m,
            "",
            &none,
            &far,
            &no_radii,
            &no_multipliers,
        ));
        let nan = BTreeMap::from([(0, vec![DrawColour::new(f64::NAN, 0.0, 0.0)])]);
        assert_value_error(d.draw_molecule_with_highlights(
            &m,
            "",
            &nan,
            &none,
            &no_radii,
            &no_multipliers,
        ));
        let atoms = BTreeMap::from([(0, vec![red()])]);
        for r in [0.0, 17.0, f64::NAN] {
            let radii = BTreeMap::from([(0, r)]);
            assert_value_error(d.draw_molecule_with_highlights(
                &m,
                "",
                &atoms,
                &none,
                &radii,
                &no_multipliers,
            ));
        }
        let multipliers = BTreeMap::from([(0, 100)]);
        assert_value_error(d.draw_molecule_with_highlights(
            &m,
            "",
            &atoms,
            &none,
            &no_radii,
            &multipliers,
        ));
        assert_value_error(d.draw_molecule_with_highlights(
            &m,
            "<sub>x",
            &atoms,
            &none,
            &no_radii,
            &no_multipliers,
        ));
        let opts = MolDrawOptions {
            highlight_radius: 0.0,
            ..opts
        };
        d.set_draw_options(&opts).unwrap();
        assert_value_error(d.draw_molecule_with_highlights(
            &m,
            "",
            &atoms,
            &none,
            &no_radii,
            &no_multipliers,
        ));
        // Nothing to highlight needs no radius.
        d.draw_molecule_with_highlights(&m, "", &none, &none, &no_radii, &no_multipliers)
            .unwrap();
    }
}

#[test]
fn grids() {
    let mols: Vec<Mol> = ["CCO", "c1ccccc1", "CN", "O=C=O"]
        .into_iter()
        .map(mol)
        .collect();
    let panels = MolDraw2DSVGParams {
        panel_width: 300,
        panel_height: 300,
        ..Default::default()
    };
    let mut d = MolDraw2DSVG::new_with_params(600, 600, &panels).unwrap();
    let params = DrawMoleculesParams {
        legends: Some(["a", "b", "c", "d"].map(String::from).to_vec()),
        highlight_atoms: Some(vec![vec![0], vec![], vec![1], vec![]]),
        ..Default::default()
    };
    d.draw_molecules_with_params(&mols, &params).unwrap();
    d.finish_drawing();
    let text = d.get_drawing_text().unwrap();
    // Python's counts.
    assert_eq!(text.matches("class='bond-0").count(), 7);
    assert_eq!(text.matches("#FF7F7F").count(), 4);
    // A drawer that has drawn holds molecules RDKit would lay out again.
    assert_value_error(d.draw_molecules(&mols));

    let fresh = || MolDraw2DSVG::new_with_params(600, 600, &panels).unwrap();
    fresh().draw_molecules(&[]).unwrap();
    let short = DrawMoleculesParams {
        legends: Some(vec!["a".into()]),
        ..Default::default()
    };
    assert_value_error(fresh().draw_molecules_with_params(&mols, &short));
    let bad_index = DrawMoleculesParams {
        highlight_atoms: Some(vec![vec![9], vec![], vec![], vec![]]),
        ..Default::default()
    };
    assert_index_error(fresh().draw_molecules_with_params(&mols, &bad_index));
    let bad_legend = DrawMoleculesParams {
        legends: Some(["a", "<sub>b", "c", "d"].map(String::from).to_vec()),
        ..Default::default()
    };
    assert_value_error(fresh().draw_molecules_with_params(&mols, &bad_legend));
    // RDKit puts molecule i in row i / (width / panelWidth).
    let wide = MolDraw2DSVGParams {
        panel_width: 700,
        panel_height: 300,
        ..Default::default()
    };
    assert_value_error(
        MolDraw2DSVG::new_with_params(600, 600, &wide)
            .unwrap()
            .draw_molecules(&mols),
    );
    assert_value_error(MolDraw2DSVG::new(-1, -1).unwrap().draw_molecules(&mols));
    // Empty lists are not given, as in Python.
    let empty = DrawMoleculesParams {
        legends: Some(vec![]),
        ..Default::default()
    };
    fresh().draw_molecules_with_params(&mols, &empty).unwrap();
}

#[test]
fn reactions() {
    use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
    let rxn = reaction_from_smarts("[C:1]=[O:2].[N:3]>>[C:1](-[N:3])-[O:2]").unwrap();
    let mut d = MolDraw2DSVG::new(600, 200).unwrap();
    let params = DrawReactionParams {
        highlight_by_reactant: true,
        highlight_colors_reactants: Some(vec![red(), blue()]),
        ..Default::default()
    };
    d.draw_reaction_with_params(&rxn, &params).unwrap();
    d.finish_drawing();
    let text = d.get_drawing_text().unwrap();
    assert_eq!(
        (
            text.matches("#FF0000").count(),
            text.matches("#0000FF").count()
        ),
        (12, 4)
    );
    assert_value_error(d.draw_reaction(&rxn));

    let mut flexi = MolDraw2DSVG::new(-1, -1).unwrap();
    flexi.draw_reaction(&rxn).unwrap();
    let (w, h) = (flexi.width(), flexi.height());
    assert!((w - 231).abs() <= 1 && (h - 25).abs() <= 1, "{w} {h}");

    let fresh = || MolDraw2DSVG::new(600, 200).unwrap();
    fresh()
        .draw_reaction(&reaction_from_smarts("CC>O>C=C").unwrap())
        .unwrap();
    // highlightByReactant with no colors: RDKit takes them modulo zero.
    let no_colours = DrawReactionParams {
        highlight_by_reactant: true,
        highlight_colors_reactants: Some(vec![]),
        ..Default::default()
    };
    assert_value_error(fresh().draw_reaction_with_params(&rxn, &no_colours));
    let bad_colour = DrawReactionParams {
        highlight_by_reactant: true,
        highlight_colors_reactants: Some(vec![DrawColour::new(2.0, 0.0, 0.0)]),
        ..Default::default()
    };
    assert_value_error(fresh().draw_reaction_with_params(&rxn, &bad_colour));
    // RDKit indexes confIds by position in each group: two reactants here.
    let short = DrawReactionParams {
        conf_ids: Some(vec![-1]),
        ..Default::default()
    };
    assert_value_error(fresh().draw_reaction_with_params(&rxn, &short));
    // Without reactants or products RDKit converts DBL_MAX to int.
    for smarts in [">>", ">O>"] {
        assert_value_error(fresh().draw_reaction(&reaction_from_smarts(smarts).unwrap()));
    }
}

#[test]
fn multi_colour_bond_highlights_need_bond_length() {
    // RDKit grows the atoms' highlights forever to meet the NaN lines of a
    // bond with no length.
    let mut m = depicted("CCO");
    let p = m.get_conformer().unwrap().get_atom_position(0).unwrap();
    m.get_conformer_mut()
        .unwrap()
        .set_atom_position(1, p)
        .unwrap();
    let atoms = BTreeMap::from([(0, vec![red()]), (1, vec![blue()])]);
    let bonds = BTreeMap::from([(0, vec![red()])]);
    let mut d = svg();
    assert_value_error(d.draw_molecule_with_highlights(
        &m,
        "",
        &atoms,
        &bonds,
        &BTreeMap::new(),
        &BTreeMap::new(),
    ));
}

#[test]
fn acs1996_drawings() {
    let size = |svg: &str| {
        let start = svg.find("width='").unwrap() + 7;
        let w: i32 = svg[start..].split("px").next().unwrap().parse().unwrap();
        let start = svg.find("height='").unwrap() + 8;
        let h: i32 = svg[start..].split("px").next().unwrap().parse().unwrap();
        (w, h)
    };
    let m = mol("c1ccccc1O");
    // Python's sizes and counts.
    let text = mol_to_acs1996_svg(&m).unwrap();
    assert_eq!((size(&text), text.matches("<path").count()), ((58, 27), 17));
    let params = MolToACS1996SVGParams {
        legend: "phenol".into(),
        highlight_atoms: vec![0, 1],
        ..Default::default()
    };
    assert_eq!(params.conf_id, -1);
    let text = mol_to_acs1996_svg_with_params(&m, &params).unwrap();
    assert_eq!(
        (
            size(&text),
            text.matches("<path").count(),
            text.matches("#FF7F7F").count()
        ),
        ((58, 30), 23, 4)
    );
    let mut d = MolDraw2DSVG::new(-1, -1).unwrap();
    draw_molecule_acs1996(&mut d, &m).unwrap();
    assert_eq!((d.width(), d.height()), (58, 27));
    let opts = d.draw_options();
    assert_close_to(opts.scaling_factor, 9.6, 1e-9);
    assert_eq!(opts.fixed_font_size, 10);
    if std::env::var_os("RDBASE").is_none() {
        assert_eq!(opts.font_file, "BuiltinRobotoRegular");
    }

    // RDKit's check of the mean bond length, and the conformer.
    assert_value_error(mol_to_acs1996_svg(&mol("C")));
    let bad_conf = MolToACS1996SVGParams {
        conf_id: 3,
        ..Default::default()
    };
    assert_value_error(mol_to_acs1996_svg_with_params(&m, &bad_conf));
    // A bond of 0.01 makes RDKit's scaling factor 1440, beyond the bounds
    // the options are checked to (1e-12 makes it 1.44e13).
    for len in [0.01, 1e-12] {
        let mut tiny = depicted("CC");
        let mut conf = tiny.get_conformer_mut().unwrap();
        conf.set_atom_position(0, Point3D::new(0.0, 0.0, 0.0))
            .unwrap();
        conf.set_atom_position(1, Point3D::new(len, 0.0, 0.0))
            .unwrap();
        assert_value_error(mol_to_acs1996_svg(&tiny));
    }
}
