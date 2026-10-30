// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// MolDrawOptions and the functions that change it. Expected values are from
// Python's rdMolDraw2D; Python reads color tuples as 32-bit floats, so
// colors are compared with a tolerance.

mod common;

use std::collections::BTreeMap;

use common::assert_close_to;
use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::*;

fn assert_colour(actual: DrawColour, expected: (f64, f64, f64, f64)) {
    for (a, e) in [
        (actual.r, expected.0),
        (actual.g, expected.1),
        (actual.b, expected.2),
        (actual.a, expected.3),
    ] {
        assert_close_to(a, e, 1e-6);
    }
}

#[test]
fn defaults_match_python() {
    let o = MolDrawOptions::default();
    assert!(!o.add_atom_indices);
    assert_eq!(o.annotation_font_scale, 0.5);
    assert_colour(o.background_colour, (1.0, 1.0, 1.0, 1.0));
    assert_eq!(o.base_font_size, -1.0);
    assert_eq!(o.bond_line_width, 2.0);
    assert_colour(o.bond_note_colour, (0.5, 0.5, 1.0, 1.0));
    assert!(o.circle_atoms);
    assert!(o.continuous_highlight);
    assert_eq!(o.drawing_extents_include, DrawElement::ALL);
    assert_eq!(o.drawing_extents_include.bits(), 2147483647);
    assert_eq!(o.fixed_bond_length, -1.0);
    assert_eq!(o.fixed_font_size, -1);
    assert_eq!(o.fixed_scale, -1.0);
    assert_eq!(o.flag_close_contacts_dist, 3);
    assert_eq!(o.font_file, "");
    assert_eq!(o.highlight_bond_width_multiplier, 8);
    assert_colour(o.highlight_colour, (1.0, 0.5, 0.5, 1.0));
    assert_eq!(o.highlight_radius, 0.3);
    assert!(o.include_metadata);
    assert_eq!(o.legend_font_size, 16);
    assert_eq!(o.legend_fraction, 0.1);
    assert_eq!(o.max_font_size, 40);
    assert_eq!(o.min_font_size, 6);
    assert_eq!(
        o.multi_colour_highlight_style,
        MultiColourHighlightStyle::CircleAndLine
    );
    assert_eq!(o.multiple_bond_offset, 0.15);
    assert_eq!(o.padding, 0.05);
    assert!(o.prepare_mols_before_drawing);
    assert_colour(o.query_colour, (0.5, 0.5, 0.5, 1.0));
    assert_eq!(o.reagent_padding, 0.0);
    assert!(o.scale_highlight_bond_width);
    assert_eq!(o.scaling_factor, 20.0);
    assert!(o.use_complex_query_atom_symbols);
    assert_eq!(o.variable_atom_radius, 0.4);
    assert_colour(o.variable_attachment_colour, (0.8, 0.8, 0.8, 1.0));
    assert_eq!(o.variable_bond_width_multiplier, 16);
    assert!(o.atom_labels.is_empty());
    assert!(o.atom_regions.is_empty());
    let palette = o.get_atom_palette();
    assert_eq!(
        palette.keys().copied().collect::<Vec<_>>(),
        [-1, 0, 1, 6, 7, 8, 9, 15, 16, 17, 35, 53, 201]
    );
    assert_colour(palette[&7], (0.0, 0.0, 1.0, 1.0));
    assert_colour(palette[&17], (0.0, 0.802, 0.0, 1.0));
}

#[test]
fn draw_element_flags() {
    assert_eq!(DrawElement::NONE.bits(), 0);
    assert_eq!(DrawElement::BONDS.bits(), 2);
    assert_eq!(DrawElement::ATOMLABELS.bits(), 4);
    assert_eq!(DrawElement::POSTSHAPES.bits(), 64);
    let both = DrawElement::BONDS | DrawElement::ATOMLABELS;
    assert_eq!(both.bits(), 6);
    assert_eq!(format!("{:?}", DrawElement::BONDS), "DrawElement::BONDS");
}

#[test]
fn palettes() {
    let mut o = MolDrawOptions::default();
    o.use_bw_atom_palette();
    assert_eq!(o.get_atom_palette().len(), 1);
    assert_colour(o.get_atom_palette()[&-1], (0.0, 0.0, 0.0, 1.0));
    o.use_cdk_atom_palette();
    assert_colour(o.get_atom_palette()[&5], (1.0, 0.71, 0.71, 1.0));
    o.use_avalon_atom_palette();
    assert_colour(o.get_atom_palette()[&9], (0.0, 0.498, 0.0, 1.0));
    o.use_default_atom_palette();
    assert_eq!(o, MolDrawOptions::default());

    let mut update = BTreeMap::new();
    update.insert(6, DrawColour::new(0.5, 0.5, 0.5));
    o.update_atom_palette(&update);
    assert_eq!(o.get_atom_palette().len(), 13);
    assert_colour(o.get_atom_palette()[&6], (0.5, 0.5, 0.5, 1.0));
    o.set_atom_palette(&update);
    assert_eq!(o.get_atom_palette(), update);
}

#[test]
fn colour_setters_check_their_range() {
    let mut o = MolDrawOptions::default();
    let err = o
        .set_background_colour(DrawColour::new(1.5, 0.0, 0.0))
        .unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    assert_colour(o.get_background_colour(), (1.0, 1.0, 1.0, 1.0));
    assert!(
        o.set_highlight_colour(DrawColour::rgba(0.1, 0.2, f64::NAN, 1.0))
            .is_err()
    );
    o.set_background_colour((0.5, 0.25, 0.0).into()).unwrap();
    assert_colour(o.get_background_colour(), (0.5, 0.25, 0.0, 1.0));
    o.set_legend_colour((0.1, 0.2, 0.3, 0.4).into()).unwrap();
    assert_colour(o.legend_colour, (0.1, 0.2, 0.3, 0.4));
}

#[test]
fn dark_and_monochrome_modes() {
    let mut o = MolDrawOptions::default();
    set_dark_mode(&mut o);
    assert_colour(o.get_background_colour(), (0.0, 0.0, 0.0, 1.0));
    assert_colour(o.get_symbol_colour(), (0.9, 0.9, 0.9, 1.0));
    assert_colour(o.get_atom_palette()[&7], (0.33, 0.41, 0.92, 1.0));

    let mut o = MolDrawOptions::default();
    let fg = DrawColour::new(0.2, 0.3, 0.4);
    set_monochrome_mode(&mut o, fg, DrawColour::new(0.9, 0.9, 0.9)).unwrap();
    assert_eq!(o.get_atom_palette().len(), 1);
    assert_eq!(o.get_atom_palette()[&-1], fg);
    assert_eq!(o.get_legend_colour(), fg);
    assert_eq!(o.get_bond_note_colour(), fg);
    assert_colour(o.get_background_colour(), (0.9, 0.9, 0.9, 1.0));
    let before = o.clone();
    let err = set_monochrome_mode(&mut o, DrawColour::new(2.0, 0.0, 0.0), fg).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    assert_eq!(o, before);
}

#[test]
fn acs1996_mode() {
    let mut o = MolDrawOptions::default();
    set_acs1996_mode(&mut o, 1.5).unwrap();
    assert_close_to(o.bond_line_width, 0.6, 1e-12);
    assert_close_to(o.scaling_factor, 9.6, 1e-12);
    assert_close_to(o.fixed_bond_length, 9.6, 1e-12);
    assert_close_to(o.multiple_bond_offset, 0.18, 1e-12);
    assert_eq!(o.highlight_bond_width_multiplier, 32);
    assert_eq!(o.fixed_font_size, 10);
    assert_close_to(o.additional_atom_label_padding, 0.066, 1e-12);
    assert!(!o.scale_bond_width);
    // FreeSans.ttf under $RDBASE, when that is set and has it.
    assert!(!o.font_file.is_empty());
    let err = set_acs1996_mode(&mut o, 0.0).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}

#[test]
fn update_from_json() {
    let mut o = MolDrawOptions::default();
    update_mol_draw_options_from_json(
        &mut o,
        r#"{"padding": 0.2, "legendFontSize": 20, "atomLabels": {"0": "X", "3": "Y"},
            "highlightColour": [0.1, 0.2, 0.3], "atomColourPalette": {"6": [0.5, 0.5, 0.5, 0.5]},
            "multiColourHighlightStyle": "Lasso",
            "drawingExtentsInclude": {"BONDS": true, "ATOMLABELS": true},
            "colourScheme": "dark", "fontFile": "foo", "padding_bad": 1, "rotate": "abc"}"#,
    )
    .unwrap();
    assert_eq!(o.padding, 0.2);
    assert_eq!(o.legend_font_size, 20);
    assert_eq!(
        o.atom_labels,
        BTreeMap::from([(0, "X".to_string()), (3, "Y".to_string())])
    );
    assert_colour(o.highlight_colour, (0.1, 0.2, 0.3, 1.0));
    assert_colour(o.get_atom_palette()[&6], (0.5, 0.5, 0.5, 0.5));
    // RDKit 2025.09 reads no colourScheme, so it is ignored like any unknown
    // key.
    assert_colour(o.get_atom_palette()[&7], (0.0, 0.0, 1.0, 1.0));
    assert_eq!(
        o.multi_colour_highlight_style,
        MultiColourHighlightStyle::Lasso
    );
    assert_eq!(
        o.drawing_extents_include,
        DrawElement::BONDS | DrawElement::ATOMLABELS
    );
    assert_colour(o.background_colour, (1.0, 1.0, 1.0, 1.0));
    assert_eq!(o.font_file, "foo");
    // A value of the wrong type is ignored.
    assert_eq!(o.rotate, 0.0);

    let mut o = MolDrawOptions::default();
    update_mol_draw_options_from_json(&mut o, r#"{"atomColourPalette": "cdk"}"#).unwrap();
    assert_colour(o.get_atom_palette()[&7], (0.188, 0.314, 0.972, 1.0));
    update_mol_draw_options_from_json(&mut o, "").unwrap();
}

// RDKit's get_rgba steps past the end of a color with fewer than three
// components (Python segfaults); the wrapper refuses it.
#[test]
fn json_colours_need_three_components() {
    for json in [
        r#"{"highlightColour": [1]}"#,
        r#"{"highlightColour": 5}"#,
        r#"{"backgroundColour": [1, 0]}"#,
        r#"{"legendColour": []}"#,
        r#"{"atomColourPalette": {"6": [1]}}"#,
        r#"{"queryColour": [0, 0, 0], "queryColour": [0]}"#,
    ] {
        let mut o = MolDrawOptions::default();
        let err = update_mol_draw_options_from_json(&mut o, json).unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{json}: {err:?}");
        assert_eq!(o, MolDrawOptions::default());
    }
}

#[test]
fn json_errors() {
    let mut o = MolDrawOptions::default();
    let err = update_mol_draw_options_from_json(&mut o, "{not json").unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
    // Atom label keys are atom indices.
    let err =
        update_mol_draw_options_from_json(&mut o, r#"{"atomLabels": {"-1": "X"}}"#).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    // Boost's parser refuses lone surrogates.
    let err = update_mol_draw_options_from_json(&mut o, r#"{"fontFile": "\udc80"}"#).unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
    assert_eq!(o, MolDrawOptions::default());
}
