// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Free drawing on MolDraw2D: lines, shapes, text, GetDrawCoords, and
// SetScale. Expected values are from Python's rdMolDraw2D. Each refused
// input has a test.

mod common;

use common::assert_close_to;
use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::*;
use rdkit_cxx::chem::{Mol, mol_from_smiles};
use rdkit_cxx::geometry::Point2D;

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn p(x: f64, y: f64) -> Point2D {
    Point2D::new(x, y)
}

fn svg() -> MolDraw2DSVG {
    MolDraw2DSVG::new(300, 300).unwrap()
}

fn assert_value_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    let err = result.unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}

fn assert_runtime_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    let err = result.unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
}

fn assert_point(actual: Point2D, x: f64, y: f64) {
    assert_close_to(actual.x, x, 1e-6);
    assert_close_to(actual.y, y, 1e-6);
}

const RED: DrawColour = DrawColour::new(1.0, 0.0, 0.0);
const GREEN: DrawColour = DrawColour::new(0.0, 1.0, 0.0);
const BLUE: DrawColour = DrawColour::new(0.0, 0.0, 1.0);

#[test]
fn shapes_and_text() {
    let mut d = svg();
    d.draw_line_with_raw_coords(&p(10.0, 10.0), &p(100.0, 100.0), true)
        .unwrap();
    d.draw_rect_with_raw_coords(&p(10.0, 10.0), &p(100.0, 50.0), true)
        .unwrap();
    d.draw_ellipse_with_raw_coords(&p(10.0, 10.0), &p(100.0, 50.0), true)
        .unwrap();
    d.draw_arc_with_raw_coords(&p(150.0, 150.0), 20.0, 0.0, 90.0, true)
        .unwrap();
    d.draw_triangle_with_raw_coords(&p(10.0, 10.0), &p(100.0, 50.0), &p(30.0, 80.0), true)
        .unwrap();
    d.draw_polygon_with_raw_coords(
        &[p(10.0, 10.0), p(100.0, 50.0), p(30.0, 80.0), p(5.0, 60.0)],
        true,
    )
    .unwrap();
    let raw = DrawArrowParams {
        raw_coords: true,
        ..Default::default()
    };
    d.draw_arrow_with_params(&p(10.0, 200.0), &p(200.0, 200.0), &raw)
        .unwrap();
    let red_arrow = DrawArrowParams {
        as_polygon: true,
        frac: 0.1,
        angle: std::f64::consts::FRAC_PI_4,
        color: RED,
        raw_coords: true,
    };
    d.draw_arrow_with_params(&p(10.0, 250.0), &p(200.0, 250.0), &red_arrow)
        .unwrap();
    let raw_wavy = DrawWavyLineParams {
        raw_coords: true,
        ..Default::default()
    };
    d.draw_wavy_line_with_params(&p(10.0, 280.0), &p(200.0, 280.0), BLUE, GREEN, &raw_wavy)
        .unwrap();
    // The text takes the current color, which the wavy line left blue.
    d.draw_string_with_raw_coords("CH<sub>3</sub>", &p(250.0, 30.0), true)
        .unwrap();
    d.draw_string_aligned_with_raw_coords("label", &p(250.0, 60.0), 1, true)
        .unwrap();
    d.finish_drawing();
    let text = d.get_drawing_text().unwrap();
    // Python's counts.
    assert_eq!(
        (
            text.matches("<path").count(),
            text.matches("<ellipse").count(),
            text.matches("#FF0000").count(),
            text.matches("#0000FF").count()
        ),
        (18, 1, 3, 9)
    );
    assert!(text.contains("d='M 10.0,10.0 L 100.0,100.0'"));
    assert!(text.contains("<path d='M10,280 C10,279.95 21.875,279.95 21.875,280 C21.875"));

    let mut d = MolDraw2DSVG::new_with_params(
        300,
        300,
        &MolDraw2DSVGParams {
            no_freetype: true,
            ..Default::default()
        },
    )
    .unwrap();
    d.draw_string_with_raw_coords("CH<sub>3</sub>", &p(250.0, 30.0), true)
        .unwrap();
    d.draw_string_aligned_with_raw_coords("label", &p(250.0, 60.0), 1, true)
        .unwrap();
    assert_eq!(d.get_drawing_text().unwrap().matches("<text").count(), 8);
}

#[test]
fn params_defaults_match_python() {
    let arrow = DrawArrowParams::default();
    assert!(!arrow.as_polygon && !arrow.raw_coords);
    assert_eq!(arrow.frac, 0.05);
    assert_eq!(arrow.angle, std::f64::consts::PI / 6.0);
    assert_eq!(arrow.color, DrawColour::new(0.0, 0.0, 0.0));
    let attachment = DrawAttachmentLineParams::default();
    assert_eq!(
        (attachment.len, attachment.n_segments, attachment.raw_coords),
        (1.0, 16, false)
    );
    let wavy = DrawWavyLineParams::default();
    assert_eq!(
        (wavy.n_segments, wavy.vert_offset, wavy.raw_coords),
        (16, 0.05, false)
    );
}

#[test]
fn draw_coords() {
    let m = mol("c1ccccc1O");
    let mut d = svg();
    // Nothing drawn: Python's RuntimeErrors (no scaling info, bad active
    // mol index).
    assert_runtime_error(d.get_draw_coords(&p(0.0, 0.0)));
    assert_runtime_error(d.get_draw_coords_int(0));
    assert_runtime_error(d.draw_line(&p(0.0, 0.0), &p(1.0, 1.0)));
    d.draw_molecule(&m).unwrap();
    assert_point(
        d.get_draw_coords(&p(0.0, 0.0)).unwrap(),
        117.23073350340938,
        150.0,
    );
    assert_point(
        d.get_draw_coords_int(0).unwrap(),
        134.26918908731096,
        218.86009842558735,
    );
    assert_runtime_error(d.get_draw_coords_int(7));
    assert_runtime_error(d.get_draw_coords_int(u32::MAX));
    d.draw_line(&p(0.0, 0.0), &p(1.0, 1.0)).unwrap();

    // The attachment line is in the molecule's frame whatever raw_coords
    // is, which RDKit passes as its waves' height (Python's paths).
    d.draw_attachment_line(&p(0.0, 0.0), &p(1.0, 0.0), RED)
        .unwrap();
    let waving = DrawAttachmentLineParams {
        raw_coords: true,
        ..Default::default()
    };
    d.draw_attachment_line_with_params(&p(0.0, 0.0), &p(1.0, 0.0), BLUE, &waving)
        .unwrap();
    let text = d.get_drawing_text().unwrap();
    assert!(
        text.contains("<path d='M170.239,123.496 C170.239,123.496 170.239,126.809 170.239,126.809")
    );
    assert!(
        text.contains("<path d='M170.239,123.496 C117.231,123.496 117.231,126.809 170.239,126.809")
    );
}

#[test]
fn set_scale() {
    let m = mol("c1ccccc1O");
    let mut d = svg();
    d.set_scale(300, 300, &p(0.0, 0.0), &p(10.0, 10.0)).unwrap();
    assert_point(d.get_draw_coords(&p(5.0, 5.0)).unwrap(), 150.0, 150.0);
    assert_point(
        d.get_draw_coords(&p(0.0, 0.0)).unwrap(),
        27.272727272727273,
        272.72727272727275,
    );
    // Free drawing in the molecule's frame needs no molecule now.
    d.draw_line(&p(0.0, 0.0), &p(1.0, 1.0)).unwrap();
    assert_runtime_error(d.get_draw_coords_int(0));
    d.set_scale_with_mol(200, 100, &p(-1.0, -2.0), &p(3.0, 4.0), &m)
        .unwrap();
    assert_point(
        d.get_draw_coords(&p(0.0, 0.0)).unwrap(),
        146.57015584415583,
        173.63636363636363,
    );
    d.draw_molecule(&m).unwrap();
    assert_point(
        d.get_draw_coords_int(0).unwrap(),
        145.95327272727272,
        181.35051962286352,
    );
    // A grid or reaction lays out every molecule the drawer holds anew.
    assert_value_error(d.draw_molecules(&[mol("CC")]));

    let mut d = svg();
    d.set_scale_with_mol(300, 300, &p(0.0, 0.0), &p(10.0, 10.0), &m)
        .unwrap();
    assert_point(
        d.get_draw_coords(&p(0.0, 0.0)).unwrap(),
        66.95699510070769,
        239.5199672306539,
    );

    // The largest scale, about 2.5e6 pixels an angstrom on this area.
    let mut d = svg();
    d.set_scale(300, 300, &p(0.0, 0.0), &p(1e-4, 1e-4)).unwrap();
    assert_point(
        d.get_draw_coords(&p(1.0, 1.0)).unwrap(),
        2454572.727272727,
        -2454272.7272727266,
    );
    d.draw_molecule(&m).unwrap();
    assert_point(
        d.get_draw_coords_int(0).unwrap(),
        788988.3116883128,
        3188820.804842706,
    );
    // Molecules drawn after take the scale: here, beyond 1e9 pixels.
    let mut d = svg();
    d.set_scale(32767, 32767, &p(0.0, 0.0), &p(1e-4, 1e-4))
        .unwrap();
    assert_value_error(d.draw_molecule(&mol(&"C".repeat(20))));
    d.draw_molecule(&mol("C")).unwrap();

    // RDKit's checks of the size (Python's RuntimeError), and the wrappers'.
    let mut d = svg();
    assert_runtime_error(d.set_scale(0, 300, &p(0.0, 0.0), &p(1.0, 1.0)));
    assert_runtime_error(d.set_scale(300, -1, &p(0.0, 0.0), &p(1.0, 1.0)));
    assert_value_error(d.set_scale(40000, 300, &p(0.0, 0.0), &p(1.0, 1.0)));
    for bad in [f64::NAN, f64::INFINITY, 1e6] {
        assert_value_error(d.set_scale(300, 300, &p(bad, 0.0), &p(1.0, 1.0)));
        assert_value_error(d.set_scale(300, 300, &p(0.0, 0.0), &p(1.0, bad)));
    }
    // The molecule is checked as for drawing.
    let opts = MolDrawOptions {
        atom_regions: vec![vec![0, 99]],
        ..Default::default()
    };
    d.set_draw_options(&opts).unwrap();
    assert!(
        d.set_scale_with_mol(300, 300, &p(0.0, 0.0), &p(1.0, 1.0), &m)
            .is_err()
    );
    // A flexible canvas has no panels to lay the transformation out in.
    let mut flexi = MolDraw2DSVG::new(-1, -1).unwrap();
    assert_value_error(flexi.set_scale(300, 300, &p(0.0, 0.0), &p(1.0, 1.0)));
}

#[test]
fn free_drawing_is_checked() {
    let mut d = svg();
    for bad in [f64::NAN, f64::INFINITY, 1e10] {
        assert_value_error(d.draw_line_with_raw_coords(&p(bad, 0.0), &p(1.0, 1.0), true));
        assert_value_error(d.draw_rect_with_raw_coords(&p(0.0, 0.0), &p(1.0, bad), true));
        assert_value_error(d.draw_string_with_raw_coords("x", &p(bad, 0.0), true));
    }
    // In the molecule's frame the bound is also on the pixels.
    d.set_scale(300, 300, &p(0.0, 0.0), &p(1e-4, 1e-4)).unwrap();
    assert_value_error(d.draw_line(&p(0.0, 0.0), &p(1e4, 0.0)));

    // RDKit's own checks, as Python raises them.
    let mut d = svg();
    assert_runtime_error(d.draw_polygon_with_raw_coords(&[p(0.0, 0.0), p(1.0, 1.0)], true));
    assert_value_error(d.draw_polygon_with_raw_coords(&[], true));
    let flat = DrawArrowParams {
        angle: 0.0,
        raw_coords: true,
        ..Default::default()
    };
    assert_runtime_error(d.draw_arrow_with_params(&p(0.0, 0.0), &p(10.0, 0.0), &flat));
    let one_segment = DrawWavyLineParams {
        n_segments: 1,
        raw_coords: true,
        ..Default::default()
    };
    assert_runtime_error(d.draw_wavy_line_with_params(
        &p(0.0, 0.0),
        &p(10.0, 0.0),
        RED,
        RED,
        &one_segment,
    ));
    // A width scaled by the last molecule drawn needs one.
    let scaled = MolDrawOptions {
        scale_bond_width: true,
        ..Default::default()
    };
    d.set_draw_options(&scaled).unwrap();
    assert_runtime_error(d.draw_line_with_raw_coords(&p(0.0, 0.0), &p(1.0, 1.0), true));
    d.set_draw_options(&MolDrawOptions::default()).unwrap();

    // RDKit rounds an odd count up: u32::MAX segments become none, and its
    // drawer then reads the first (a crash in Python).
    for n_segments in [u32::MAX, 100001] {
        let many = DrawWavyLineParams {
            n_segments,
            raw_coords: true,
            ..Default::default()
        };
        assert_value_error(d.draw_wavy_line_with_params(
            &p(0.0, 0.0),
            &p(10.0, 0.0),
            RED,
            RED,
            &many,
        ));
    }
    let raw_wavy = DrawWavyLineParams {
        raw_coords: true,
        ..Default::default()
    };
    assert_value_error(d.draw_wavy_line_with_params(
        &p(5.0, 5.0),
        &p(5.0, 5.0),
        RED,
        RED,
        &raw_wavy,
    ));
    // drawArc converts its step count to int (a NaN angle is undefined).
    for (a1, a2) in [(0.0, f64::NAN), (f64::NAN, 90.0), (0.0, 1e5), (90.0, 0.0)] {
        assert_value_error(d.draw_arc_with_raw_coords(&p(50.0, 50.0), 10.0, a1, a2, true));
    }
    assert_value_error(d.draw_arc_with_raw_coords(&p(50.0, 50.0), f64::NAN, 0.0, 90.0, true));
    // The attachment line is perpendicular to a line of no length.
    d.draw_molecule(&mol("CC")).unwrap();
    assert_value_error(d.draw_attachment_line(&p(1.0, 1.0), &p(1.0, 1.0), RED));
    let no_len = DrawAttachmentLineParams {
        len: 0.0,
        ..Default::default()
    };
    assert_value_error(d.draw_attachment_line_with_params(
        &p(0.0, 0.0),
        &p(1.0, 0.0),
        RED,
        &no_len,
    ));
    assert_value_error(d.set_colour(DrawColour::new(2.0, 0.0, 0.0)));
    assert_value_error(d.draw_attachment_line(
        &p(0.0, 0.0),
        &p(1.0, 0.0),
        DrawColour::new(-1.0, 0.0, 0.0),
    ));

    // Text: alignments Python takes, and text RDKit can lay out.
    assert_value_error(d.draw_string_aligned_with_raw_coords("x", &p(10.0, 10.0), 3, true));
    for text in ["<sub>x</sub>", "<sup>2</sup>", "", "a:b>c<d"] {
        assert_value_error(d.draw_string_with_raw_coords(text, &p(10.0, 10.0), true));
    }
    assert_value_error(d.draw_string_aligned_with_raw_coords("", &p(10.0, 10.0), 0, true));
    let mut ascii = MolDraw2DSVG::new_with_params(
        300,
        300,
        &MolDraw2DSVGParams {
            no_freetype: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_value_error(ascii.draw_string_aligned_with_raw_coords(
        "caf\u{e9}",
        &p(10.0, 10.0),
        0,
        true,
    ));
    d.draw_string_aligned_with_raw_coords("caf\u{e9}", &p(10.0, 10.0), 0, true)
        .unwrap();
    // Text laid out as a label is split by isupper, which takes ASCII.
    assert_value_error(d.draw_string_with_raw_coords("caf\u{e9}", &p(10.0, 10.0), true));
}

#[test]
fn flexible_canvases_need_a_drawing_first() {
    // Before its first drawing a flexible canvas has nothing to draw on:
    // Python's Cairo drawer crashes drawing text to its null context.
    let mut cairo = MolDraw2DCairo::new(-1, -1).unwrap();
    assert_value_error(cairo.draw_string_with_raw_coords("CO", &p(10.0, 10.0), true));
    assert_value_error(cairo.draw_line_with_raw_coords(&p(0.0, 0.0), &p(1.0, 1.0), true));
    let mut svg = MolDraw2DSVG::new(-1, -1).unwrap();
    assert_value_error(svg.draw_rect_with_raw_coords(&p(0.0, 0.0), &p(1.0, 1.0), true));
    cairo.draw_molecule(&mol("CO")).unwrap();
    cairo
        .draw_string_with_raw_coords("CO", &p(10.0, 10.0), true)
        .unwrap();
    cairo.get_drawing_text().unwrap();
}

#[test]
fn comic_mode_lines_need_length() {
    // Comic mode draws arrows and triangles by hand, which loops forever on
    // a line it shifts to no length (Python hangs on this triangle).
    let mut d = svg();
    let comic = MolDrawOptions {
        comic_mode: true,
        ..Default::default()
    };
    d.set_draw_options(&comic).unwrap();
    assert_value_error(d.draw_triangle_with_raw_coords(
        &p(10.0, 10.0),
        &p(10.0, 10.0),
        &p(50.0, 50.0),
        true,
    ));
    assert_value_error(d.draw_triangle_with_raw_coords(
        &p(10.0, 10.0),
        &p(10.04, 10.0),
        &p(50.0, 50.0),
        true,
    ));
    let raw = DrawArrowParams {
        raw_coords: true,
        ..Default::default()
    };
    assert_value_error(d.draw_arrow_with_params(&p(10.0, 10.0), &p(10.0, 10.0), &raw));
    d.draw_triangle_with_raw_coords(&p(10.0, 10.0), &p(100.0, 10.0), &p(50.0, 50.0), true)
        .unwrap();
    d.draw_arrow_with_params(&p(10.0, 200.0), &p(200.0, 200.0), &raw)
        .unwrap();
    // Without comic mode such lines are drawn.
    let mut plain = svg();
    plain
        .draw_triangle_with_raw_coords(&p(10.0, 10.0), &p(10.0, 10.0), &p(50.0, 50.0), true)
        .unwrap();
}
