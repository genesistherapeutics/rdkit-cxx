// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// ContourParams, ContourAndDrawGrid, and ContourAndDrawGaussians. Expected
// values are from Python's rdMolDraw2D. Each refused input has a test.

mod common;

use common::assert_close_to;
use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::*;
use rdkit_cxx::chem::rd_depictor::compute_2d_coords;
use rdkit_cxx::chem::{Mol, mol_from_smiles};
use rdkit_cxx::geometry::Point2D;

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

// Python's xs = [0.1 * i for i in range(21)], ys for range(11), and
// data[i][j] = sin(xs[i]) * cos(2 * ys[j]).
fn grid() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let xs: Vec<f64> = (0..21).map(|i| 0.1 * i as f64).collect();
    let ys: Vec<f64> = (0..11).map(|i| 0.1 * i as f64).collect();
    let mut data = Vec::new();
    for x in &xs {
        for y in &ys {
            data.push(x.sin() * (2.0 * y).cos());
        }
    }
    (xs, ys, data)
}

fn counts(d: &mut MolDraw2DSVG) -> (usize, usize) {
    d.finish_drawing();
    let text = d.get_drawing_text().unwrap();
    (
        text.matches("<path").count(),
        text.matches("stroke-dasharray").count(),
    )
}

/// Checks a drawing's paths and dashed paths against Python's counts on
/// x86-64. The contour lines are traced through floating-point crossings,
/// which other platforms split differently (arm64 and macOS draw a piece or
/// two more or fewer), so the counts may differ by a few.
fn assert_counts(d: &mut MolDraw2DSVG, (paths, dashed): (usize, usize)) {
    let (actual_paths, actual_dashed) = counts(d);
    assert!(
        actual_paths.abs_diff(paths) <= 2 + paths / 10,
        "{actual_paths} paths, not about {paths}"
    );
    assert!(
        actual_dashed.abs_diff(dashed) <= 2,
        "{actual_dashed} dashed paths, not about {dashed}"
    );
    assert_eq!(actual_dashed == 0, dashed == 0);
}

// Phenol, depicted, and gaussians on its atoms.
fn phenol() -> (Mol, Vec<Point2D>, Vec<f64>, Vec<f64>) {
    let mut m = mol_from_smiles("c1ccccc1O").unwrap();
    compute_2d_coords(&mut m).unwrap();
    let n = m.get_num_atoms();
    let conf = m.get_conformer().unwrap();
    let locs = (0..n)
        .map(|i| {
            let p = conf.get_atom_position(i).unwrap();
            Point2D::new(p.x, p.y)
        })
        .collect();
    let heights = (0..n).map(|i| 0.1 * i as f64 - 0.3).collect();
    let widths = vec![0.4; n as usize];
    (m, locs, heights, widths)
}

#[test]
fn contour_params_defaults_match_python() {
    let ps = ContourParams::default();
    assert!(ps.set_scale && ps.dash_negative && !ps.fill_grid && ps.draw_as_lines);
    assert!(!ps.use_fill_threshold && ps.fill_threshold_is_fraction);
    assert_eq!(
        (
            ps.grid_resolution,
            ps.contour_width,
            ps.extra_grid_padding,
            ps.coord_scale_for_quantization,
            ps.isoval_scale_for_quantization,
            ps.fill_threshold
        ),
        (0.15, 1.0, 0.0, 1000.0, 1e6, 0.01)
    );
    assert_eq!(ps.contour_colour, DrawColour::rgba(0.5, 0.5, 0.5, 0.5));
    assert_eq!(
        ps.colour_map,
        [
            DrawColour::rgba(0.557, 0.004, 0.322, 0.5),
            DrawColour::rgba(1.0, 1.0, 1.0, 0.5),
            DrawColour::rgba(0.153, 0.392, 0.098, 0.5)
        ]
    );
    let mut ps = ps;
    ps.set_contour_colour(DrawColour::new(1.0, 0.0, 0.0));
    ps.set_colour_map(&[
        DrawColour::new(0.0, 0.0, 0.0),
        DrawColour::new(1.0, 1.0, 1.0),
    ]);
    assert_eq!(ps.contour_colour, DrawColour::new(1.0, 0.0, 0.0));
    assert_eq!(ps.colour_map.len(), 2);
}

#[test]
fn contours_of_a_grid() {
    let (xs, ys, data) = grid();
    // Python's counts of paths and dashed paths, and the scale set.
    let mut d = svg();
    contour_and_draw_grid(&mut d, &data, &xs, &ys).unwrap();
    assert_point(
        d.get_draw_coords(&Point2D::new(0.0, 0.0)).unwrap(),
        27.272727272727273,
        211.36363636363637,
    );
    assert_counts(&mut d, (17, 4));

    let mut d = svg();
    let fill = ContourAndDrawGridParams {
        n_contours: 5,
        params: ContourParams {
            fill_grid: true,
            ..Default::default()
        },
        ..Default::default()
    };
    contour_and_draw_grid_with_params(&mut d, &data, &xs, &ys, &fill).unwrap();
    assert_counts(&mut d, (192, 1));

    // Given levels and no fill, RDKit finds no range and draws nothing.
    let mut d = svg();
    let levels = ContourAndDrawGridParams {
        n_contours: 3,
        levels: vec![-0.5, 0.0, 0.5],
        params: ContourParams {
            draw_as_lines: false,
            ..Default::default()
        },
        ..Default::default()
    };
    contour_and_draw_grid_with_params(&mut d, &data, &xs, &ys, &levels).unwrap();
    assert_eq!(counts(&mut d), (0, 0));
    // Nor does it for a constant grid, which needs no scale.
    let unscaled = ContourAndDrawGridParams {
        params: ContourParams {
            set_scale: false,
            ..Default::default()
        },
        ..Default::default()
    };
    contour_and_draw_grid_with_params(
        &mut svg(),
        &[0.0; 9],
        &[0.0, 1.0, 2.0],
        &[0.0, 1.0, 2.0],
        &unscaled,
    )
    .unwrap();
}

#[test]
fn contours_of_gaussians() {
    let (m, locs, heights, widths) = phenol();
    let mut d = svg();
    contour_and_draw_gaussians(&mut d, &locs, &heights, &widths).unwrap();
    assert_point(
        d.get_draw_coords(&Point2D::new(0.0, 0.0)).unwrap(),
        165.93860684769777,
        150.0,
    );
    assert_counts(&mut d, (24, 12));

    let mut d = svg();
    let params = ContourAndDrawGaussiansParams {
        n_contours: 6,
        params: ContourParams {
            fill_grid: true,
            grid_resolution: 0.3,
            ..Default::default()
        },
        mol: Some(&m),
        ..Default::default()
    };
    contour_and_draw_gaussians_with_params(&mut d, &locs, &heights, &widths, &params).unwrap();
    d.draw_molecule(&m).unwrap();
    assert_point(
        d.get_draw_coords(&Point2D::new(0.0, 0.0)).unwrap(),
        173.72888034295804,
        150.0,
    );
    assert_counts(&mut d, (199, 6));
}

#[test]
fn grids_are_checked() {
    let (xs, ys, data) = grid();
    // Python's checks, and RDKit's (RuntimeError in Python).
    assert_value_error(contour_and_draw_grid(&mut svg(), &data, &[], &ys));
    assert_value_error(contour_and_draw_grid(&mut svg(), &data, &xs, &[]));
    assert_value_error(contour_and_draw_grid(&mut svg(), &data[1..], &xs, &ys));
    assert_runtime_error(contour_and_draw_grid(
        &mut svg(),
        &[0.0, 1.0, 2.0],
        &[0.0],
        &[0.0, 1.0, 2.0],
    ));
    let one_colour = ContourAndDrawGridParams {
        params: ContourParams {
            colour_map: vec![DrawColour::new(1.0, 0.0, 0.0)],
            ..Default::default()
        },
        ..Default::default()
    };
    assert_runtime_error(contour_and_draw_grid_with_params(
        &mut svg(),
        &data,
        &xs,
        &ys,
        &one_colour,
    ));
    let few_levels = ContourAndDrawGridParams {
        n_contours: 20,
        levels: vec![0.1, 0.2],
        params: ContourParams {
            fill_grid: true,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_value_error(contour_and_draw_grid_with_params(
        &mut svg(),
        &data,
        &xs,
        &ys,
        &few_levels,
    ));
    // Without set_scale, the grid is placed by what was drawn.
    let unscaled = ContourAndDrawGridParams {
        params: ContourParams {
            set_scale: false,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_runtime_error(contour_and_draw_grid_with_params(
        &mut svg(),
        &data,
        &xs,
        &ys,
        &unscaled,
    ));
    let mut drawn = svg();
    drawn
        .draw_molecule(&mol_from_smiles("CO").unwrap())
        .unwrap();
    contour_and_draw_grid_with_params(&mut drawn, &data, &xs, &ys, &unscaled).unwrap();

    // Values RDKit's fill would make NaN colours of, which the SVG drawer
    // converts to int.
    for bad in [f64::NAN, f64::INFINITY, 1e301] {
        let mut values = data.clone();
        values[3] = bad;
        assert_value_error(contour_and_draw_grid(&mut svg(), &values, &xs, &ys));
        let mut bad_xs = xs.clone();
        bad_xs[3] = bad;
        assert_value_error(contour_and_draw_grid(&mut svg(), &data, &bad_xs, &ys));
    }
    let nan_level = ContourAndDrawGridParams {
        n_contours: 2,
        levels: vec![0.0, f64::NAN],
        ..Default::default()
    };
    assert_value_error(contour_and_draw_grid_with_params(
        &mut svg(),
        &data,
        &xs,
        &ys,
        &nan_level,
    ));
    let bad_params = [
        ContourParams {
            colour_map: vec![
                DrawColour::new(1.0, 0.0, 0.0),
                DrawColour::new(2.0, 0.0, 0.0),
            ],
            ..Default::default()
        },
        ContourParams {
            contour_colour: DrawColour::new(f64::NAN, 0.0, 0.0),
            ..Default::default()
        },
        ContourParams {
            contour_width: f64::NAN,
            ..Default::default()
        },
        ContourParams {
            coord_scale_for_quantization: f64::INFINITY,
            ..Default::default()
        },
    ];
    for params in bad_params {
        let p = ContourAndDrawGridParams {
            params,
            ..Default::default()
        };
        assert_value_error(contour_and_draw_grid_with_params(
            &mut svg(),
            &data,
            &xs,
            &ys,
            &p,
        ));
    }
    let many = ContourAndDrawGridParams {
        n_contours: 100001,
        ..Default::default()
    };
    assert_value_error(contour_and_draw_grid_with_params(
        &mut svg(),
        &data,
        &xs,
        &ys,
        &many,
    ));
    // The grid's corners are SetScale's.
    let far: Vec<f64> = xs.iter().map(|x| x * 1e6).collect();
    assert_value_error(contour_and_draw_grid(&mut svg(), &data, &far, &ys));
    // A flexible canvas has nothing to draw on before its first drawing.
    let mut flexi = MolDraw2DSVG::new(-1, -1).unwrap();
    assert_value_error(contour_and_draw_grid(&mut flexi, &data, &xs, &ys));
}

#[test]
fn gaussians_are_checked() {
    let (_, locs, heights, widths) = phenol();
    assert_value_error(contour_and_draw_gaussians(
        &mut svg(),
        &[],
        &heights,
        &widths,
    ));
    assert_value_error(contour_and_draw_gaussians(&mut svg(), &locs, &[], &widths));
    assert_value_error(contour_and_draw_gaussians(&mut svg(), &locs, &heights, &[]));
    assert_runtime_error(contour_and_draw_gaussians(
        &mut svg(),
        &locs,
        &heights[..3],
        &widths,
    ));
    // RDKit divides by each width and its square.
    for bad in [0.0, -1.0, 1e-101, f64::NAN] {
        let mut w = widths.clone();
        w[0] = bad;
        assert_value_error(contour_and_draw_gaussians(&mut svg(), &locs, &heights, &w));
    }
    let mut h = heights.clone();
    h[0] = f64::INFINITY;
    assert_value_error(contour_and_draw_gaussians(&mut svg(), &locs, &h, &widths));
    let mut far = locs.clone();
    far[0] = Point2D::new(1e6, 0.0);
    assert_value_error(contour_and_draw_gaussians(
        &mut svg(),
        &far,
        &heights,
        &widths,
    ));
    // RDKit converts the grid's size from double and allocates it.
    for resolution in [0.0, -0.1, f64::NAN, 1e-6] {
        let p = ContourAndDrawGaussiansParams {
            params: ContourParams {
                grid_resolution: resolution,
                ..Default::default()
            },
            ..Default::default()
        };
        assert_value_error(contour_and_draw_gaussians_with_params(
            &mut svg(),
            &locs,
            &heights,
            &widths,
            &p,
        ));
    }
    let padded = ContourAndDrawGaussiansParams {
        params: ContourParams {
            extra_grid_padding: f64::NAN,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_value_error(contour_and_draw_gaussians_with_params(
        &mut svg(),
        &locs,
        &heights,
        &widths,
        &padded,
    ));
    // Without set_scale it needs something drawn to size the grid by.
    let unscaled = ContourAndDrawGaussiansParams {
        params: ContourParams {
            set_scale: false,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_runtime_error(contour_and_draw_gaussians_with_params(
        &mut svg(),
        &locs,
        &heights,
        &widths,
        &unscaled,
    ));
}
