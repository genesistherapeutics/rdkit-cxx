// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Contour plots: `ContourParams`, `contour_and_draw_grid`, and
//! `contour_and_draw_gaussians`.

use super::drawer::MolDraw2D;
use super::options::DrawColour;
use crate::Result;
use crate::chem::Mol;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::error::expect_ok;
use crate::ffi::draw::ffi;
use crate::geometry::Point2D;

/// Parameters for drawing contours. Mirrors
/// `rdMolDraw2D.ContourParams`; `Default` is RDKit's.
#[derive(Clone, Debug, PartialEq)]
pub struct ContourParams {
    /// Sets the drawer's scale to the grid (with `set_scale`), as for
    /// drawing the grid first.
    pub set_scale: bool,
    /// Draws negative contours dashed.
    pub dash_negative: bool,
    /// Fills the grid's boxes with colors from `colour_map`.
    pub fill_grid: bool,
    /// The spacing of the grid `contour_and_draw_gaussians` computes, in
    /// the molecule's frame.
    pub grid_resolution: f64,
    /// The contours' line width.
    pub contour_width: f64,
    /// Extra space around the gaussians' grid, in the molecule's frame.
    pub extra_grid_padding: f64,
    /// The contours' color.
    pub contour_colour: DrawColour,
    /// The fill's colors: two, from the lowest value to the highest, or
    /// three, for negative values, zero, and positive values.
    pub colour_map: Vec<DrawColour>,
    /// Draws the contours as continuous lines rather than segments.
    pub draw_as_lines: bool,
    /// The scale by which RDKit rounds coordinates to join segments into
    /// lines.
    pub coord_scale_for_quantization: f64,
    /// The scale by which RDKit rounds levels to join segments into lines.
    pub isoval_scale_for_quantization: f64,
    /// Fills only boxes whose value's magnitude reaches `fill_threshold`.
    pub use_fill_threshold: bool,
    /// The magnitude a box's value must reach to be filled.
    pub fill_threshold: f64,
    /// Takes `fill_threshold` as a fraction of the data's range.
    pub fill_threshold_is_fraction: bool,
}

impl Default for ContourParams {
    fn default() -> ContourParams {
        ContourParams::from_ffi(expect_ok(ffi::contour_params_default()))
    }
}

impl ContourParams {
    /// Sets `contour_colour`. Mirrors `ContourParams.setContourColour`.
    pub fn set_contour_colour(&mut self, colour: DrawColour) {
        self.contour_colour = colour;
    }

    /// Sets `colour_map`. Mirrors `ContourParams.setColourMap`.
    pub fn set_colour_map(&mut self, colours: &[DrawColour]) {
        self.colour_map = colours.to_vec();
    }

    fn from_ffi(data: ffi::ContourParamsData) -> ContourParams {
        ContourParams {
            set_scale: data.set_scale,
            dash_negative: data.dash_negative,
            fill_grid: data.fill_grid,
            grid_resolution: data.grid_resolution,
            contour_width: data.contour_width,
            extra_grid_padding: data.extra_grid_padding,
            contour_colour: data.contour_colour,
            colour_map: data.colour_map,
            draw_as_lines: data.draw_as_lines,
            coord_scale_for_quantization: data.coord_scale_for_quantization,
            isoval_scale_for_quantization: data.isoval_scale_for_quantization,
            use_fill_threshold: data.use_fill_threshold,
            fill_threshold: data.fill_threshold,
            fill_threshold_is_fraction: data.fill_threshold_is_fraction,
        }
    }

    fn to_ffi(&self) -> ffi::ContourParamsData {
        ffi::ContourParamsData {
            set_scale: self.set_scale,
            dash_negative: self.dash_negative,
            fill_grid: self.fill_grid,
            grid_resolution: self.grid_resolution,
            contour_width: self.contour_width,
            extra_grid_padding: self.extra_grid_padding,
            contour_colour: self.contour_colour,
            colour_map: self.colour_map.clone(),
            draw_as_lines: self.draw_as_lines,
            coord_scale_for_quantization: self.coord_scale_for_quantization,
            isoval_scale_for_quantization: self.isoval_scale_for_quantization,
            use_fill_threshold: self.use_fill_threshold,
            fill_threshold: self.fill_threshold,
            fill_threshold_is_fraction: self.fill_threshold_is_fraction,
        }
    }
}

/// Parameters for `contour_and_draw_grid_with_params`, with Python's
/// defaults (A2, A3: `params` is the `ContourParams` Python takes).
#[derive(Clone, Debug)]
pub struct ContourAndDrawGridParams<'a> {
    /// The number of contours.
    pub n_contours: u32,
    /// The contours' levels, in increasing order; empty for `n_contours`
    /// levels spread over the data's range.
    pub levels: Vec<f64>,
    /// How the contours are drawn.
    pub params: ContourParams,
    /// A molecule to fit the scale to as well, when `params.set_scale` is
    /// set, as for drawing it over the plot.
    pub mol: Option<&'a Mol>,
}

impl Default for ContourAndDrawGridParams<'_> {
    fn default() -> Self {
        ContourAndDrawGridParams {
            n_contours: 10,
            levels: Vec::new(),
            params: ContourParams::default(),
            mol: None,
        }
    }
}

/// Parameters for `contour_and_draw_gaussians_with_params`, with Python's
/// defaults (A2, A3).
#[derive(Clone, Debug)]
pub struct ContourAndDrawGaussiansParams<'a> {
    /// The number of contours.
    pub n_contours: u32,
    /// The contours' levels, in increasing order; empty for `n_contours`
    /// levels spread over the grid's range.
    pub levels: Vec<f64>,
    /// How the grid is computed and the contours drawn.
    pub params: ContourParams,
    /// A molecule to fit the scale to as well, when `params.set_scale` is
    /// set.
    pub mol: Option<&'a Mol>,
}

impl Default for ContourAndDrawGaussiansParams<'_> {
    fn default() -> Self {
        ContourAndDrawGaussiansParams {
            n_contours: 10,
            levels: Vec::new(),
            params: ContourParams::default(),
            mol: None,
        }
    }
}

/// Draws contours of values on a grid, in the molecule's frame: `data`
/// holds a value for each x and y coordinate, x-major (`data[i *
/// ycoords.len() + j]` at `xcoords[i]`, `ycoords[j]`), as Python's array of
/// shape `(len(xcoords), len(ycoords))`. Mirrors
/// `rdMolDraw2D.ContourAndDrawGrid(drawer, data, xcoords, ycoords)`; see
/// `contour_and_draw_grid_with_params` for its errors.
pub fn contour_and_draw_grid(
    drawer: &mut MolDraw2D,
    data: &[f64],
    xcoords: &[f64],
    ycoords: &[f64],
) -> Result<()> {
    contour_and_draw_grid_with_params(
        drawer,
        data,
        xcoords,
        ycoords,
        &ContourAndDrawGridParams::default(),
    )
}

/// As `contour_and_draw_grid`, with levels, parameters, or a molecule.
/// Mirrors `rdMolDraw2D.ContourAndDrawGrid`. As in RDKit, given levels and
/// no fill it draws nothing (it finds the data's range only to fill or to
/// make levels).
///
/// As in Python, `Error::Value` for no coordinates, data of another size,
/// or more contours than levels, and `Error::Runtime` (RDKit's checks) for
/// fewer than two colors in `colour_map`, fewer than two coordinates
/// along either axis to contour, or nothing drawn to place the grid by
/// without `set_scale`. Also `Error::Value` for values or levels that are
/// not finite or beyond 1e300 (an infinite one makes RDKit's fill colors
/// NaN), colors outside 0 to 1, a contour width outside the line width
/// bounds (RDKit sets it as the drawer's), more than 100000 contours,
/// coordinates beyond the free-drawing bounds, and as `MolDraw2D::set_scale`
/// when `params.set_scale` is set.
pub fn contour_and_draw_grid_with_params(
    drawer: &mut MolDraw2D,
    data: &[f64],
    xcoords: &[f64],
    ycoords: &[f64],
    params: &ContourAndDrawGridParams,
) -> Result<()> {
    let mols: Vec<(&str, &Mol)> = params.mol.map(|m| ("", m)).into_iter().collect();
    with_borrowed_mols_shared(mols, |list| {
        Ok(ffi::contour_and_draw_grid(
            drawer.pin_mut(),
            data,
            xcoords,
            ycoords,
            params.n_contours,
            &params.levels,
            &params.params.to_ffi(),
            list,
        )?)
    })
}

/// Draws contours of a sum of gaussians, on a grid of
/// `params.grid_resolution` over the drawing (after setting its scale to
/// the gaussians, by default). Mirrors
/// `rdMolDraw2D.ContourAndDrawGaussians(drawer, locs, heights, widths)`;
/// see `contour_and_draw_gaussians_with_params` for its errors.
pub fn contour_and_draw_gaussians(
    drawer: &mut MolDraw2D,
    locs: &[Point2D],
    heights: &[f64],
    widths: &[f64],
) -> Result<()> {
    contour_and_draw_gaussians_with_params(
        drawer,
        locs,
        heights,
        widths,
        &ContourAndDrawGaussiansParams::default(),
    )
}

/// As `contour_and_draw_gaussians`, with levels, parameters, or a
/// molecule. Mirrors `rdMolDraw2D.ContourAndDrawGaussians`.
///
/// As in Python, `Error::Value` for no gaussians, and `Error::Runtime`
/// (RDKit's check) for lists of different lengths. Also `Error::Value` for
/// locations beyond 1e5, heights beyond 1e100, widths outside 1e-100 to
/// 1e100 (RDKit divides by them, and by their squares), a grid resolution
/// that is not positive, a grid of more than 1e8 points (RDKit converts the
/// count from double and allocates it), and as
/// `contour_and_draw_grid_with_params`.
pub fn contour_and_draw_gaussians_with_params(
    drawer: &mut MolDraw2D,
    locs: &[Point2D],
    heights: &[f64],
    widths: &[f64],
    params: &ContourAndDrawGaussiansParams,
) -> Result<()> {
    let mols: Vec<(&str, &Mol)> = params.mol.map(|m| ("", m)).into_iter().collect();
    with_borrowed_mols_shared(mols, |list| {
        Ok(ffi::contour_and_draw_gaussians(
            drawer.pin_mut(),
            locs,
            heights,
            widths,
            params.n_contours,
            &params.levels,
            &params.params.to_ffi(),
            list,
        )?)
    })
}
