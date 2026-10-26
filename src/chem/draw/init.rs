// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/Chem/Draw/__init__.py` that draws
//! molecules to images, re-exported from `rdkit_cxx::chem::draw` as Python's
//! is from `rdkit.Chem.Draw`. These are the functions Jupyter's display of
//! molecules (`IPythonConsole`) and the RDKit Cookbook's drawing recipes
//! call.
//!
//! Python returns PIL images, or PNG data with `returnPNG=True`; these
//! return the PNG data. The `MolDrawOptions` Python's `**kwargs` set by name
//! are set on the `draw_options` given instead.

use std::collections::BTreeMap;

use crate::chem::Mol;
use crate::chem::draw::rd_mol_draw_2d::{
    DrawColour, DrawMoleculeParams, DrawMoleculesParams, MolDraw2D, MolDraw2DCairo,
    MolDraw2DCairoParams, MolDraw2DSVG, MolDraw2DSVGParams, MolDrawOptions,
    PrepareMolForDrawingParams, prepare_mol_for_drawing_with_params,
};
use crate::rd_base::BlockLogs;
use crate::{Error, Result};

/// Whether to kekulize for drawing: not when an aromatic bond is a query,
/// which kekulization would lose. Mirrors `Draw.shouldKekulize`.
pub fn should_kekulize(mol: &Mol, kekulize: bool) -> bool {
    kekulize
        && !mol
            .get_bonds()
            .any(|bond| bond.get_is_aromatic() && bond.has_query())
}

/// The molecule as Python's `_prepareMol` leaves it for drawing: prepared
/// (kekulized where it can be, wedged, with coordinates) unless the options
/// say molecules are drawn as they are.
fn prepare_mol(
    mol: &Mol,
    draw_options: Option<&MolDrawOptions>,
    kekulize: bool,
    wedge_bonds: bool,
) -> Result<Mol> {
    // Python updates the caller's molecule; this updates a copy.
    let mut mol = mol.clone();
    if mol.needs_update_property_cache() {
        mol.update_property_cache_with_strict(false)?;
    }
    let kekulize = should_kekulize(&mol, kekulize);
    if draw_options.is_some_and(|options| !options.prepare_mols_before_drawing) {
        return Ok(mol);
    }
    let prepared = {
        let _block = BlockLogs::default();
        prepare_mol_for_drawing_with_params(
            &mol,
            &PrepareMolForDrawingParams {
                kekulize,
                wedge_bonds,
                ..Default::default()
            },
        )
    };
    match prepared {
        // Python catches the ValueError of a molecule that will not
        // kekulize, and draws it aromatic.
        Err(
            Error::Kekulize { .. }
            | Error::AtomKekulize { .. }
            | Error::MolSanitize { .. }
            | Error::AtomValence { .. }
            | Error::Value { .. },
        ) => prepare_mol_for_drawing_with_params(
            &mol,
            &PrepareMolForDrawingParams {
                kekulize: false,
                wedge_bonds,
                ..Default::default()
            },
        ),
        prepared => prepared,
    }
}

/// Sets the options on a drawer, with molecules drawn as given (they have
/// been prepared already).
fn set_options(
    drawer: &mut MolDraw2D,
    draw_options: Option<&MolDrawOptions>,
    highlight_color: Option<DrawColour>,
) -> Result<()> {
    let mut options = match draw_options {
        Some(options) => options.clone(),
        None => drawer.draw_options(),
    };
    if let Some(colour) = highlight_color {
        options.set_highlight_colour(colour)?;
    }
    options.prepare_mols_before_drawing = false;
    drawer.set_draw_options(&options)
}

/// Parameters for `mol_to_image_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolToImageParams {
    /// The image's width and height, in pixels.
    pub size: (i32, i32),
    /// Kekulizes the molecule before drawing.
    pub kekulize: bool,
    /// Draws wedge (stereo) bonds.
    pub wedge_bonds: bool,
    /// The drawing options; `None` for the defaults.
    pub options: Option<MolDrawOptions>,
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight.
    pub highlight_bonds: Vec<u32>,
    /// The highlight colour; `None` for the options'.
    pub highlight_color: Option<DrawColour>,
    /// Text drawn under (or by) the molecule.
    pub legend: String,
}

impl Default for MolToImageParams {
    fn default() -> MolToImageParams {
        MolToImageParams {
            size: (300, 300),
            kekulize: true,
            wedge_bonds: true,
            options: None,
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            highlight_color: None,
            legend: String::new(),
        }
    }
}

/// A PNG drawing of the molecule. Mirrors `Draw.MolToImage`, which returns
/// a PIL image of it.
pub fn mol_to_image(mol: &Mol) -> Result<Vec<u8>> {
    mol_to_image_with_params(mol, &MolToImageParams::default())
}

/// Mirrors `Draw.MolToImage(mol, size, kekulize, wedgeBonds, fitImage,
/// options, highlightAtoms=..., highlightBonds=..., highlightColor=...,
/// legend=...)`. Python ignores `fitImage`, which is not ported.
pub fn mol_to_image_with_params(mol: &Mol, params: &MolToImageParams) -> Result<Vec<u8>> {
    let mol = prepare_mol(
        mol,
        params.options.as_ref(),
        params.kekulize,
        params.wedge_bonds,
    )?;
    let mut drawer = MolDraw2DCairo::new(params.size.0, params.size.1)?;
    set_options(&mut drawer, params.options.as_ref(), params.highlight_color)?;
    drawer.draw_molecule_with_params(
        &mol,
        &DrawMoleculeParams {
            highlight_atoms: params.highlight_atoms.clone(),
            highlight_bonds: Some(params.highlight_bonds.clone()),
            legend: params.legend.clone(),
            ..Default::default()
        },
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Parameters for `mol_to_svg_with_params`, with Python's defaults (A2).
/// Python's `wedgeBonds` is not ported: its `MolToSVG` ignores it, and
/// always wedges.
#[derive(Clone, Debug, PartialEq)]
pub struct MolToSVGParams {
    /// The image's width and height, in pixels.
    pub size: (i32, i32),
    /// Kekulizes the molecule before drawing.
    pub kekulize: bool,
    /// The drawing options; `None` for the defaults.
    pub draw_options: Option<MolDrawOptions>,
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight.
    pub highlight_bonds: Vec<u32>,
    /// Text drawn under (or by) the molecule.
    pub legend: String,
}

impl Default for MolToSVGParams {
    fn default() -> MolToSVGParams {
        MolToSVGParams {
            size: (300, 300),
            kekulize: true,
            draw_options: None,
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            legend: String::new(),
        }
    }
}

/// An SVG drawing of the molecule. Mirrors `Draw.MolToSVG`. (Jupyter's
/// display of molecules, `IPythonConsole`, draws them this way, at its
/// `molSize`, with its `drawOptions`.)
pub fn mol_to_svg(mol: &Mol) -> Result<String> {
    mol_to_svg_with_params(mol, &MolToSVGParams::default())
}

/// Mirrors `Draw.MolToSVG(mol, size, kekulize, wedgeBonds, drawOptions,
/// highlightAtoms=..., highlightBonds=..., legend=...)`.
pub fn mol_to_svg_with_params(mol: &Mol, params: &MolToSVGParams) -> Result<String> {
    let mol = prepare_mol(mol, params.draw_options.as_ref(), params.kekulize, true)?;
    let mut drawer = MolDraw2DSVG::new(params.size.0, params.size.1)?;
    set_options(&mut drawer, params.draw_options.as_ref(), None)?;
    drawer.draw_molecule_with_params(
        &mol,
        &DrawMoleculeParams {
            highlight_atoms: params.highlight_atoms.clone(),
            highlight_bonds: Some(params.highlight_bonds.clone()),
            legend: params.legend.clone(),
            ..Default::default()
        },
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Parameters for `mols_to_grid_image_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolsToGridImageParams {
    /// The molecules in each row.
    pub mols_per_row: u32,
    /// Each molecule's panel width and height, in pixels.
    pub sub_img_size: (i32, i32),
    /// Each molecule's legend; `None` for none.
    pub legends: Option<Vec<String>>,
    /// The atoms to highlight in each molecule.
    pub highlight_atom_lists: Option<Vec<Vec<u32>>>,
    /// The bonds to highlight in each molecule.
    pub highlight_bond_lists: Option<Vec<Vec<u32>>>,
    /// The drawing options; `None` for the defaults.
    pub draw_options: Option<MolDrawOptions>,
}

impl Default for MolsToGridImageParams {
    fn default() -> MolsToGridImageParams {
        MolsToGridImageParams {
            mols_per_row: 3,
            sub_img_size: (200, 200),
            legends: None,
            highlight_atom_lists: None,
            highlight_bond_lists: None,
            draw_options: None,
        }
    }
}

/// A PNG drawing of the molecules in a grid. Mirrors
/// `Draw.MolsToGridImage`, which returns a PIL image of it (or its PNG data,
/// with `returnPNG=True`).
pub fn mols_to_grid_image(mols: &[Mol]) -> Result<Vec<u8>> {
    mols_to_grid_image_with_params(mols, &MolsToGridImageParams::default())
}

/// Mirrors `Draw.MolsToGridImage(mols, molsPerRow, subImgSize, legends,
/// highlightAtomLists, highlightBondLists, drawOptions=...)`. Lists longer
/// than `mols` are cut to its length, as in Python; a shorter one is
/// `Error::Value`, as `MolDraw2D::draw_molecules_with_params` says, and so
/// are no molecules or no molecules per row (Python divides by zero).
pub fn mols_to_grid_image_with_params(
    mols: &[Mol],
    params: &MolsToGridImageParams,
) -> Result<Vec<u8>> {
    let (width, height, panel_width, panel_height) = grid_size(mols, params)?;
    let mut drawer = MolDraw2DCairo::new_with_params(
        width,
        height,
        &MolDraw2DCairoParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )?;
    draw_grid(&mut drawer, mols, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// An SVG drawing of the molecules in a grid. Mirrors
/// `Draw.MolsToGridImage(mols, useSVG=True)`.
pub fn mols_to_grid_image_use_svg(mols: &[Mol]) -> Result<String> {
    mols_to_grid_image_use_svg_with_params(mols, &MolsToGridImageParams::default())
}

/// Mirrors `Draw.MolsToGridImage(mols, ..., useSVG=True)`, with the errors
/// of `mols_to_grid_image_with_params`.
pub fn mols_to_grid_image_use_svg_with_params(
    mols: &[Mol],
    params: &MolsToGridImageParams,
) -> Result<String> {
    let (width, height, panel_width, panel_height) = grid_size(mols, params)?;
    let mut drawer = MolDraw2DSVG::new_with_params(
        width,
        height,
        &MolDraw2DSVGParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )?;
    draw_grid(&mut drawer, mols, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// The grid's canvas and panel sizes.
fn grid_size(mols: &[Mol], params: &MolsToGridImageParams) -> Result<(i32, i32, i32, i32)> {
    if params.mols_per_row == 0 {
        return Err(Error::Value {
            message: "molsPerRow must be positive".to_owned(),
        });
    }
    let n_rows = mols.len().div_ceil(params.mols_per_row as usize);
    let (sub_width, sub_height) = params.sub_img_size;
    let size = |count: usize, sub: i32| {
        i32::try_from(count)
            .ok()
            .and_then(|count| count.checked_mul(sub))
            .ok_or_else(|| Error::Value {
                message: "the grid is too large to draw".to_owned(),
            })
    };
    Ok((
        size(params.mols_per_row as usize, sub_width)?,
        size(n_rows, sub_height)?,
        sub_width,
        sub_height,
    ))
}

fn draw_grid(drawer: &mut MolDraw2D, mols: &[Mol], params: &MolsToGridImageParams) -> Result<()> {
    if let Some(options) = &params.draw_options {
        drawer.set_draw_options(options)?;
    }
    fn truncated<T: Clone>(list: &Option<Vec<T>>, len: usize) -> Option<Vec<T>> {
        list.as_ref()
            .map(|items| items.iter().take(len).cloned().collect())
    }
    let legends =
        truncated(&params.legends, mols.len()).unwrap_or_else(|| vec![String::new(); mols.len()]);
    drawer.draw_molecules_with_params(
        mols,
        &DrawMoleculesParams {
            legends: Some(legends),
            highlight_atoms: truncated(&params.highlight_atom_lists, mols.len()),
            highlight_bonds: truncated(&params.highlight_bond_lists, mols.len()),
            ..Default::default()
        },
    )
}

/// Parameters for `draw_mol_with_matches_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug)]
pub struct DrawMolWithMatchesParams {
    /// The image's width and height, in pixels.
    pub mol_size: (i32, i32),
    /// The query that matched; with it, the bonds it matched are
    /// highlighted too.
    pub qry: Option<Mol>,
    /// Text drawn under the molecule.
    pub label: String,
    /// The drawing options.
    pub options: MolDrawOptions,
    /// A colour for each match, in turn; `None` for the options' highlight
    /// colour.
    pub colors: Option<Vec<DrawColour>>,
    /// Kekulizes the molecule before drawing.
    pub kekulize: bool,
    /// The conformer to draw; -1 for the first.
    pub conf_id: i32,
}

impl Default for DrawMolWithMatchesParams {
    fn default() -> DrawMolWithMatchesParams {
        DrawMolWithMatchesParams {
            mol_size: (350, 300),
            qry: None,
            label: String::new(),
            options: MolDrawOptions::default(),
            colors: None,
            kekulize: true,
            conf_id: -1,
        }
    }
}

/// A PNG drawing of the molecule with each substructure match (a list of
/// atom indices) highlighted. Mirrors `Draw.DrawMolWithMatches`. (Jupyter's
/// display of molecules draws a molecule that holds the matches of its last
/// substructure search this way.)
pub fn draw_mol_with_matches(mol: &Mol, matches: &[Vec<u32>]) -> Result<Vec<u8>> {
    draw_mol_with_matches_with_params(mol, matches, &DrawMolWithMatchesParams::default())
}

/// Mirrors `Draw.DrawMolWithMatches(mol, matches, molSize, qry, label,
/// options, colors, kekulize=..., confId=...)`. `Error::Value` for an empty
/// `colors` (Python divides by zero) or a match that does not fit the query.
pub fn draw_mol_with_matches_with_params(
    mol: &Mol,
    matches: &[Vec<u32>],
    params: &DrawMolWithMatchesParams,
) -> Result<Vec<u8>> {
    let mut drawer = MolDraw2DCairo::new(params.mol_size.0, params.mol_size.1)?;
    draw_with_matches(&mut drawer, mol, matches, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// An SVG drawing of the molecule with each match highlighted. Mirrors
/// `Draw.DrawMolWithMatches(mol, matches, doPNG=False)`.
pub fn draw_mol_with_matches_as_svg(mol: &Mol, matches: &[Vec<u32>]) -> Result<String> {
    draw_mol_with_matches_as_svg_with_params(mol, matches, &DrawMolWithMatchesParams::default())
}

/// Mirrors `Draw.DrawMolWithMatches(mol, matches, ..., doPNG=False)`, with
/// the errors of `draw_mol_with_matches_with_params`.
pub fn draw_mol_with_matches_as_svg_with_params(
    mol: &Mol,
    matches: &[Vec<u32>],
    params: &DrawMolWithMatchesParams,
) -> Result<String> {
    let mut drawer = MolDraw2DSVG::new(params.mol_size.0, params.mol_size.1)?;
    draw_with_matches(&mut drawer, mol, matches, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Adds `colour` to each element's colours, once.
fn add_colors_to_map(els: &[u32], cols: &mut BTreeMap<u32, Vec<DrawColour>>, colour: DrawColour) {
    for &el in els {
        let colours = cols.entry(el).or_default();
        if !colours.contains(&colour) {
            colours.push(colour);
        }
    }
}

fn draw_with_matches(
    drawer: &mut MolDraw2D,
    mol: &Mol,
    matches: &[Vec<u32>],
    params: &DrawMolWithMatchesParams,
) -> Result<()> {
    let options = &params.options;
    let mol = prepare_mol(mol, Some(options), params.kekulize, true)?;
    let mut acols = BTreeMap::new();
    let mut bcols = BTreeMap::new();
    let mut h_rads = BTreeMap::new();
    let mut h_lw_mult = BTreeMap::new();
    for (i, alist) in matches.iter().enumerate() {
        let mut blist = Vec::new();
        if let Some(qry) = &params.qry {
            for bnd in qry.get_bonds() {
                let atom = |idx: u32| {
                    alist.get(idx as usize).copied().ok_or(Error::Index {
                        index: idx,
                        len: alist.len() as u32,
                    })
                };
                let b = mol
                    .get_bond_between_atoms(
                        atom(bnd.get_begin_atom_idx())?,
                        atom(bnd.get_end_atom_idx())?,
                    )?
                    .ok_or_else(|| Error::Value {
                        message: "a match lacks a bond of the query".to_owned(),
                    })?;
                blist.push(b.get_idx());
            }
        }
        for &ai in alist {
            h_rads.insert(ai, options.highlight_radius);
        }
        for &bi in &blist {
            h_lw_mult.insert(bi, options.highlight_bond_width_multiplier);
        }
        let colour = match &params.colors {
            Some(colors) if colors.is_empty() => {
                return Err(Error::Value {
                    message: "colors must not be empty".to_owned(),
                });
            }
            Some(colors) => colors[i % colors.len()],
            None => options.get_highlight_colour(),
        };
        add_colors_to_map(alist, &mut acols, colour);
        add_colors_to_map(&blist, &mut bcols, colour);
    }
    let mut options = options.clone();
    options.prepare_mols_before_drawing = false;
    drawer.set_draw_options(&options)?;
    drawer.draw_molecule_with_highlights_with_conf_id(
        &mol,
        &params.label,
        &acols,
        &bcols,
        &h_rads,
        &h_lw_mult,
        params.conf_id,
    )
}
