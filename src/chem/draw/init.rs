// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/Chem/Draw/__init__.py`, re-exported
//! from `rdkit_cxx::chem::draw` as Python's is from `rdkit.Chem.Draw`:
//! drawings of molecules, grids of them, and reactions to images or files,
//! the drawings of fingerprint bits' environments, and `DebugDraw`.
//!
//! Python returns PIL images, or PNG data with `returnPNG=True`; these
//! return the PNG data. The `MolDrawOptions` Python's `**kwargs` set by name
//! are set on the `draw_options` given instead. Where Python changes the
//! caller's molecule for a drawing (computing coordinates, noting atoms),
//! these change a copy.

use std::collections::BTreeMap;
use std::path::Path;

use crate::chem::draw::rd_mol_draw_2d::{
    DrawColour, DrawMoleculeParams, DrawMoleculesParams, DrawReactionParams, MolDraw2D,
    MolDraw2DCairo, MolDraw2DCairoParams, MolDraw2DSVG, MolDraw2DSVGParams, MolDrawOptions,
    PrepareMolForDrawingParams, prepare_mol_for_drawing_with_params,
};
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::rd_depictor::compute_2d_coords;
use crate::chem::{
    Conformer, Mol, fast_find_rings, find_atom_environment_of_radius_n, mol_fragment_to_smiles,
    mol_from_smiles, path_to_submol_with_outputs,
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
    /// The highlight color; `None` for the options'.
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
    let (width, height, panel_width, panel_height) = grid_size(mols.len(), params)?;
    let mut drawer = MolDraw2DCairo::new_with_params(
        width,
        height,
        &MolDraw2DCairoParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )?;
    let mols: Vec<Option<&Mol>> = mols.iter().map(Some).collect();
    draw_grid(&mut drawer, &mols, params)?;
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
    let (width, height, panel_width, panel_height) = grid_size(mols.len(), params)?;
    let mut drawer = MolDraw2DSVG::new_with_params(
        width,
        height,
        &MolDraw2DSVGParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )?;
    let mols: Vec<Option<&Mol>> = mols.iter().map(Some).collect();
    draw_grid(&mut drawer, &mols, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// The grid's canvas and panel sizes.
fn grid_size(num_mols: usize, params: &MolsToGridImageParams) -> Result<(i32, i32, i32, i32)> {
    if params.mols_per_row == 0 {
        return Err(Error::Value {
            message: "molsPerRow must be positive".to_owned(),
        });
    }
    let n_rows = num_mols.div_ceil(params.mols_per_row as usize);
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

fn draw_grid(
    drawer: &mut MolDraw2D,
    mols: &[Option<&Mol>],
    params: &MolsToGridImageParams,
) -> Result<()> {
    if let Some(options) = &params.draw_options {
        drawer.set_draw_options(options)?;
    }
    fn truncated<T: Clone>(list: &Option<Vec<T>>, len: usize) -> Option<Vec<T>> {
        list.as_ref()
            .map(|items| items.iter().take(len).cloned().collect())
    }
    let legends =
        truncated(&params.legends, mols.len()).unwrap_or_else(|| vec![String::new(); mols.len()]);
    drawer.draw_optional_molecules_with_params(
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
    /// A color for each match, in turn; `None` for the options' highlight
    /// color.
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

/// Adds `colour` to each element's colors, once.
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
                    alist
                        .get(idx as usize)
                        .copied()
                        .ok_or(Error::index(idx, alist.len()))
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

/// Parameters for `mol_to_file_with_params`, with Python's defaults (A2).
/// Python ignores `fitImage`, which is not ported.
#[derive(Clone, Debug, PartialEq)]
pub struct MolToFileParams {
    /// The image's width and height, in pixels.
    pub size: (i32, i32),
    /// Kekulizes the molecule before drawing.
    pub kekulize: bool,
    /// Draws wedge (stereo) bonds, in PNG files (Python's SVG ignores it).
    pub wedge_bonds: bool,
    /// `"png"` or `"svg"`; `None` for the file name's extension.
    pub image_type: Option<String>,
    /// The drawing options; `None` for the defaults.
    pub options: Option<MolDrawOptions>,
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight.
    pub highlight_bonds: Vec<u32>,
    /// Text drawn under (or by) the molecule.
    pub legend: String,
}

impl Default for MolToFileParams {
    fn default() -> MolToFileParams {
        MolToFileParams {
            size: (300, 300),
            kekulize: true,
            wedge_bonds: true,
            image_type: None,
            options: None,
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            legend: String::new(),
        }
    }
}

/// Writes a drawing of the molecule to a PNG or SVG file, by the file
/// name's extension. Mirrors `Draw.MolToFile`.
pub fn mol_to_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_file_with_params(mol, filename, &MolToFileParams::default())
}

/// Mirrors `Draw.MolToFile(mol, filename, size, kekulize, wedgeBonds,
/// imageType, fitImage, options, highlightAtoms=..., highlightBonds=...,
/// legend=...)`. `Error::Value` for an empty file name or an image type
/// other than `png` and `svg` (the extension is matched as Python's
/// `os.path.splitext` finds it, case and all), and `Error::Os` when the
/// file cannot be written.
pub fn mol_to_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolToFileParams,
) -> Result<()> {
    let filename = filename.as_ref();
    if filename.as_os_str().is_empty() {
        return Err(Error::Value {
            message: "no filename provided".to_owned(),
        });
    }
    let image_type = match &params.image_type {
        Some(image_type) => image_type.clone(),
        None => filename
            .extension()
            .map(|ext| ext.to_string_lossy().into_owned())
            .unwrap_or_default(),
    };
    let data = match image_type.as_str() {
        "png" => mol_to_image_with_params(
            mol,
            &MolToImageParams {
                size: params.size,
                kekulize: params.kekulize,
                wedge_bonds: params.wedge_bonds,
                options: params.options.clone(),
                highlight_atoms: params.highlight_atoms.clone(),
                highlight_bonds: params.highlight_bonds.clone(),
                highlight_color: None,
                legend: params.legend.clone(),
            },
        )?,
        "svg" => mol_to_svg_with_params(
            mol,
            &MolToSVGParams {
                size: params.size,
                kekulize: params.kekulize,
                draw_options: params.options.clone(),
                highlight_atoms: params.highlight_atoms.clone(),
                highlight_bonds: params.highlight_bonds.clone(),
                legend: params.legend.clone(),
            },
        )?
        .into_bytes(),
        _ => {
            return Err(Error::Value {
                message: "only PNG and SVG image types are supported".to_owned(),
            });
        }
    };
    std::fs::write(filename, data).map_err(|err| Error::Os {
        message: format!("{}: {err}", filename.display()),
    })
}

/// Parameters for `mols_matrix_to_grid_image_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolsMatrixToGridImageParams {
    /// Each molecule's panel width and height, in pixels.
    pub sub_img_size: (i32, i32),
    /// Each molecule's legend, row by row; `None` for none.
    pub legends_matrix: Option<Vec<Vec<String>>>,
    /// The atoms to highlight in each molecule, row by row.
    pub highlight_atom_lists_matrix: Option<Vec<Vec<Vec<u32>>>>,
    /// The bonds to highlight in each molecule, row by row.
    pub highlight_bond_lists_matrix: Option<Vec<Vec<Vec<u32>>>>,
    /// The drawing options; `None` for the defaults.
    pub draw_options: Option<MolDrawOptions>,
}

impl Default for MolsMatrixToGridImageParams {
    fn default() -> MolsMatrixToGridImageParams {
        MolsMatrixToGridImageParams {
            sub_img_size: (200, 200),
            legends_matrix: None,
            highlight_atom_lists_matrix: None,
            highlight_bond_lists_matrix: None,
            draw_options: None,
        }
    }
}

/// A PNG drawing of the molecules in a grid, a row of the grid for each row
/// of `mols_matrix`, padded to the longest; `None` leaves a panel empty.
/// Mirrors `Draw.MolsMatrixToGridImage`, which returns a PIL image of it
/// (or its PNG data, with `returnPNG=True`).
pub fn mols_matrix_to_grid_image(mols_matrix: &[Vec<Option<&Mol>>]) -> Result<Vec<u8>> {
    mols_matrix_to_grid_image_with_params(mols_matrix, &MolsMatrixToGridImageParams::default())
}

/// Mirrors `Draw.MolsMatrixToGridImage(molsMatrix, subImgSize,
/// legendsMatrix, highlightAtomListsMatrix, highlightBondListsMatrix,
/// drawOptions=...)`. As in Python, `Error::Value` for a matrix whose rows
/// differ from `mols_matrix`'s in number or length, and for no rows; also
/// for rows that are all empty (Python divides by zero), and as
/// `MolDraw2D::draw_molecules_with_params`.
pub fn mols_matrix_to_grid_image_with_params(
    mols_matrix: &[Vec<Option<&Mol>>],
    params: &MolsMatrixToGridImageParams,
) -> Result<Vec<u8>> {
    let (mols, grid) = mols_nested_to_linear(mols_matrix, params)?;
    let (width, height, panel_width, panel_height) = grid_size(mols.len(), &grid)?;
    let mut drawer = MolDraw2DCairo::new_with_params(
        width,
        height,
        &MolDraw2DCairoParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )?;
    draw_grid(&mut drawer, &mols, &grid)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// An SVG drawing of the molecules in a grid. Mirrors
/// `Draw.MolsMatrixToGridImage(molsMatrix, useSVG=True)`.
pub fn mols_matrix_to_grid_image_use_svg(mols_matrix: &[Vec<Option<&Mol>>]) -> Result<String> {
    mols_matrix_to_grid_image_use_svg_with_params(
        mols_matrix,
        &MolsMatrixToGridImageParams::default(),
    )
}

/// Mirrors `Draw.MolsMatrixToGridImage(molsMatrix, ..., useSVG=True)`, with
/// the errors of `mols_matrix_to_grid_image_with_params`.
pub fn mols_matrix_to_grid_image_use_svg_with_params(
    mols_matrix: &[Vec<Option<&Mol>>],
    params: &MolsMatrixToGridImageParams,
) -> Result<String> {
    let (mols, grid) = mols_nested_to_linear(mols_matrix, params)?;
    let (width, height, panel_width, panel_height) = grid_size(mols.len(), &grid)?;
    let mut drawer = MolDraw2DSVG::new_with_params(
        width,
        height,
        &MolDraw2DSVGParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )?;
    draw_grid(&mut drawer, &mols, &grid)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Checks that a matrix has `mols_matrix`'s shape, with Python's messages.
fn check_matrix_shape<T>(
    mols_matrix: &[Vec<Option<&Mol>>],
    matrix: &Option<Vec<Vec<T>>>,
    name: &str,
) -> Result<()> {
    let Some(matrix) = matrix else {
        return Ok(());
    };
    let n_mols_rows = mols_matrix.len();
    if matrix.len() != n_mols_rows {
        return Err(Error::Value {
            message: format!(
                "If {name} is provided it must be the same length (have the same number of \
                 sub-iterables) as molsMatrix, {n_mols_rows}; its length is {}.",
                matrix.len()
            ),
        });
    }
    for (row_index, (row, mols_row)) in matrix.iter().zip(mols_matrix).enumerate() {
        if row.len() != mols_row.len() {
            return Err(Error::Value {
                message: format!(
                    "If {name} is provided each of its sub-iterables must be the same length as \
                     the corresponding sub-iterable of molsMatrix. For sub-iterable of index \
                     {row_index}, its length in molsMatrix is {} while its length in {name} is {}.",
                    mols_row.len(),
                    row.len()
                ),
            });
        }
    }
    Ok(())
}

/// Python's `_MolsNestedToLinear`: the molecules, padded row by row with
/// `None`, and the grid parameters of `MolsToGridImage`.
fn mols_nested_to_linear<'a>(
    mols_matrix: &[Vec<Option<&'a Mol>>],
    params: &MolsMatrixToGridImageParams,
) -> Result<(Vec<Option<&'a Mol>>, MolsToGridImageParams)> {
    check_matrix_shape(mols_matrix, &params.legends_matrix, "legendsMatrix")?;
    check_matrix_shape(
        mols_matrix,
        &params.highlight_atom_lists_matrix,
        "highlightAtomListsMatrix",
    )?;
    check_matrix_shape(
        mols_matrix,
        &params.highlight_bond_lists_matrix,
        "highlightBondListsMatrix",
    )?;
    let mols_per_row = mols_matrix
        .iter()
        .map(Vec::len)
        .max()
        .ok_or_else(|| Error::Value {
            message: "max() arg is an empty sequence".to_owned(),
        })?;
    fn pad<T: Clone>(matrix: &[Vec<T>], len: usize, with: T) -> Vec<T> {
        matrix
            .iter()
            .flat_map(|row| {
                row.iter()
                    .cloned()
                    .chain(std::iter::repeat_n(with.clone(), len - row.len()))
            })
            .collect()
    }
    let mols = pad(mols_matrix, mols_per_row, None);
    let grid = MolsToGridImageParams {
        mols_per_row: u32::try_from(mols_per_row).map_err(|_| Error::Value {
            message: "the grid is too large to draw".to_owned(),
        })?,
        sub_img_size: params.sub_img_size,
        legends: params
            .legends_matrix
            .as_ref()
            .map(|m| pad(m, mols_per_row, String::new())),
        highlight_atom_lists: params
            .highlight_atom_lists_matrix
            .as_ref()
            .map(|m| pad(m, mols_per_row, Vec::new())),
        highlight_bond_lists: params
            .highlight_bond_lists_matrix
            .as_ref()
            .map(|m| pad(m, mols_per_row, Vec::new())),
        draw_options: params.draw_options.clone(),
    };
    Ok((mols, grid))
}

/// Parameters for `reaction_to_image_with_params`, with Python's defaults
/// (A2). The last three are `MolDraw2D.DrawReaction`'s, which Python's
/// `**kwargs` pass on.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionToImageParams {
    /// Each template's panel width, and the drawing's height, in pixels.
    pub sub_img_size: (i32, i32),
    /// The drawing options; `None` for the defaults.
    pub draw_options: Option<MolDrawOptions>,
    /// Highlights each reactant's mapped atoms, and the product atoms they
    /// map to, in a color of its own.
    pub highlight_by_reactant: bool,
    /// The reactants' colors, in turn; `None` for the options'
    /// `highlight_colour_palette`.
    pub highlight_colors_reactants: Option<Vec<DrawColour>>,
    /// The conformer of each reactant, product, and agent, by its position
    /// in its group.
    pub conf_ids: Option<Vec<i32>>,
}

impl Default for ReactionToImageParams {
    fn default() -> ReactionToImageParams {
        ReactionToImageParams {
            sub_img_size: (200, 200),
            draw_options: None,
            highlight_by_reactant: false,
            highlight_colors_reactants: None,
            conf_ids: None,
        }
    }
}

/// A PNG drawing of the reaction. Mirrors `Draw.ReactionToImage`, which
/// returns a PIL image of it (or its PNG data, with `returnPNG=True`).
pub fn reaction_to_image(rxn: &ChemicalReaction) -> Result<Vec<u8>> {
    reaction_to_image_with_params(rxn, &ReactionToImageParams::default())
}

/// Mirrors `Draw.ReactionToImage(rxn, subImgSize, drawOptions=...,
/// highlightByReactant=..., highlightColorsReactants=..., confIds=...)`,
/// with the errors of `MolDraw2D::draw_reaction_with_params`.
pub fn reaction_to_image_with_params(
    rxn: &ChemicalReaction,
    params: &ReactionToImageParams,
) -> Result<Vec<u8>> {
    let mut drawer = MolDraw2DCairo::new(reaction_width(rxn, params)?, params.sub_img_size.1)?;
    draw_reaction(&mut drawer, rxn, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// An SVG drawing of the reaction. Mirrors
/// `Draw.ReactionToImage(rxn, useSVG=True)`.
pub fn reaction_to_image_use_svg(rxn: &ChemicalReaction) -> Result<String> {
    reaction_to_image_use_svg_with_params(rxn, &ReactionToImageParams::default())
}

/// Mirrors `Draw.ReactionToImage(rxn, ..., useSVG=True)`, with the errors
/// of `reaction_to_image_with_params`.
pub fn reaction_to_image_use_svg_with_params(
    rxn: &ChemicalReaction,
    params: &ReactionToImageParams,
) -> Result<String> {
    let mut drawer = MolDraw2DSVG::new(reaction_width(rxn, params)?, params.sub_img_size.1)?;
    draw_reaction(&mut drawer, rxn, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// A panel for each template, and one for the arrow.
fn reaction_width(rxn: &ChemicalReaction, params: &ReactionToImageParams) -> Result<i32> {
    let panels = u64::from(rxn.get_num_reactant_templates())
        + u64::from(rxn.get_num_product_templates())
        + 1;
    i32::try_from(panels)
        .ok()
        .and_then(|panels| panels.checked_mul(params.sub_img_size.0))
        .ok_or_else(|| Error::Value {
            message: "the reaction is too large to draw".to_owned(),
        })
}

fn draw_reaction(
    drawer: &mut MolDraw2D,
    rxn: &ChemicalReaction,
    params: &ReactionToImageParams,
) -> Result<()> {
    if let Some(options) = &params.draw_options {
        drawer.set_draw_options(options)?;
    }
    drawer.draw_reaction_with_params(
        rxn,
        &DrawReactionParams {
            highlight_by_reactant: params.highlight_by_reactant,
            highlight_colors_reactants: params.highlight_colors_reactants.clone(),
            conf_ids: params.conf_ids.clone(),
        },
    )
}

/// Sets the options to draw in comic mode, in RDKit's Comic Neue font.
/// Mirrors `Draw.SetComicMode`; `Error::Value` when the font's path (in
/// RDKit's data directory, `$RDBASE/Data` when `RDBASE` is set) is not
/// UTF-8.
pub fn set_comic_mode(opts: &mut MolDrawOptions) -> Result<()> {
    let font = crate::rd_paths::rd_data_dir()
        .join("Fonts")
        .join("ComicNeue-Regular.ttf");
    opts.font_file = font
        .into_os_string()
        .into_string()
        .map_err(|_| Error::Value {
            message: "the font's path is not UTF-8".to_owned(),
        })?;
    opts.comic_mode = true;
    Ok(())
}

/// Parameters for `debug_draw_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DebugDrawParams {
    /// The image's width and height, in pixels (unused with a drawer).
    pub size: (i32, i32),
    /// Draws atoms in black and white.
    pub use_bw: bool,
    /// Notes each atom's hydrogen count.
    pub include_h_labels: bool,
    /// Notes each atom's index.
    pub add_atom_indices: bool,
    /// Notes each bond's index.
    pub add_bond_indices: bool,
}

impl Default for DebugDrawParams {
    fn default() -> DebugDrawParams {
        DebugDrawParams {
            size: (350, 350),
            use_bw: true,
            include_h_labels: true,
            add_atom_indices: true,
            add_bond_indices: false,
        }
    }
}

/// An SVG drawing for debugging: atom indices and hydrogen counts noted,
/// aromatic atoms and bonds highlighted, and the molecule drawn as it is,
/// not kekulized. Mirrors `Draw.DebugDraw`.
pub fn debug_draw(mol: &Mol) -> Result<String> {
    debug_draw_with_params(mol, &DebugDrawParams::default())
}

/// Mirrors `Draw.DebugDraw(mol, size, asSVG=True, useBW, includeHLabels,
/// addAtomIndices, addBondIndices)`. Python notes the hydrogen counts on
/// the caller's molecule (its atoms' `atomNote`); this notes them on a
/// copy.
pub fn debug_draw_with_params(mol: &Mol, params: &DebugDrawParams) -> Result<String> {
    let mut drawer = MolDraw2DSVG::new(params.size.0, params.size.1)?;
    debug_draw_with_drawer(mol, &mut drawer, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// A PNG drawing for debugging. Mirrors `Draw.DebugDraw(mol, asSVG=False)`.
pub fn debug_draw_as_png(mol: &Mol) -> Result<Vec<u8>> {
    debug_draw_as_png_with_params(mol, &DebugDrawParams::default())
}

/// Mirrors `Draw.DebugDraw(mol, size, asSVG=False, ...)`, as
/// `debug_draw_with_params`.
pub fn debug_draw_as_png_with_params(mol: &Mol, params: &DebugDrawParams) -> Result<Vec<u8>> {
    let mut drawer = MolDraw2DCairo::new(params.size.0, params.size.1)?;
    debug_draw_with_drawer(mol, &mut drawer, params)?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Draws for debugging on `drawer`; `size` is unused. Mirrors
/// `Draw.DebugDraw(mol, drawer=drawer, ...)`, which also finishes the
/// drawing and returns its text, as the drawer's `finish_drawing` and
/// `get_drawing_text` do here.
pub fn debug_draw_with_drawer(
    mol: &Mol,
    drawer: &mut MolDraw2D,
    params: &DebugDrawParams,
) -> Result<()> {
    let mut options = drawer.draw_options();
    if params.use_bw {
        options.use_bw_atom_palette();
    }
    options.add_atom_indices = params.add_atom_indices;
    options.add_bond_indices = params.add_bond_indices;
    options.annotation_font_scale = 0.75;
    let mut mol = mol.clone();
    if params.include_h_labels {
        for idx in 0..mol.get_num_atoms() {
            let num_hs = mol.get_atom_with_idx(idx)?.get_total_num_hs()?;
            if num_hs != 0 {
                mol.get_atom_with_idx_mut(idx)?
                    .set_prop("atomNote", &format!("H{num_hs}"));
            }
        }
    }
    let arom_atoms: Vec<u32> = mol
        .get_atoms()
        .filter(|atom| atom.get_is_aromatic())
        .map(|atom| atom.get_idx())
        .collect();
    let clrs = arom_atoms
        .iter()
        .map(|&idx| (idx, DrawColour::new(0.9, 0.9, 0.2)))
        .collect();
    let arom_bonds: Vec<u32> = mol
        .get_bonds()
        .filter(|bond| bond.get_is_aromatic())
        .map(|bond| bond.get_idx())
        .collect();
    // Python prepares a copy, which it does not draw; the call still fails
    // where Python's would.
    prepare_mol_for_drawing_with_params(
        &mol,
        &PrepareMolForDrawingParams {
            kekulize: false,
            add_chiral_hs: false,
            ..Default::default()
        },
    )?;
    options.prepare_mols_before_drawing = false;
    drawer.set_draw_options(&options)?;
    drawer.draw_molecule_with_params(
        &mol,
        &DrawMoleculeParams {
            highlight_atoms: arom_atoms,
            highlight_bonds: Some(arom_bonds),
            highlight_atom_colors: clrs,
            ..Default::default()
        },
    )
}

/// A fingerprint bit's environment, drawn as a molecule of its own with
/// highlights. Mirrors `Draw.FingerprintEnv`, the named tuple the bit
/// drawings build.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct FingerprintEnv {
    /// The environment, with the molecule's coordinates. Mirrors
    /// `FingerprintEnv.submol`.
    pub submol: Mol,
    /// The atoms to highlight. Mirrors `FingerprintEnv.highlightAtoms`.
    pub highlight_atoms: Vec<u32>,
    /// The highlighted atoms' colors. Mirrors `FingerprintEnv.atomColors`.
    pub atom_colors: BTreeMap<u32, DrawColour>,
    /// The bonds to highlight. Mirrors `FingerprintEnv.highlightBonds`.
    pub highlight_bonds: Vec<u32>,
    /// The highlighted bonds' colors. Mirrors `FingerprintEnv.bondColors`.
    pub bond_colors: BTreeMap<u32, DrawColour>,
    /// The highlighted atoms' radii. Mirrors `FingerprintEnv.highlightRadii`.
    pub highlight_radii: BTreeMap<u32, f64>,
}

/// The molecule, with 2D coordinates when it has no conformer. Python
/// computes them on the caller's molecule; this computes them on a copy.
fn with_coords(mol: &Mol) -> Result<Mol> {
    let mut mol = mol.clone();
    if mol.get_num_conformers() == 0 {
        compute_2d_coords(&mut mol)?;
    }
    Ok(mol)
}

/// The environment's submolecule (`Chem.PathToSubmol`) with the
/// molecule's coordinates, its map from the molecule's atoms, and the
/// submolecule's bonds that are in `bit_path`. A Morgan environment of no
/// bonds is the SMILES of its atoms, `atoms_to_use`.
fn env_submol(
    mol: &Mol,
    path: &[u32],
    atoms_to_use: Option<&[u32]>,
    bit_path: &[u32],
) -> Result<(Mol, BTreeMap<u32, u32>, Vec<u32>)> {
    let (mut submol, amap) = if let Some(atoms) = atoms_to_use.filter(|_| path.is_empty()) {
        (
            mol_from_smiles(&mol_fragment_to_smiles(mol, atoms)?)?,
            BTreeMap::new(),
        )
    } else {
        let (submol, outputs) = path_to_submol_with_outputs(mol, path, false)?;
        (submol, outputs.atom_map)
    };
    fast_find_rings(&mut submol)?;
    let mut conf = Conformer::new(submol.get_num_atoms());
    {
        let conf_ori = mol.get_conformer_with_id(0)?;
        for (&i1, &i2) in &amap {
            conf.set_atom_position(i2, conf_ori.get_atom_position(i1)?)?;
        }
    }
    submol.add_conformer(&conf)?;
    let mut env_submol = Vec::new();
    for &b in bit_path {
        let bond = mol.get_bond_with_idx(b)?;
        let begin = amap[&bond.get_begin_atom_idx()];
        let end = amap[&bond.get_end_atom_idx()];
        let sub_bond = submol
            .get_bond_between_atoms(begin, end)?
            .expect("the submolecule has the path's bonds");
        env_submol.push(sub_bond.get_idx());
    }
    Ok((submol, amap, env_submol))
}

/// Colors the bonds outside the bit's path.
fn extra_bonds(submol: &Mol, env_submol: &[u32], color: DrawColour, env: &mut FingerprintEnv) {
    for bond in submol.get_bonds() {
        let bidx = bond.get_idx();
        if !env_submol.contains(&bidx) {
            env.bond_colors.insert(bidx, color);
            env.highlight_bonds.push(bidx);
        }
    }
}

fn get_morgan_env(
    mol: &Mol,
    atom_id: u32,
    radius: u32,
    params: &DrawMorganEnvParams,
) -> Result<FingerprintEnv> {
    let mol = with_coords(mol)?;
    let bit_path = find_atom_environment_of_radius_n(&mol, radius, atom_id)?;
    // The atoms to highlight.
    let mut atoms_to_use = std::collections::BTreeSet::from([atom_id]);
    for &b in &bit_path {
        let bond = mol.get_bond_with_idx(b)?;
        atoms_to_use.insert(bond.get_begin_atom_idx());
        atoms_to_use.insert(bond.get_end_atom_idx());
    }
    // The environment, enlarged by one more bond. (The submolecule's atoms
    // keep the molecule's order, so the order of the bonds does not matter.)
    let mut enlarged_env = std::collections::BTreeSet::new();
    for &atom in &atoms_to_use {
        for bond in mol.get_atom_with_idx(atom)?.get_bonds() {
            if !bit_path.contains(&bond.get_idx()) {
                enlarged_env.insert(bond.get_idx());
            }
        }
    }
    let mut enlarged_env: Vec<u32> = enlarged_env.into_iter().collect();
    enlarged_env.extend_from_slice(&bit_path);
    let atoms: Vec<u32> = atoms_to_use.iter().copied().collect();
    let (submol, amap, env_submol) = env_submol(&mol, &enlarged_env, Some(&atoms), &bit_path)?;
    let mut env = FingerprintEnv {
        submol,
        highlight_atoms: Vec::new(),
        atom_colors: BTreeMap::new(),
        highlight_bonds: Vec::new(),
        bond_colors: BTreeMap::new(),
        highlight_radii: BTreeMap::new(),
    };
    for (&aidx, &sub_idx) in &amap {
        if atoms_to_use.contains(&aidx) {
            let atom = mol.get_atom_with_idx(aidx)?;
            let color = match () {
                _ if params.center_color.is_some() && aidx == atom_id => params.center_color,
                _ if params.aromatic_color.is_some() && atom.get_is_aromatic() => {
                    params.aromatic_color
                }
                _ if params.ring_color.is_some() && atom.is_in_ring() => params.ring_color,
                _ => None,
            };
            if let Some(color) = color {
                env.atom_colors.insert(sub_idx, color);
                env.highlight_atoms.push(sub_idx);
                env.highlight_radii.insert(sub_idx, params.base_rad);
            }
        } else {
            let mut atom = env.submol.get_atom_with_idx_mut(sub_idx)?;
            atom.set_atomic_num(0);
            atom.update_property_cache()?;
        }
    }
    let submol = env.submol.clone();
    extra_bonds(&submol, &env_submol, params.extra_color, &mut env);
    Ok(env)
}

fn get_rdkit_env(
    mol: &Mol,
    bond_path: &[u32],
    params: &DrawRDKitEnvParams,
) -> Result<FingerprintEnv> {
    let mol = with_coords(mol)?;
    let mut atoms_to_use = std::collections::BTreeSet::new();
    for &b in bond_path {
        let bond = mol.get_bond_with_idx(b)?;
        atoms_to_use.insert(bond.get_begin_atom_idx());
        atoms_to_use.insert(bond.get_end_atom_idx());
    }
    let (submol, amap, env_submol) = env_submol(&mol, bond_path, None, bond_path)?;
    let mut env = FingerprintEnv {
        submol,
        highlight_atoms: Vec::new(),
        atom_colors: BTreeMap::new(),
        highlight_bonds: Vec::new(),
        bond_colors: BTreeMap::new(),
        highlight_radii: BTreeMap::new(),
    };
    for (&aidx, &sub_idx) in &amap {
        if atoms_to_use.contains(&aidx) {
            let aromatic = mol.get_atom_with_idx(aidx)?.get_is_aromatic();
            let color = match () {
                _ if params.aromatic_color.is_some() && aromatic => params.aromatic_color,
                _ if params.non_aromatic_color.is_some() && !aromatic => params.non_aromatic_color,
                _ => None,
            };
            if let Some(color) = color {
                env.atom_colors.insert(sub_idx, color);
                env.highlight_atoms.push(sub_idx);
                env.highlight_radii.insert(sub_idx, params.base_rad);
            }
        }
    }
    let submol = env.submol.clone();
    extra_bonds(&submol, &env_submol, params.extra_color, &mut env);
    Ok(env)
}

/// Sets the options for drawing environments: Python's, else the
/// drawer's, without continuous highlights or metadata.
fn set_env_options(drawer: &mut MolDraw2D, draw_options: Option<&MolDrawOptions>) -> Result<()> {
    let mut options = match draw_options {
        Some(options) => options.clone(),
        None => drawer.draw_options(),
    };
    options.continuous_highlight = false;
    options.include_metadata = false;
    drawer.set_draw_options(&options)
}

fn draw_env(
    drawer: &mut MolDraw2D,
    env: &FingerprintEnv,
    draw_options: Option<&MolDrawOptions>,
) -> Result<()> {
    set_env_options(drawer, draw_options)?;
    drawer.draw_molecule_with_params(
        &env.submol,
        &DrawMoleculeParams {
            highlight_atoms: env.highlight_atoms.clone(),
            highlight_bonds: Some(env.highlight_bonds.clone()),
            highlight_atom_colors: env.atom_colors.clone(),
            highlight_bond_colors: env.bond_colors.clone(),
            highlight_atom_radii: env.highlight_radii.clone(),
            ..Default::default()
        },
    )
}

fn draw_envs(
    drawer: &mut MolDraw2D,
    envs: &[FingerprintEnv],
    legends: &Option<Vec<String>>,
    draw_options: Option<&MolDrawOptions>,
) -> Result<()> {
    set_env_options(drawer, draw_options)?;
    let mols: Vec<Option<&Mol>> = envs.iter().map(|env| Some(&env.submol)).collect();
    drawer.draw_optional_molecules_with_params(
        &mols,
        &DrawMoleculesParams {
            legends: Some(
                legends
                    .clone()
                    .unwrap_or_else(|| vec![String::new(); envs.len()]),
            ),
            highlight_atoms: Some(envs.iter().map(|e| e.highlight_atoms.clone()).collect()),
            highlight_atom_colors: Some(envs.iter().map(|e| e.atom_colors.clone()).collect()),
            highlight_bonds: Some(envs.iter().map(|e| e.highlight_bonds.clone()).collect()),
            highlight_bond_colors: Some(envs.iter().map(|e| e.bond_colors.clone()).collect()),
            highlight_atom_radii: Some(envs.iter().map(|e| e.highlight_radii.clone()).collect()),
            ..Default::default()
        },
    )
}

/// A grid's canvas and panel sizes, for `count` panels.
fn envs_grid(
    count: usize,
    mols_per_row: u32,
    sub_img_size: (i32, i32),
) -> Result<MolsToGridImageParams> {
    let params = MolsToGridImageParams {
        mols_per_row,
        sub_img_size,
        ..Default::default()
    };
    grid_size(count, &params)?;
    Ok(params)
}

fn new_grid_svg(count: usize, mols_per_row: u32, sub_img_size: (i32, i32)) -> Result<MolDraw2DSVG> {
    let (width, height, panel_width, panel_height) =
        grid_size(count, &envs_grid(count, mols_per_row, sub_img_size)?)?;
    MolDraw2DSVG::new_with_params(
        width,
        height,
        &MolDraw2DSVGParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )
}

fn new_grid_cairo(
    count: usize,
    mols_per_row: u32,
    sub_img_size: (i32, i32),
) -> Result<MolDraw2DCairo> {
    let (width, height, panel_width, panel_height) =
        grid_size(count, &envs_grid(count, mols_per_row, sub_img_size)?)?;
    MolDraw2DCairo::new_with_params(
        width,
        height,
        &MolDraw2DCairoParams {
            panel_width,
            panel_height,
            ..Default::default()
        },
    )
}

/// Parameters for `draw_morgan_env_with_params`, with Python's defaults
/// (A2). A color of `None` leaves those atoms unhighlighted.
#[derive(Clone, Debug, PartialEq)]
pub struct DrawMorganEnvParams {
    /// The image's width and height, in pixels.
    pub mol_size: (i32, i32),
    /// The highlight radius of the environment's atoms.
    pub base_rad: f64,
    /// The color of its aromatic atoms.
    pub aromatic_color: Option<DrawColour>,
    /// The color of its other atoms in rings.
    pub ring_color: Option<DrawColour>,
    /// The color of its central atom.
    pub center_color: Option<DrawColour>,
    /// The color of the bonds around it.
    pub extra_color: DrawColour,
    /// The drawing options; `None` for the drawer's.
    pub draw_options: Option<MolDrawOptions>,
}

impl Default for DrawMorganEnvParams {
    fn default() -> DrawMorganEnvParams {
        DrawMorganEnvParams {
            mol_size: (150, 150),
            base_rad: 0.3,
            aromatic_color: Some(DrawColour::new(0.9, 0.9, 0.2)),
            ring_color: Some(DrawColour::new(0.8, 0.8, 0.8)),
            center_color: Some(DrawColour::new(0.6, 0.6, 0.9)),
            extra_color: DrawColour::new(0.9, 0.9, 0.9),
            draw_options: None,
        }
    }
}

/// An SVG drawing of the Morgan environment of `radius` around atom
/// `atom_id`. Mirrors `Draw.DrawMorganEnv`.
pub fn draw_morgan_env(mol: &Mol, atom_id: u32, radius: u32) -> Result<String> {
    draw_morgan_env_with_params(mol, atom_id, radius, &DrawMorganEnvParams::default())
}

/// Mirrors `Draw.DrawMorganEnv(mol, atomId, radius, molSize, baseRad,
/// useSVG=True, aromaticColor, ringColor, centerColor, extraColor,
/// drawOptions)`. Python gives a molecule without conformers 2D
/// coordinates; this gives them to a copy.
pub fn draw_morgan_env_with_params(
    mol: &Mol,
    atom_id: u32,
    radius: u32,
    params: &DrawMorganEnvParams,
) -> Result<String> {
    let env = get_morgan_env(mol, atom_id, radius, params)?;
    let mut drawer = MolDraw2DSVG::new(params.mol_size.0, params.mol_size.1)?;
    draw_env(&mut drawer, &env, params.draw_options.as_ref())?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// A PNG drawing of the Morgan environment. Mirrors
/// `Draw.DrawMorganEnv(mol, atomId, radius, useSVG=False)`.
pub fn draw_morgan_env_as_png(mol: &Mol, atom_id: u32, radius: u32) -> Result<Vec<u8>> {
    draw_morgan_env_as_png_with_params(mol, atom_id, radius, &DrawMorganEnvParams::default())
}

/// Mirrors `Draw.DrawMorganEnv(mol, atomId, radius, ..., useSVG=False)`.
pub fn draw_morgan_env_as_png_with_params(
    mol: &Mol,
    atom_id: u32,
    radius: u32,
    params: &DrawMorganEnvParams,
) -> Result<Vec<u8>> {
    let env = get_morgan_env(mol, atom_id, radius, params)?;
    let mut drawer = MolDraw2DCairo::new(params.mol_size.0, params.mol_size.1)?;
    draw_env(&mut drawer, &env, params.draw_options.as_ref())?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Parameters for `draw_morgan_envs_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DrawMorganEnvsParams {
    /// The environments in each row.
    pub mols_per_row: u32,
    /// Each environment's panel width and height, in pixels.
    pub sub_img_size: (i32, i32),
    /// The highlight radius of the environments' atoms.
    pub base_rad: f64,
    /// The color of their aromatic atoms.
    pub aromatic_color: Option<DrawColour>,
    /// The color of their other atoms in rings.
    pub ring_color: Option<DrawColour>,
    /// The color of their central atoms.
    pub center_color: Option<DrawColour>,
    /// The color of the bonds around them.
    pub extra_color: DrawColour,
    /// Each environment's legend; `None` for none.
    pub legends: Option<Vec<String>>,
    /// The drawing options; `None` for the drawer's.
    pub draw_options: Option<MolDrawOptions>,
}

impl Default for DrawMorganEnvsParams {
    fn default() -> DrawMorganEnvsParams {
        let env = DrawMorganEnvParams::default();
        DrawMorganEnvsParams {
            mols_per_row: 3,
            sub_img_size: (150, 150),
            base_rad: env.base_rad,
            aromatic_color: env.aromatic_color,
            ring_color: env.ring_color,
            center_color: env.center_color,
            extra_color: env.extra_color,
            legends: None,
            draw_options: None,
        }
    }
}

impl DrawMorganEnvsParams {
    fn env(&self) -> DrawMorganEnvParams {
        DrawMorganEnvParams {
            mol_size: self.sub_img_size,
            base_rad: self.base_rad,
            aromatic_color: self.aromatic_color,
            ring_color: self.ring_color,
            center_color: self.center_color,
            extra_color: self.extra_color,
            draw_options: None,
        }
    }
}

fn morgan_envs(
    envs: &[(&Mol, u32, u32)],
    params: &DrawMorganEnvsParams,
) -> Result<Vec<FingerprintEnv>> {
    let env_params = params.env();
    envs.iter()
        .map(|&(mol, atom_id, radius)| get_morgan_env(mol, atom_id, radius, &env_params))
        .collect()
}

/// An SVG drawing of Morgan environments, each a `(mol, atom_id, radius)`,
/// in a grid. Mirrors `Draw.DrawMorganEnvs`.
pub fn draw_morgan_envs(envs: &[(&Mol, u32, u32)]) -> Result<String> {
    draw_morgan_envs_with_params(envs, &DrawMorganEnvsParams::default())
}

/// Mirrors `Draw.DrawMorganEnvs(envs, molsPerRow, subImgSize, baseRad,
/// useSVG=True, aromaticColor, ringColor, centerColor, extraColor, legends,
/// drawOptions)`, with the errors of
/// `MolDraw2D::draw_molecules_with_params`, and `Error::Value` for no
/// molecules per row (Python divides by zero).
pub fn draw_morgan_envs_with_params(
    envs: &[(&Mol, u32, u32)],
    params: &DrawMorganEnvsParams,
) -> Result<String> {
    let fp_envs = morgan_envs(envs, params)?;
    let mut drawer = new_grid_svg(envs.len(), params.mols_per_row, params.sub_img_size)?;
    draw_envs(
        &mut drawer,
        &fp_envs,
        &params.legends,
        params.draw_options.as_ref(),
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// A PNG drawing of Morgan environments in a grid. Mirrors
/// `Draw.DrawMorganEnvs(envs, useSVG=False)`.
pub fn draw_morgan_envs_as_png(envs: &[(&Mol, u32, u32)]) -> Result<Vec<u8>> {
    draw_morgan_envs_as_png_with_params(envs, &DrawMorganEnvsParams::default())
}

/// Mirrors `Draw.DrawMorganEnvs(envs, ..., useSVG=False)`, with the errors
/// of `draw_morgan_envs_with_params`.
pub fn draw_morgan_envs_as_png_with_params(
    envs: &[(&Mol, u32, u32)],
    params: &DrawMorganEnvsParams,
) -> Result<Vec<u8>> {
    let fp_envs = morgan_envs(envs, params)?;
    let mut drawer = new_grid_cairo(envs.len(), params.mols_per_row, params.sub_img_size)?;
    draw_envs(
        &mut drawer,
        &fp_envs,
        &params.legends,
        params.draw_options.as_ref(),
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Morgan bit information, as the fingerprint generators'
/// `AdditionalOutput::get_bit_info_map` gives it: each bit's `(atom,
/// radius)` environments.
/// A Rust addition, naming Python's `bitInfo` dictionary.
pub type MorganBitInfo = BTreeMap<u64, Vec<(u32, u32)>>;

/// RDKit-fingerprint bit information, as `rdk_fingerprint_with_outputs`
/// gives it (`bit_info`): each bit's paths, as bond indices.
/// A Rust addition, naming Python's `bitInfo` dictionary.
pub type RDKitBitInfo = BTreeMap<u64, Vec<Vec<u32>>>;

/// `bit_info[bit_id][which_example]`, with Python's `KeyError` and
/// `IndexError`.
fn bit_example<T>(
    bit_info: &BTreeMap<u64, Vec<T>>,
    bit_id: u64,
    which_example: usize,
) -> Result<&T> {
    let examples = bit_info.get(&bit_id).ok_or_else(|| Error::Key {
        key: bit_id.to_string(),
    })?;
    examples
        .get(which_example)
        .ok_or(Error::index(which_example, examples.len()))
}

/// Parameters for `draw_morgan_bit_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DrawMorganBitParams {
    /// Which of the bit's environments to draw.
    pub which_example: usize,
    /// How to draw it, as `draw_morgan_env_with_params`.
    pub env: DrawMorganEnvParams,
}

/// An SVG drawing of the first environment that sets a Morgan bit.
/// `Error::Key` for a bit not in `bit_info`. Mirrors `Draw.DrawMorganBit`.
pub fn draw_morgan_bit(mol: &Mol, bit_id: u64, bit_info: &MorganBitInfo) -> Result<String> {
    draw_morgan_bit_with_params(mol, bit_id, bit_info, &DrawMorganBitParams::default())
}

/// Mirrors `Draw.DrawMorganBit(mol, bitId, bitInfo, whichExample,
/// **kwargs)`; `Error::Index` for an example past the bit's last.
pub fn draw_morgan_bit_with_params(
    mol: &Mol,
    bit_id: u64,
    bit_info: &MorganBitInfo,
    params: &DrawMorganBitParams,
) -> Result<String> {
    let &(atom_id, radius) = bit_example(bit_info, bit_id, params.which_example)?;
    draw_morgan_env_with_params(mol, atom_id, radius, &params.env)
}

/// A PNG drawing of a Morgan bit. Mirrors
/// `Draw.DrawMorganBit(mol, bitId, bitInfo, useSVG=False)`.
pub fn draw_morgan_bit_as_png(mol: &Mol, bit_id: u64, bit_info: &MorganBitInfo) -> Result<Vec<u8>> {
    draw_morgan_bit_as_png_with_params(mol, bit_id, bit_info, &DrawMorganBitParams::default())
}

/// Mirrors `Draw.DrawMorganBit(mol, bitId, bitInfo, ..., useSVG=False)`.
pub fn draw_morgan_bit_as_png_with_params(
    mol: &Mol,
    bit_id: u64,
    bit_info: &MorganBitInfo,
    params: &DrawMorganBitParams,
) -> Result<Vec<u8>> {
    let &(atom_id, radius) = bit_example(bit_info, bit_id, params.which_example)?;
    draw_morgan_env_as_png_with_params(mol, atom_id, radius, &params.env)
}

fn morgan_bit_envs<'a>(
    tpls: &[(&'a Mol, u64, &MorganBitInfo, usize)],
) -> Result<Vec<(&'a Mol, u32, u32)>> {
    tpls.iter()
        .map(|&(mol, bit_id, bit_info, which_example)| {
            let &(atom_id, radius) = bit_example(bit_info, bit_id, which_example)?;
            Ok((mol, atom_id, radius))
        })
        .collect()
}

/// An SVG drawing of Morgan bits, each a `(mol, bit_id, bit_info,
/// which_example)`, in a grid. Mirrors `Draw.DrawMorganBits` (whose tuples
/// may leave out `whichExample`, 0).
pub fn draw_morgan_bits(tpls: &[(&Mol, u64, &MorganBitInfo, usize)]) -> Result<String> {
    draw_morgan_bits_with_params(tpls, &DrawMorganEnvsParams::default())
}

/// Mirrors `Draw.DrawMorganBits(tpls, **kwargs)`, which passes the
/// keywords to `DrawMorganEnvs`.
pub fn draw_morgan_bits_with_params(
    tpls: &[(&Mol, u64, &MorganBitInfo, usize)],
    params: &DrawMorganEnvsParams,
) -> Result<String> {
    draw_morgan_envs_with_params(&morgan_bit_envs(tpls)?, params)
}

/// A PNG drawing of Morgan bits in a grid. Mirrors
/// `Draw.DrawMorganBits(tpls, useSVG=False)`.
pub fn draw_morgan_bits_as_png(tpls: &[(&Mol, u64, &MorganBitInfo, usize)]) -> Result<Vec<u8>> {
    draw_morgan_bits_as_png_with_params(tpls, &DrawMorganEnvsParams::default())
}

/// Mirrors `Draw.DrawMorganBits(tpls, ..., useSVG=False)`.
pub fn draw_morgan_bits_as_png_with_params(
    tpls: &[(&Mol, u64, &MorganBitInfo, usize)],
    params: &DrawMorganEnvsParams,
) -> Result<Vec<u8>> {
    draw_morgan_envs_as_png_with_params(&morgan_bit_envs(tpls)?, params)
}

/// Parameters for `draw_rdkit_env_with_params`, with Python's defaults
/// (A2). A color of `None` leaves those atoms unhighlighted.
#[derive(Clone, Debug, PartialEq)]
pub struct DrawRDKitEnvParams {
    /// The image's width and height, in pixels.
    pub mol_size: (i32, i32),
    /// The highlight radius of the path's atoms.
    pub base_rad: f64,
    /// The color of its aromatic atoms.
    pub aromatic_color: Option<DrawColour>,
    /// The color of the bonds around it.
    pub extra_color: DrawColour,
    /// The color of its other atoms.
    pub non_aromatic_color: Option<DrawColour>,
    /// The drawing options; `None` for the drawer's.
    pub draw_options: Option<MolDrawOptions>,
}

impl Default for DrawRDKitEnvParams {
    fn default() -> DrawRDKitEnvParams {
        DrawRDKitEnvParams {
            mol_size: (150, 150),
            base_rad: 0.3,
            aromatic_color: Some(DrawColour::new(0.9, 0.9, 0.2)),
            extra_color: DrawColour::new(0.9, 0.9, 0.9),
            non_aromatic_color: None,
            draw_options: None,
        }
    }
}

/// An SVG drawing of the path of bonds `bond_path`, one that sets a bit of
/// the RDKit fingerprint. Mirrors `Draw.DrawRDKitEnv`.
pub fn draw_rdkit_env(mol: &Mol, bond_path: &[u32]) -> Result<String> {
    draw_rdkit_env_with_params(mol, bond_path, &DrawRDKitEnvParams::default())
}

/// Mirrors `Draw.DrawRDKitEnv(mol, bondPath, molSize, baseRad, useSVG=True,
/// aromaticColor, extraColor, nonAromaticColor, drawOptions)`. Python gives
/// a molecule without conformers 2D coordinates; this gives them to a copy.
pub fn draw_rdkit_env_with_params(
    mol: &Mol,
    bond_path: &[u32],
    params: &DrawRDKitEnvParams,
) -> Result<String> {
    let env = get_rdkit_env(mol, bond_path, params)?;
    let mut drawer = MolDraw2DSVG::new(params.mol_size.0, params.mol_size.1)?;
    draw_env(&mut drawer, &env, params.draw_options.as_ref())?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// A PNG drawing of the path. Mirrors
/// `Draw.DrawRDKitEnv(mol, bondPath, useSVG=False)`.
pub fn draw_rdkit_env_as_png(mol: &Mol, bond_path: &[u32]) -> Result<Vec<u8>> {
    draw_rdkit_env_as_png_with_params(mol, bond_path, &DrawRDKitEnvParams::default())
}

/// Mirrors `Draw.DrawRDKitEnv(mol, bondPath, ..., useSVG=False)`.
pub fn draw_rdkit_env_as_png_with_params(
    mol: &Mol,
    bond_path: &[u32],
    params: &DrawRDKitEnvParams,
) -> Result<Vec<u8>> {
    let env = get_rdkit_env(mol, bond_path, params)?;
    let mut drawer = MolDraw2DCairo::new(params.mol_size.0, params.mol_size.1)?;
    draw_env(&mut drawer, &env, params.draw_options.as_ref())?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Parameters for `draw_rdkit_envs_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DrawRDKitEnvsParams {
    /// The paths in each row.
    pub mols_per_row: u32,
    /// Each path's panel width and height, in pixels.
    pub sub_img_size: (i32, i32),
    /// The highlight radius of the paths' atoms.
    pub base_rad: f64,
    /// The color of their aromatic atoms.
    pub aromatic_color: Option<DrawColour>,
    /// The color of the bonds around them.
    pub extra_color: DrawColour,
    /// The color of their other atoms.
    pub non_aromatic_color: Option<DrawColour>,
    /// Each path's legend; `None` for none.
    pub legends: Option<Vec<String>>,
    /// The drawing options; `None` for the drawer's.
    pub draw_options: Option<MolDrawOptions>,
}

impl Default for DrawRDKitEnvsParams {
    fn default() -> DrawRDKitEnvsParams {
        let env = DrawRDKitEnvParams::default();
        DrawRDKitEnvsParams {
            mols_per_row: 3,
            sub_img_size: (150, 150),
            base_rad: env.base_rad,
            aromatic_color: env.aromatic_color,
            extra_color: env.extra_color,
            non_aromatic_color: env.non_aromatic_color,
            legends: None,
            draw_options: None,
        }
    }
}

fn rdkit_envs(
    envs: &[(&Mol, &[u32])],
    params: &DrawRDKitEnvsParams,
) -> Result<Vec<FingerprintEnv>> {
    let env_params = DrawRDKitEnvParams {
        mol_size: params.sub_img_size,
        base_rad: params.base_rad,
        aromatic_color: params.aromatic_color,
        extra_color: params.extra_color,
        non_aromatic_color: params.non_aromatic_color,
        draw_options: None,
    };
    envs.iter()
        .map(|&(mol, bond_path)| get_rdkit_env(mol, bond_path, &env_params))
        .collect()
}

/// An SVG drawing of paths, each a `(mol, bond_path)`, in a grid. Mirrors
/// `Draw.DrawRDKitEnvs`.
pub fn draw_rdkit_envs(envs: &[(&Mol, &[u32])]) -> Result<String> {
    draw_rdkit_envs_with_params(envs, &DrawRDKitEnvsParams::default())
}

/// Mirrors `Draw.DrawRDKitEnvs(envs, molsPerRow, subImgSize, baseRad,
/// useSVG=True, aromaticColor, extraColor, nonAromaticColor, legends,
/// drawOptions)`, with the errors of `draw_morgan_envs_with_params`.
pub fn draw_rdkit_envs_with_params(
    envs: &[(&Mol, &[u32])],
    params: &DrawRDKitEnvsParams,
) -> Result<String> {
    let fp_envs = rdkit_envs(envs, params)?;
    let mut drawer = new_grid_svg(envs.len(), params.mols_per_row, params.sub_img_size)?;
    draw_envs(
        &mut drawer,
        &fp_envs,
        &params.legends,
        params.draw_options.as_ref(),
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// A PNG drawing of paths in a grid. Mirrors
/// `Draw.DrawRDKitEnvs(envs, useSVG=False)`.
pub fn draw_rdkit_envs_as_png(envs: &[(&Mol, &[u32])]) -> Result<Vec<u8>> {
    draw_rdkit_envs_as_png_with_params(envs, &DrawRDKitEnvsParams::default())
}

/// Mirrors `Draw.DrawRDKitEnvs(envs, ..., useSVG=False)`.
pub fn draw_rdkit_envs_as_png_with_params(
    envs: &[(&Mol, &[u32])],
    params: &DrawRDKitEnvsParams,
) -> Result<Vec<u8>> {
    let fp_envs = rdkit_envs(envs, params)?;
    let mut drawer = new_grid_cairo(envs.len(), params.mols_per_row, params.sub_img_size)?;
    draw_envs(
        &mut drawer,
        &fp_envs,
        &params.legends,
        params.draw_options.as_ref(),
    )?;
    drawer.finish_drawing();
    drawer.get_drawing_text()
}

/// Parameters for `draw_rdkit_bit_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DrawRDKitBitParams {
    /// Which of the bit's paths to draw.
    pub which_example: usize,
    /// How to draw it, as `draw_rdkit_env_with_params`.
    pub env: DrawRDKitEnvParams,
}

/// An SVG drawing of the first path that sets a bit of the RDKit
/// fingerprint. `Error::Key` for a bit not in `bit_info`. Mirrors
/// `Draw.DrawRDKitBit`.
pub fn draw_rdkit_bit(mol: &Mol, bit_id: u64, bit_info: &RDKitBitInfo) -> Result<String> {
    draw_rdkit_bit_with_params(mol, bit_id, bit_info, &DrawRDKitBitParams::default())
}

/// Mirrors `Draw.DrawRDKitBit(mol, bitId, bitInfo, whichExample,
/// **kwargs)`; `Error::Index` for an example past the bit's last.
pub fn draw_rdkit_bit_with_params(
    mol: &Mol,
    bit_id: u64,
    bit_info: &RDKitBitInfo,
    params: &DrawRDKitBitParams,
) -> Result<String> {
    let path = bit_example(bit_info, bit_id, params.which_example)?;
    draw_rdkit_env_with_params(mol, path, &params.env)
}

/// A PNG drawing of a bit of the RDKit fingerprint. Mirrors
/// `Draw.DrawRDKitBit(mol, bitId, bitInfo, useSVG=False)`.
pub fn draw_rdkit_bit_as_png(mol: &Mol, bit_id: u64, bit_info: &RDKitBitInfo) -> Result<Vec<u8>> {
    draw_rdkit_bit_as_png_with_params(mol, bit_id, bit_info, &DrawRDKitBitParams::default())
}

/// Mirrors `Draw.DrawRDKitBit(mol, bitId, bitInfo, ..., useSVG=False)`.
pub fn draw_rdkit_bit_as_png_with_params(
    mol: &Mol,
    bit_id: u64,
    bit_info: &RDKitBitInfo,
    params: &DrawRDKitBitParams,
) -> Result<Vec<u8>> {
    let path = bit_example(bit_info, bit_id, params.which_example)?;
    draw_rdkit_env_as_png_with_params(mol, path, &params.env)
}

fn rdkit_bit_envs<'a>(
    tpls: &[(&'a Mol, u64, &'a RDKitBitInfo, usize)],
) -> Result<Vec<(&'a Mol, &'a [u32])>> {
    tpls.iter()
        .map(|&(mol, bit_id, bit_info, which_example)| {
            Ok((
                mol,
                bit_example(bit_info, bit_id, which_example)?.as_slice(),
            ))
        })
        .collect()
}

/// An SVG drawing of bits of the RDKit fingerprint, each a `(mol, bit_id,
/// bit_info, which_example)`, in a grid. Mirrors `Draw.DrawRDKitBits`
/// (whose tuples may leave out `whichExample`, 0).
pub fn draw_rdkit_bits(tpls: &[(&Mol, u64, &RDKitBitInfo, usize)]) -> Result<String> {
    draw_rdkit_bits_with_params(tpls, &DrawRDKitEnvsParams::default())
}

/// Mirrors `Draw.DrawRDKitBits(tpls, **kwargs)`, which passes the keywords
/// to `DrawRDKitEnvs`.
pub fn draw_rdkit_bits_with_params(
    tpls: &[(&Mol, u64, &RDKitBitInfo, usize)],
    params: &DrawRDKitEnvsParams,
) -> Result<String> {
    draw_rdkit_envs_with_params(&rdkit_bit_envs(tpls)?, params)
}

/// A PNG drawing of bits of the RDKit fingerprint in a grid. Mirrors
/// `Draw.DrawRDKitBits(tpls, useSVG=False)`.
pub fn draw_rdkit_bits_as_png(tpls: &[(&Mol, u64, &RDKitBitInfo, usize)]) -> Result<Vec<u8>> {
    draw_rdkit_bits_as_png_with_params(tpls, &DrawRDKitEnvsParams::default())
}

/// Mirrors `Draw.DrawRDKitBits(tpls, ..., useSVG=False)`.
pub fn draw_rdkit_bits_as_png_with_params(
    tpls: &[(&Mol, u64, &RDKitBitInfo, usize)],
    params: &DrawRDKitEnvsParams,
) -> Result<Vec<u8>> {
    draw_rdkit_envs_as_png_with_params(&rdkit_bit_envs(tpls)?, params)
}
