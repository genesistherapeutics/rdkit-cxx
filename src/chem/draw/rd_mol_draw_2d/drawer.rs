// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The drawers, and the functions that draw molecules with them.

use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Deref, DerefMut};

use cxx::UniquePtr;

use super::options::{DrawColour, MolDrawOptions};
use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::draw::ffi;
use crate::geometry::Point2D;

/// A molecule drawer. Mirrors `rdMolDraw2D.MolDraw2D`, the base class of
/// `MolDraw2DSVG` and `MolDraw2DCairo`, which dereference to it, so its
/// methods work on both, and functions that take a drawer take
/// `&mut MolDraw2D`.
///
/// It is `Send` and not `Sync`. A drawing call that fails partway leaves
/// the drawer unable to draw: later drawing calls return `Error::Value`,
/// since RDKit may have kept a half-built molecule (see the module
/// documentation).
pub struct MolDraw2D {
    inner: UniquePtr<ffi::Drawer>,
}

/// Parameters for `MolDraw2D::draw_molecule_with_params`, with Python's
/// defaults (A2). Python's two `DrawMolecule` overloads differ in whether
/// they take `highlightBonds`: without it (`None` here) RDKit highlights the
/// bonds between highlighted atoms.
#[derive(Clone, Debug, PartialEq)]
pub struct DrawMoleculeParams {
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight; `None` highlights those between highlighted
    /// atoms (Python's first overload).
    pub highlight_bonds: Option<Vec<u32>>,
    /// Highlight colours of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colours of bonds, by index.
    pub highlight_bond_colors: BTreeMap<u32, DrawColour>,
    /// Highlight radii of atoms, by index, in molecule coordinates.
    pub highlight_atom_radii: BTreeMap<u32, f64>,
    /// The conformer to draw; -1 for the first.
    pub conf_id: i32,
    /// Text drawn under (or by) the molecule.
    pub legend: String,
}

impl Default for DrawMoleculeParams {
    fn default() -> DrawMoleculeParams {
        DrawMoleculeParams {
            highlight_atoms: Vec::new(),
            highlight_bonds: None,
            highlight_atom_colors: BTreeMap::new(),
            highlight_bond_colors: BTreeMap::new(),
            highlight_atom_radii: BTreeMap::new(),
            conf_id: -1,
            legend: String::new(),
        }
    }
}

/// Parameters for `MolDraw2D::get_mol_size_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetMolSizeParams {
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight.
    pub highlight_bonds: Vec<u32>,
    /// Highlight colours of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colours of bonds, by index.
    pub highlight_bond_colors: BTreeMap<u32, DrawColour>,
    /// Highlight radii of atoms, by index, in molecule coordinates.
    pub highlight_atom_radii: BTreeMap<u32, f64>,
    /// The conformer to draw; -1 for the first.
    pub conf_id: i32,
    /// Text drawn under (or by) the molecule.
    pub legend: String,
}

impl Default for GetMolSizeParams {
    fn default() -> GetMolSizeParams {
        GetMolSizeParams {
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            highlight_atom_colors: BTreeMap::new(),
            highlight_bond_colors: BTreeMap::new(),
            highlight_atom_radii: BTreeMap::new(),
            conf_id: -1,
            legend: String::new(),
        }
    }
}

/// The highlights and text of one molecule, as the bridge carries them.
pub(crate) struct MoleculeArgs<'a> {
    pub(crate) highlight_atoms: &'a [u32],
    pub(crate) highlight_bonds: Option<&'a [u32]>,
    pub(crate) highlight_atom_colors: &'a BTreeMap<u32, DrawColour>,
    pub(crate) highlight_bond_colors: &'a BTreeMap<u32, DrawColour>,
    pub(crate) highlight_atom_radii: &'a BTreeMap<u32, f64>,
    pub(crate) conf_id: i32,
    pub(crate) legend: &'a str,
}

impl MoleculeArgs<'_> {
    pub(crate) fn to_ffi(&self) -> ffi::DrawMoleculeData {
        ffi::DrawMoleculeData {
            highlight_atoms: self.highlight_atoms.to_vec(),
            has_highlight_atoms: true,
            highlight_bonds: self.highlight_bonds.unwrap_or_default().to_vec(),
            has_highlight_bonds: self.highlight_bonds.is_some(),
            highlight_atom_colour_keys: self.highlight_atom_colors.keys().copied().collect(),
            highlight_atom_colour_values: self.highlight_atom_colors.values().copied().collect(),
            highlight_bond_colour_keys: self.highlight_bond_colors.keys().copied().collect(),
            highlight_bond_colour_values: self.highlight_bond_colors.values().copied().collect(),
            highlight_atom_radii_keys: self.highlight_atom_radii.keys().copied().collect(),
            highlight_atom_radii_values: self.highlight_atom_radii.values().copied().collect(),
            conf_id: self.conf_id,
            legend: self.legend.to_owned(),
        }
    }
}

impl MolDraw2D {
    pub(crate) fn from_ffi(inner: UniquePtr<ffi::Drawer>) -> MolDraw2D {
        MolDraw2D { inner }
    }

    pub(crate) fn as_ffi(&self) -> &ffi::Drawer {
        self.inner.as_ref().expect("a drawer is never null")
    }

    pub(crate) fn pin_mut(&mut self) -> std::pin::Pin<&mut ffi::Drawer> {
        self.inner.pin_mut()
    }

    /// The canvas width in pixels. Mirrors `MolDraw2D.Width`; a flexible
    /// canvas (-1) takes its size from its first drawing.
    pub fn width(&self) -> i32 {
        ffi::drawer_width(self.as_ffi())
    }

    /// The canvas height in pixels. Mirrors `MolDraw2D.Height`.
    pub fn height(&self) -> i32 {
        ffi::drawer_height(self.as_ffi())
    }

    /// The font size, roughly in pixels. Mirrors `MolDraw2D.FontSize`.
    pub fn font_size(&self) -> f64 {
        ffi::drawer_font_size(self.as_ffi())
    }

    /// Sets the font size, roughly in pixels, which RDKit keeps within the
    /// options' minimum and maximum font sizes. Mirrors
    /// `MolDraw2D.SetFontSize`; `Error::Value` unless it is from 1 to 512.
    pub fn set_font_size(&mut self, new_size: f64) -> Result<()> {
        Ok(ffi::drawer_set_font_size(self.pin_mut(), new_size)?)
    }

    /// The bond line width (the `bond_line_width` option). Mirrors
    /// `MolDraw2D.LineWidth`.
    pub fn line_width(&self) -> f64 {
        ffi::drawer_line_width(self.as_ffi())
    }

    /// Sets the bond line width. Mirrors `MolDraw2D.SetLineWidth`;
    /// `Error::Value` unless it is from 0 to 64.
    pub fn set_line_width(&mut self, width: f64) -> Result<()> {
        Ok(ffi::drawer_set_line_width(self.pin_mut(), width)?)
    }

    /// Whether polygons are filled. Mirrors `MolDraw2D.FillPolys`.
    pub fn fill_polys(&self) -> bool {
        ffi::drawer_fill_polys(self.as_ffi())
    }

    /// Sets whether polygons are filled. Mirrors `MolDraw2D.SetFillPolys`.
    pub fn set_fill_polys(&mut self, val: bool) {
        ffi::drawer_set_fill_polys(self.pin_mut(), val);
    }

    /// Whether flexible mode is on. Mirrors `MolDraw2D.FlexiMode`.
    pub fn flexi_mode(&self) -> bool {
        ffi::drawer_flexi_mode(self.as_ffi())
    }

    /// In flexible mode, molecules are always drawn with the default bond
    /// length, font size, and so on. Mirrors `MolDraw2D.SetFlexiMode`.
    pub fn set_flexi_mode(&mut self, mode: bool) {
        ffi::drawer_set_flexi_mode(self.pin_mut(), mode);
    }

    /// The drawing offset in pixels. Mirrors `MolDraw2D.Offset`.
    pub fn offset(&self) -> Point2D {
        ffi::drawer_offset(self.as_ffi())
    }

    /// Sets the drawing offset in pixels. Mirrors `MolDraw2D.SetOffset`;
    /// `Error::Value` for an offset beyond 32767.
    pub fn set_offset(&mut self, x: i32, y: i32) -> Result<()> {
        Ok(ffi::drawer_set_offset(self.pin_mut(), x, y)?)
    }

    /// Sets the colour of drawing and filling. Mirrors `MolDraw2D.SetColour`;
    /// `Error::Value` unless every component is from 0 to 1.
    pub fn set_colour(&mut self, tpl: DrawColour) -> Result<()> {
        Ok(ffi::drawer_set_colour(self.pin_mut(), &tpl)?)
    }

    /// Fills the drawing with the background colour. Mirrors
    /// `MolDraw2D.ClearDrawing`.
    pub fn clear_drawing(&mut self) {
        expect_ok(ffi::drawer_clear_drawing(self.pin_mut()));
    }

    /// A copy of the drawing options. Mirrors `MolDraw2D.drawOptions`, which
    /// returns the options to edit in place; set them back with
    /// `set_draw_options`.
    pub fn draw_options(&self) -> MolDrawOptions {
        let data = expect_ok(ffi::drawer_draw_options(self.as_ffi()));
        expect_ok_options(MolDrawOptions::from_ffi(data))
    }

    /// Replaces the drawing options. Mirrors `MolDraw2D.SetDrawOptions`;
    /// `Error::Value` for options RDKit's layout cannot use (see the module
    /// documentation), leaving the options unchanged.
    pub fn set_draw_options(&mut self, opts: &MolDrawOptions) -> Result<()> {
        Ok(ffi::drawer_set_draw_options(
            self.pin_mut(),
            &opts.to_ffi(),
        )?)
    }

    /// Draws a molecule. Mirrors `MolDraw2D.DrawMolecule(mol)`.
    ///
    /// A molecule without conformers is drawn with 2D coordinates RDKit's
    /// depictor computes, as in Python. `Error::Value` for a molecule,
    /// legend, or canvas RDKit cannot draw safely (see the module
    /// documentation).
    pub fn draw_molecule(&mut self, mol: &Mol) -> Result<()> {
        self.draw_molecule_with_params(mol, &DrawMoleculeParams::default())
    }

    /// Draws a molecule with highlights, a conformer, or a legend. Mirrors
    /// `MolDraw2D.DrawMolecule`, both overloads (see `DrawMoleculeParams`).
    /// `Error::Index` for a highlighted atom or bond past the end, and
    /// `Error::Value` for a missing conformer, a colour outside 0 to 1, or
    /// a radius outside 0 to 16, and as `draw_molecule`.
    pub fn draw_molecule_with_params(
        &mut self,
        mol: &Mol,
        params: &DrawMoleculeParams,
    ) -> Result<()> {
        let args = MoleculeArgs {
            highlight_atoms: &params.highlight_atoms,
            highlight_bonds: params.highlight_bonds.as_deref(),
            highlight_atom_colors: &params.highlight_atom_colors,
            highlight_bond_colors: &params.highlight_bond_colors,
            highlight_atom_radii: &params.highlight_atom_radii,
            conf_id: params.conf_id,
            legend: &params.legend,
        };
        Ok(ffi::drawer_draw_molecule(
            self.pin_mut(),
            mol.as_ffi(),
            &args.to_ffi(),
        )?)
    }

    /// The width and height a drawing of the molecule would take at the
    /// drawer's settings. Mirrors `MolDraw2D.GetMolSize(mol)`; errors as
    /// `draw_molecule`.
    pub fn get_mol_size(&mut self, mol: &Mol) -> Result<(i32, i32)> {
        self.get_mol_size_with_params(mol, &GetMolSizeParams::default())
    }

    /// As `get_mol_size`, with highlights, a conformer, or a legend.
    /// Mirrors `MolDraw2D.GetMolSize`; errors as
    /// `draw_molecule_with_params`.
    pub fn get_mol_size_with_params(
        &mut self,
        mol: &Mol,
        params: &GetMolSizeParams,
    ) -> Result<(i32, i32)> {
        let args = MoleculeArgs {
            highlight_atoms: &params.highlight_atoms,
            highlight_bonds: Some(params.highlight_bonds.as_slice()),
            highlight_atom_colors: &params.highlight_atom_colors,
            highlight_bond_colors: &params.highlight_bond_colors,
            highlight_atom_radii: &params.highlight_atom_radii,
            conf_id: params.conf_id,
            legend: &params.legend,
        };
        let size = ffi::drawer_get_mol_size(self.pin_mut(), mol.as_ffi(), &args.to_ffi())?;
        Ok((size.width, size.height))
    }
}

impl fmt::Debug for MolDraw2D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MolDraw2D")
            .field("width", &self.width())
            .field("height", &self.height())
            .finish_non_exhaustive()
    }
}

fn expect_ok_options(result: Result<MolDrawOptions>) -> MolDrawOptions {
    result.unwrap_or_else(|err| panic!("unexpected error reading drawing options: {err}"))
}

/// An SVG drawer. Mirrors `rdMolDraw2D.MolDraw2DSVG`; dereferences to
/// `MolDraw2D` for the drawing methods.
#[derive(Debug)]
pub struct MolDraw2DSVG(MolDraw2D);

/// Parameters for `MolDraw2DSVG::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolDraw2DSVGParams {
    /// The width of each panel of a grid of molecules; -1 for the canvas
    /// width.
    pub panel_width: i32,
    /// The height of each panel; -1 for the canvas height.
    pub panel_height: i32,
    /// Draws text with RDKit's own metrics, as SVG text, rather than as
    /// FreeType glyph outlines. Its text must then be ASCII.
    pub no_freetype: bool,
}

impl Default for MolDraw2DSVGParams {
    fn default() -> MolDraw2DSVGParams {
        MolDraw2DSVGParams {
            panel_width: -1,
            panel_height: -1,
            no_freetype: false,
        }
    }
}

/// Parameters for `MolDraw2DSVG::tag_atoms_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct TagAtomsParams {
    /// The radius of each atom's selection circle, in molecule coordinates.
    pub radius: f64,
    /// Event attributes for each atom's circle: each `name` is set to
    /// `handler(index);`.
    pub events: BTreeMap<String, String>,
}

impl Default for TagAtomsParams {
    fn default() -> TagAtomsParams {
        TagAtomsParams {
            radius: 0.2,
            events: BTreeMap::new(),
        }
    }
}

impl MolDraw2DSVG {
    /// A drawer for a `width` by `height` canvas, in pixels; -1 for both
    /// sizes the canvas to its first drawing. Mirrors
    /// `MolDraw2DSVG(width, height)`; `Error::Value` for a size other than
    /// -1 or 1 to 32767, or only one of them -1.
    pub fn new(width: i32, height: i32) -> Result<MolDraw2DSVG> {
        MolDraw2DSVG::new_with_params(width, height, &MolDraw2DSVGParams::default())
    }

    /// As `new`, with panels for grids or without FreeType. Mirrors
    /// `MolDraw2DSVG(width, height, panelWidth, panelHeight, noFreetype)`;
    /// `Error::Value` for a panel size other than -1 or 1 to 32767, or
    /// fixed panels on a flexible canvas, and as `new`.
    pub fn new_with_params(
        width: i32,
        height: i32,
        params: &MolDraw2DSVGParams,
    ) -> Result<MolDraw2DSVG> {
        let inner = ffi::mol_draw_2d_svg_new(
            width,
            height,
            params.panel_width,
            params.panel_height,
            params.no_freetype,
        )?;
        Ok(MolDraw2DSVG(MolDraw2D::from_ffi(inner)))
    }

    /// Ends the SVG. Mirrors `MolDraw2DSVG.FinishDrawing`.
    pub fn finish_drawing(&mut self) {
        expect_ok(ffi::svg_finish_drawing(self.0.pin_mut()));
    }

    /// The SVG so far. Mirrors `MolDraw2DSVG.GetDrawingText`;
    /// `Error::Value` if it is not UTF-8 (text RDKit wrote byte by byte,
    /// without FreeType), as Python raises `UnicodeDecodeError`.
    pub fn get_drawing_text(&self) -> Result<String> {
        to_string(expect_ok(ffi::svg_drawing_text(self.0.as_ffi())))
    }

    /// Adds RDKit's metadata for the molecule (its atoms and bonds as
    /// SMILES, and where they were drawn) to the SVG. Mirrors
    /// `MolDraw2DSVG.AddMoleculeMetadata(mol)`; `Error::Value` for a
    /// molecule without conformers, or a drawer that has not drawn the
    /// molecule it places the atoms by.
    pub fn add_molecule_metadata(&mut self, mol: &Mol) -> Result<()> {
        self.add_molecule_metadata_with_conf_id(mol, -1)
    }

    /// As `add_molecule_metadata`, for the conformer `conf_id`. Mirrors
    /// `MolDraw2DSVG.AddMoleculeMetadata(mol, confId)`.
    pub fn add_molecule_metadata_with_conf_id(&mut self, mol: &Mol, conf_id: i32) -> Result<()> {
        Ok(ffi::svg_add_molecule_metadata(
            self.0.pin_mut(),
            mol.as_ffi(),
            conf_id,
        )?)
    }

    /// Adds invisible circles and lines over the last molecule drawn, for
    /// selecting its atoms and bonds. Mirrors `MolDraw2DSVG.TagAtoms(mol)`;
    /// `Error::Runtime` (RDKit's own check) when no molecule was drawn or
    /// `mol` has atoms the drawing lacks.
    pub fn tag_atoms(&mut self, mol: &Mol) -> Result<()> {
        self.tag_atoms_with_params(mol, &TagAtomsParams::default())
    }

    /// As `tag_atoms`, with a radius and event handlers. Mirrors
    /// `MolDraw2DSVG.TagAtoms`.
    pub fn tag_atoms_with_params(&mut self, mol: &Mol, params: &TagAtomsParams) -> Result<()> {
        let keys: Vec<String> = params.events.keys().cloned().collect();
        let values: Vec<String> = params.events.values().cloned().collect();
        Ok(ffi::svg_tag_atoms(
            self.0.pin_mut(),
            mol.as_ffi(),
            params.radius,
            &keys,
            &values,
        )?)
    }
}

impl Deref for MolDraw2DSVG {
    type Target = MolDraw2D;

    fn deref(&self) -> &MolDraw2D {
        &self.0
    }
}

impl DerefMut for MolDraw2DSVG {
    fn deref_mut(&mut self) -> &mut MolDraw2D {
        &mut self.0
    }
}

/// Parameters for `prepare_mol_for_drawing_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct PrepareMolForDrawingParams {
    /// Kekulizes the molecule, keeping its aromatic flags.
    pub kekulize: bool,
    /// Adds hydrogens to chiral atoms in two or more rings, so wedges can
    /// be drawn to them.
    pub add_chiral_hs: bool,
    /// Wedges the bonds at chiral atoms.
    pub wedge_bonds: bool,
    /// Computes 2D coordinates even when the molecule has a conformer.
    pub force_coords: bool,
    /// Draws double bonds with unknown stereo as wavy bonds.
    pub wavy_bonds: bool,
}

impl Default for PrepareMolForDrawingParams {
    fn default() -> PrepareMolForDrawingParams {
        PrepareMolForDrawingParams {
            kekulize: true,
            add_chiral_hs: true,
            wedge_bonds: true,
            force_coords: false,
            wavy_bonds: false,
        }
    }
}

/// A copy of the molecule prepared for drawing: kekulized, with hydrogens
/// on chiral ring atoms, wedged bonds, and 2D coordinates if it has no
/// conformer. Mirrors `rdMolDraw2D.PrepareMolForDrawing(mol)`;
/// `Error::Value` for a molecule RDKit's depictor cannot lay out
/// (`rd_depictor::compute_2d_coords` lists what it checks).
pub fn prepare_mol_for_drawing(mol: &Mol) -> Result<Mol> {
    prepare_mol_for_drawing_with_params(mol, &PrepareMolForDrawingParams::default())
}

/// As `prepare_mol_for_drawing`, choosing the steps. Mirrors
/// `rdMolDraw2D.PrepareMolForDrawing`.
pub fn prepare_mol_for_drawing_with_params(
    mol: &Mol,
    params: &PrepareMolForDrawingParams,
) -> Result<Mol> {
    let res = ffi::prepare_mol_for_drawing(
        mol.as_ffi(),
        params.kekulize,
        params.add_chiral_hs,
        params.wedge_bonds,
        params.force_coords,
        params.wavy_bonds,
    )?;
    Ok(Mol::from_ffi(res))
}

/// Parameters for `prepare_and_draw_molecule_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct PrepareAndDrawMoleculeParams {
    /// Text drawn under (or by) the molecule.
    pub legend: String,
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight.
    pub highlight_bonds: Vec<u32>,
    /// Highlight colours of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colours of bonds, by index.
    pub highlight_bond_colors: BTreeMap<u32, DrawColour>,
    /// Highlight radii of atoms, by index, in molecule coordinates.
    pub highlight_atom_radii: BTreeMap<u32, f64>,
    /// The conformer to draw; -1 for the first.
    pub conf_id: i32,
    /// Kekulizes the molecule, keeping its aromatic flags.
    pub kekulize: bool,
}

impl Default for PrepareAndDrawMoleculeParams {
    fn default() -> PrepareAndDrawMoleculeParams {
        PrepareAndDrawMoleculeParams {
            legend: String::new(),
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            highlight_atom_colors: BTreeMap::new(),
            highlight_bond_colors: BTreeMap::new(),
            highlight_atom_radii: BTreeMap::new(),
            conf_id: -1,
            kekulize: true,
        }
    }
}

/// Prepares a copy of the molecule (as `prepare_mol_for_drawing`) and draws
/// it. Mirrors `rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)`; errors as
/// `prepare_mol_for_drawing` and `MolDraw2D::draw_molecule`.
pub fn prepare_and_draw_molecule(drawer: &mut MolDraw2D, mol: &Mol) -> Result<()> {
    prepare_and_draw_molecule_with_params(drawer, mol, &PrepareAndDrawMoleculeParams::default())
}

/// As `prepare_and_draw_molecule`, with highlights, a conformer, or a
/// legend. Mirrors `rdMolDraw2D.PrepareAndDrawMolecule`; errors as
/// `MolDraw2D::draw_molecule_with_params`.
pub fn prepare_and_draw_molecule_with_params(
    drawer: &mut MolDraw2D,
    mol: &Mol,
    params: &PrepareAndDrawMoleculeParams,
) -> Result<()> {
    let args = MoleculeArgs {
        highlight_atoms: &params.highlight_atoms,
        highlight_bonds: Some(params.highlight_bonds.as_slice()),
        highlight_atom_colors: &params.highlight_atom_colors,
        highlight_bond_colors: &params.highlight_bond_colors,
        highlight_atom_radii: &params.highlight_atom_radii,
        conf_id: params.conf_id,
        legend: &params.legend,
    };
    Ok(ffi::drawer_prepare_and_draw_molecule(
        drawer.pin_mut(),
        mol.as_ffi(),
        &args.to_ffi(),
        params.kekulize,
    )?)
}

/// Parameters for `mol_to_svg_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolToSVGParams {
    /// The width in pixels.
    pub width: u32,
    /// The height in pixels.
    pub height: u32,
    /// Atoms to highlight (and the bonds between them).
    pub highlight_atoms: Vec<u32>,
    /// Unused, as in RDKit.
    pub kekulize: bool,
    /// Multiplies the bond line width.
    pub line_width_mult: u32,
    /// Circles highlighted atoms.
    pub include_atom_circles: bool,
    /// The conformer to draw; -1 for the first.
    pub conf_id: i32,
}

impl Default for MolToSVGParams {
    fn default() -> MolToSVGParams {
        MolToSVGParams {
            width: 300,
            height: 300,
            highlight_atoms: Vec::new(),
            kekulize: true,
            line_width_mult: 1,
            include_atom_circles: true,
            conf_id: -1,
        }
    }
}

/// An SVG drawing of the molecule, 300 pixels square. Mirrors
/// `rdMolDraw2D.MolToSVG(mol)`; errors as `MolDraw2D::draw_molecule`.
pub fn mol_to_svg(mol: &Mol) -> Result<String> {
    mol_to_svg_with_params(mol, &MolToSVGParams::default())
}

/// As `mol_to_svg`, with a size, highlights, or a conformer. Mirrors
/// `rdMolDraw2D.MolToSVG`; `Error::Value` for a size beyond 32767 (or 0),
/// or a line width multiplier that makes the line wider than 64, and as
/// `MolDraw2D::draw_molecule_with_params`.
pub fn mol_to_svg_with_params(mol: &Mol, params: &MolToSVGParams) -> Result<String> {
    // Sizes past i32::MAX are refused as any size past 32767.
    let size = |v: u32| i32::try_from(v).unwrap_or(i32::MAX);
    let text = ffi::mol_to_svg(
        mol.as_ffi(),
        size(params.width),
        size(params.height),
        &params.highlight_atoms,
        true,
        params.line_width_mult,
        params.include_atom_circles,
        params.conf_id,
    )?;
    to_string(text)
}

/// The mean bond length of the molecule's first conformer; 0 without bonds.
/// Mirrors `rdMolDraw2D.MeanBondLength(mol)`; `Error::Value` for a molecule
/// with bonds and no conformers.
pub fn mean_bond_length(mol: &Mol) -> Result<f64> {
    mean_bond_length_with_conf_id(mol, -1)
}

/// As `mean_bond_length`, for the conformer `conf_id`. Mirrors
/// `rdMolDraw2D.MeanBondLength(mol, confId)`.
pub fn mean_bond_length_with_conf_id(mol: &Mol, conf_id: i32) -> Result<f64> {
    Ok(ffi::mean_bond_length(mol.as_ffi(), conf_id)?)
}

/// Updates the drawer's options from JSON, as
/// `update_mol_draw_options_from_json` does. Mirrors
/// `rdMolDraw2D.UpdateDrawerParamsFromJSON`; errors as that function and
/// `MolDraw2D::set_draw_options`, leaving the options unchanged.
pub fn update_drawer_params_from_json(drawer: &mut MolDraw2D, json: &str) -> Result<()> {
    Ok(ffi::drawer_update_params_from_json(drawer.pin_mut(), json)?)
}

/// Sets dark-mode colours on the drawer's options. Mirrors
/// `rdMolDraw2D.SetDarkMode(MolDraw2D)`.
pub fn set_dark_mode_mol_draw_2d(d2d: &mut MolDraw2D) {
    ffi::drawer_set_dark_mode(d2d.pin_mut());
}

/// Draws everything in `fg_colour` on `bg_colour`. Mirrors
/// `rdMolDraw2D.SetMonochromeMode(MolDraw2D, ...)`; `Error::Value` unless
/// every component is from 0 to 1.
pub fn set_monochrome_mode_mol_draw_2d(
    drawer: &mut MolDraw2D,
    fg_colour: DrawColour,
    bg_colour: DrawColour,
) -> Result<()> {
    Ok(ffi::drawer_set_monochrome_mode(
        drawer.pin_mut(),
        &fg_colour,
        &bg_colour,
    )?)
}
