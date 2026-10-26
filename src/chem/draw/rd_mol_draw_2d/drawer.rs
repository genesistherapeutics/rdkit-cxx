// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The drawers, and the functions that draw molecules with them.

use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Deref, DerefMut};
use std::path::Path;

use cxx::UniquePtr;

use super::options::{DrawColour, MolDrawOptions};
use crate::Result;
use crate::chem::Mol;
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::error::{expect_ok, to_string};
use crate::ffi::draw::ffi;
use crate::ffi::geometry::ffi::IndexList;
use crate::geometry::Point2D;
use crate::util::path_bytes;

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
    /// Highlight colors of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colors of bonds, by index.
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
    /// Highlight colors of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colors of bonds, by index.
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

/// Parameters for `MolDraw2D::draw_molecules_with_params`, with Python's
/// defaults (A2). Each list, when given, has an entry for each molecule;
/// an empty list is not given, as in Python.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DrawMoleculesParams {
    /// Atoms to highlight in each molecule.
    pub highlight_atoms: Option<Vec<Vec<u32>>>,
    /// Bonds to highlight in each molecule; without them, and with
    /// `continuous_highlight`, the bonds between highlighted atoms.
    pub highlight_bonds: Option<Vec<Vec<u32>>>,
    /// Highlight colors of each molecule's atoms.
    pub highlight_atom_colors: Option<Vec<BTreeMap<u32, DrawColour>>>,
    /// Highlight colors of each molecule's bonds.
    pub highlight_bond_colors: Option<Vec<BTreeMap<u32, DrawColour>>>,
    /// Highlight radii of each molecule's atoms.
    pub highlight_atom_radii: Option<Vec<BTreeMap<u32, f64>>>,
    /// The conformer of each molecule to draw.
    pub conf_ids: Option<Vec<i32>>,
    /// The legend of each molecule.
    pub legends: Option<Vec<String>>,
}

/// Parameters for `MolDraw2D::draw_reaction_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DrawReactionParams {
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

    /// Sets the color of drawing and filling. Mirrors `MolDraw2D.SetColour`;
    /// `Error::Value` unless every component is from 0 to 1, and
    /// `Error::Runtime` for a Cairo drawer on a flexible canvas before
    /// anything is drawn on it, as in Python.
    pub fn set_colour(&mut self, tpl: DrawColour) -> Result<()> {
        Ok(ffi::drawer_set_colour(self.pin_mut(), &tpl)?)
    }

    /// Fills the drawing with the background color. Mirrors
    /// `MolDraw2D.ClearDrawing`; `Error::Value` on a flexible canvas before
    /// anything is drawn on it, which RDKit would set up at -1 by -1
    /// pixels, for good.
    pub fn clear_drawing(&mut self) -> Result<()> {
        Ok(ffi::drawer_clear_drawing(self.pin_mut())?)
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
    /// `Error::Value` for a missing conformer, a color outside 0 to 1, or
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

    /// Draws a molecule with highlights of several colors: each atom in
    /// `highlight_atom_map` and bond in `highlight_bond_map` with its
    /// colors, in the `multi_colour_highlight_style` of the options, atoms
    /// with the radii in `highlight_radii` (else `highlight_radius`) and
    /// bonds with their line widths multiplied as in
    /// `highlight_linewidth_multipliers`. Mirrors
    /// `MolDraw2D.DrawMoleculeWithHighlights`.
    ///
    /// `Error::Index` for an atom or bond past the end, and `Error::Value`
    /// for a color outside 0 to 1, a radius not above 0 or beyond 16 (or a
    /// `highlight_radius` of 0: RDKit grows a highlight of no radius until
    /// it meets its bonds, which never happens), a multiplier outside 0 to
    /// 64, and as `draw_molecule`.
    pub fn draw_molecule_with_highlights(
        &mut self,
        mol: &Mol,
        legend: &str,
        highlight_atom_map: &BTreeMap<u32, Vec<DrawColour>>,
        highlight_bond_map: &BTreeMap<u32, Vec<DrawColour>>,
        highlight_radii: &BTreeMap<u32, f64>,
        highlight_linewidth_multipliers: &BTreeMap<u32, i32>,
    ) -> Result<()> {
        self.draw_molecule_with_highlights_with_conf_id(
            mol,
            legend,
            highlight_atom_map,
            highlight_bond_map,
            highlight_radii,
            highlight_linewidth_multipliers,
            -1,
        )
    }

    /// As `draw_molecule_with_highlights`, with Python's `confId`. Mirrors
    /// `MolDraw2D.DrawMoleculeWithHighlights(..., confId)`. RDKit passes
    /// `confId` where its layout takes whether to draw annotations, so it
    /// draws the first conformer whatever `conf_id` is, and no annotations
    /// for a `conf_id` of 0, as in Python.
    #[allow(clippy::too_many_arguments)]
    pub fn draw_molecule_with_highlights_with_conf_id(
        &mut self,
        mol: &Mol,
        legend: &str,
        highlight_atom_map: &BTreeMap<u32, Vec<DrawColour>>,
        highlight_bond_map: &BTreeMap<u32, Vec<DrawColour>>,
        highlight_radii: &BTreeMap<u32, f64>,
        highlight_linewidth_multipliers: &BTreeMap<u32, i32>,
        conf_id: i32,
    ) -> Result<()> {
        let lists = |map: &BTreeMap<u32, Vec<DrawColour>>| -> Vec<ffi::ColourList> {
            map.values()
                .map(|colours| ffi::ColourList {
                    colours: colours.clone(),
                })
                .collect()
        };
        let data = ffi::DrawMoleculeWithHighlightsData {
            legend: legend.to_owned(),
            atom_keys: highlight_atom_map.keys().copied().collect(),
            atom_colours: lists(highlight_atom_map),
            bond_keys: highlight_bond_map.keys().copied().collect(),
            bond_colours: lists(highlight_bond_map),
            radii_keys: highlight_radii.keys().copied().collect(),
            radii_values: highlight_radii.values().copied().collect(),
            multiplier_keys: highlight_linewidth_multipliers.keys().copied().collect(),
            multiplier_values: highlight_linewidth_multipliers.values().copied().collect(),
            conf_id,
        };
        Ok(ffi::drawer_draw_molecule_with_highlights(
            self.pin_mut(),
            mol.as_ffi(),
            &data,
        )?)
    }

    /// Draws the molecules in a grid of panels (the panel size the drawer
    /// was made with). Mirrors `MolDraw2D.DrawMolecules(mols)`.
    ///
    /// The drawer needs a canvas of fixed size, panels that fit on it, and
    /// nothing drawn yet: RDKit lays out again every molecule a drawer holds,
    /// with the options of the moment, when it draws a grid. `Error::Value`
    /// otherwise, and as `draw_molecule` for each molecule.
    pub fn draw_molecules(&mut self, mols: &[Mol]) -> Result<()> {
        self.draw_molecules_with_params(mols, &DrawMoleculesParams::default())
    }

    /// As `draw_molecules`, with highlights, conformers, or legends.
    /// Mirrors `MolDraw2D.DrawMolecules`; `Error::Value` for a list whose
    /// length is not the number of molecules, and as
    /// `draw_molecule_with_params`.
    pub fn draw_molecules_with_params(
        &mut self,
        mols: &[Mol],
        params: &DrawMoleculesParams,
    ) -> Result<()> {
        fn given<T>(list: &Option<Vec<T>>) -> Option<&Vec<T>> {
            list.as_ref().filter(|v| !v.is_empty())
        }
        let colour_maps =
            |list: &Option<Vec<BTreeMap<u32, DrawColour>>>| -> Vec<ffi::IndexColourMap> {
                given(list)
                    .map(|maps| {
                        maps.iter()
                            .map(|m| ffi::IndexColourMap {
                                keys: m.keys().copied().collect(),
                                values: m.values().copied().collect(),
                            })
                            .collect()
                    })
                    .unwrap_or_default()
            };
        let index_lists = |list: &Option<Vec<Vec<u32>>>| -> Vec<IndexList> {
            given(list)
                .map(|lists| {
                    lists
                        .iter()
                        .map(|v| IndexList { indices: v.clone() })
                        .collect()
                })
                .unwrap_or_default()
        };
        let data = ffi::DrawMoleculesData {
            has_highlight_atoms: given(&params.highlight_atoms).is_some(),
            highlight_atoms: index_lists(&params.highlight_atoms),
            has_highlight_bonds: given(&params.highlight_bonds).is_some(),
            highlight_bonds: index_lists(&params.highlight_bonds),
            has_atom_colours: given(&params.highlight_atom_colors).is_some(),
            atom_colours: colour_maps(&params.highlight_atom_colors),
            has_bond_colours: given(&params.highlight_bond_colors).is_some(),
            bond_colours: colour_maps(&params.highlight_bond_colors),
            has_radii: given(&params.highlight_atom_radii).is_some(),
            radii: given(&params.highlight_atom_radii)
                .map(|maps| {
                    maps.iter()
                        .map(|m| ffi::IndexRadiusMap {
                            keys: m.keys().copied().collect(),
                            values: m.values().copied().collect(),
                        })
                        .collect()
                })
                .unwrap_or_default(),
            has_conf_ids: given(&params.conf_ids).is_some(),
            conf_ids: given(&params.conf_ids).cloned().unwrap_or_default(),
            has_legends: given(&params.legends).is_some(),
            legends: given(&params.legends).cloned().unwrap_or_default(),
        };
        let inner = self.pin_mut();
        with_borrowed_mols_shared(mols.iter().map(|m| ("", m)), |list| {
            Ok(ffi::drawer_draw_molecules(inner, list, &data)?)
        })
    }

    /// Draws a reaction: its reactants, an arrow with its agents, and its
    /// products. Mirrors `MolDraw2D.DrawReaction(rxn)`.
    ///
    /// The drawer needs nothing drawn yet: RDKit lays the reaction out
    /// against every molecule a drawer holds. `Error::Value` otherwise, for
    /// a reaction with neither reactants nor products, and as
    /// `draw_molecule` for each template.
    pub fn draw_reaction(&mut self, rxn: &ChemicalReaction) -> Result<()> {
        self.draw_reaction_with_params(rxn, &DrawReactionParams::default())
    }

    /// As `draw_reaction`, with highlights or conformers. Mirrors
    /// `MolDraw2D.DrawReaction`; `Error::Value` for an empty list of
    /// colors with `highlight_by_reactant` (RDKit takes each reactant's
    /// color modulo their number), fewer conformer ids than the largest
    /// group of templates, or a color outside 0 to 1.
    pub fn draw_reaction_with_params(
        &mut self,
        rxn: &ChemicalReaction,
        params: &DrawReactionParams,
    ) -> Result<()> {
        let empty_colours = Vec::new();
        let empty_ids = Vec::new();
        Ok(ffi::drawer_draw_reaction(
            self.pin_mut(),
            rxn.as_ffi(),
            params.highlight_by_reactant,
            params
                .highlight_colors_reactants
                .as_ref()
                .unwrap_or(&empty_colours),
            params.highlight_colors_reactants.is_some(),
            params.conf_ids.as_ref().unwrap_or(&empty_ids),
            params.conf_ids.is_some(),
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

    /// Draws a line in the current color, from points in the molecule's
    /// frame (by the last molecule drawn, or `set_scale`). Mirrors
    /// `MolDraw2D.DrawLine(cds1, cds2)`; see the free-drawing errors in the
    /// module documentation.
    pub fn draw_line(&mut self, cds1: &Point2D, cds2: &Point2D) -> Result<()> {
        self.draw_line_with_raw_coords(cds1, cds2, false)
    }

    /// As `draw_line`, from points in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawLine(cds1, cds2, rawCoords)`.
    pub fn draw_line_with_raw_coords(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_line(
            self.pin_mut(),
            cds1,
            cds2,
            raw_coords,
        )?)
    }

    /// Draws an arrow from `cds1` to `cds2`. Mirrors
    /// `MolDraw2D.DrawArrow(cds1, cds2)`.
    pub fn draw_arrow(&mut self, cds1: &Point2D, cds2: &Point2D) -> Result<()> {
        self.draw_arrow_with_params(cds1, cds2, &DrawArrowParams::default())
    }

    /// As `draw_arrow`, with the head's shape and color. Mirrors
    /// `MolDraw2D.DrawArrow`; `Error::Runtime` for an angle not above 1e-6
    /// (RDKit's check), and `Error::Value` for a head whose points come out
    /// beyond the coordinate bounds.
    pub fn draw_arrow_with_params(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        params: &DrawArrowParams,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_arrow(
            self.pin_mut(),
            cds1,
            cds2,
            params.as_polygon,
            params.frac,
            params.angle,
            &params.color,
            params.raw_coords,
        )?)
    }

    /// Draws a triangle in the current color, filled when `fill_polys`
    /// is. Mirrors `MolDraw2D.DrawTriangle(cds1, cds2, cds3)`.
    pub fn draw_triangle(&mut self, cds1: &Point2D, cds2: &Point2D, cds3: &Point2D) -> Result<()> {
        self.draw_triangle_with_raw_coords(cds1, cds2, cds3, false)
    }

    /// As `draw_triangle`, in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawTriangle(cds1, cds2, cds3, rawCoords)`.
    pub fn draw_triangle_with_raw_coords(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        cds3: &Point2D,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_triangle(
            self.pin_mut(),
            cds1,
            cds2,
            cds3,
            raw_coords,
        )?)
    }

    /// Draws a polygon through the points. Mirrors
    /// `MolDraw2D.DrawPolygon(cds)`; as in Python, `Error::Value` for no
    /// points and `Error::Runtime` (RDKit's check) for one or two.
    pub fn draw_polygon(&mut self, cds: &[Point2D]) -> Result<()> {
        self.draw_polygon_with_raw_coords(cds, false)
    }

    /// As `draw_polygon`, in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawPolygon(cds, rawCoords)`.
    pub fn draw_polygon_with_raw_coords(
        &mut self,
        cds: &[Point2D],
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_polygon(self.pin_mut(), cds, raw_coords)?)
    }

    /// Draws the ellipse in the rectangle with corners `cds1` and `cds2`.
    /// Mirrors `MolDraw2D.DrawEllipse(cds1, cds2)`.
    pub fn draw_ellipse(&mut self, cds1: &Point2D, cds2: &Point2D) -> Result<()> {
        self.draw_ellipse_with_raw_coords(cds1, cds2, false)
    }

    /// As `draw_ellipse`, in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawEllipse(cds1, cds2, rawCoords)`.
    pub fn draw_ellipse_with_raw_coords(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_ellipse(
            self.pin_mut(),
            cds1,
            cds2,
            raw_coords,
        )?)
    }

    /// Draws the rectangle with corners `cds1` and `cds2`. Mirrors
    /// `MolDraw2D.DrawRect(cds1, cds2)`.
    pub fn draw_rect(&mut self, cds1: &Point2D, cds2: &Point2D) -> Result<()> {
        self.draw_rect_with_raw_coords(cds1, cds2, false)
    }

    /// As `draw_rect`, in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawRect(cds1, cds2, rawCoords)`.
    pub fn draw_rect_with_raw_coords(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_rect(
            self.pin_mut(),
            cds1,
            cds2,
            raw_coords,
        )?)
    }

    /// Draws an arc of a circle from `angle1` to `angle2`, in degrees.
    /// Mirrors `MolDraw2D.DrawArc(center, radius, angle1, angle2)`;
    /// `Error::Value` for angles beyond 10000 degrees or `angle2` below
    /// `angle1`, which RDKit steps through in too few steps (dividing by
    /// zero, or with too few points to draw).
    pub fn draw_arc(
        &mut self,
        center: &Point2D,
        radius: f64,
        angle1: f64,
        angle2: f64,
    ) -> Result<()> {
        self.draw_arc_with_raw_coords(center, radius, angle1, angle2, false)
    }

    /// As `draw_arc`, in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawArc(center, radius, angle1, angle2, rawCoords)`.
    pub fn draw_arc_with_raw_coords(
        &mut self,
        center: &Point2D,
        radius: f64,
        angle1: f64,
        angle2: f64,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_arc(
            self.pin_mut(),
            center,
            radius,
            angle1,
            angle2,
            raw_coords,
        )?)
    }

    /// Draws a line across the end of the line from `cds1` to `cds2`, as
    /// for an attachment point. Mirrors
    /// `MolDraw2D.DrawAttachmentLine(cds1, cds2, color)`.
    pub fn draw_attachment_line(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        color: DrawColour,
    ) -> Result<()> {
        self.draw_attachment_line_with_params(
            cds1,
            cds2,
            color,
            &DrawAttachmentLineParams::default(),
        )
    }

    /// As `draw_attachment_line`, with its length and segments. Mirrors
    /// `MolDraw2D.DrawAttachmentLine`. As in RDKit, `raw_coords` is passed
    /// where the line takes its waves' height: the line is always drawn in
    /// the molecule's frame (so it needs a molecule drawn, or `set_scale`),
    /// waving by 1 with `raw_coords` and not at all without.
    pub fn draw_attachment_line_with_params(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        color: DrawColour,
        params: &DrawAttachmentLineParams,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_attachment_line(
            self.pin_mut(),
            cds1,
            cds2,
            &color,
            params.len,
            params.n_segments,
            params.raw_coords,
        )?)
    }

    /// Draws a wavy line, in `color1` (`color2` is unused, as in RDKit's
    /// SVG and Cairo drawers). Mirrors
    /// `MolDraw2D.DrawWavyLine(cds1, cds2, color1, color2)`.
    pub fn draw_wavy_line(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        color1: DrawColour,
        color2: DrawColour,
    ) -> Result<()> {
        self.draw_wavy_line_with_params(cds1, cds2, color1, color2, &DrawWavyLineParams::default())
    }

    /// As `draw_wavy_line`, with its segments and waves' height. Mirrors
    /// `MolDraw2D.DrawWavyLine`; `Error::Runtime` for fewer than 2 segments
    /// (RDKit's check), and `Error::Value` for more than 100000 (RDKit
    /// rounds an odd count up, which makes `u32::MAX` none, and then reads
    /// the first) or a line of no length.
    pub fn draw_wavy_line_with_params(
        &mut self,
        cds1: &Point2D,
        cds2: &Point2D,
        color1: DrawColour,
        color2: DrawColour,
        params: &DrawWavyLineParams,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_wavy_line(
            self.pin_mut(),
            cds1,
            cds2,
            &color1,
            &color2,
            params.n_segments,
            params.vert_offset,
            params.raw_coords,
        )?)
    }

    /// Draws text centered on `pos`, laid out as an atom label is (split at
    /// its atoms and charges, `<sub>` and `<sup>` markup read). Mirrors
    /// `MolDraw2D.DrawString(string, pos)`; `Error::Value` for text atom
    /// labels may not have (see the module documentation).
    pub fn draw_string(&mut self, string: &str, pos: &Point2D) -> Result<()> {
        self.draw_string_with_raw_coords(string, pos, false)
    }

    /// As `draw_string`, at a position in pixels with `raw_coords`. Mirrors
    /// `MolDraw2D.DrawString(string, pos, rawCoords)`.
    pub fn draw_string_with_raw_coords(
        &mut self,
        string: &str,
        pos: &Point2D,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_string(
            self.pin_mut(),
            string,
            pos,
            0,
            false,
            raw_coords,
        )?)
    }

    /// Draws text aligned on `pos`: `align` 0 for its middle, 1 for its
    /// start, 2 for its end. Mirrors `MolDraw2D.DrawString(string, pos,
    /// align)`; `Error::Value` for another `align`, as in Python, or text
    /// notes may not have.
    pub fn draw_string_aligned(&mut self, string: &str, pos: &Point2D, align: i32) -> Result<()> {
        self.draw_string_aligned_with_raw_coords(string, pos, align, false)
    }

    /// As `draw_string_aligned`, at a position in pixels with `raw_coords`.
    /// Mirrors `MolDraw2D.DrawString(string, pos, align, rawCoords)`.
    pub fn draw_string_aligned_with_raw_coords(
        &mut self,
        string: &str,
        pos: &Point2D,
        align: i32,
        raw_coords: bool,
    ) -> Result<()> {
        Ok(ffi::drawer_draw_string(
            self.pin_mut(),
            string,
            pos,
            align,
            true,
            raw_coords,
        )?)
    }

    /// Where a point in the molecule's frame is drawn, in pixels: by the
    /// transformation `set_scale` set, or else the last molecule drawn.
    /// Mirrors `MolDraw2D.GetDrawCoords(point)`; `Error::Runtime` with
    /// neither, as in Python.
    pub fn get_draw_coords(&self, point: &Point2D) -> Result<Point2D> {
        Ok(ffi::drawer_get_draw_coords(self.as_ffi(), point)?)
    }

    /// Where an atom of the last molecule drawn is drawn, in pixels.
    /// Mirrors `MolDraw2D.GetDrawCoords(atomIndex)`; `Error::Runtime`
    /// without a molecule drawn, or for an atom it lacks, as in Python.
    pub fn get_draw_coords_int(&self, atom_index: u32) -> Result<Point2D> {
        Ok(ffi::drawer_get_atom_draw_coords(self.as_ffi(), atom_index)?)
    }

    /// Fixes the transformation from the molecule's frame to pixels, for a
    /// `width` by `height` area to show `minv` to `maxv`: free drawing, and
    /// molecules drawn after, use it. Mirrors
    /// `MolDraw2D.SetScale(width, height, minv, maxv)`; `Error::Runtime`
    /// for a size below 1 (RDKit's check), and `Error::Value` for one above
    /// 32767, corners beyond 1e5, or a flexible canvas (or flexi mode),
    /// whose panels RDKit lays the transformation out in.
    pub fn set_scale(
        &mut self,
        width: i32,
        height: i32,
        minv: &Point2D,
        maxv: &Point2D,
    ) -> Result<()> {
        Ok(ffi::drawer_set_scale(
            self.pin_mut(),
            width,
            height,
            minv,
            maxv,
        )?)
    }

    /// As `set_scale`, widening the area to the molecule as RDKit lays it
    /// out. Mirrors `MolDraw2D.SetScale(width, height, minv, maxv, mol)`;
    /// errors as `set_scale` and `draw_molecule`.
    pub fn set_scale_with_mol(
        &mut self,
        width: i32,
        height: i32,
        minv: &Point2D,
        maxv: &Point2D,
        mol: &Mol,
    ) -> Result<()> {
        Ok(ffi::drawer_set_scale_with_mol(
            self.pin_mut(),
            width,
            height,
            minv,
            maxv,
            mol.as_ffi(),
        )?)
    }
}

/// Parameters for `MolDraw2D::draw_arrow_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DrawArrowParams {
    /// Draws the head as a triangle, rather than two lines.
    pub as_polygon: bool,
    /// The head's length, as a fraction of the arrow's.
    pub frac: f64,
    /// The angle between the line and each side of the head, in radians.
    pub angle: f64,
    /// The arrow's color (Python's `None` is black).
    pub color: DrawColour,
    /// Takes the points in pixels.
    pub raw_coords: bool,
}

impl Default for DrawArrowParams {
    fn default() -> DrawArrowParams {
        DrawArrowParams {
            as_polygon: false,
            frac: 0.05,
            angle: std::f64::consts::PI / 6.0,
            color: DrawColour::new(0.0, 0.0, 0.0),
            raw_coords: false,
        }
    }
}

/// Parameters for `MolDraw2D::draw_attachment_line_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DrawAttachmentLineParams {
    /// The line's length, in the molecule's frame.
    pub len: f64,
    /// The number of waves' segments.
    pub n_segments: u32,
    /// Passed where the line takes its waves' height (see
    /// `draw_attachment_line_with_params`).
    pub raw_coords: bool,
}

impl Default for DrawAttachmentLineParams {
    fn default() -> DrawAttachmentLineParams {
        DrawAttachmentLineParams {
            len: 1.0,
            n_segments: 16,
            raw_coords: false,
        }
    }
}

/// Parameters for `MolDraw2D::draw_wavy_line_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DrawWavyLineParams {
    /// The number of segments (rounded up to even).
    pub n_segments: u32,
    /// The waves' height, in the units of the points.
    pub vert_offset: f64,
    /// Takes the points in pixels.
    pub raw_coords: bool,
}

impl Default for DrawWavyLineParams {
    fn default() -> DrawWavyLineParams {
        DrawWavyLineParams {
            n_segments: 16,
            vert_offset: 0.05,
            raw_coords: false,
        }
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

/// A PNG drawer, drawing with Cairo. Mirrors `rdMolDraw2D.MolDraw2DCairo`;
/// dereferences to `MolDraw2D` for the drawing methods. RDKit built without
/// Cairo has none: `new` then returns `Error::Runtime`.
#[derive(Debug)]
pub struct MolDraw2DCairo(MolDraw2D);

/// Parameters for `MolDraw2DCairo::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolDraw2DCairoParams {
    /// The width of each panel of a grid of molecules; -1 for the canvas
    /// width.
    pub panel_width: i32,
    /// The height of each panel; -1 for the canvas height.
    pub panel_height: i32,
    /// Draws text with Cairo's own fonts rather than FreeType's glyph
    /// outlines. Its text must then be ASCII.
    pub no_freetype: bool,
}

impl Default for MolDraw2DCairoParams {
    fn default() -> MolDraw2DCairoParams {
        MolDraw2DCairoParams {
            panel_width: -1,
            panel_height: -1,
            no_freetype: false,
        }
    }
}

impl MolDraw2DCairo {
    /// A drawer for a `width` by `height` canvas, in pixels; -1 for both
    /// sizes the canvas to its first drawing. Mirrors
    /// `MolDraw2DCairo(width, height)`; `Error::Value` for a size other
    /// than -1 or 1 to 32767, only one of them -1, or a canvas of more than
    /// 2^31 bytes (4 a pixel), which Cairo cannot make.
    pub fn new(width: i32, height: i32) -> Result<MolDraw2DCairo> {
        MolDraw2DCairo::new_with_params(width, height, &MolDraw2DCairoParams::default())
    }

    /// As `new`, with panels for grids or without FreeType. Mirrors
    /// `MolDraw2DCairo(width, height, panelWidth, panelHeight, noFreetype)`;
    /// `Error::Value` for a panel size other than -1 or 1 to 32767, or
    /// fixed panels on a flexible canvas, and as `new`.
    pub fn new_with_params(
        width: i32,
        height: i32,
        params: &MolDraw2DCairoParams,
    ) -> Result<MolDraw2DCairo> {
        let inner = ffi::mol_draw_2d_cairo_new(
            width,
            height,
            params.panel_width,
            params.panel_height,
            params.no_freetype,
        )?;
        Ok(MolDraw2DCairo(MolDraw2D::from_ffi(inner)))
    }

    /// Does nothing, as in RDKit: a Cairo drawing needs no finishing.
    /// Mirrors `MolDraw2DCairo.FinishDrawing`.
    pub fn finish_drawing(&mut self) {
        expect_ok(ffi::cairo_finish_drawing(self.0.pin_mut()));
    }

    /// The drawing as PNG data, with RDKit's metadata for what was drawn
    /// (when the `include_metadata` option is on). Mirrors
    /// `MolDraw2DCairo.GetDrawingText`; `Error::Runtime` on a flexible
    /// canvas before anything is drawn on it, as in Python, or when Cairo
    /// could not make the canvas (no memory), where Python returns no PNG.
    pub fn get_drawing_text(&self) -> Result<Vec<u8>> {
        Ok(ffi::cairo_drawing_text(self.0.as_ffi())?)
    }

    /// Writes the PNG data of `get_drawing_text` to a file. Mirrors
    /// `MolDraw2DCairo.WriteDrawingText`; `Error::Os` for a file that
    /// cannot be written, where RDKit logs an error, and as
    /// `get_drawing_text`.
    pub fn write_drawing_text<P: AsRef<Path>>(&self, f_name: P) -> Result<()> {
        Ok(ffi::cairo_write_drawing_text(
            self.0.as_ffi(),
            path_bytes(f_name.as_ref()),
        )?)
    }
}

impl Deref for MolDraw2DCairo {
    type Target = MolDraw2D;

    fn deref(&self) -> &MolDraw2D {
        &self.0
    }
}

impl DerefMut for MolDraw2DCairo {
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
    /// Highlight colors of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colors of bonds, by index.
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

/// Parameters for `draw_molecule_acs1996_with_params` and
/// `mol_to_acs1996_svg_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DrawMoleculeACS1996Params {
    /// Text drawn under (or by) the molecule.
    pub legend: String,
    /// Atoms to highlight.
    pub highlight_atoms: Vec<u32>,
    /// Bonds to highlight (none are highlighted for the atoms alone).
    pub highlight_bonds: Vec<u32>,
    /// Highlight colors of atoms, by index.
    pub highlight_atom_colors: BTreeMap<u32, DrawColour>,
    /// Highlight colors of bonds, by index.
    pub highlight_bond_colors: BTreeMap<u32, DrawColour>,
    /// Highlight radii of atoms, by index, in molecule coordinates.
    pub highlight_atom_radii: BTreeMap<u32, f64>,
    /// The conformer to draw; -1 for the first.
    pub conf_id: i32,
}

/// Parameters for `mol_to_acs1996_svg_with_params`, with Python's defaults
/// (A2): the same as `draw_molecule_acs1996_with_params`'s.
pub type MolToACS1996SVGParams = DrawMoleculeACS1996Params;

impl Default for DrawMoleculeACS1996Params {
    fn default() -> DrawMoleculeACS1996Params {
        DrawMoleculeACS1996Params {
            legend: String::new(),
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            highlight_atom_colors: BTreeMap::new(),
            highlight_bond_colors: BTreeMap::new(),
            highlight_atom_radii: BTreeMap::new(),
            conf_id: -1,
        }
    }
}

impl DrawMoleculeACS1996Params {
    fn args(&self) -> MoleculeArgs<'_> {
        MoleculeArgs {
            highlight_atoms: &self.highlight_atoms,
            highlight_bonds: Some(self.highlight_bonds.as_slice()),
            highlight_atom_colors: &self.highlight_atom_colors,
            highlight_bond_colors: &self.highlight_bond_colors,
            highlight_atom_radii: &self.highlight_atom_radii,
            conf_id: self.conf_id,
            legend: &self.legend,
        }
    }
}

/// Draws the molecule as the ACS 1996 guidelines ask, setting the drawer's
/// options as `set_acs1996_mode` does for its mean bond length (best on a
/// flexible canvas). Mirrors `rdMolDraw2D.DrawMoleculeACS1996(drawer,
/// mol)`; see `draw_molecule_acs1996_with_params` for its errors.
pub fn draw_molecule_acs1996(drawer: &mut MolDraw2D, mol: &Mol) -> Result<()> {
    draw_molecule_acs1996_with_params(drawer, mol, &DrawMoleculeACS1996Params::default())
}

/// As `draw_molecule_acs1996`, with a legend, highlights, or a conformer.
/// Mirrors `rdMolDraw2D.DrawMoleculeACS1996`. `Error::Value` for a
/// molecule without bonds (RDKit's check of the mean bond length), or one
/// whose options come out beyond the drawer's bounds (a mean bond length
/// below 0.0144 makes the scaling factor over 1000), and as
/// `MolDraw2D::draw_molecule_with_params`.
pub fn draw_molecule_acs1996_with_params(
    drawer: &mut MolDraw2D,
    mol: &Mol,
    params: &DrawMoleculeACS1996Params,
) -> Result<()> {
    Ok(ffi::drawer_draw_molecule_acs1996(
        drawer.pin_mut(),
        mol.as_ffi(),
        &params.args().to_ffi(),
    )?)
}

/// An SVG drawing of the molecule as the ACS 1996 guidelines ask, on a
/// flexible canvas. Mirrors `rdMolDraw2D.MolToACS1996SVG(mol)`; errors as
/// `draw_molecule_acs1996_with_params`.
pub fn mol_to_acs1996_svg(mol: &Mol) -> Result<String> {
    mol_to_acs1996_svg_with_params(mol, &MolToACS1996SVGParams::default())
}

/// As `mol_to_acs1996_svg`, with a legend, highlights, or a conformer.
/// Mirrors `rdMolDraw2D.MolToACS1996SVG`.
pub fn mol_to_acs1996_svg_with_params(mol: &Mol, params: &MolToACS1996SVGParams) -> Result<String> {
    to_string(ffi::mol_to_acs1996_svg(
        mol.as_ffi(),
        &params.args().to_ffi(),
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

/// Sets dark-mode colors on the drawer's options. Mirrors
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
