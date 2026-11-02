// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! `MolDrawOptions`, colors, and the enums the options use.

use std::collections::BTreeMap;

use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::expect_ok;
use crate::ffi::draw::ffi;
use crate::ffi::geometry::ffi::IndexList;
use crate::{Error, Result};

pub use crate::ffi::draw::ffi::DrawColour;

/// A color component as Python's bindings convert a tuple's: through C's
/// `float`.
const fn as_float(x: f64) -> f64 {
    x as f32 as f64
}

impl DrawColour {
    /// An opaque color. Mirrors a 3-tuple `(r, g, b)`, whose components
    /// Python's bindings round to `f32`, as this does; a struct literal keeps
    /// them as written.
    pub const fn new(r: f64, g: f64, b: f64) -> DrawColour {
        DrawColour::rgba(r, g, b, 1.0)
    }

    /// A color with an alpha. Mirrors a 4-tuple `(r, g, b, a)`, rounded as
    /// for `new`.
    pub const fn rgba(r: f64, g: f64, b: f64, a: f64) -> DrawColour {
        DrawColour {
            r: as_float(r),
            g: as_float(g),
            b: as_float(b),
            a: as_float(a),
        }
    }

    /// `Error::Value` unless every component is from 0 to 1, as Python's
    /// color setters require.
    pub(crate) fn check(&self) -> Result<()> {
        Ok(ffi::check_draw_colour(self)?)
    }
}

impl From<(f64, f64, f64)> for DrawColour {
    fn from((r, g, b): (f64, f64, f64)) -> DrawColour {
        DrawColour::new(r, g, b)
    }
}

impl From<(f64, f64, f64, f64)> for DrawColour {
    fn from((r, g, b, a): (f64, f64, f64, f64)) -> DrawColour {
        DrawColour::rgba(r, g, b, a)
    }
}

rdkit_enum! {
    /// Where a molecule's legend goes. Mirrors
    /// `rdMolDraw2D.LegendPosition`.
    LegendPosition: ffi::LegendPosition {
        Bottom = Bottom,
        Top = Top,
        Left = Left,
        Right = Right,
    }
}

rdkit_enum! {
    /// How `draw_molecule_with_highlights` draws several colors on one
    /// atom or bond. Mirrors `rdMolDraw2D.MultiColourHighlightStyle`.
    MultiColourHighlightStyle: ffi::MultiColourHighlightStyle {
        CircleAndLine = CIRCLEANDLINE,
        Lasso = LASSO,
    }
}

/// Parts of a drawing, as bits, for `MolDrawOptions::drawing_extents_include`.
/// Mirrors `rdMolDraw2D.DrawElement` (CONVENTIONS.md F1); the constants keep
/// Python's spelling.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct DrawElement(u32);

macro_rules! draw_elements {
    ($($name:ident),* $(,)?) => {
        impl DrawElement {
            $(
                #[doc = concat!("Mirrors `DrawElement.", stringify!($name), "`.")]
                pub const $name: DrawElement = DrawElement(ffi::DrawElement::$name.repr);
            )*

            const NAMES: &[(&str, DrawElement)] = &[$((stringify!($name), DrawElement::$name)),*];
        }
    };
}

draw_elements!(
    NONE,
    PRESHAPES,
    BONDS,
    ATOMLABELS,
    HIGHLIGHTS,
    ANNOTATIONS,
    RADICALS,
    POSTSHAPES,
    ALL,
);

crate::flags::flag_ops!(DrawElement);

/// Drawing options. Mirrors `rdMolDraw2D.MolDrawOptions` (A7): a drawer's
/// options are read with `MolDraw2D::draw_options` and replaced with
/// `MolDraw2D::set_draw_options`, where Python edits them in place.
///
/// Any values may be stored; a drawing call checks them, and refuses with
/// `Error::Value` values RDKit's layout cannot use (see the module
/// documentation).
#[derive(Clone, Debug, PartialEq)]
pub struct MolDrawOptions {
    /// Adds atom indices to drawings.
    pub add_atom_indices: bool,
    /// Adds bond indices to drawings.
    pub add_bond_indices: bool,
    /// Adds R/S and E/Z labels to drawings.
    pub add_stereo_annotation: bool,
    /// Adds the enhanced stereo group labels.
    pub add_stereo_group_annotation: bool,
    /// Additional padding around atom labels, as a fraction of the font
    /// size.
    pub additional_atom_label_padding: f64,
    /// The color of annotations.
    pub annotation_colour: DrawColour,
    /// The size of annotation fonts relative to atom labels.
    pub annotation_font_scale: f64,
    /// Always highlights atoms with circles rather than ellipses around
    /// longer labels.
    pub atom_highlights_are_circles: bool,
    /// Labels deuterium as D and tritium as T.
    pub atom_label_deuterium_tritium: bool,
    /// Replacement labels for atoms, by index.
    pub atom_labels: BTreeMap<u32, String>,
    /// The color of atom indices and notes.
    pub atom_note_colour: DrawColour,
    /// Groups of atoms to outline.
    pub atom_regions: Vec<Vec<u32>>,
    /// The background color.
    pub background_colour: DrawColour,
    /// The relative size of fonts; nonpositive uses RDKit's 0.6.
    pub base_font_size: f64,
    /// The line width for bonds, in pixels.
    pub bond_line_width: f64,
    /// The color of bond indices and notes.
    pub bond_note_colour: DrawColour,
    /// Puts brackets around atom lists in query atoms.
    pub brackets_around_atom_lists: bool,
    /// Moves each molecule's center to the origin before drawing.
    pub centre_molecules_before_drawing: bool,
    /// Draws circles under highlighted atoms.
    pub circle_atoms: bool,
    /// Clears the background before drawing a molecule.
    pub clear_background: bool,
    /// Draws bonds as hand-drawn lines.
    pub comic_mode: bool,
    /// Highlights by drawing an outline under the molecule.
    pub continuous_highlight: bool,
    /// Draws every molecule of a grid at one scale.
    pub draw_mols_same_scale: bool,
    /// The parts of a drawing its extent includes.
    pub drawing_extents_include: DrawElement,
    /// Draws breaks at dummy atoms.
    pub dummies_are_attachments: bool,
    /// Adds isotope labels to dummy atoms.
    pub dummy_isotope_labels: bool,
    /// Draws terminal methyls explicitly.
    pub explicit_methyl: bool,
    /// Fills the areas that highlight atoms and atom regions.
    pub fill_highlights: bool,
    /// When positive, fixes the bond length, in pixels, unless the drawing
    /// would not fit.
    pub fixed_bond_length: f64,
    /// When not -1, the font size in pixels, whatever the scale.
    pub fixed_font_size: i32,
    /// When positive, fixes the scale to this fraction of the drawing's
    /// width, unless the drawing would not fit.
    pub fixed_scale: f64,
    /// When positive, the distance in pixels under which close contacts are
    /// highlighted.
    pub flag_close_contacts_dist: i32,
    /// The font file for FreeType, or `BuiltinTelexRegular` (the default
    /// when empty) or `BuiltinRobotoRegular`.
    pub font_file: String,
    /// What the default bond line width is multiplied by for highlights.
    pub highlight_bond_width_multiplier: i32,
    /// The default highlight color.
    pub highlight_colour: DrawColour,
    /// The default radius of highlight circles, in molecule coordinates.
    pub highlight_radius: f64,
    /// Includes atom tags in the output.
    pub include_atom_tags: bool,
    /// Adds an "ABS" annotation when the chiral flag is set.
    pub include_chiral_flag_label: bool,
    /// Includes metadata that can rebuild the molecules and reactions drawn.
    /// A Cairo drawer's metadata holds their pickles, so it then refuses what
    /// `Mol::to_binary` and `ChemicalReaction::to_binary` refuse.
    pub include_metadata: bool,
    /// Includes radicals in the drawing.
    pub include_radicals: bool,
    /// Adds isotope labels to atoms that are not dummies.
    pub isotope_labels: bool,
    /// The color of legends.
    pub legend_colour: DrawColour,
    /// The legend's font size, in pixels.
    pub legend_font_size: i32,
    /// The fraction of a panel the legend takes.
    pub legend_fraction: f64,
    /// Where the legend goes.
    pub legend_position: LegendPosition,
    /// Draws a legend on the left or right one character per line.
    pub legend_vertical_text: bool,
    /// The largest font size, in pixels; -1 for no limit.
    pub max_font_size: i32,
    /// The smallest font size, in pixels; -1 for no limit.
    pub min_font_size: i32,
    /// How several highlight colors on one atom or bond are drawn.
    pub multi_colour_highlight_style: MultiColourHighlightStyle,
    /// The offset of the extra lines of a multiple bond, as a fraction of
    /// the mean bond length.
    pub multiple_bond_offset: f64,
    /// Leaves out atom labels.
    pub no_atom_labels: bool,
    /// The fraction of empty space around the molecule.
    pub padding: f64,
    /// Prepares molecules (kekulizing, wedging, adding chiral hydrogens)
    /// before drawing them.
    pub prepare_mols_before_drawing: bool,
    /// The color of query bonds.
    pub query_colour: DrawColour,
    /// The fraction of empty space around each component of a reaction
    /// (`componentPadding` in C++).
    pub reagent_padding: f64,
    /// Rotates molecules about their centers by this many degrees.
    pub rotate: f64,
    /// Scales bond widths with the drawing.
    pub scale_bond_width: bool,
    /// Scales highlighted bond widths with the drawing.
    pub scale_highlight_bond_width: bool,
    /// Pixels per unit of molecule coordinates when the canvas is sized to
    /// the drawing.
    pub scaling_factor: f64,
    /// Shows every CIP code, including those of OR and AND stereo groups.
    pub show_all_cip_codes: bool,
    /// Shows one molecule-level label when every stereocenter is in one
    /// stereo group.
    pub simplified_stereo_group_label: bool,
    /// Draws every bond in the symbol color.
    pub single_colour_bonds: bool,
    /// Draws wedged and dashed bonds in the symbol color.
    pub single_colour_wedge_bonds: bool,
    /// Splits bonds into a segment per atom.
    pub split_bonds: bool,
    /// Draws highlighted heteroatoms in their usual colors.
    pub standard_colours_for_highlighted_atoms: bool,
    /// The label of ABS stereo groups.
    pub stereo_group_abs_label: String,
    /// The label of AND stereo groups.
    pub stereo_group_and_label: String,
    /// The label of OR stereo groups.
    pub stereo_group_or_label: String,
    /// The color of symbols and arrows in reactions.
    pub symbol_colour: DrawColour,
    /// Draws double bonds with unspecified stereo crossed and possible
    /// stereocenters with unspecified stereo with wavy bonds.
    pub unspecified_stereo_is_unknown: bool,
    /// Draws any-atom, any-heteroatom, and any-halogen queries as A, Q,
    /// and X.
    pub use_complex_query_atom_symbols: bool,
    /// Prefers the wedging of a molecule read from a mol block.
    pub use_mol_block_wedging: bool,
    /// The radius of atoms in variable attachment points.
    pub variable_atom_radius: f64,
    /// The color of variable attachment points.
    pub variable_attachment_colour: DrawColour,
    /// What the bond width is multiplied by for variable attachment points.
    pub variable_bond_width_multiplier: i32,
    /// The colors of atoms and bonds by atomic number, -1 for the default
    /// (`atomColourPalette` in C++; Python reads and writes it with
    /// `getAtomPalette` and `setAtomPalette`).
    pub atom_colour_palette: BTreeMap<i32, DrawColour>,
}

/// The defaults of a new C++ `MolDrawOptions` (A7).
impl Default for MolDrawOptions {
    fn default() -> MolDrawOptions {
        MolDrawOptions::from_ffi(expect_ok(ffi::mol_draw_options_default()))
            .expect("RDKit's default drawing options are UTF-8")
    }
}

fn text(bytes: Vec<u8>, what: &str) -> Result<String> {
    String::from_utf8(bytes).map_err(|err| Error::Value {
        message: format!("invalid UTF-8 in {what}: {err}"),
    })
}

macro_rules! colour_accessors {
    ($($field:ident, $get:ident, $set:ident, $python:literal;)*) => {
        impl MolDrawOptions {
            $(
                #[doc = concat!("The `", stringify!($field), "` field. Mirrors `MolDrawOptions.get", $python, "`.")]
                pub fn $get(&self) -> DrawColour {
                    self.$field
                }

                #[doc = concat!("Sets the `", stringify!($field), "` field. Mirrors `MolDrawOptions.set", $python, "`; `Error::Value` unless every component is from 0 to 1.")]
                pub fn $set(&mut self, tpl: DrawColour) -> Result<()> {
                    tpl.check()?;
                    self.$field = tpl;
                    Ok(())
                }
            )*
        }
    };
}

colour_accessors! {
    annotation_colour, get_annotation_colour, set_annotation_colour, "AnnotationColour";
    atom_note_colour, get_atom_note_colour, set_atom_note_colour, "AtomNoteColour";
    background_colour, get_background_colour, set_background_colour, "BackgroundColour";
    bond_note_colour, get_bond_note_colour, set_bond_note_colour, "BondNoteColour";
    highlight_colour, get_highlight_colour, set_highlight_colour, "HighlightColour";
    legend_colour, get_legend_colour, set_legend_colour, "LegendColour";
    query_colour, get_query_colour, set_query_colour, "QueryColour";
    symbol_colour, get_symbol_colour, set_symbol_colour, "SymbolColour";
    variable_attachment_colour, get_variable_attachment_colour, set_variable_attachment_colour, "VariableAttachmentColour";
}

impl MolDrawOptions {
    /// Uses RDKit's default atom palette. Mirrors
    /// `MolDrawOptions.useDefaultAtomPalette`.
    pub fn use_default_atom_palette(&mut self) {
        self.preset(0);
    }

    /// Uses a black and white atom palette. Mirrors
    /// `MolDrawOptions.useBWAtomPalette`.
    pub fn use_bw_atom_palette(&mut self) {
        self.preset(1);
    }

    /// Uses the Avalon renderer's atom palette. Mirrors
    /// `MolDrawOptions.useAvalonAtomPalette`.
    pub fn use_avalon_atom_palette(&mut self) {
        self.preset(2);
    }

    /// Uses (part of) CDK's atom palette. Mirrors
    /// `MolDrawOptions.useCDKAtomPalette`.
    pub fn use_cdk_atom_palette(&mut self) {
        self.preset(3);
    }

    /// The atom palette, by atomic number (-1 for the default). Mirrors
    /// `MolDrawOptions.getAtomPalette`.
    pub fn get_atom_palette(&self) -> BTreeMap<i32, DrawColour> {
        self.atom_colour_palette.clone()
    }

    /// Replaces the atom palette. Mirrors `MolDrawOptions.setAtomPalette`.
    pub fn set_atom_palette(&mut self, cmap: &BTreeMap<i32, DrawColour>) {
        self.atom_colour_palette = cmap.clone();
    }

    /// Adds or replaces entries of the atom palette. Mirrors
    /// `MolDrawOptions.updateAtomPalette`.
    pub fn update_atom_palette(&mut self, cmap: &BTreeMap<i32, DrawColour>) {
        self.atom_colour_palette
            .extend(cmap.iter().map(|(k, v)| (*k, *v)));
    }

    fn preset(&mut self, preset: u8) {
        *self = self.applied(|data| ffi::mol_draw_options_preset(data, preset));
    }

    /// The options after a C++ operation that cannot fail on options that
    /// came from Rust.
    fn applied(
        &self,
        op: impl FnOnce(
            &ffi::MolDrawOptionsData,
        ) -> std::result::Result<ffi::MolDrawOptionsData, cxx::Exception>,
    ) -> MolDrawOptions {
        MolDrawOptions::from_ffi(expect_ok(op(&self.to_ffi())))
            .expect("RDKit kept the drawing options' strings UTF-8")
    }

    pub(crate) fn from_ffi(data: ffi::MolDrawOptionsData) -> Result<MolDrawOptions> {
        let atom_labels = data
            .atom_label_keys
            .into_iter()
            .zip(data.atom_label_values)
            .map(|(k, v)| Ok((k, text(v.bytes, "atomLabels")?)))
            .collect::<Result<_>>()?;
        Ok(MolDrawOptions {
            add_atom_indices: data.add_atom_indices,
            add_bond_indices: data.add_bond_indices,
            add_stereo_annotation: data.add_stereo_annotation,
            add_stereo_group_annotation: data.add_stereo_group_annotation,
            additional_atom_label_padding: data.additional_atom_label_padding,
            annotation_colour: data.annotation_colour,
            annotation_font_scale: data.annotation_font_scale,
            atom_highlights_are_circles: data.atom_highlights_are_circles,
            atom_label_deuterium_tritium: data.atom_label_deuterium_tritium,
            atom_labels,
            atom_note_colour: data.atom_note_colour,
            atom_regions: data.atom_regions.into_iter().map(|r| r.indices).collect(),
            background_colour: data.background_colour,
            base_font_size: data.base_font_size,
            bond_line_width: data.bond_line_width,
            bond_note_colour: data.bond_note_colour,
            brackets_around_atom_lists: data.brackets_around_atom_lists,
            centre_molecules_before_drawing: data.centre_molecules_before_drawing,
            circle_atoms: data.circle_atoms,
            clear_background: data.clear_background,
            comic_mode: data.comic_mode,
            continuous_highlight: data.continuous_highlight,
            draw_mols_same_scale: data.draw_mols_same_scale,
            drawing_extents_include: DrawElement(data.drawing_extents_include),
            dummies_are_attachments: data.dummies_are_attachments,
            dummy_isotope_labels: data.dummy_isotope_labels,
            explicit_methyl: data.explicit_methyl,
            fill_highlights: data.fill_highlights,
            fixed_bond_length: data.fixed_bond_length,
            fixed_font_size: data.fixed_font_size,
            fixed_scale: data.fixed_scale,
            flag_close_contacts_dist: data.flag_close_contacts_dist,
            font_file: text(data.font_file, "fontFile")?,
            highlight_bond_width_multiplier: data.highlight_bond_width_multiplier,
            highlight_colour: data.highlight_colour,
            highlight_radius: data.highlight_radius,
            include_atom_tags: data.include_atom_tags,
            include_chiral_flag_label: data.include_chiral_flag_label,
            include_metadata: data.include_metadata,
            include_radicals: data.include_radicals,
            isotope_labels: data.isotope_labels,
            legend_colour: data.legend_colour,
            legend_font_size: data.legend_font_size,
            legend_fraction: data.legend_fraction,
            legend_position: LegendPosition::from_ffi(data.legend_position),
            legend_vertical_text: data.legend_vertical_text,
            max_font_size: data.max_font_size,
            min_font_size: data.min_font_size,
            multi_colour_highlight_style: MultiColourHighlightStyle::from_ffi(
                data.multi_colour_highlight_style,
            ),
            multiple_bond_offset: data.multiple_bond_offset,
            no_atom_labels: data.no_atom_labels,
            padding: data.padding,
            prepare_mols_before_drawing: data.prepare_mols_before_drawing,
            query_colour: data.query_colour,
            reagent_padding: data.reagent_padding,
            rotate: data.rotate,
            scale_bond_width: data.scale_bond_width,
            scale_highlight_bond_width: data.scale_highlight_bond_width,
            scaling_factor: data.scaling_factor,
            show_all_cip_codes: data.show_all_cip_codes,
            simplified_stereo_group_label: data.simplified_stereo_group_label,
            single_colour_bonds: data.single_colour_bonds,
            single_colour_wedge_bonds: data.single_colour_wedge_bonds,
            split_bonds: data.split_bonds,
            standard_colours_for_highlighted_atoms: data.standard_colours_for_highlighted_atoms,
            stereo_group_abs_label: text(data.stereo_group_abs_label, "stereoGroupAbsLabel")?,
            stereo_group_and_label: text(data.stereo_group_and_label, "stereoGroupAndLabel")?,
            stereo_group_or_label: text(data.stereo_group_or_label, "stereoGroupOrLabel")?,
            symbol_colour: data.symbol_colour,
            unspecified_stereo_is_unknown: data.unspecified_stereo_is_unknown,
            use_complex_query_atom_symbols: data.use_complex_query_atom_symbols,
            use_mol_block_wedging: data.use_mol_block_wedging,
            variable_atom_radius: data.variable_atom_radius,
            variable_attachment_colour: data.variable_attachment_colour,
            variable_bond_width_multiplier: data.variable_bond_width_multiplier,
            atom_colour_palette: data
                .atom_colour_palette_keys
                .into_iter()
                .zip(data.atom_colour_palette_values)
                .collect(),
        })
    }

    pub(crate) fn to_ffi(&self) -> ffi::MolDrawOptionsData {
        ffi::MolDrawOptionsData {
            atom_label_deuterium_tritium: self.atom_label_deuterium_tritium,
            dummies_are_attachments: self.dummies_are_attachments,
            circle_atoms: self.circle_atoms,
            split_bonds: self.split_bonds,
            highlight_colour: self.highlight_colour,
            continuous_highlight: self.continuous_highlight,
            fill_highlights: self.fill_highlights,
            highlight_radius: self.highlight_radius,
            flag_close_contacts_dist: self.flag_close_contacts_dist,
            include_atom_tags: self.include_atom_tags,
            clear_background: self.clear_background,
            background_colour: self.background_colour,
            query_colour: self.query_colour,
            legend_font_size: self.legend_font_size,
            legend_fraction: self.legend_fraction,
            max_font_size: self.max_font_size,
            min_font_size: self.min_font_size,
            fixed_font_size: self.fixed_font_size,
            annotation_font_scale: self.annotation_font_scale,
            font_file: self.font_file.as_bytes().to_vec(),
            legend_colour: self.legend_colour,
            legend_position: self.legend_position.to_ffi(),
            legend_vertical_text: self.legend_vertical_text,
            multiple_bond_offset: self.multiple_bond_offset,
            padding: self.padding,
            reagent_padding: self.reagent_padding,
            additional_atom_label_padding: self.additional_atom_label_padding,
            atom_label_keys: self.atom_labels.keys().copied().collect(),
            atom_label_values: self
                .atom_labels
                .values()
                .map(|v| ffi::DrawText {
                    bytes: v.as_bytes().to_vec(),
                })
                .collect(),
            no_atom_labels: self.no_atom_labels,
            atom_regions: self
                .atom_regions
                .iter()
                .map(|r| IndexList { indices: r.clone() })
                .collect(),
            symbol_colour: self.symbol_colour,
            annotation_colour: self.annotation_colour,
            atom_note_colour: self.atom_note_colour,
            bond_note_colour: self.bond_note_colour,
            bond_line_width: self.bond_line_width,
            scale_bond_width: self.scale_bond_width,
            scale_highlight_bond_width: self.scale_highlight_bond_width,
            highlight_bond_width_multiplier: self.highlight_bond_width_multiplier,
            prepare_mols_before_drawing: self.prepare_mols_before_drawing,
            atom_colour_palette_keys: self.atom_colour_palette.keys().copied().collect(),
            atom_colour_palette_values: self.atom_colour_palette.values().copied().collect(),
            fixed_scale: self.fixed_scale,
            fixed_bond_length: self.fixed_bond_length,
            rotate: self.rotate,
            add_atom_indices: self.add_atom_indices,
            add_bond_indices: self.add_bond_indices,
            isotope_labels: self.isotope_labels,
            dummy_isotope_labels: self.dummy_isotope_labels,
            add_stereo_annotation: self.add_stereo_annotation,
            show_all_cip_codes: self.show_all_cip_codes,
            atom_highlights_are_circles: self.atom_highlights_are_circles,
            multi_colour_highlight_style: self.multi_colour_highlight_style.to_ffi(),
            centre_molecules_before_drawing: self.centre_molecules_before_drawing,
            drawing_extents_include: self.drawing_extents_include.bits(),
            explicit_methyl: self.explicit_methyl,
            include_radicals: self.include_radicals,
            include_metadata: self.include_metadata,
            comic_mode: self.comic_mode,
            variable_bond_width_multiplier: self.variable_bond_width_multiplier,
            variable_atom_radius: self.variable_atom_radius,
            variable_attachment_colour: self.variable_attachment_colour,
            include_chiral_flag_label: self.include_chiral_flag_label,
            simplified_stereo_group_label: self.simplified_stereo_group_label,
            stereo_group_and_label: self.stereo_group_and_label.as_bytes().to_vec(),
            stereo_group_or_label: self.stereo_group_or_label.as_bytes().to_vec(),
            stereo_group_abs_label: self.stereo_group_abs_label.as_bytes().to_vec(),
            add_stereo_group_annotation: self.add_stereo_group_annotation,
            unspecified_stereo_is_unknown: self.unspecified_stereo_is_unknown,
            single_colour_wedge_bonds: self.single_colour_wedge_bonds,
            single_colour_bonds: self.single_colour_bonds,
            use_mol_block_wedging: self.use_mol_block_wedging,
            scaling_factor: self.scaling_factor,
            base_font_size: self.base_font_size,
            draw_mols_same_scale: self.draw_mols_same_scale,
            use_complex_query_atom_symbols: self.use_complex_query_atom_symbols,
            brackets_around_atom_lists: self.brackets_around_atom_lists,
            standard_colours_for_highlighted_atoms: self.standard_colours_for_highlighted_atoms,
        }
    }
}

/// Updates options from JSON, whose keys are the options' C++ names
/// (`componentPadding`, `atomColourPalette`, ...) plus `colourScheme`
/// (`"dark"` or `"monochrome"`). Mirrors
/// `rdMolDraw2D.UpdateMolDrawOptionsFromJSON`; keys it does not know, and
/// values of the wrong type, are ignored, as in RDKit.
///
/// `Error::Runtime` for text that is not JSON, and `Error::Value` for a
/// color with fewer than three components (RDKit reads past the end of
/// it), atom labels whose keys are negative or not integers, or text that
/// is not UTF-8. The options are unchanged on error.
pub fn update_mol_draw_options_from_json(opts: &mut MolDrawOptions, json: &str) -> Result<()> {
    let data = ffi::update_mol_draw_options_from_json(&opts.to_ffi(), json)?;
    *opts = MolDrawOptions::from_ffi(data)?;
    Ok(())
}

/// Sets options for drawings close to the ACS 1996 guidelines, for
/// molecules whose mean bond length is `mean_bond_length` (from
/// `mean_bond_length`). Mirrors `rdMolDraw2D.SetACS1996Mode`.
///
/// It sets thin bonds, a scaling factor of 14.4 pixels per mean bond
/// length, monochrome colors, a fixed font size of 10, and, unless a font
/// file is set, `$RDBASE/Data/Fonts/FreeSans.ttf` when that exists, else
/// `BuiltinRobotoRegular`. `Error::Value` for a mean bond length that is
/// not positive, or a font path that is not UTF-8.
pub fn set_acs1996_mode(draw_options: &mut MolDrawOptions, mean_bond_length: f64) -> Result<()> {
    let data = ffi::mol_draw_options_acs1996(&draw_options.to_ffi(), mean_bond_length)?;
    *draw_options = MolDrawOptions::from_ffi(data)?;
    Ok(())
}

/// Sets dark-mode colors. Mirrors `rdMolDraw2D.SetDarkMode(MolDrawOptions)`.
pub fn set_dark_mode(d2d: &mut MolDrawOptions) {
    d2d.preset(4);
}

/// Draws everything in `fg_colour` on `bg_colour`. Mirrors
/// `rdMolDraw2D.SetMonochromeMode(MolDrawOptions, ...)`; `Error::Value`
/// unless every component is from 0 to 1.
pub fn set_monochrome_mode(
    options: &mut MolDrawOptions,
    fg_colour: DrawColour,
    bg_colour: DrawColour,
) -> Result<()> {
    let data = ffi::mol_draw_options_monochrome(&options.to_ffi(), &fg_colour, &bg_colour)?;
    *options = MolDrawOptions::from_ffi(data)?;
    Ok(())
}
