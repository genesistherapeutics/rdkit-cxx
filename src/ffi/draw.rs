// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A color as red, green, blue, and alpha components, each from 0 to 1.
    /// Mirrors the `(r, g, b, a)` tuples of `rdMolDraw2D` (RDKit's
    /// `DrawColour`).
    #[derive(Clone, Copy, Debug, PartialEq)]
    struct DrawColour {
        /// Red.
        pub r: f64,
        /// Green.
        pub g: f64,
        /// Blue.
        pub b: f64,
        /// Alpha (opacity).
        pub a: f64,
    }

    /// `RDKit::MultiColourHighlightStyle`; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum MultiColourHighlightStyle {
        CIRCLEANDLINE,
        LASSO,
    }

    /// Bytes RDKit holds as a `std::string`, checked as UTF-8 in Rust (T2).
    struct DrawText {
        bytes: Vec<u8>,
    }

    /// `RDKit::MolDrawOptions` (A7): maps as parallel keys and values, and
    /// strings as bytes. The C++ names are in `draw.cpp`.
    struct MolDrawOptionsData {
        atom_label_deuterium_tritium: bool,
        dummies_are_attachments: bool,
        circle_atoms: bool,
        split_bonds: bool,
        highlight_colour: DrawColour,
        continuous_highlight: bool,
        fill_highlights: bool,
        highlight_radius: f64,
        flag_close_contacts_dist: i32,
        include_atom_tags: bool,
        clear_background: bool,
        background_colour: DrawColour,
        query_colour: DrawColour,
        legend_font_size: i32,
        legend_fraction: f64,
        max_font_size: i32,
        min_font_size: i32,
        fixed_font_size: i32,
        annotation_font_scale: f64,
        font_file: Vec<u8>,
        legend_colour: DrawColour,
        multiple_bond_offset: f64,
        padding: f64,
        reagent_padding: f64,
        additional_atom_label_padding: f64,
        atom_label_keys: Vec<u32>,
        atom_label_values: Vec<DrawText>,
        no_atom_labels: bool,
        atom_regions: Vec<IndexList>,
        symbol_colour: DrawColour,
        annotation_colour: DrawColour,
        atom_note_colour: DrawColour,
        bond_note_colour: DrawColour,
        bond_line_width: f64,
        scale_bond_width: bool,
        scale_highlight_bond_width: bool,
        highlight_bond_width_multiplier: i32,
        prepare_mols_before_drawing: bool,
        atom_colour_palette_keys: Vec<i32>,
        atom_colour_palette_values: Vec<DrawColour>,
        fixed_scale: f64,
        fixed_bond_length: f64,
        rotate: f64,
        add_atom_indices: bool,
        add_bond_indices: bool,
        isotope_labels: bool,
        dummy_isotope_labels: bool,
        add_stereo_annotation: bool,
        atom_highlights_are_circles: bool,
        multi_colour_highlight_style: MultiColourHighlightStyle,
        centre_molecules_before_drawing: bool,
        explicit_methyl: bool,
        include_radicals: bool,
        include_metadata: bool,
        comic_mode: bool,
        variable_bond_width_multiplier: i32,
        variable_atom_radius: f64,
        variable_attachment_colour: DrawColour,
        include_chiral_flag_label: bool,
        simplified_stereo_group_label: bool,
        unspecified_stereo_is_unknown: bool,
        single_colour_wedge_bonds: bool,
        use_mol_block_wedging: bool,
        scaling_factor: f64,
        base_font_size: f64,
        draw_mols_same_scale: bool,
        use_complex_query_atom_symbols: bool,
        brackets_around_atom_lists: bool,
        standard_colours_for_highlighted_atoms: bool,
    }

    /// A drawing call's molecule-specific arguments: Python's
    /// `highlightAtoms`, `highlightBonds`, `highlightAtomColors`,
    /// `highlightBondColors`, `highlightAtomRadii`, `confId`, and `legend`,
    /// maps as parallel keys and values.
    struct DrawMoleculeData {
        highlight_atoms: Vec<u32>,
        has_highlight_atoms: bool,
        highlight_bonds: Vec<u32>,
        has_highlight_bonds: bool,
        highlight_atom_colour_keys: Vec<u32>,
        highlight_atom_colour_values: Vec<DrawColour>,
        highlight_bond_colour_keys: Vec<u32>,
        highlight_bond_colour_values: Vec<DrawColour>,
        highlight_atom_radii_keys: Vec<u32>,
        highlight_atom_radii_values: Vec<f64>,
        conf_id: i32,
        legend: String,
    }

    /// The colors of one atom's or bond's multi-color highlight.
    struct ColourList {
        colours: Vec<DrawColour>,
    }

    /// `MolDraw2D.DrawMoleculeWithHighlights`'s arguments, maps as parallel
    /// keys and values.
    struct DrawMoleculeWithHighlightsData {
        legend: String,
        atom_keys: Vec<u32>,
        atom_colours: Vec<ColourList>,
        bond_keys: Vec<u32>,
        bond_colours: Vec<ColourList>,
        radii_keys: Vec<u32>,
        radii_values: Vec<f64>,
        multiplier_keys: Vec<u32>,
        multiplier_values: Vec<i32>,
        conf_id: i32,
    }

    /// Highlight colors by atom or bond index.
    struct IndexColourMap {
        keys: Vec<u32>,
        values: Vec<DrawColour>,
    }

    /// Highlight radii by atom index.
    struct IndexRadiusMap {
        keys: Vec<u32>,
        values: Vec<f64>,
    }

    /// `MolDraw2D.DrawMolecules`'s optional lists, one entry per molecule,
    /// each with whether it was given.
    struct DrawMoleculesData {
        has_highlight_atoms: bool,
        highlight_atoms: Vec<IndexList>,
        has_highlight_bonds: bool,
        highlight_bonds: Vec<IndexList>,
        has_atom_colours: bool,
        atom_colours: Vec<IndexColourMap>,
        has_bond_colours: bool,
        bond_colours: Vec<IndexColourMap>,
        has_radii: bool,
        radii: Vec<IndexRadiusMap>,
        has_conf_ids: bool,
        conf_ids: Vec<i32>,
        has_legends: bool,
        legends: Vec<String>,
    }

    /// A width and height in pixels.
    struct DrawSize {
        width: i32,
        height: i32,
    }

    /// `MolDraw2DUtils::ContourParams`, as `ContourParams` holds it.
    struct ContourParamsData {
        set_scale: bool,
        dash_negative: bool,
        fill_grid: bool,
        grid_resolution: f64,
        contour_width: f64,
        extra_grid_padding: f64,
        contour_colour: DrawColour,
        colour_map: Vec<DrawColour>,
        draw_as_lines: bool,
        coord_scale_for_quantization: f64,
        isoval_scale_for_quantization: f64,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/draw.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type Point2D = crate::ffi::geometry::ffi::Point2D;
        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type ChemicalReaction = crate::ffi::chem_reactions::ffi::ChemicalReaction;
        type MultiColourHighlightStyle;

        fn mol_draw_options_default() -> Result<MolDrawOptionsData>;
        /// Applies one of RDKit's color presets: 0 to 3 the default, black
        /// and white, Avalon, and CDK atom palettes, 4 dark mode.
        fn mol_draw_options_preset(
            data: &MolDrawOptionsData,
            preset: u8,
        ) -> Result<MolDrawOptionsData>;
        fn mol_draw_options_monochrome(
            data: &MolDrawOptionsData,
            fg_colour: &DrawColour,
            bg_colour: &DrawColour,
        ) -> Result<MolDrawOptionsData>;
        fn mol_draw_options_acs1996(
            data: &MolDrawOptionsData,
            mean_bond_length: f64,
        ) -> Result<MolDrawOptionsData>;
        fn update_mol_draw_options_from_json(
            data: &MolDrawOptionsData,
            json: &str,
        ) -> Result<MolDrawOptionsData>;
        fn check_draw_colour(colour: &DrawColour) -> Result<()>;

        /// An RDKit MolDraw2DSVG or MolDraw2DCairo, with what the wrappers
        /// track about it (`draw.h`).
        type Drawer;

        fn mol_draw_2d_svg_new(
            width: i32,
            height: i32,
            panel_width: i32,
            panel_height: i32,
            no_freetype: bool,
        ) -> Result<UniquePtr<Drawer>>;
        fn drawer_width(d: &Drawer) -> Result<i32>;
        fn drawer_height(d: &Drawer) -> Result<i32>;
        fn drawer_font_size(d: &Drawer) -> Result<f64>;
        fn drawer_set_font_size(d: Pin<&mut Drawer>, size: f64) -> Result<()>;
        fn drawer_line_width(d: &Drawer) -> Result<f64>;
        fn drawer_set_line_width(d: Pin<&mut Drawer>, width: f64) -> Result<()>;
        fn drawer_fill_polys(d: &Drawer) -> Result<bool>;
        fn drawer_set_fill_polys(d: Pin<&mut Drawer>, val: bool) -> Result<()>;
        fn drawer_flexi_mode(d: &Drawer) -> Result<bool>;
        fn drawer_set_flexi_mode(d: Pin<&mut Drawer>, mode: bool) -> Result<()>;
        fn drawer_offset(d: &Drawer) -> Result<Point2D>;
        fn drawer_set_offset(d: Pin<&mut Drawer>, x: i32, y: i32) -> Result<()>;
        fn drawer_set_colour(d: Pin<&mut Drawer>, colour: &DrawColour) -> Result<()>;
        fn drawer_clear_drawing(d: Pin<&mut Drawer>) -> Result<()>;
        fn drawer_draw_options(d: &Drawer) -> Result<MolDrawOptionsData>;
        fn drawer_set_draw_options(d: Pin<&mut Drawer>, data: &MolDrawOptionsData) -> Result<()>;
        fn drawer_update_params_from_json(d: Pin<&mut Drawer>, json: &str) -> Result<()>;
        fn drawer_set_dark_mode(d: Pin<&mut Drawer>) -> Result<()>;
        fn drawer_set_monochrome_mode(
            d: Pin<&mut Drawer>,
            fg_colour: &DrawColour,
            bg_colour: &DrawColour,
        ) -> Result<()>;
        fn drawer_draw_molecule(
            d: Pin<&mut Drawer>,
            mol: &RWMol,
            data: &DrawMoleculeData,
        ) -> Result<()>;
        fn drawer_get_mol_size(
            d: Pin<&mut Drawer>,
            mol: &RWMol,
            data: &DrawMoleculeData,
        ) -> Result<DrawSize>;
        fn drawer_draw_molecule_with_highlights(
            d: Pin<&mut Drawer>,
            mol: &RWMol,
            data: &DrawMoleculeWithHighlightsData,
        ) -> Result<()>;
        fn drawer_draw_molecules(
            d: Pin<&mut Drawer>,
            mols: Pin<&mut BorrowedMols>,
            data: &DrawMoleculesData,
        ) -> Result<()>;
        fn borrowed_mols_push_none(list: Pin<&mut BorrowedMols>) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn drawer_draw_reaction(
            d: Pin<&mut Drawer>,
            rxn: &ChemicalReaction,
            highlight_by_reactant: bool,
            colours: &Vec<DrawColour>,
            has_colours: bool,
            conf_ids: &Vec<i32>,
            has_conf_ids: bool,
        ) -> Result<()>;
        fn drawer_prepare_and_draw_molecule(
            d: Pin<&mut Drawer>,
            mol: &RWMol,
            data: &DrawMoleculeData,
            kekulize: bool,
        ) -> Result<()>;
        fn drawer_draw_line(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            raw_coords: bool,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn drawer_draw_arrow(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            as_polygon: bool,
            frac: f64,
            angle: f64,
            colour: &DrawColour,
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_draw_triangle(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            cds3: &Point2D,
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_draw_polygon(
            d: Pin<&mut Drawer>,
            cds: &[Point2D],
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_draw_ellipse(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_draw_rect(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_draw_arc(
            d: Pin<&mut Drawer>,
            centre: &Point2D,
            radius: f64,
            angle1: f64,
            angle2: f64,
            raw_coords: bool,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn drawer_draw_attachment_line(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            colour: &DrawColour,
            len: f64,
            n_segments: u32,
            raw_coords: bool,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn drawer_draw_wavy_line(
            d: Pin<&mut Drawer>,
            cds1: &Point2D,
            cds2: &Point2D,
            colour1: &DrawColour,
            colour2: &DrawColour,
            n_segments: u32,
            vert_offset: f64,
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_draw_string(
            d: Pin<&mut Drawer>,
            text: &str,
            pos: &Point2D,
            align: i32,
            has_align: bool,
            raw_coords: bool,
        ) -> Result<()>;
        fn drawer_get_draw_coords(d: &Drawer, point: &Point2D) -> Result<Point2D>;
        fn drawer_get_atom_draw_coords(d: &Drawer, atom_index: u32) -> Result<Point2D>;
        fn drawer_set_scale(
            d: Pin<&mut Drawer>,
            width: i32,
            height: i32,
            minv: &Point2D,
            maxv: &Point2D,
        ) -> Result<()>;
        fn drawer_set_scale_with_mol(
            d: Pin<&mut Drawer>,
            width: i32,
            height: i32,
            minv: &Point2D,
            maxv: &Point2D,
            mol: &RWMol,
        ) -> Result<()>;
        fn drawer_draw_molecule_acs1996(
            d: Pin<&mut Drawer>,
            mol: &RWMol,
            data: &DrawMoleculeData,
        ) -> Result<()>;
        fn mol_to_acs1996_svg(mol: &RWMol, data: &DrawMoleculeData)
        -> Result<UniquePtr<CxxString>>;
        fn svg_finish_drawing(d: Pin<&mut Drawer>) -> Result<()>;
        fn svg_drawing_text(d: &Drawer) -> Result<UniquePtr<CxxString>>;
        fn svg_add_molecule_metadata(d: Pin<&mut Drawer>, mol: &RWMol, conf_id: i32) -> Result<()>;
        fn svg_tag_atoms(
            d: Pin<&mut Drawer>,
            mol: &RWMol,
            radius: f64,
            event_keys: &Vec<String>,
            event_values: &Vec<String>,
        ) -> Result<()>;
        fn mol_draw_2d_cairo_new(
            width: i32,
            height: i32,
            panel_width: i32,
            panel_height: i32,
            no_freetype: bool,
        ) -> Result<UniquePtr<Drawer>>;
        fn cairo_finish_drawing(d: Pin<&mut Drawer>) -> Result<()>;
        fn cairo_drawing_text(d: &Drawer) -> Result<Vec<u8>>;
        fn cairo_write_drawing_text(d: &Drawer, path: &[u8]) -> Result<()>;

        fn prepare_mol_for_drawing(
            mol: &RWMol,
            kekulize: bool,
            add_chiral_hs: bool,
            wedge_bonds: bool,
            force_coords: bool,
            wavy_bonds: bool,
        ) -> Result<UniquePtr<RWMol>>;
        #[allow(clippy::too_many_arguments)]
        fn mol_to_svg(
            mol: &RWMol,
            width: i32,
            height: i32,
            atoms: &Vec<u32>,
            has_atoms: bool,
            line_width_mult: u32,
            include_atom_circles: bool,
            conf_id: i32,
        ) -> Result<UniquePtr<CxxString>>;
        fn mean_bond_length(mol: &RWMol, conf_id: i32) -> Result<f64>;

        fn contour_params_default() -> Result<ContourParamsData>;
        #[allow(clippy::too_many_arguments)]
        fn contour_and_draw_grid(
            d: Pin<&mut Drawer>,
            data: &[f64],
            xcoords: &[f64],
            ycoords: &[f64],
            n_contours: u32,
            levels: &[f64],
            params: &ContourParamsData,
            mols: Pin<&mut BorrowedMols>,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn contour_and_draw_gaussians(
            d: Pin<&mut Drawer>,
            locs: &[Point2D],
            heights: &[f64],
            widths: &[f64],
            n_contours: u32,
            levels: &[f64],
            params: &ContourParamsData,
            mols: Pin<&mut BorrowedMols>,
        ) -> Result<()>;
    }

    impl Vec<DrawText> {}
    impl Vec<DrawColour> {}
}

// SAFETY: a Drawer owns its RDKit drawer, which holds copies of what it drew
// (DrawMol copies each molecule and its highlights) and its own text drawer
// and output (an SVG stream, or a Cairo surface and context, which Cairo
// lets any one thread at a time use), with no thread affinity: moving it to
// another thread moves all it owns. It is not Sync: RDKit's const methods write to it
// (MolDraw2DSVG::addMoleculeMetadata appends to its output), so the Rust
// methods that call them take it mutably (CONVENTIONS.md S4).
unsafe impl Send for ffi::Drawer {}
