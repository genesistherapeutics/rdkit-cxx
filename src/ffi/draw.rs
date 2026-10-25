// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A colour as red, green, blue, and alpha components, each from 0 to 1.
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

    /// `RDKit::MolDrawOptions::LegendPosition`; cxx checks the values.
    #[repr(i32)]
    enum LegendPosition {
        Bottom,
        Top,
        Left,
        Right,
    }

    /// `RDKit::MultiColourHighlightStyle`; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum MultiColourHighlightStyle {
        CIRCLEANDLINE,
        LASSO,
    }

    /// `RDKit::DrawElement`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum DrawElement {
        NONE = 0,
        PRESHAPES = 1,
        BONDS = 2,
        ATOMLABELS = 4,
        HIGHLIGHTS = 8,
        ANNOTATIONS = 16,
        RADICALS = 32,
        POSTSHAPES = 64,
        ALL = 0x7fffffff,
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
        legend_position: LegendPosition,
        legend_vertical_text: bool,
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
        show_all_cip_codes: bool,
        atom_highlights_are_circles: bool,
        multi_colour_highlight_style: MultiColourHighlightStyle,
        centre_molecules_before_drawing: bool,
        drawing_extents_include: u32,
        explicit_methyl: bool,
        include_radicals: bool,
        include_metadata: bool,
        comic_mode: bool,
        variable_bond_width_multiplier: i32,
        variable_atom_radius: f64,
        variable_attachment_colour: DrawColour,
        include_chiral_flag_label: bool,
        simplified_stereo_group_label: bool,
        stereo_group_and_label: Vec<u8>,
        stereo_group_or_label: Vec<u8>,
        stereo_group_abs_label: Vec<u8>,
        add_stereo_group_annotation: bool,
        unspecified_stereo_is_unknown: bool,
        single_colour_wedge_bonds: bool,
        single_colour_bonds: bool,
        use_mol_block_wedging: bool,
        scaling_factor: f64,
        base_font_size: f64,
        draw_mols_same_scale: bool,
        use_complex_query_atom_symbols: bool,
        brackets_around_atom_lists: bool,
        standard_colours_for_highlighted_atoms: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/draw.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type LegendPosition;
        type MultiColourHighlightStyle;
        type DrawElement;

        fn mol_draw_options_default() -> Result<MolDrawOptionsData>;
        /// Applies one of RDKit's colour presets: 0 to 3 the default, black
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
    }

    impl Vec<DrawText> {}
    impl Vec<DrawColour> {}
}
