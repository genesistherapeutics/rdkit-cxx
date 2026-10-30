// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/draw.h"

#include <cstdint>
#include <string>
#include <vector>

#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/draw_checks.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/draw.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

std::string from_bytes(const rust::Vec<std::uint8_t> &bytes) {
  return std::string(bytes.begin(), bytes.end());
}

rust::Vec<std::uint8_t> to_bytes(const std::string &text) {
  rust::Vec<std::uint8_t> res;
  res.reserve(text.size());
  for (char c : text) {
    res.push_back(static_cast<std::uint8_t>(c));
  }
  return res;
}

}  // namespace

RDKit::DrawColour colour_from_data(const DrawColour &colour) {
  return RDKit::DrawColour(colour.r, colour.g, colour.b, colour.a);
}

DrawColour colour_to_data(const RDKit::DrawColour &colour) {
  return DrawColour{colour.r, colour.g, colour.b, colour.a};
}

RDKit::MolDrawOptions options_from_data(const MolDrawOptionsData &data) {
  RDKit::MolDrawOptions opts;
  opts.atomLabelDeuteriumTritium = data.atom_label_deuterium_tritium;
  opts.dummiesAreAttachments = data.dummies_are_attachments;
  opts.circleAtoms = data.circle_atoms;
  opts.splitBonds = data.split_bonds;
  opts.highlightColour = colour_from_data(data.highlight_colour);
  opts.continuousHighlight = data.continuous_highlight;
  opts.fillHighlights = data.fill_highlights;
  opts.highlightRadius = data.highlight_radius;
  opts.flagCloseContactsDist = data.flag_close_contacts_dist;
  opts.includeAtomTags = data.include_atom_tags;
  opts.clearBackground = data.clear_background;
  opts.backgroundColour = colour_from_data(data.background_colour);
  opts.queryColour = colour_from_data(data.query_colour);
  opts.legendFontSize = data.legend_font_size;
  opts.legendFraction = data.legend_fraction;
  opts.maxFontSize = data.max_font_size;
  opts.minFontSize = data.min_font_size;
  opts.fixedFontSize = data.fixed_font_size;
  opts.annotationFontScale = data.annotation_font_scale;
  opts.fontFile = from_bytes(data.font_file);
  opts.legendColour = colour_from_data(data.legend_colour);
  opts.multipleBondOffset = data.multiple_bond_offset;
  opts.padding = data.padding;
  opts.componentPadding = data.reagent_padding;
  opts.additionalAtomLabelPadding = data.additional_atom_label_padding;
  if (data.atom_label_keys.size() != data.atom_label_values.size()) {
    throw ValueErrorException("atom label keys and values differ in number");
  }
  for (std::size_t i = 0; i < data.atom_label_keys.size(); ++i) {
    // Labels are found by atom index, so a key past INT_MAX, which no atom
    // has, may become any int.
    opts.atomLabels[static_cast<int>(data.atom_label_keys[i])] =
        from_bytes(data.atom_label_values[i].bytes);
  }
  opts.noAtomLabels = data.no_atom_labels;
  for (const auto &region : data.atom_regions) {
    std::vector<int> atoms;
    for (auto idx : region.indices) {
      atoms.push_back(static_cast<int>(idx));
    }
    opts.atomRegions.push_back(std::move(atoms));
  }
  opts.symbolColour = colour_from_data(data.symbol_colour);
  opts.annotationColour = colour_from_data(data.annotation_colour);
  opts.atomNoteColour = colour_from_data(data.atom_note_colour);
  opts.bondNoteColour = colour_from_data(data.bond_note_colour);
  opts.bondLineWidth = data.bond_line_width;
  opts.scaleBondWidth = data.scale_bond_width;
  opts.scaleHighlightBondWidth = data.scale_highlight_bond_width;
  opts.highlightBondWidthMultiplier = data.highlight_bond_width_multiplier;
  opts.prepareMolsBeforeDrawing = data.prepare_mols_before_drawing;
  if (data.atom_colour_palette_keys.size() !=
      data.atom_colour_palette_values.size()) {
    throw ValueErrorException("atom palette keys and values differ in number");
  }
  opts.atomColourPalette.clear();
  for (std::size_t i = 0; i < data.atom_colour_palette_keys.size(); ++i) {
    opts.atomColourPalette[data.atom_colour_palette_keys[i]] =
        colour_from_data(data.atom_colour_palette_values[i]);
  }
  opts.fixedScale = data.fixed_scale;
  opts.fixedBondLength = data.fixed_bond_length;
  opts.rotate = data.rotate;
  opts.addAtomIndices = data.add_atom_indices;
  opts.addBondIndices = data.add_bond_indices;
  opts.isotopeLabels = data.isotope_labels;
  opts.dummyIsotopeLabels = data.dummy_isotope_labels;
  opts.addStereoAnnotation = data.add_stereo_annotation;
  opts.showAllCIPCodes = data.show_all_cip_codes;
  opts.atomHighlightsAreCircles = data.atom_highlights_are_circles;
  opts.multiColourHighlightStyle = data.multi_colour_highlight_style;
  opts.centreMoleculesBeforeDrawing = data.centre_molecules_before_drawing;
  opts.explicitMethyl = data.explicit_methyl;
  opts.includeRadicals = data.include_radicals;
  opts.includeMetadata = data.include_metadata;
  opts.comicMode = data.comic_mode;
  opts.variableBondWidthMultiplier = data.variable_bond_width_multiplier;
  opts.variableAtomRadius = data.variable_atom_radius;
  opts.variableAttachmentColour =
      colour_from_data(data.variable_attachment_colour);
  opts.includeChiralFlagLabel = data.include_chiral_flag_label;
  opts.simplifiedStereoGroupLabel = data.simplified_stereo_group_label;
  opts.unspecifiedStereoIsUnknown = data.unspecified_stereo_is_unknown;
  opts.singleColourWedgeBonds = data.single_colour_wedge_bonds;
  opts.useMolBlockWedging = data.use_mol_block_wedging;
  opts.scalingFactor = data.scaling_factor;
  opts.baseFontSize = data.base_font_size;
  opts.drawMolsSameScale = data.draw_mols_same_scale;
  opts.useComplexQueryAtomSymbols = data.use_complex_query_atom_symbols;
  opts.bracketsAroundAtomLists = data.brackets_around_atom_lists;
  opts.standardColoursForHighlightedAtoms =
      data.standard_colours_for_highlighted_atoms;
  return opts;
}

MolDrawOptionsData options_to_data(const RDKit::MolDrawOptions &opts) {
  MolDrawOptionsData data;
  data.atom_label_deuterium_tritium = opts.atomLabelDeuteriumTritium;
  data.dummies_are_attachments = opts.dummiesAreAttachments;
  data.circle_atoms = opts.circleAtoms;
  data.split_bonds = opts.splitBonds;
  data.highlight_colour = colour_to_data(opts.highlightColour);
  data.continuous_highlight = opts.continuousHighlight;
  data.fill_highlights = opts.fillHighlights;
  data.highlight_radius = opts.highlightRadius;
  data.flag_close_contacts_dist = opts.flagCloseContactsDist;
  data.include_atom_tags = opts.includeAtomTags;
  data.clear_background = opts.clearBackground;
  data.background_colour = colour_to_data(opts.backgroundColour);
  data.query_colour = colour_to_data(opts.queryColour);
  data.legend_font_size = opts.legendFontSize;
  data.legend_fraction = opts.legendFraction;
  data.max_font_size = opts.maxFontSize;
  data.min_font_size = opts.minFontSize;
  data.fixed_font_size = opts.fixedFontSize;
  data.annotation_font_scale = opts.annotationFontScale;
  data.font_file = to_bytes(opts.fontFile);
  data.legend_colour = colour_to_data(opts.legendColour);
  data.multiple_bond_offset = opts.multipleBondOffset;
  data.padding = opts.padding;
  data.reagent_padding = opts.componentPadding;
  data.additional_atom_label_padding = opts.additionalAtomLabelPadding;
  for (const auto &[idx, label] : opts.atomLabels) {
    // Only the JSON parser can give a label a negative key, which
    // update_mol_draw_options_from_json refuses.
    if (idx < 0) {
      throw ValueErrorException("atom label keys must be atom indices");
    }
    data.atom_label_keys.push_back(static_cast<std::uint32_t>(idx));
    data.atom_label_values.push_back(DrawText{to_bytes(label)});
  }
  data.no_atom_labels = opts.noAtomLabels;
  for (const auto &region : opts.atomRegions) {
    IndexList atoms;
    for (auto idx : region) {
      atoms.indices.push_back(static_cast<std::uint32_t>(idx));
    }
    data.atom_regions.push_back(std::move(atoms));
  }
  data.symbol_colour = colour_to_data(opts.symbolColour);
  data.annotation_colour = colour_to_data(opts.annotationColour);
  data.atom_note_colour = colour_to_data(opts.atomNoteColour);
  data.bond_note_colour = colour_to_data(opts.bondNoteColour);
  data.bond_line_width = opts.bondLineWidth;
  data.scale_bond_width = opts.scaleBondWidth;
  data.scale_highlight_bond_width = opts.scaleHighlightBondWidth;
  data.highlight_bond_width_multiplier = opts.highlightBondWidthMultiplier;
  data.prepare_mols_before_drawing = opts.prepareMolsBeforeDrawing;
  for (const auto &[atomic_num, colour] : opts.atomColourPalette) {
    data.atom_colour_palette_keys.push_back(atomic_num);
    data.atom_colour_palette_values.push_back(colour_to_data(colour));
  }
  data.fixed_scale = opts.fixedScale;
  data.fixed_bond_length = opts.fixedBondLength;
  data.rotate = opts.rotate;
  data.add_atom_indices = opts.addAtomIndices;
  data.add_bond_indices = opts.addBondIndices;
  data.isotope_labels = opts.isotopeLabels;
  data.dummy_isotope_labels = opts.dummyIsotopeLabels;
  data.add_stereo_annotation = opts.addStereoAnnotation;
  data.show_all_cip_codes = opts.showAllCIPCodes;
  data.atom_highlights_are_circles = opts.atomHighlightsAreCircles;
  data.multi_colour_highlight_style = opts.multiColourHighlightStyle;
  data.centre_molecules_before_drawing = opts.centreMoleculesBeforeDrawing;
  data.explicit_methyl = opts.explicitMethyl;
  data.include_radicals = opts.includeRadicals;
  data.include_metadata = opts.includeMetadata;
  data.comic_mode = opts.comicMode;
  data.variable_bond_width_multiplier = opts.variableBondWidthMultiplier;
  data.variable_atom_radius = opts.variableAtomRadius;
  data.variable_attachment_colour =
      colour_to_data(opts.variableAttachmentColour);
  data.include_chiral_flag_label = opts.includeChiralFlagLabel;
  data.simplified_stereo_group_label = opts.simplifiedStereoGroupLabel;
  data.unspecified_stereo_is_unknown = opts.unspecifiedStereoIsUnknown;
  data.single_colour_wedge_bonds = opts.singleColourWedgeBonds;
  data.use_mol_block_wedging = opts.useMolBlockWedging;
  data.scaling_factor = opts.scalingFactor;
  data.base_font_size = opts.baseFontSize;
  data.draw_mols_same_scale = opts.drawMolsSameScale;
  data.use_complex_query_atom_symbols = opts.useComplexQueryAtomSymbols;
  data.brackets_around_atom_lists = opts.bracketsAroundAtomLists;
  data.standard_colours_for_highlighted_atoms =
      opts.standardColoursForHighlightedAtoms;
  return data;
}

MolDrawOptionsData mol_draw_options_default() {
  return options_to_data(RDKit::MolDrawOptions());
}

MolDrawOptionsData mol_draw_options_preset(const MolDrawOptionsData &data,
                                           std::uint8_t preset) {
  auto opts = options_from_data(data);
  switch (preset) {
    case 0:
      RDKit::assignDefaultPalette(opts.atomColourPalette);
      break;
    case 1:
      RDKit::assignBWPalette(opts.atomColourPalette);
      break;
    case 2:
      RDKit::assignAvalonPalette(opts.atomColourPalette);
      break;
    case 3:
      RDKit::assignCDKPalette(opts.atomColourPalette);
      break;
    case 4:
      RDKit::setDarkMode(opts);
      break;
    default:
      throw ValueErrorException("unknown drawing preset");
  }
  return options_to_data(opts);
}

MolDrawOptionsData mol_draw_options_monochrome(const MolDrawOptionsData &data,
                                               const DrawColour &fg_colour,
                                               const DrawColour &bg_colour) {
  auto opts = options_from_data(data);
  auto fg = colour_from_data(fg_colour);
  auto bg = colour_from_data(bg_colour);
  // As pyTupleToDrawColour in Wrap/rdMolDraw2D.cpp.
  check_colour(fg, "fgColour");
  check_colour(bg, "bgColour");
  RDKit::setMonochromeMode(opts, fg, bg);
  return options_to_data(opts);
}

MolDrawOptionsData mol_draw_options_acs1996(const MolDrawOptionsData &data,
                                            double mean_bond_length) {
  auto opts = options_from_data(data);
  RDKit::MolDraw2DUtils::setACS1996Options(opts, mean_bond_length);
  // The font file can come from $RDBASE, which need not be UTF-8.
  check_utf8(opts.fontFile, "fontFile");
  return options_to_data(opts);
}

MolDrawOptionsData update_mol_draw_options_from_json(
    const MolDrawOptionsData &data, rust::Str json) {
  auto opts = options_from_data(data);
  std::string text(json);
  check_json_options(text);
  RDKit::MolDraw2DUtils::updateMolDrawOptionsFromJSON(opts, text);
  // Boost's parser refuses lone surrogates, so the text stays UTF-8; this
  // keeps the conversion's guarantee should a parser write them.
  check_utf8(opts.fontFile, "fontFile");
  for (const auto &[idx, label] : opts.atomLabels) {
    check_utf8(label, "atomLabels");
  }
  return options_to_data(opts);
}

void check_draw_colour(const DrawColour &colour) {
  check_colour(colour_from_data(colour), "color");
}

}  // namespace rdkit_cxx
