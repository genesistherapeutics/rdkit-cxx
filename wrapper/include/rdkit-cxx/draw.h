/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDraw2D: molecule drawing. What RDKit's drawing code trusts, and how
// these wrappers check it, is in draw_checks.h.

#pragma once

// MolDraw2D's headers do not change with RDK_BUILD_CAIRO_SUPPORT or
// RDK_BUILD_FREETYPE_SUPPORT (but for a Windows-only member), and the build
// defines both as the library was built (wrapper/CMakeLists.txt); RDConfig.h
// comes first so every source sees the same configuration.
#include <cstdint>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <RDGeneral/RDConfig.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared types, defined by the bridge (src/ffi/draw.rs).
struct DrawColour;
struct MolDrawOptionsData;

// Enums cxx checks against RDKit's (CONVENTIONS.md N6). DrawElement is a
// plain enum unless USE_BETTER_ENUMS is defined, which the wrappers never do.
using LegendPosition = RDKit::MolDrawOptions::LegendPosition;
using MultiColourHighlightStyle = RDKit::MultiColourHighlightStyle;
using DrawElement = RDKit::DrawElement;

MolDrawOptionsData mol_draw_options_default();
MolDrawOptionsData mol_draw_options_preset(const MolDrawOptionsData &data,
                                           std::uint8_t preset);
MolDrawOptionsData mol_draw_options_monochrome(const MolDrawOptionsData &data,
                                               const DrawColour &fg_colour,
                                               const DrawColour &bg_colour);
MolDrawOptionsData mol_draw_options_acs1996(const MolDrawOptionsData &data,
                                            double mean_bond_length);
MolDrawOptionsData update_mol_draw_options_from_json(
    const MolDrawOptionsData &data, rust::Str json);
void check_draw_colour(const DrawColour &colour);

// Conversions for the other drawing wrappers.
RDKit::MolDrawOptions options_from_data(const MolDrawOptionsData &data);
MolDrawOptionsData options_to_data(const RDKit::MolDrawOptions &opts);
RDKit::DrawColour colour_from_data(const DrawColour &colour);
DrawColour colour_to_data(const RDKit::DrawColour &colour);

}  // namespace rdkit_cxx
