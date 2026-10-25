/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of what RDKit's drawing code (MolDraw2D, DrawMol, the text drawers,
// and MolDrawOptionsJSONParsers.cpp) trusts (CONVENTIONS.md U1). Each throws
// ValueErrorException.

#pragma once

#include <string>

#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace rdkit_cxx {

// The largest canvas or panel dimension, in pixels: Cairo's limit for an
// image surface, and a bound under which RDKit's pixel arithmetic, in int,
// cannot overflow.
constexpr int kMaxCanvasSize = 32767;

// Throws unless every component is finite and from 0 to 1. RDKit's SVG
// writer converts 255 times each component to int (DrawColourToSVG), which a
// value outside int's range, or NaN, makes undefined.
void check_colour(const RDKit::DrawColour &colour, const char *what);

// Throws unless every option is within the bounds RDKit's layout needs:
// finite values, and sizes, scales, and fractions under which the int
// arithmetic of DrawMol and MolDraw2D stays in range and no scale or font
// size comes out zero, negative, or infinite (draw_checks.cpp lists them).
void check_draw_options(const RDKit::MolDrawOptions &opts);

// RDKit's JSON option parser reads each colour's components by stepping an
// iterator three times without checking for the end of the node
// (get_rgba in MolDrawOptionsJSONParsers.cpp). Throws unless every colour
// option, and every atom palette entry, has at least three components.
void check_json_options(const std::string &json);

// Throws unless the bytes are UTF-8 (CONVENTIONS.md T2).
void check_utf8(const std::string &text, const char *what);

}  // namespace rdkit_cxx
