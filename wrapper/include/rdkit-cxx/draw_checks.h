/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of what RDKit's drawing code (MolDraw2D, DrawMol, the text drawers,
// and MolDrawOptionsJSONParsers.cpp) trusts (CONVENTIONS.md U1). Each throws
// ValueErrorException, or rdkit_cxx::IndexError for an index.

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <GraphMol/ROMol.h>

namespace rdkit_cxx {

// The largest canvas or panel dimension, in pixels: Cairo's limit for an
// image surface, and a bound under which RDKit's pixel arithmetic, in int,
// cannot overflow.
constexpr int kMaxCanvasSize = 32767;

// The fewest pixels a canvas or panel of fixed size may leave for the
// drawing once its padding is taken. DrawMol sizes the drawing in whole
// pixels (int), and divides by the size it takes for the scale; a drawing
// area under a pixel rounds to zero, and the zero scale RDKit derives from
// it gives NaN coordinates that it converts to int.
constexpr int kMinDrawSize = 1;

// The longest text, in bytes, a drawing takes: a legend, a note, a label.
// DrawMol measures text in molecule coordinates and converts the drawing's
// size in pixels to int; within this length (and the font bounds of
// check_draw_options) the size stays in int's range.
constexpr std::size_t kMaxDrawTextSize = 10000;

// The largest highlight radius, in molecule coordinates, as for the
// highlightRadius option.
constexpr double kMaxHighlightRadius = 16.0;

// Throws unless every component is finite and from 0 to 1. RDKit's SVG
// writer converts 255 times each component to int (DrawColourToSVG), which a
// value outside int's range, or NaN, makes undefined.
void check_colour(const RDKit::DrawColour &colour, const char *what);

// Throws unless every option is within the bounds RDKit's layout needs:
// finite values, and sizes, scales, and fractions under which the int
// arithmetic of DrawMol and MolDraw2D stays in range and no scale or font
// size comes out zero, negative, or infinite (draw_checks.cpp lists them),
// and its text as check_annotation and check_label check it.
void check_draw_options(const RDKit::MolDrawOptions &opts);

// RDKit's JSON option parser reads each color's components by stepping an
// iterator three times without checking for the end of the node
// (get_rgba in MolDrawOptionsJSONParsers.cpp). Throws unless every color
// option, and every atom palette entry, has at least three components.
void check_json_options(const std::string &json);

// Throws unless the bytes are UTF-8 (CONVENTIONS.md T2).
void check_utf8(const std::string &text, const char *what);

// A drawer's canvas and panel sizes, as MolDraw2D's constructors take them:
// each -1 (flexible: the drawing sets it) or from 1 to kMaxCanvasSize, with
// the width and height both flexible or both fixed (DrawMol sizes one
// flexible dimension from the other in int, which a drawing much wider than
// it is tall overflows), and panels flexible on a flexible canvas.
void check_canvas_size(int width, int height, int panel_width,
                       int panel_height);

// Throws unless each fixed dimension leaves kMinDrawSize pixels once
// padding is taken from both sides.
void check_draw_area(int width, int height, const RDKit::MolDrawOptions &opts);

// What RDKit's text drawers trust about text. DrawText lays text out as
// runs of characters, switching to subscripts and superscripts at <sub>,
// <sup>, </sub>, and </sup>, and positions a subscript from the last normal
// character before it and a superscript from the last before it or the
// first after it (adjustStringRectsForSuperSubScript), and aligns a run from
// its first character (alignString): a run with no normal character, or a
// subscript before its first one, is read at index -1. The text drawer that
// does not use FreeType sizes each character from a table indexed by char
// (char_widths in MolDraw2DDetails.h), which a byte above 0x7f indexes before
// its start, and divides by the widest character's width, which a run of
// control characters makes zero. `freetype` says whether the drawer's text
// drawer uses FreeType.

// Text drawn as one run (a note, a CIP label, a legend line, a string drawn
// with an alignment): non-empty, at most kMaxDrawTextSize bytes, with a
// normal character and no subscript before the first one, and without
// FreeType, ASCII with a printable character.
void check_annotation(const std::string &text, bool freetype, const char *what);

// An atom label as DrawMol::getAtomSymbol builds it, which DrawText splits
// into pieces (atomLabelToPieces) and draws as runs that depend on the
// label's orientation, any of which DrawMol may choose: each piece, and the
// pieces joined forwards and backwards, must pass check_annotation. Unless
// the label starts with <lit>, atomLabelToPieces passes its bytes to
// isupper, which a byte above 0x7f (a negative char) makes undefined, so it
// must be ASCII. AtomSymbol::adjustColons strips <...> markup from a label
// with a colon in a loop that never ends when a '>' comes before the next
// '<'; such a label is refused too.
void check_label(const std::string &label, bool freetype, const char *what);

// A legend as DrawMol draws it: split at newlines into lines, each checked
// as check_annotation checks it, at least one of them not empty.
void check_legend(const std::string &legend, bool freetype);

// Everything DrawMol reads from a molecule it draws with these options,
// beyond what RDKit checks itself: its text (labels, notes, CIP labels,
// substance group text, stereo group labels), as check_annotation and
// check_label check it; its coordinates, and those of its substance groups'
// brackets and data fields, finite and within the depiction bound
// (kMaxDepictionCoordinate in depiction_checks.h); its ring information,
// whose rings DrawMol walks bond by bond (DrawMol::bondInsideRing reads the
// coordinates of atom -1 for a ring that is not a cycle); the atom regions
// of the options; and the link-node labels DrawMol places by a direction it
// normalizes. `conf_id` is the conformer to be drawn; a molecule without
// conformers is checked as RDKit's depiction will lay it out, which the
// caller checks with check_compute_2d_coords_mol. `annotations` says whether
// the drawer draws annotations (MolDraw2D::supportsAnnotations), and
// `symmetrized` whether RDKit symmetrizes the molecule's rings before it
// draws (MolDraw2D::drawMolecule does; drawMoleculeWithHighlights does not,
// and draws whatever rings the molecule has).
void check_drawn_mol(const RDKit::ROMol &mol, int conf_id,
                     const RDKit::MolDrawOptions &opts, bool freetype,
                     bool annotations, bool symmetrized);

// Highlight indices (RDKit keeps them as int and draws them by index, as
// Python's pythonObjectToVect checks), colors, and radii.
void check_highlight_atoms(const RDKit::ROMol &mol,
                           const std::vector<int> &atoms);
void check_highlight_bonds(const RDKit::ROMol &mol,
                           const std::vector<int> &bonds);
void check_highlight_colours(const std::map<int, RDKit::DrawColour> &colours);
void check_highlight_radii(const std::map<int, double> &radii);

}  // namespace rdkit_cxx
