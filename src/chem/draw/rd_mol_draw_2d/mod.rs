// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.Draw.rdMolDraw2D`: 2D drawings of molecules and
//! reactions as SVG or PNG.
//!
//! A drawer (`MolDraw2DSVG`, `MolDraw2DCairo`) dereferences to `MolDraw2D`,
//! which has the drawing methods. Its options are a plain struct
//! (`MolDrawOptions`): `draw_options` returns a copy, and `set_draw_options`
//! replaces them, where Python edits `drawOptions()` in place.
//!
//! # What drawing refuses
//!
//! RDKit's drawing code trusts much of its input, and some of what Python
//! accepts reads memory out of bounds, divides by zero into coordinates and
//! sizes it then converts to `int`, or loops forever. The wrappers check
//! that input first and refuse it with `Error::Value` (`Error::Index` for an
//! index past the end), leaving the drawer as it was:
//!
//! - Canvas and panel sizes from 1 to 32767 pixels, or -1 (flexible: the
//!   drawing sets them), with width and height both flexible or both fixed,
//!   and flexible panels on a flexible canvas. A fixed canvas or panel must
//!   leave a pixel to draw in once padding is taken.
//! - Options within the bounds RDKit's layout needs (`set_draw_options`
//!   checks them): padding up to 0.45, font sizes up to 512 pixels, scaling
//!   factors from 0.1 to 1000, line widths up to 64, highlight radii up to
//!   16, colour components from 0 to 1, and so on.
//! - Text RDKit can lay out: legends, notes, labels, and other text of at
//!   most 10000 bytes, each line or label with a normal (not super- or
//!   subscript) character, and no subscript before the first; atom labels
//!   RDKit splits into pieces in ASCII; and with `no_freetype`, ASCII text
//!   with a printable character.
//! - Molecules with coordinates (and substance group brackets and data
//!   positions) finite and within 1e5, ring information whose rings are
//!   cycles, data substance groups with their fields, link nodes whose bonds
//!   have length, and, with attachment points, a mean bond length of at
//!   least 0.001.
//! - Fonts (`font_file`) RDKit can size: not bitmap fonts, whose glyphs it
//!   divides by zero to size.
//! - Drawings RDKit cannot scale: each is first laid out as RDKit will, and
//!   refused when nothing of it is left to size it by (bonds shorter than
//!   0.01, which RDKit does not draw, or `drawing_extents_include` leaving
//!   out what is drawn), when, on a flexible canvas, it comes out too small
//!   to scale, or when its layout has shapes RDKit cannot draw (coordinates
//!   or arc angles that are not finite, wavy lines with no room for their
//!   waves).
//!
//! A drawing call that fails inside RDKit after it has started laying out
//! the molecule leaves the drawer unable to draw (later drawing calls
//! return `Error::Value`), since RDKit may keep the half-built molecule and
//! index it later.
//!
//! # Threads
//!
//! Drawers are `Send` and not `Sync`. Drawing a molecule without conformers
//! computes its 2D coordinates first, holding the depiction lock (as
//! `rd_depictor` does), and comic mode holds the lock on random numbers,
//! since it draws its lines with `std::rand`. Comic mode also loops forever
//! on a line of no length, as in Python; the wrappers refuse molecules with
//! two bonded atoms at one position in comic mode.

mod drawer;
mod options;

pub use drawer::*;
pub use options::*;
