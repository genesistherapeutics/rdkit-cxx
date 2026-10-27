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
// defines both as the library was built (wrapper/CMakeLists.txt).
#include <cstdint>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <GraphMol/RWMol.h>
#include <RDGeneral/RDConfig.h>

#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace RDKit::MolDraw2D_detail {
class DrawMol;
}

namespace rdkit_cxx {

// Shared types, defined by the bridges (src/ffi/draw.rs, geometry.rs).
struct DrawColour;
struct MolDrawOptionsData;
struct DrawMoleculeData;
struct DrawMoleculeWithHighlightsData;
struct DrawMoleculesData;
struct DrawSize;
struct ContourParamsData;
struct Point2D;

using RDKit::RWMol;
using ChemicalReaction = RDKit::ChemicalReaction;

// Enums cxx checks against RDKit's (CONVENTIONS.md N6).
using MultiColourHighlightStyle = RDKit::MultiColourHighlightStyle;

// The access a drawer's RDKit subclass gives to its text drawer
// (draw_drawer.cpp).
class TextDrawerAccess;

// A drawer Rust owns: an RDKit MolDraw2DSVG or MolDraw2DCairo, and what the
// wrappers track of the state RDKit keeps private (draw_drawer.cpp).
struct Drawer {
  Drawer(std::unique_ptr<RDKit::MolDraw2D> drawer, TextDrawerAccess &text,
         bool freetype, bool svg);
  ~Drawer();
  Drawer(const Drawer &) = delete;
  Drawer &operator=(const Drawer &) = delete;

  std::unique_ptr<RDKit::MolDraw2D> drawer;
  TextDrawerAccess *text;
  // Whether the text drawer uses FreeType (draw_checks.h).
  bool freetype;
  // MolDraw2DSVG, or else MolDraw2DCairo.
  bool svg;
  // MolDraw2D's drawMols_.size(), activeMolIdx_, and whether it has a
  // globalDrawTrans_, which its drawing methods read without checking
  // (MolDraw2D::getDrawCoords reads drawMols_[activeMolIdx_]), as the
  // wrappers' successful calls leave them.
  int num_mols = 0;
  int active_mol = -1;
  bool global_transform = false;
  // A copy of the globalDrawTrans_ SetScale made, built as RDKit built it,
  // by which the wrappers lay out later drawings as RDKit will
  // (DrawMol::setTransformation). It refers to the drawer's options and
  // text drawer, so it is declared after the drawer, to be destroyed first.
  std::unique_ptr<RDKit::MolDraw2D_detail::DrawMol> global;
  // Set when a drawing call throws: RDKit may then have kept a half-built
  // molecule, or not advanced activeMolIdx_ past it, so the drawer draws
  // no more.
  bool failed = false;
  // The highlights of each DrawMoleculeWithHighlights call. RDKit's
  // DrawMolMCH keeps references to them (mcHighlightAtomMap_ and the
  // others), and reads them again whenever it lays the molecule out anew,
  // as later calls can (drawMolecules redraws every molecule the drawer
  // holds), so they live as long as the drawer.
  struct MultiHighlights {
    std::map<int, std::vector<RDKit::DrawColour>> atoms;
    std::map<int, std::vector<RDKit::DrawColour>> bonds;
    std::map<int, double> radii;
    std::map<int, int> multipliers;
  };
  std::list<MultiHighlights> multi_highlights;
};

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

// Drawers.
std::unique_ptr<Drawer> mol_draw_2d_svg_new(int width, int height,
                                            int panel_width, int panel_height,
                                            bool no_freetype);
int drawer_width(const Drawer &d);
int drawer_height(const Drawer &d);
double drawer_font_size(const Drawer &d);
void drawer_set_font_size(Drawer &d, double size);
double drawer_line_width(const Drawer &d);
void drawer_set_line_width(Drawer &d, double width);
bool drawer_fill_polys(const Drawer &d);
void drawer_set_fill_polys(Drawer &d, bool val);
bool drawer_flexi_mode(const Drawer &d);
void drawer_set_flexi_mode(Drawer &d, bool mode);
Point2D drawer_offset(const Drawer &d);
void drawer_set_offset(Drawer &d, int x, int y);
void drawer_set_colour(Drawer &d, const DrawColour &colour);
void drawer_clear_drawing(Drawer &d);
MolDrawOptionsData drawer_draw_options(const Drawer &d);
void drawer_set_draw_options(Drawer &d, const MolDrawOptionsData &data);
void drawer_update_params_from_json(Drawer &d, rust::Str json);
void drawer_set_dark_mode(Drawer &d);
void drawer_set_monochrome_mode(Drawer &d, const DrawColour &fg_colour,
                                const DrawColour &bg_colour);
void drawer_draw_molecule(Drawer &d, const RWMol &mol,
                          const DrawMoleculeData &data);
DrawSize drawer_get_mol_size(Drawer &d, const RWMol &mol,
                             const DrawMoleculeData &data);
void drawer_draw_molecule_with_highlights(
    Drawer &d, const RWMol &mol, const DrawMoleculeWithHighlightsData &data);
void drawer_draw_molecules(Drawer &d, BorrowedMols &mols,
                           const DrawMoleculesData &data);
void borrowed_mols_push_none(BorrowedMols &list);
void drawer_draw_reaction(Drawer &d, const ChemicalReaction &rxn,
                          bool highlight_by_reactant,
                          const rust::Vec<DrawColour> &colours,
                          bool has_colours,
                          const rust::Vec<std::int32_t> &conf_ids,
                          bool has_conf_ids);
void drawer_prepare_and_draw_molecule(Drawer &d, const RWMol &mol,
                                      const DrawMoleculeData &data,
                                      bool kekulize);
void drawer_draw_line(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                      bool raw_coords);
void drawer_draw_arrow(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                       bool as_polygon, double frac, double angle,
                       const DrawColour &colour, bool raw_coords);
void drawer_draw_triangle(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                          const Point2D &cds3, bool raw_coords);
void drawer_draw_polygon(Drawer &d, rust::Slice<const Point2D> cds,
                         bool raw_coords);
void drawer_draw_ellipse(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                         bool raw_coords);
void drawer_draw_rect(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                      bool raw_coords);
void drawer_draw_arc(Drawer &d, const Point2D &centre, double radius,
                     double angle1, double angle2, bool raw_coords);
void drawer_draw_attachment_line(Drawer &d, const Point2D &cds1,
                                 const Point2D &cds2, const DrawColour &colour,
                                 double len, std::uint32_t n_segments,
                                 bool raw_coords);
void drawer_draw_wavy_line(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                           const DrawColour &colour1, const DrawColour &colour2,
                           std::uint32_t n_segments, double vert_offset,
                           bool raw_coords);
void drawer_draw_string(Drawer &d, rust::Str text, const Point2D &pos,
                        std::int32_t align, bool has_align, bool raw_coords);
Point2D drawer_get_draw_coords(const Drawer &d, const Point2D &point);
Point2D drawer_get_atom_draw_coords(const Drawer &d, std::uint32_t atom_index);
void drawer_set_scale(Drawer &d, int width, int height, const Point2D &minv,
                      const Point2D &maxv);
void drawer_set_scale_with_mol(Drawer &d, int width, int height,
                               const Point2D &minv, const Point2D &maxv,
                               const RWMol &mol);
void drawer_draw_molecule_acs1996(Drawer &d, const RWMol &mol,
                                  const DrawMoleculeData &data);
std::unique_ptr<std::string> mol_to_acs1996_svg(const RWMol &mol,
                                                const DrawMoleculeData &data);
void svg_finish_drawing(Drawer &d);
std::unique_ptr<std::string> svg_drawing_text(const Drawer &d);
void svg_add_molecule_metadata(Drawer &d, const RWMol &mol, int conf_id);
void svg_tag_atoms(Drawer &d, const RWMol &mol, double radius,
                   const rust::Vec<rust::String> &event_keys,
                   const rust::Vec<rust::String> &event_values);
std::unique_ptr<Drawer> mol_draw_2d_cairo_new(int width, int height,
                                              int panel_width, int panel_height,
                                              bool no_freetype);
void cairo_finish_drawing(Drawer &d);
rust::Vec<std::uint8_t> cairo_drawing_text(const Drawer &d);
void cairo_write_drawing_text(const Drawer &d,
                              rust::Slice<const std::uint8_t> path);

// Module functions.
std::unique_ptr<RWMol> prepare_mol_for_drawing(const RWMol &mol, bool kekulize,
                                               bool add_chiral_hs,
                                               bool wedge_bonds,
                                               bool force_coords,
                                               bool wavy_bonds);
std::unique_ptr<std::string> mol_to_svg(const RWMol &mol, int width, int height,
                                        const rust::Vec<std::uint32_t> &atoms,
                                        bool has_atoms,
                                        std::uint32_t line_width_mult,
                                        bool include_atom_circles, int conf_id);
double mean_bond_length(const RWMol &mol, int conf_id);
ContourParamsData contour_params_default();
void contour_and_draw_grid(Drawer &d, rust::Slice<const double> data,
                           rust::Slice<const double> xcoords,
                           rust::Slice<const double> ycoords,
                           std::uint32_t n_contours,
                           rust::Slice<const double> levels,
                           const ContourParamsData &params, BorrowedMols &mols);
void contour_and_draw_gaussians(Drawer &d, rust::Slice<const Point2D> locs,
                                rust::Slice<const double> heights,
                                rust::Slice<const double> widths,
                                std::uint32_t n_contours,
                                rust::Slice<const double> levels,
                                const ContourParamsData &params,
                                BorrowedMols &mols);

// Conversions for the other drawing wrappers.
RDKit::MolDrawOptions options_from_data(const MolDrawOptionsData &data);
MolDrawOptionsData options_to_data(const RDKit::MolDrawOptions &opts);
RDKit::DrawColour colour_from_data(const DrawColour &colour);
DrawColour colour_to_data(const RDKit::DrawColour &colour);

}  // namespace rdkit_cxx
