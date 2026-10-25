// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The drawers: MolDraw2DSVG (and MolDraw2DCairo), and the functions that
// draw with them.

#include <climits>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>
// The text drawer's and DrawMol's interfaces, which RDKit does not install;
// wrapper/vendor holds RDKit's own copies, from the release the crate
// supports.
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawText.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/draw.h"
#include "rdkit-cxx/draw_checks.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/draw.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

// RDKit's text drawer (the bridge has a DrawText of its own).
using TextDrawer = RDKit::MolDraw2D_detail::DrawText;

// MolDraw2D keeps its text drawer in a protected member; each drawer class
// derives from RDKit's to reach it.
class TextDrawerAccess {
 public:
  virtual TextDrawer &text_drawer() = 0;

 protected:
  ~TextDrawerAccess() = default;
};

namespace {

class SvgDrawer : public RDKit::MolDraw2DSVG, public TextDrawerAccess {
 public:
  SvgDrawer(int width, int height, int panel_width, int panel_height,
            bool no_freetype)
      : RDKit::MolDraw2DSVG(width, height, panel_width, panel_height,
                            no_freetype) {}
  TextDrawer &text_drawer() override { return *text_drawer_; }
};

// MolDraw2DSVG and MolDraw2DCairo fall back to RDKit's own text metrics when
// FreeType fails to start, so whether a drawer uses FreeType is read from its
// text drawer: only FreeType's keeps a font file name.
bool uses_freetype(TextDrawer &text) {
  text.setFontFile("BuiltinTelexRegular");
  const bool freetype = text.getFontFile() == "BuiltinTelexRegular";
  text.setFontFile("");
  return freetype;
}

// The largest glyph, in multiples of the font size, a font may have. RDKit
// sizes each glyph from the font's units per em, which a bitmap font has
// none of (FreeType gives 0, and DrawTextFT divides by it, making every
// coordinate NaN), and a font's glyphs set the size of every label.
constexpr double kMaxGlyphSize = 8.0;

// Loads the font as MolDraw2D::setupTextDrawer will (it then finds the font
// loaded), and measures two of each byte's glyph: their width and height
// must be finite and within kMaxGlyphSize font sizes (each).
void check_font(Drawer &d) {
  if (!d.freetype) {
    return;
  }
  auto &text = d.text->text_drawer();
  try {
    text.setFontFile(d.drawer->drawOptions().fontFile);
  } catch (const std::runtime_error &) {
    text.setFontFile("");
    // setupTextDrawer then logs the default font's path, which throws
    // without $RDBASE; throwing it here leaves the drawer usable.
    static_cast<void>(text.getFontFile());
  }
  const double scale = text.fontScale();
  text.setFontScale(1.0, true);
  const double size = text.fontSize();
  bool ok = std::isfinite(size) && size > 0;
  for (int b = 0; ok && b < 256; ++b) {
    const std::string pair(2, static_cast<char>(b));
    double x_min, y_min, x_max, y_max;
    const bool dont_split = true;
    text.getStringExtremes(pair, RDKit::MolDraw2D_detail::OrientType::E, x_min,
                           y_min, x_max, y_max, dont_split);
    ok = (x_max - x_min) <= 2 * kMaxGlyphSize * size &&
         (y_max - y_min) <= kMaxGlyphSize * size;
  }
  text.setFontScale(scale, true);
  if (!ok) {
    throw ValueErrorException(
        "the font (fontFile) has glyphs RDKit cannot size: a bitmap font, or "
        "glyphs more than 8 times the font size");
  }
}

void check_usable(const Drawer &d) {
  if (d.failed) {
    throw ValueErrorException(
        "the drawer cannot draw after a drawing call that failed partway; "
        "make a new one");
  }
}

// Runs a drawing call, marking the drawer failed if it throws.
template <class F>
void drawing_call(Drawer &d, F &&f) {
  check_usable(d);
  try {
    f();
  } catch (...) {
    d.failed = true;
    throw;
  }
}

std::vector<int> to_ints(const rust::Vec<std::uint32_t> &v) {
  return std::vector<int>(v.begin(), v.end());
}

// A molecule as DrawMol will draw it, and the locks drawing it needs.
// DrawMol gives a molecule without conformers 2D coordinates
// (initDrawMolecule), after MolDraw2D::drawMolecule has symmetrized its
// rings; this does the same on a copy, under the depiction lock, so the
// coordinates can be checked before RDKit draws, and DrawMol depicts no
// more. Depiction matches ring templates to the molecule (control_c_lock.h).
// Comic mode draws its lines with std::rand (handdrawnLine in
// MolDraw2DDetails.cpp), a global with no lock on some platforms, so the
// drawing holds random_mutex() too.
class DrawnMol {
 public:
  DrawnMol(const RWMol &mol, const RDKit::MolDrawOptions &opts) : d_mol(&mol) {
    if (!mol.getNumConformers()) {
      auto depict = lock_depiction();
      std::lock_guard<std::mutex> control_c(control_c_mutex());
      d_copy = std::make_unique<RWMol>(mol);
      MatrixCacheGuard guard(*d_copy);
      if (!d_copy->getRingInfo()->isSymmSssr()) {
        RDKit::MolOps::symmetrizeSSSR(*d_copy);
      }
      check_compute_2d_coords_mol(*d_copy, false);
      const bool canon_orient = true;
      RDDepict::compute2DCoords(*d_copy, nullptr, canon_orient);
      d_mol = d_copy.get();
    }
    if (opts.comicMode) {
      d_random = std::unique_lock<std::mutex>(random_mutex());
    }
  }
  const RWMol &mol() const { return *d_mol; }

 private:
  const RWMol *d_mol;
  std::unique_ptr<RWMol> d_copy;
  std::unique_lock<std::mutex> d_random;
};

// DrawMoleculeData as RDKit takes it.
struct Highlights {
  std::vector<int> atoms;
  bool has_atoms;
  std::vector<int> bonds;
  bool has_bonds;
  std::map<int, RDKit::DrawColour> atom_colours;
  std::map<int, RDKit::DrawColour> bond_colours;
  std::map<int, double> radii;
};

Highlights read_highlights(const RWMol &mol, const DrawMoleculeData &data) {
  Highlights h;
  h.atoms = to_ints(data.highlight_atoms);
  h.has_atoms = data.has_highlight_atoms;
  h.bonds = to_ints(data.highlight_bonds);
  h.has_bonds = data.has_highlight_bonds;
  check_highlight_atoms(mol, h.atoms);
  check_highlight_bonds(mol, h.bonds);
  for (std::size_t i = 0; i < data.highlight_atom_colour_keys.size(); ++i) {
    h.atom_colours[static_cast<int>(data.highlight_atom_colour_keys[i])] =
        colour_from_data(data.highlight_atom_colour_values[i]);
  }
  for (std::size_t i = 0; i < data.highlight_bond_colour_keys.size(); ++i) {
    h.bond_colours[static_cast<int>(data.highlight_bond_colour_keys[i])] =
        colour_from_data(data.highlight_bond_colour_values[i]);
  }
  for (std::size_t i = 0; i < data.highlight_atom_radii_keys.size(); ++i) {
    h.radii[static_cast<int>(data.highlight_atom_radii_keys[i])] =
        data.highlight_atom_radii_values[i];
  }
  check_highlight_colours(h.atom_colours);
  check_highlight_colours(h.bond_colours);
  check_highlight_radii(h.radii);
  return h;
}

// RDKit lays a drawing out from its extents, which DrawMol::findExtremes
// takes from the shapes, labels, and annotations drawingExtentsInclude
// names; DrawMol leaves out bonds shorter than 0.01, so a molecule that
// small, or extents that leave out everything drawn, have none, and
// findExtremes leaves the extremes at +-DBL_MAX / 2, from which every drawn
// coordinate is measured. A flexible canvas takes its size from the
// extents: DrawMol::calculateScale sets the width to scalingFactor times
// the extent (plus padding) in whole pixels, and the scale to the pixels
// left over the extent, so a drawing under a pixel or so gets a scale of
// zero, and RDKit then divides by it (createDrawObjects's relative font
// scale with a fixed font size, DrawMol::setScale's font scale, the wavy
// lines' segment counts). Either way it converts the NaN and infinite sizes
// and coordinates that come of it to int. So each drawing is first laid out
// on a DrawMol of the wrappers' own, as MolDraw2D will lay it out
// (createDrawObjects, repeated step by step), checking the extremes and
// each scale before RDKit divides by them and each size before RDKit
// converts it to int; the drawing that follows takes the same steps with
// the same numbers.
namespace detail = RDKit::MolDraw2D_detail;

[[noreturn]] void refuse_layout() {
  throw ValueErrorException(
      "the drawing is too small (or too large) for a flexible canvas: "
      "RDKit's scale for it comes out zero or not finite (increase "
      "scalingFactor, or draw on a canvas of fixed size)");
}

void checked_scale(detail::DrawMol &dm) {
  dm.findExtremes();
  for (double v : {dm.xMin_, dm.xMax_, dm.yMin_, dm.yMax_}) {
    if (!(std::abs(v) <= 1e12)) {
      throw ValueErrorException(
          "nothing of the drawing is left to size it by: its bonds are "
          "shorter than 0.01 (which RDKit does not draw), or "
          "drawingExtentsInclude leaves out what is drawn");
    }
  }
  // calculateScale converts scalingFactor times each extent, with padding,
  // to int.
  const double pad = 1 + 2 * dm.marginPadding_;
  for (double range : {dm.xRange_, dm.yRange_}) {
    const double size = dm.drawOptions_.scalingFactor * range * pad;
    if (!(size >= 0 && size <= INT_MAX / 4)) {
      refuse_layout();
    }
  }
  dm.calculateScale();
  if (!(dm.scale_ > 0 && std::isfinite(dm.scale_) && dm.fontScale_ > 0 &&
        std::isfinite(dm.fontScale_))) {
    refuse_layout();
  }
}

// The shapes a layout leaves, as their myDraw methods will read them:
// DrawShapeArc::myDraw has drawArc step its angles in int, and
// DrawShapeWavyLine::myDraw converts its length over its waves' width to
// int, which MolDraw2D::drawWavyLine takes as an unsigned count of at least
// two (a PRECONDITION).
void check_shapes(
    const std::vector<std::unique_ptr<detail::DrawShape>> &shapes) {
  for (const auto &shape : shapes) {
    bool ok = true;
    for (const auto &p : shape->points_) {
      ok = ok && std::abs(p.x) <= 1e9 && std::abs(p.y) <= 1e9;
    }
    if (const auto *arc =
            dynamic_cast<const detail::DrawShapeArc *>(shape.get())) {
      ok = ok && std::abs(arc->ang1_) <= 1e4 && std::abs(arc->ang2_) <= 1e4;
    }
    if (const auto *wavy =
            dynamic_cast<const detail::DrawShapeWavyLine *>(shape.get())) {
      const double segments = (wavy->points_[0] - wavy->points_[1]).length() /
                              (wavy->offset_ * 2.0 / 3.0);
      ok = ok && segments >= 1.5 && segments <= 1e6;
    }
    if (!ok) {
      throw ValueErrorException(
          "the drawing's geometry is degenerate: RDKit lays it out with "
          "shapes it cannot draw (coordinates or angles that are not finite, "
          "or wavy lines with no room for their waves)");
    }
  }
}

// MolDraw2D::createDrawObjects, with the checks.
void checked_layout(detail::DrawMol &dm) {
  auto &text = dm.textDrawer_;
  text.setFontScale(dm.fontScale_, true);
  dm.partitionForLegend();
  dm.extractAll(dm.scale_);
  checked_scale(dm);
  const bool ignore_font_limits = dm.drawOptions_.fixedFontSize != -1;
  if (!text.setFontScale(dm.fontScale_, ignore_font_limits) ||
      ignore_font_limits) {
    if (dm.flexiCanvasX_ && dm.flexiCanvasY_ &&
        (dm.fontScale_ - dm.scale_) > 1e-4) {
      dm.width_ = -1;
      dm.height_ = -1;
      const auto current = text.fontScale();
      const auto relative = dm.fontScale_ / dm.scale_;
      dm.resetEverything();
      dm.fontScale_ = relative;
      text.setFontScale(relative, true);
      dm.extractAll(dm.scale_);
      checked_scale(dm);
      text.setFontScale(current, true);
    }
    dm.setScale(dm.scale_, text.fontScale(), ignore_font_limits);
  } else {
    dm.finishCreateDrawObjects();
  }
  if (!(dm.width_ >= 1 && dm.height_ >= 1)) {
    refuse_layout();
  }
  check_shapes(dm.bonds_);
  check_shapes(dm.preShapes_);
  check_shapes(dm.postShapes_);
  check_shapes(dm.highlights_);
}

// MolDraw2D::setupTextDrawer's settings, which the drawing will make too.
void set_up_text_drawer(Drawer &d) {
  auto &text = d.text->text_drawer();
  const auto &opts = d.drawer->drawOptions();
  text.setMaxFontSize(opts.maxFontSize);
  text.setMinFontSize(opts.minFontSize);
  if (opts.baseFontSize > 0.0) {
    text.setBaseFontSize(opts.baseFontSize);
  }
}

// Lays out a molecule as MolDraw2D::drawMolecule and getMolSize will, at
// the panel size.
void check_layout(Drawer &d, const RWMol &mol, const std::string &legend,
                  const Highlights &h, const std::vector<int> *bonds,
                  int conf_id) {
  set_up_text_drawer(d);
  RDKit::ROMol copy(mol);
  if (!mol.getRingInfo()->isSymmSssr()) {
    RDKit::MolOps::symmetrizeSSSR(copy);
  }
  detail::DrawMol probe(copy, legend, d.drawer->panelWidth(),
                        d.drawer->panelHeight(), d.drawer->drawOptions(),
                        d.text->text_drawer(), h.has_atoms ? &h.atoms : nullptr,
                        bonds, &h.atom_colours, &h.bond_colours, nullptr,
                        &h.radii, d.drawer->supportsAnnotations(), conf_id);
  checked_layout(probe);
}

// What drawing a molecule on the drawer trusts, beyond the options, which
// the setters check: the drawing area, the font, the legend, and the
// molecule as DrawMol sees it.
void check_drawing(Drawer &d, const RWMol &mol, int conf_id,
                   const std::string &legend) {
  const auto &opts = d.drawer->drawOptions();
  check_draw_area(d.drawer->panelWidth(), d.drawer->panelHeight(), opts);
  check_font(d);
  check_legend(legend, opts, d.freetype);
  check_drawn_mol(mol, conf_id, opts, d.freetype,
                  d.drawer->supportsAnnotations());
}

// Without highlight bonds, Python's first DrawMolecule overload highlights
// the bonds between highlighted atoms (bonds_from_atoms); its other drawing
// functions highlight none.
void draw_molecule(Drawer &d, const RWMol &mol, const DrawMoleculeData &data,
                   bool bonds_from_atoms) {
  check_usable(d);
  const auto h = read_highlights(mol, data);
  const std::string legend(data.legend);
  DrawnMol drawn(mol, d.drawer->drawOptions());
  check_drawing(d, drawn.mol(), data.conf_id, legend);
  // The bonds RDKit highlights: as given, or, for Python's first overload,
  // those between highlighted atoms, as MolDraw2D::drawMolecule finds them.
  std::vector<int> atom_bonds;
  const std::vector<int> *bonds = h.has_bonds ? &h.bonds : nullptr;
  if (!h.has_bonds && bonds_from_atoms) {
    if (h.has_atoms) {
      detail::getBondHighlightsForAtoms(drawn.mol(), h.atoms, atom_bonds);
    }
    bonds = &atom_bonds;
  }
  check_layout(d, drawn.mol(), legend, h, bonds, data.conf_id);
  drawing_call(d, [&] {
    d.drawer->drawMolecule(
        drawn.mol(), legend, h.has_atoms ? &h.atoms : nullptr, bonds,
        &h.atom_colours, &h.bond_colours, &h.radii, data.conf_id);
  });
  ++d.num_mols;
  ++d.active_mol;
}

}  // namespace

Drawer::Drawer(std::unique_ptr<RDKit::MolDraw2D> drawer, TextDrawerAccess &text,
               bool freetype, bool svg)
    : drawer(std::move(drawer)), text(&text), freetype(freetype), svg(svg) {}

Drawer::~Drawer() = default;

std::unique_ptr<Drawer> mol_draw_2d_svg_new(int width, int height,
                                            int panel_width, int panel_height,
                                            bool no_freetype) {
  check_canvas_size(width, height, panel_width, panel_height);
  auto svg = std::make_unique<SvgDrawer>(width, height, panel_width,
                                         panel_height, no_freetype);
  TextDrawerAccess &text = *svg;
  const bool freetype = uses_freetype(text.text_drawer());
  return std::make_unique<Drawer>(std::move(svg), text, freetype, true);
}

int drawer_width(const Drawer &d) { return d.drawer->width(); }

int drawer_height(const Drawer &d) { return d.drawer->height(); }

double drawer_font_size(const Drawer &d) { return d.drawer->fontSize(); }

void drawer_set_font_size(Drawer &d, double size) {
  // DrawText::setFontSize divides by the base font size and keeps the scale
  // within the font size limits when they are set; the SVG text drawer
  // converts the font size to unsigned (DrawTextSVG::drawChar).
  if (!(size >= 1 && size <= 512)) {
    throw ValueErrorException("the font size must be from 1 to 512");
  }
  d.drawer->setFontSize(size);
}

double drawer_line_width(const Drawer &d) { return d.drawer->lineWidth(); }

void drawer_set_line_width(Drawer &d, double width) {
  // MolDraw2D::setLineWidth sets the bondLineWidth option.
  auto opts = d.drawer->drawOptions();
  opts.bondLineWidth = width;
  check_draw_options(opts);
  d.drawer->setLineWidth(width);
}

bool drawer_fill_polys(const Drawer &d) { return d.drawer->fillPolys(); }

void drawer_set_fill_polys(Drawer &d, bool val) { d.drawer->setFillPolys(val); }

bool drawer_flexi_mode(const Drawer &d) { return d.drawer->flexiMode(); }

void drawer_set_flexi_mode(Drawer &d, bool mode) {
  d.drawer->setFlexiMode(mode);
}

Point2D drawer_offset(const Drawer &d) {
  const auto p = d.drawer->offset();
  return Point2D{p.x, p.y};
}

void drawer_set_offset(Drawer &d, int x, int y) {
  // DrawMol and drawMolecules add the offsets to pixel positions in int.
  if (std::abs(static_cast<long long>(x)) > kMaxCanvasSize ||
      std::abs(static_cast<long long>(y)) > kMaxCanvasSize) {
    throw ValueErrorException("the offsets must be within 32767 pixels");
  }
  d.drawer->setOffset(x, y);
}

void drawer_set_colour(Drawer &d, const DrawColour &colour) {
  const auto c = colour_from_data(colour);
  check_colour(c, "colour");
  d.drawer->setColour(c);
}

void drawer_clear_drawing(Drawer &d) { d.drawer->clearDrawing(); }

MolDrawOptionsData drawer_draw_options(const Drawer &d) {
  return options_to_data(d.drawer->drawOptions());
}

void drawer_set_draw_options(Drawer &d, const MolDrawOptionsData &data) {
  auto opts = options_from_data(data);
  check_draw_options(opts);
  d.drawer->drawOptions() = opts;
}

void drawer_update_params_from_json(Drawer &d, rust::Str json) {
  auto opts = d.drawer->drawOptions();
  const std::string text(json);
  check_json_options(text);
  RDKit::MolDraw2DUtils::updateMolDrawOptionsFromJSON(opts, text);
  check_draw_options(opts);
  check_utf8(opts.fontFile, "fontFile");
  d.drawer->drawOptions() = opts;
}

void drawer_set_dark_mode(Drawer &d) { RDKit::setDarkMode(*d.drawer); }

void drawer_set_monochrome_mode(Drawer &d, const DrawColour &fg_colour,
                                const DrawColour &bg_colour) {
  const auto fg = colour_from_data(fg_colour);
  const auto bg = colour_from_data(bg_colour);
  check_colour(fg, "fgColour");
  check_colour(bg, "bgColour");
  RDKit::setMonochromeMode(*d.drawer, fg, bg);
}

void drawer_draw_molecule(Drawer &d, const RWMol &mol,
                          const DrawMoleculeData &data) {
  draw_molecule(d, mol, data, true);
}

DrawSize drawer_get_mol_size(Drawer &d, const RWMol &mol,
                             const DrawMoleculeData &data) {
  check_usable(d);
  const auto h = read_highlights(mol, data);
  const std::string legend(data.legend);
  DrawnMol drawn(mol, d.drawer->drawOptions());
  check_drawing(d, drawn.mol(), data.conf_id, legend);
  check_layout(d, drawn.mol(), legend, h, h.has_bonds ? &h.bonds : nullptr,
               data.conf_id);
  // GetMolSize lays the molecule out in a DrawMol of its own, and leaves the
  // drawer's molecules as they were. Python's passes the bonds it is given,
  // or none.
  const auto size = d.drawer->getMolSize(
      drawn.mol(), legend, h.has_atoms ? &h.atoms : nullptr,
      h.has_bonds ? &h.bonds : nullptr, &h.atom_colours, &h.bond_colours,
      &h.radii, data.conf_id);
  return DrawSize{size.first, size.second};
}

void drawer_prepare_and_draw_molecule(Drawer &d, const RWMol &mol,
                                      const DrawMoleculeData &data,
                                      bool kekulize) {
  // As MolDraw2DUtils::prepareAndDrawMolecule: prepares a copy, then draws
  // it with prepareMolsBeforeDrawing off.
  check_usable(d);
  auto copy = prepare_mol_for_drawing(mol, kekulize, true, true, false, false);
  auto &opts = d.drawer->drawOptions();
  const bool prepare = opts.prepareMolsBeforeDrawing;
  opts.prepareMolsBeforeDrawing = false;
  try {
    draw_molecule(d, *copy, data, false);
  } catch (...) {
    opts.prepareMolsBeforeDrawing = prepare;
    throw;
  }
  opts.prepareMolsBeforeDrawing = prepare;
}

void svg_finish_drawing(Drawer &d) {
  static_cast<RDKit::MolDraw2DSVG &>(*d.drawer).finishDrawing();
}

std::unique_ptr<std::string> svg_drawing_text(const Drawer &d) {
  return std::make_unique<std::string>(
      static_cast<const RDKit::MolDraw2DSVG &>(*d.drawer).getDrawingText());
}

void svg_add_molecule_metadata(Drawer &d, const RWMol &mol, int conf_id) {
  check_usable(d);
  // It places each atom by getDrawCoords, which reads
  // drawMols_[activeMolIdx_] (index -1 after an empty reaction), unless the
  // atom records its drawn position.
  if (mol.getNumAtoms() && !d.global_transform &&
      !(d.active_mol >= 0 && d.active_mol < d.num_mols)) {
    throw ValueErrorException(
        "AddMoleculeMetadata needs a drawn molecule to place the atoms by");
  }
  // It reads the conformer once per atom, after writing the atom's start;
  // reading it first fails before anything is written.
  if (mol.getNumAtoms()) {
    static_cast<void>(mol.getConformer(conf_id));
  }
  static_cast<const RDKit::MolDraw2DSVG &>(*d.drawer).addMoleculeMetadata(
      mol, conf_id);
}

void svg_tag_atoms(Drawer &d, const RWMol &mol, double radius,
                   const rust::Vec<rust::String> &event_keys,
                   const rust::Vec<rust::String> &event_values) {
  check_usable(d);
  std::map<std::string, std::string> events;
  for (std::size_t i = 0; i < event_keys.size(); ++i) {
    events[std::string(event_keys[i])] = std::string(event_values[i]);
  }
  // getDrawCoords(int) and scale() check activeMolIdx_ and the atom index.
  static_cast<RDKit::MolDraw2DSVG &>(*d.drawer).tagAtoms(mol, radius, events);
}

std::unique_ptr<RWMol> prepare_mol_for_drawing(const RWMol &mol, bool kekulize,
                                               bool add_chiral_hs,
                                               bool wedge_bonds,
                                               bool force_coords,
                                               bool wavy_bonds) {
  auto res = std::make_unique<RWMol>(mol);
  MatrixCacheGuard guard(*res);
  // It computes 2D coordinates when forced or when the molecule has none,
  // after kekulizing and adding the chiral hydrogens (without coordinates
  // then), so the depictor's checks see the molecule so prepared: chiral
  // atoms as isAtomCandForChiralH in MolDraw2DUtils.cpp picks them.
  std::unique_lock<std::mutex> depict;
  std::unique_lock<std::mutex> control_c;
  if (force_coords || !mol.getNumConformers()) {
    depict = lock_depiction();
    control_c = std::unique_lock<std::mutex>(control_c_mutex());
    RWMol with_hs(mol);
    if (kekulize) {
      RDKit::MolOps::KekulizeIfPossible(with_hs, false);
    }
    const auto *rings = with_hs.getRingInfo();
    if (add_chiral_hs && rings->isInitialized()) {
      std::vector<unsigned> chiral;
      for (const auto atom : with_hs.atoms()) {
        if (rings->numAtomRings(atom->getIdx()) > 1u &&
            (atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW ||
             atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW)) {
          chiral.push_back(atom->getIdx());
        }
      }
      if (!chiral.empty()) {
        RDKit::MolOps::addHs(with_hs, false, false, &chiral);
      }
    }
    check_compute_2d_coords_mol(with_hs, false);
  }
  RDKit::MolDraw2DUtils::prepareMolForDrawing(
      *res, kekulize, add_chiral_hs, wedge_bonds, force_coords, wavy_bonds);
  return res;
}

std::unique_ptr<std::string> mol_to_svg(const RWMol &mol, int width, int height,
                                        const rust::Vec<std::uint32_t> &atoms,
                                        bool has_atoms,
                                        std::uint32_t line_width_mult,
                                        bool include_atom_circles,
                                        int conf_id) {
  // As molToSVG in nbWrap/rdMolDraw2D.cpp, on a drawer of the wrappers'
  // own, so it is checked as any other.
  auto d = mol_draw_2d_svg_new(width, height, -1, -1, false);
  auto opts = d->drawer->drawOptions();
  opts.bondLineWidth *= line_width_mult;
  opts.circleAtoms = include_atom_circles;
  opts.prepareMolsBeforeDrawing = false;
  check_draw_options(opts);
  d->drawer->drawOptions() = opts;
  DrawMoleculeData data{};
  data.highlight_atoms = atoms;
  data.has_highlight_atoms = has_atoms;
  data.has_highlight_bonds = false;
  data.conf_id = conf_id;
  draw_molecule(*d, mol, data, true);
  svg_finish_drawing(*d);
  return svg_drawing_text(*d);
}

double mean_bond_length(const RWMol &mol, int conf_id) {
  // Coordinates are summed as doubles; any conformer is fine.
  return RDKit::MolDraw2DUtils::meanBondLength(mol, conf_id);
}

}  // namespace rdkit_cxx
