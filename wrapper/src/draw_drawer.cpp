// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The drawers: MolDraw2DSVG (and MolDraw2DCairo), and the functions that
// draw with them.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>
#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif
// The text drawer's and DrawMol's interfaces, which RDKit does not install;
// wrapper/vendor/<release> holds RDKit's own copies, from the release found,
// so the layouts the wrappers build match the library's.
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawMolMCHCircleAndLine.h>
#include <GraphMol/MolDraw2D/DrawMolMCHLasso.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

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
#include "rdkit-cxx/util.h"

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

#ifdef RDK_BUILD_CAIRO_SUPPORT
// MolDraw2DCairo makes its canvas a Cairo image surface when the canvas is
// sized (on construction, or at the first drawing on a flexible canvas).
// When Cairo cannot make it (a size Cairo or pixman refuses, or no memory),
// it gives a surface in an error state, which drawing leaves alone and from
// which getDrawingText gets no PNG; addMetadataToPNG then checks the
// missing signature (checkPNGHeader compares bytes it failed to read,
// which are uninitialized). So png() gets the image with the metadata
// (protected d_metadata) set aside, checks it, and adds the metadata as
// addMetadataToPNG would.
class CairoDrawer : public RDKit::MolDraw2DCairo, public TextDrawerAccess {
 public:
  CairoDrawer(int width, int height, int panel_width, int panel_height,
              bool no_freetype)
      : RDKit::MolDraw2DCairo(width, height, panel_width, panel_height,
                              no_freetype) {}
  TextDrawer &text_drawer() override { return *text_drawer_; }

  std::string png() {
    // A flexible canvas has no Cairo context until its first drawing sizes
    // it (getDrawingText's PRECONDITION).
    if (width() < 0) {
      throw std::runtime_error("no draw context");
    }
    std::vector<std::pair<std::string, std::string>> metadata;
    metadata.swap(d_metadata);
    std::string res;
    try {
      res = getDrawingText();
    } catch (...) {
      metadata.swap(d_metadata);
      throw;
    }
    metadata.swap(d_metadata);
    static const char signature[] = "\x89PNG\r\n\x1a\n";
    if (res.compare(0, 8, signature, 8) != 0) {
      throw std::runtime_error(
          "Cairo could not make the image: its surface failed (a canvas too "
          "large, or no memory)");
    }
    if (!d_metadata.empty()) {
      res = RDKit::addMetadataToPNGString(res, d_metadata);
    }
    return res;
  }
};
#endif

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

// The shortest line comic mode draws by hand, in the units it is given in.
// handdrawnLine (MolDraw2DDetails.cpp) shifts each end by up to 0.02 in
// each direction and then shrinks its wiggles until they are small next to
// the line's steps, which never ends for a step of no length: its scale
// factor sticks at the smallest denormal.
constexpr double kMinComicLineLength = 0.06;

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
    // Free drawing (DrawString) draws with the font the text drawer has.
    text.setFontFile("");
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

// The locks a drawing call holds for the whole call (CONVENTIONS.md U4),
// in U4's order:
// - lock_depiction(), when it depicts a molecule without conformers
//   (DrawnMol) or a reaction template without one;
// - control_c_mutex(): depiction matches ring templates to the molecule
//   (control_c_lock.h); and every DrawMol (a molecule's layout, the
//   wrappers' and RDKit's) runs prepareMolForDrawing, which, as RDKit's
//   reaction components do (makeReactionDrawMol), blocks the logs with
//   RDLog::LogStateSetter, writing every logger's enabled flag without a
//   lock. RDKit's other LogStateSetter users the wrappers reach
//   (deprotection, R-group decomposition, link nodes, the embedder's
//   torsion preferences) match substructures, so they hold this lock too,
//   and it orders them all.
// - random_mutex(), for a drawing in comic mode, which draws its lines with
//   std::rand (handdrawnLine in MolDraw2DDetails.cpp), a global with no
//   lock on some platforms.
class DrawingLocks {
 public:
  DrawingLocks(bool depicts, bool comic) {
    if (depicts) {
      d_depict = lock_depiction();
    }
    d_control_c = std::unique_lock<ControlCMutex>(control_c_mutex());
    if (comic) {
      d_random = std::unique_lock<std::mutex>(random_mutex());
    }
  }

 private:
  std::unique_lock<std::mutex> d_depict;
  std::unique_lock<ControlCMutex> d_control_c;
  std::unique_lock<std::mutex> d_random;
};

bool lacks_conformers(const RWMol &mol) { return !mol.getNumConformers(); }

// A molecule as DrawMol will draw it. DrawMol gives a molecule without
// conformers 2D coordinates (initDrawMolecule), after
// MolDraw2D::drawMolecule has symmetrized its rings (drawMoleculeWithHighlights
// does not: symmetrize); this does the same on a copy, so the coordinates
// can be checked before RDKit draws, and DrawMol depicts no more. The
// caller holds DrawingLocks, with the depiction lock for a molecule without
// conformers.
class DrawnMol {
 public:
  DrawnMol(const RWMol &mol, bool symmetrize = true) : d_mol(&mol) {
    if (lacks_conformers(mol)) {
      d_copy = std::make_unique<RWMol>(mol);
      MatrixCacheGuard guard(*d_copy);
      if (symmetrize && !d_copy->getRingInfo()->isSymmSssr()) {
        RDKit::MolOps::symmetrizeSSSR(*d_copy);
      }
      check_compute_2d_coords_mol(*d_copy, false);
      const bool canon_orient = true;
      RDDepict::compute2DCoords(*d_copy, nullptr, canon_orient);
      d_mol = d_copy.get();
    }
  }
  const RWMol &mol() const { return *d_mol; }

 private:
  const RWMol *d_mol;
  std::unique_ptr<RWMol> d_copy;
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
// takes from the shapes, labels, and annotations; DrawMol leaves out bonds
// shorter than 0.01, so a molecule that small has none, and findExtremes
// leaves the extremes at +-DBL_MAX / 2, from which every drawn coordinate is
// measured. A flexible canvas takes its size from the
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
          "shorter than 0.01 (which RDKit does not draw)");
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
          "shapes it cannot draw (coordinates or angles that are not finite "
          "or out of range, or wavy lines with no room for their waves)");
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

// After SetScale, MolDraw2D::drawTheMolecule gives each molecule the
// transformation it set (DrawMol::setTransformation), which lays the
// molecule out again at that scale.
void checked_transformation(const Drawer &d, detail::DrawMol &dm) {
  if (!d.global) {
    return;
  }
  dm.setTransformation(*d.global);
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

// A drawing's size in pixels: its DrawMol's, which sizes a flexible canvas
// (MolDraw2D::fixVariableDimensions).
using Size = std::pair<int, int>;

// Lays out a molecule as MolDraw2D::drawMolecule (drawing) and getMolSize
// will, at the panel size, giving its size.
Size check_layout(Drawer &d, const RWMol &mol, const std::string &legend,
                  const Highlights &h, const std::vector<int> *bonds,
                  int conf_id, bool drawing = true) {
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
  if (drawing) {
    checked_transformation(d, probe);
  }
  return {probe.width_, probe.height_};
}

// As check_layout, for drawMoleculeWithHighlights's DrawMolMCH, which RDKit
// builds from the molecule as given (it does not symmetrize its rings).
Size check_mch_layout(Drawer &d, const RWMol &mol, const std::string &legend,
                      const Drawer::MultiHighlights &h, int conf_id) {
  set_up_text_drawer(d);
  auto &opts = d.drawer->drawOptions();
  std::unique_ptr<detail::DrawMol> probe;
  switch (opts.multiColourHighlightStyle) {
    case RDKit::MultiColourHighlightStyle::CIRCLEANDLINE:
      probe = std::make_unique<detail::DrawMolMCHCircleAndLine>(
          mol, legend, d.drawer->panelWidth(), d.drawer->panelHeight(), opts,
          d.text->text_drawer(), h.atoms, h.bonds, h.radii, h.multipliers,
          conf_id);
      break;
    case RDKit::MultiColourHighlightStyle::LASSO:
      probe = std::make_unique<detail::DrawMolMCHLasso>(
          mol, legend, d.drawer->panelWidth(), d.drawer->panelHeight(), opts,
          d.text->text_drawer(), h.atoms, h.bonds, h.radii, h.multipliers,
          conf_id);
      break;
  }
  checked_layout(*probe);
  checked_transformation(d, *probe);
  return {probe->width_, probe->height_};
}

// The largest canvas, in bytes of its pixels (4 each), pixman makes Cairo's
// image surfaces of: it refuses larger ones, and Cairo then gives a surface
// in an error state (CairoDrawer), as it does for a side over 32767 pixels.
constexpr std::int64_t kMaxCairoCanvasBytes = INT_MAX;

void check_cairo_size(int width, int height) {
  if (!(width <= kMaxCanvasSize && height <= kMaxCanvasSize &&
        std::int64_t{width} * height * 4 <= kMaxCairoCanvasBytes)) {
    std::ostringstream msg;
    msg << "a Cairo canvas of " << width << " by " << height
        << " pixels is too large: Cairo's are at most 32767 pixels a side "
           "and 2^31 bytes (4 a pixel)";
    throw ValueErrorException(msg.str());
  }
}

// The first drawing on a flexible canvas sizes it, and MolDraw2DCairo then
// makes its surface of that size (startDrawing), which must be one Cairo
// can make.
void check_canvas_fits(const Drawer &d, Size size) {
  if (!d.svg && d.drawer->width() < 0) {
    check_cairo_size(size.first, size.second);
  }
}

// A size RDKit's layout computes in double and keeps in int: refused unless
// finite and well within int's range, so the sums RDKit makes of such sizes
// stay in it.
void check_pixels(double v) {
  if (!(std::abs(v) <= INT_MAX / 8)) {
    throw ValueErrorException(
        "the drawing comes out too large (or degenerate) for RDKit's pixel "
        "arithmetic, which it does in int");
  }
}

// DrawMol::setScale, and the checks checked_layout makes after it.
void checked_set_scale(detail::DrawMol &dm, double scale, double font_scale,
                       bool ignore_font_limits) {
  if (!(scale > 0 && std::isfinite(scale) && font_scale > 0 &&
        std::isfinite(font_scale))) {
    refuse_layout();
  }
  dm.setScale(scale, font_scale, ignore_font_limits);
  check_shapes(dm.bonds_);
  check_shapes(dm.preShapes_);
  check_shapes(dm.postShapes_);
  check_shapes(dm.highlights_);
}

// DrawMol::shrinkToFit converts its extents times its scale to int.
void checked_shrink_to_fit(detail::DrawMol &dm) {
  check_pixels(dm.xRange_ * dm.scale_);
  check_pixels(dm.yRange_ * dm.scale_);
  dm.shrinkToFit(false);
  if (!(dm.width_ >= 1 && dm.height_ >= 1)) {
    refuse_layout();
  }
}

// MolDraw2D::reactionWidth, in 64 bits.
std::int64_t reaction_width(
    const std::vector<std::unique_ptr<detail::DrawMol>> &reagents,
    const std::vector<std::unique_ptr<detail::DrawMol>> &products,
    const std::vector<std::unique_ptr<detail::DrawMol>> &agents, double padding,
    std::int64_t arrow_mult, std::int64_t gap) {
  std::int64_t tot = 0;
  if (!reagents.empty()) {
    for (const auto &dm : reagents) {
      tot += dm->width_;
    }
    tot += gap * static_cast<std::int64_t>(reagents.size() - 1);
    tot += gap / 2;
  }
  if (agents.empty()) {
    tot += arrow_mult * gap;
  } else {
    tot += gap / 2;
    for (const auto &dm : agents) {
      tot += dm->width_ + gap / 2;
    }
  }
  if (!products.empty()) {
    tot += gap / 2;
    for (const auto &dm : products) {
      tot += dm->width_;
    }
    tot += gap * static_cast<std::int64_t>(products.size() - 1);
  }
  const double total = tot + 2.0 * padding * tot;
  check_pixels(total);
  return static_cast<std::int64_t>(total);
}

// MolDraw2D::drawReaction lays its reactants, products, and agents out, each
// as a DrawMol of the panel's height and a width of its own (with no
// padding), scales them to the smallest scale of their group, shrinks them
// to fit, scales them down again to fit the canvas's width
// (calcReactionOffsets), and once more for the component padding
// (drawReactionPart), all in pixels it keeps in int. The layout is repeated
// here, on a copy of the reaction, with each scale and size checked before
// RDKit divides by it or converts it to int. It reads every DrawMol the
// drawer holds, so a reaction is drawn on a drawer that has drawn nothing
// else.
class ReactionLayout {
 public:
  ReactionLayout(Drawer &d, const RDKit::ChemicalReaction &rxn,
                 bool by_reactant,
                 const std::vector<RDKit::DrawColour> *colours,
                 const std::vector<int> *conf_ids)
      : d_drawer(d), d_rxn(rxn), d_conf_ids(conf_ids) {
    auto &opts = d.drawer->drawOptions();
    // findReactionHighlights, after ChemicalReactionToRxnMol, which sets the
    // templates' atoms' reaction roles.
    std::unique_ptr<RDKit::ROMol> rxn_mol(
        RDKit::ChemicalReactionToRxnMol(d_rxn));
    if (by_reactant) {
      const auto *cols = colours ? colours : &opts.highlightColourPalette;
      // It takes each reactant's color modulo their number.
      if (cols->empty()) {
        throw ValueErrorException(
            "highlightByReactant needs at least one color");
      }
      for (std::size_t m = 0; m < d_rxn.getReactants().size(); ++m) {
        for (const auto atom : d_rxn.getReactants()[m]->atoms()) {
          int role = -1;
          if (atom->getPropIfPresent("molRxnRole", role) && role == 1 &&
              atom->getAtomMapNum()) {
            d_colours.insert(std::make_pair(atom->getAtomMapNum(),
                                            (*cols)[m % cols->size()]));
          }
        }
      }
    }
    const int panel_height = d.drawer->panelHeight();
    const int use_height = panel_height > 0
                               ? panel_height * (1.0 - 2.0 * opts.padding)
                               : panel_height;
    double min_scale = std::numeric_limits<double>::max();
    double min_font_scale = 0;
    components(d_rxn.getReactants(), use_height, d_reagents, min_scale,
               min_font_scale);
    components(d_rxn.getProducts(), use_height, d_products, min_scale,
               min_font_scale);
    for (auto *group : {&d_reagents, &d_products}) {
      for (auto &dm : *group) {
        checked_set_scale(*dm, min_scale, min_font_scale, true);
        checked_shrink_to_fit(*dm);
      }
    }
    // plusWidth = minScale, in int: DBL_MAX for a reaction with neither
    // reactants nor products.
    if (d_reagents.empty() && d_products.empty()) {
      throw ValueErrorException(
          "a reaction needs a reactant or a product to be drawn");
    }
    check_pixels(min_scale);
    d_plus_width = static_cast<std::int64_t>(min_scale);
    const int agent_height = int(0.4 * use_height);
    min_scale = std::numeric_limits<double>::max();
    components(d_rxn.getAgents(), agent_height, d_agents, min_scale,
               min_font_scale);
    for (auto &dm : d_agents) {
      checked_set_scale(*dm, min_scale, min_font_scale, true);
      checked_shrink_to_fit(*dm);
    }
    offsets(opts);
  }

  std::size_t num_components() const {
    return d_reagents.size() + d_products.size() + d_agents.size();
  }

  // The canvas's size once drawReaction has sized it.
  Size canvas_size() const { return d_canvas; }

  // Comic mode draws the plus signs and the arrow's shaft by hand, which
  // never ends for a line of no length (kMinComicLineLength). A plus's
  // strokes are plusWidth / 3 pixels long, in int, and drawArrow pulls the
  // shaft's end back by its mitre, 0.5 * lineWidth / sin(pi / 6), with
  // the line width the last component's scale gives (getDrawLineWidth).
  void check_comic(const RDKit::MolDraw2D &drawer) const {
    const bool plus = d_final_plus_width != 0 &&
                      (d_reagents.size() > 1 || d_products.size() > 1);
    const auto stroke = d_final_plus_width > 30 ? 10 : d_final_plus_width / 3;
    if (plus && stroke == 0) {
      throw ValueErrorException(
          "comicMode cannot draw the reaction's plus signs on a canvas this "
          "small: they would be lines of no length");
    }
    double line_width = drawer.lineWidth();
    if (drawer.drawOptions().scaleBondWidth) {
      line_width = std::max(0.0, line_width * last_drawn().getScale() * 0.02);
    }
    const double shaft =
        d_arrow_length > 1.0e-6
            ? std::abs(d_arrow_length - 0.5 * line_width / std::sin(M_PI / 6))
            : 0.0;
    if (!(shaft >= kMinComicLineLength)) {
      throw ValueErrorException(
          "comicMode cannot draw the reaction's arrow on a canvas this small: "
          "its shaft would be shorter than 0.06");
    }
  }

 private:
  // makeReactionComponents and makeReactionDrawMol.
  void components(const std::vector<RDKit::ROMOL_SPTR> &bits, int height,
                  std::vector<std::unique_ptr<detail::DrawMol>> &dms,
                  double &min_scale, double &min_font_scale) {
    auto &opts = d_drawer.drawer->drawOptions();
    for (std::size_t m = 0; m < bits.size(); ++m) {
      std::vector<int> atoms, bonds;
      std::map<int, RDKit::DrawColour> atom_map, bond_map;
      const int cid = d_conf_ids ? (*d_conf_ids)[m] : -1;
      auto *frag = static_cast<RWMol *>(bits[m].get());
      for (auto atom : frag->atoms()) {
        auto ai = d_colours.find(atom->getAtomMapNum());
        for (auto bond : frag->bonds()) {
          auto begi = d_colours.find(bond->getBeginAtom()->getAtomMapNum());
          if (begi != d_colours.end()) {
            auto endi = d_colours.find(bond->getEndAtom()->getAtomMapNum());
            if (endi != d_colours.end() && begi->second == endi->second) {
              bonds.push_back(bond->getIdx());
              bond_map.insert(std::make_pair(bond->getIdx(), begi->second));
            }
          }
        }
        if (ai != d_colours.end()) {
          atoms.push_back(atom->getIdx());
          atom_map.insert(std::make_pair(atom->getIdx(), ai->second));
          atom->setAtomMapNum(0);
        }
      }
      auto &mol = *frag;
      mol.updatePropertyCache(false);
      if (opts.prepareMolsBeforeDrawing) {
        RDKit::MolOps::KekulizeIfPossible(mol, false);
        RDKit::MolOps::setHybridization(mol);
      }
      if (!mol.getNumConformers()) {
        // The caller holds the depiction locks.
        MatrixCacheGuard guard(mol);
        check_compute_2d_coords_mol(mol, false);
        const bool canon_orient = true;
        RDDepict::compute2DCoords(mol, nullptr, canon_orient);
      } else {
        if (!find_conformer(mol, cid)) {
          throw RDKit::ConformerException("Bad Conformer Id");
        }
        detail::centerMolForDrawing(mol, cid);
      }
      RDKit::MolDraw2DUtils::prepareMolForDrawing(mol, false, false);
      const bool symmetrized = false;
      check_drawn_mol(mol, cid, opts, d_drawer.freetype,
                      d_drawer.drawer->supportsAnnotations(), symmetrized);
      const double padding = opts.padding;
      opts.padding = 0.0;
      try {
        dms.push_back(std::make_unique<detail::DrawMol>(
            mol, "", -1, height, opts, d_drawer.text->text_drawer(), &atoms,
            &bonds, &atom_map, &bond_map, nullptr, nullptr,
            d_drawer.drawer->supportsAnnotations(), cid, true));
      } catch (...) {
        opts.padding = padding;
        throw;
      }
      opts.padding = padding;
      checked_layout(*dms.back());
      if (dms.back()->getScale() < min_scale) {
        min_scale = dms.back()->getScale();
        min_font_scale = dms.back()->getFontScale();
      }
    }
  }

  std::vector<detail::DrawMol *> all() const {
    std::vector<detail::DrawMol *> res;
    for (const auto *group : {&d_reagents, &d_products, &d_agents}) {
      for (const auto &dm : *group) {
        res.push_back(dm.get());
      }
    }
    return res;
  }

  // drawReaction's height for a flexible canvas, calcReactionOffsets, and
  // drawReactionPart's scales.
  void offsets(const RDKit::MolDrawOptions &opts) {
    const auto &d = *d_drawer.drawer;
    d_canvas = {d.width(), d.height()};
    if (d.height() == -1) {
      int height = -1;
      for (const auto *dm : all()) {
        height = std::max(height, dm->height_);
      }
      check_pixels(height + height * 2.0 * opts.padding);
      d_canvas.second = static_cast<int>(height + height * 2.0 * opts.padding);
    }
    const std::int64_t arrow_mult = 2;
    std::int64_t plus_width = d_plus_width;
    std::int64_t width_to_use =
        d.panelWidth() == -1 ? d.width() : d.panelWidth();
    if (width_to_use == -1) {
      width_to_use = reaction_width(d_reagents, d_products, d_agents,
                                    opts.padding, arrow_mult, plus_width);
      if (d.width() == -1) {
        d_canvas.first = static_cast<int>(width_to_use);
      }
    }
    auto tot = reaction_width(d_reagents, d_products, d_agents, opts.padding,
                              arrow_mult, plus_width);
    for (int i = 0; i < 5; ++i) {
      int max_width = std::numeric_limits<int>::lowest();
      for (const auto *dm : all()) {
        max_width = std::max(max_width, dm->width_);
      }
      plus_width = max_width / 4;
      plus_width =
          plus_width > width_to_use / 20 ? width_to_use / 20 : plus_width;
      const auto old_tot = tot;
      if (tot <= 0) {
        refuse_layout();
      }
      const double stretch =
          double(width_to_use * (1 - 2.0 * opts.padding)) / tot;
      if (!(stretch < 1.0)) {
        break;
      }
      for (auto *dm : all()) {
        checked_set_scale(*dm, stretch * dm->getScale(),
                          stretch * dm->getFontScale(), false);
        checked_shrink_to_fit(*dm);
      }
      tot = reaction_width(d_reagents, d_products, d_agents, opts.padding,
                           arrow_mult, plus_width);
      if (std::abs(double(tot - old_tot)) < 0.01 * width_to_use) {
        break;
      }
    }
    if (!d_reagents.empty() || !d_products.empty()) {
      // plusWidth = the first reactant's (or product's) scale, in int.
      check_pixels(first_scale());
      plus_width = static_cast<std::int64_t>(first_scale());
    }
    // With room for 5 more pixels at each gap, plusWidth grows by 5; the
    // arrow spans 2 plusWidths, or the agents and half a plusWidth around
    // each (calcReactionOffsets).
    d_final_plus_width = plus_width;
    std::int64_t gaps = (d_reagents.empty() ? 0 : d_reagents.size() - 1) +
                        (d_products.empty() ? 0 : d_products.size() - 1) + 1 +
                        (d_agents.empty() ? 2 : d_agents.size() - 1);
    if (width_to_use - reaction_width(d_reagents, d_products, d_agents,
                                      opts.padding, arrow_mult, plus_width) >
        gaps * 5) {
      d_final_plus_width += 5;
    }
    if (d_agents.empty()) {
      d_arrow_length = double(arrow_mult * d_final_plus_width);
    } else {
      std::int64_t length = d_final_plus_width / 2;
      for (const auto &dm : d_agents) {
        length += dm->width_ + d_final_plus_width / 2;
      }
      d_arrow_length = double(length);
    }
    // With room to spare, plusWidth grows by 5.
    tot = reaction_width(d_reagents, d_products, d_agents, opts.padding,
                         arrow_mult, plus_width + 5);
    check_pixels(double(width_to_use) + d.offset().x);
    check_pixels(tot / (1 + 2.0 * opts.padding));
    // drawReactionPart: each component's scale, less the component
    // padding.
    for (auto *dm : all()) {
      const double wp = dm->width_ * opts.componentPadding;
      const double hp = dm->height_ * opts.componentPadding;
      const double scale = std::min((dm->width_ - 2 * wp) / dm->width_,
                                    (dm->height_ - 2 * hp) / dm->height_) *
                           dm->getScale();
      checked_set_scale(*dm, scale, scale, true);
    }
  }

  // The last component drawReaction draws, whose scale is the drawer's when
  // it draws the arrow.
  const detail::DrawMol &last_drawn() const {
    if (!d_products.empty()) {
      return *d_products.back();
    }
    if (!d_agents.empty()) {
      return *d_agents.back();
    }
    return *d_reagents.back();
  }

  double first_scale() const {
    if (!d_reagents.empty()) {
      return d_reagents.front()->scale_;
    }
    if (!d_products.empty()) {
      return d_products.front()->scale_;
    }
    return 0;
  }

  Drawer &d_drawer;
  RDKit::ChemicalReaction d_rxn;
  const std::vector<int> *d_conf_ids;
  std::map<int, RDKit::DrawColour> d_colours;
  std::int64_t d_plus_width = 0;
  std::int64_t d_final_plus_width = 0;
  double d_arrow_length = 0;
  Size d_canvas;
  std::vector<std::unique_ptr<detail::DrawMol>> d_reagents, d_products,
      d_agents;
};

// What drawing a molecule on the drawer trusts, beyond the options, which
// the setters check: the drawing area, the font, the legend, and the
// molecule as DrawMol sees it.
void check_drawing(Drawer &d, const RWMol &mol, int conf_id,
                   const std::string &legend, bool annotations,
                   bool symmetrized) {
  const auto &opts = d.drawer->drawOptions();
  check_draw_area(d.drawer->panelWidth(), d.drawer->panelHeight(), opts);
  check_font(d);
  check_legend(legend, d.freetype);
  check_drawn_mol(mol, conf_id, opts, d.freetype, annotations, symmetrized);
}

// MolDraw2DCairo, with includeMetadata, pickles each molecule and reaction
// it draws into the PNG (updateMetadata in MolDraw2DCairo.cpp): the molecule
// as drawn, which only adds hydrogens to the one given, so passes
// check_picklable whenever that one does.
void check_metadata(const Drawer &d, const RDKit::ROMol &mol) {
  if (!d.svg && d.drawer->drawOptions().includeMetadata) {
    check_picklable(mol);
  }
}

bool comic(const Drawer &d) { return d.drawer->drawOptions().comicMode; }

// Without highlight bonds, Python's first DrawMolecule overload highlights
// the bonds between highlighted atoms (bonds_from_atoms); its other drawing
// functions highlight none. The caller holds DrawingLocks.
void draw_molecule(Drawer &d, const RWMol &mol, const DrawMoleculeData &data,
                   bool bonds_from_atoms) {
  check_usable(d);
  const auto h = read_highlights(mol, data);
  const std::string legend(data.legend);
  DrawnMol drawn(mol);
  check_drawing(d, drawn.mol(), data.conf_id, legend,
                d.drawer->supportsAnnotations(), true);
  check_metadata(d, drawn.mol());
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
  check_canvas_fits(
      d, check_layout(d, drawn.mol(), legend, h, bonds, data.conf_id));
  drawing_call(d, [&] {
    d.drawer->drawMolecule(
        drawn.mol(), legend, h.has_atoms ? &h.atoms : nullptr, bonds,
        &h.atom_colours, &h.bond_colours, &h.radii, data.conf_id);
  });
  ++d.num_mols;
  ++d.active_mol;
}

// DrawMoleculeWithHighlights's maps as RDKit takes them, checked: indices
// (DrawMolMCHLasso sets bits of an atom's index without checking it),
// colors, radii, and line width multipliers, which getHighlightBondWidth
// multiplies the line width by and DrawMolMCHCircleAndLine converts to int.
// A highlight of no radius makes
// DrawMolMCHCircleAndLine::fixHighlightJoinProblems grow it, by multiples,
// forever, so radii (and highlightRadius, for atoms without one) must be
// positive.
Drawer::MultiHighlights read_multi_highlights(
    const RWMol &mol, const DrawMoleculeWithHighlightsData &data,
    const RDKit::MolDrawOptions &opts) {
  Drawer::MultiHighlights h;
  auto colours = [](const ColourList &list) {
    std::vector<RDKit::DrawColour> res;
    for (const auto &c : list.colours) {
      res.push_back(colour_from_data(c));
      check_colour(res.back(), "highlight color");
    }
    return res;
  };
  for (std::size_t i = 0; i < data.atom_keys.size(); ++i) {
    check_index(data.atom_keys[i], mol.getNumAtoms());
    h.atoms[static_cast<int>(data.atom_keys[i])] =
        colours(data.atom_colours[i]);
  }
  for (std::size_t i = 0; i < data.bond_keys.size(); ++i) {
    check_index(data.bond_keys[i], mol.getNumBonds());
    h.bonds[static_cast<int>(data.bond_keys[i])] =
        colours(data.bond_colours[i]);
  }
  for (std::size_t i = 0; i < data.radii_keys.size(); ++i) {
    const double r = data.radii_values[i];
    if (!(r > 0 && r <= kMaxHighlightRadius)) {
      throw ValueErrorException(
          "multi-color highlight radii must be above 0 and at most 16");
    }
    h.radii[static_cast<int>(data.radii_keys[i])] = r;
  }
  for (std::size_t i = 0; i < data.multiplier_keys.size(); ++i) {
    const int m = data.multiplier_values[i];
    if (!(m >= 0 && m <= 64)) {
      throw ValueErrorException(
          "highlight line width multipliers must be from 0 to 64");
    }
    h.multipliers[static_cast<int>(data.multiplier_keys[i])] = m;
  }
  if ((!h.atoms.empty() || !h.bonds.empty()) && !(opts.highlightRadius > 0)) {
    throw ValueErrorException(
        "multi-color highlights need a positive highlightRadius");
  }
  return h;
}

// DrawMolMCHCircleAndLine draws a highlighted bond's colors along the
// bond's perpendicular, which a bond whose atoms are at one position makes
// NaN, and fixHighlightJoinProblems then grows the atoms' highlights,
// forever, to meet those lines.
void check_mch_bonds(const RWMol &mol, const Drawer::MultiHighlights &h,
                     const RDKit::MolDrawOptions &opts) {
  if (opts.multiColourHighlightStyle !=
          RDKit::MultiColourHighlightStyle::CIRCLEANDLINE ||
      h.bonds.empty() || !mol.getNumConformers()) {
    return;
  }
  const auto &conf = mol.getConformer();
  for (const auto &[idx, colours] : h.bonds) {
    const auto *bond = mol.getBondWithIdx(idx);
    const auto &a = conf.getAtomPos(bond->getBeginAtomIdx());
    const auto &b = conf.getAtomPos(bond->getEndAtomIdx());
    if (a.x == b.x && a.y == b.y) {
      throw ValueErrorException(
          "a highlighted bond's atoms are at one position, so RDKit's "
          "multi-color highlight of it never ends");
    }
  }
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
  check_colour(c, "color");
  // MolDraw2DCairo::setColour sets its Cairo context's color, which it has
  // none of before a flexible canvas is sized (RDKit's PRECONDITION).
  if (!d.svg && d.drawer->width() < 0) {
    throw std::runtime_error("no draw context");
  }
  d.drawer->setColour(c);
}

void drawer_clear_drawing(Drawer &d) {
  // MolDraw2D::clearDrawing sets the canvas up (initDrawing) if it is not
  // yet, which a flexible canvas is not until its first drawing sizes it:
  // RDKit would make it -1 by -1 pixels, and not again (an SVG header of
  // that size, or a Cairo surface in an error state, which draws nothing).
  if (d.drawer->width() < 0) {
    throw ValueErrorException(
        "a flexible canvas has no size to clear until something is drawn on "
        "it");
  }
  d.drawer->clearDrawing();
}

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
  const DrawingLocks locks(lacks_conformers(mol), comic(d));
  draw_molecule(d, mol, data, true);
}

DrawSize drawer_get_mol_size(Drawer &d, const RWMol &mol,
                             const DrawMoleculeData &data) {
  check_usable(d);
  const auto h = read_highlights(mol, data);
  const std::string legend(data.legend);
  // It draws nothing, so no std::rand.
  const DrawingLocks locks(lacks_conformers(mol), false);
  DrawnMol drawn(mol);
  check_drawing(d, drawn.mol(), data.conf_id, legend,
                d.drawer->supportsAnnotations(), true);
  const bool drawing = false;
  check_layout(d, drawn.mol(), legend, h, h.has_bonds ? &h.bonds : nullptr,
               data.conf_id, drawing);
  // GetMolSize lays the molecule out in a DrawMol of its own, and leaves the
  // drawer's molecules as they were. Python's passes the bonds it is given,
  // or none.
  const auto size = d.drawer->getMolSize(
      drawn.mol(), legend, h.has_atoms ? &h.atoms : nullptr,
      h.has_bonds ? &h.bonds : nullptr, &h.atom_colours, &h.bond_colours,
      &h.radii, data.conf_id);
  return DrawSize{size.first, size.second};
}

void drawer_draw_molecule_with_highlights(
    Drawer &d, const RWMol &mol, const DrawMoleculeWithHighlightsData &data) {
  check_usable(d);
  const auto &opts = d.drawer->drawOptions();
  auto h = read_multi_highlights(mol, data, opts);
  const std::string legend(data.legend);
  const bool symmetrize = false;
  const DrawingLocks locks(lacks_conformers(mol), comic(d));
  DrawnMol drawn(mol, symmetrize);
  check_mch_bonds(drawn.mol(), h, opts);
  // DrawMolMCH passes confId where DrawMol takes includeAnnotations, so
  // RDKit draws the first conformer, with annotations unless confId is 0.
  check_drawing(d, drawn.mol(), -1, legend, data.conf_id != 0, symmetrize);
  check_metadata(d, drawn.mol());
  check_canvas_fits(d,
                    check_mch_layout(d, drawn.mol(), legend, h, data.conf_id));
  d.multi_highlights.push_back(std::move(h));
  const auto &kept = d.multi_highlights.back();
  drawing_call(d, [&] {
    d.drawer->drawMoleculeWithHighlights(drawn.mol(), legend, kept.atoms,
                                         kept.bonds, kept.radii,
                                         kept.multipliers, data.conf_id);
  });
  ++d.num_mols;
  ++d.active_mol;
}

// MolDraw2D::drawMolecules lays out each molecule as drawMolecule does (but
// with its rings as given), at the panel size, offset to its panel, then
// gives them all the smallest scale (drawMolsSameScale). It does that for
// every molecule the drawer holds, with the options of the moment, so a
// grid is drawn on a drawer that has drawn nothing else; and it places
// molecule i in row i / (width / panelWidth), which divides by zero for a
// panel wider than the canvas. A missing molecule (a null pointer, Python's
// None) leaves its panel empty, as RDKit skips it.
void drawer_draw_molecules(Drawer &d, BorrowedMols &mols,
                           const DrawMoleculesData &data) {
  check_usable(d);
  const auto &list = mols.mols();
  const std::size_t n = list.size();
  if (!n) {
    return;
  }
  auto check_length = [n](bool has, std::size_t len, const char *what) {
    if (has && len != n) {
      throw ValueErrorException(std::string("If ") + what +
                                " is provided it must be the same length as "
                                "the molecule list.");
    }
  };
  check_length(data.has_highlight_atoms, data.highlight_atoms.size(),
               "highlightAtoms");
  check_length(data.has_highlight_bonds, data.highlight_bonds.size(),
               "highlightBonds");
  check_length(data.has_atom_colours, data.atom_colours.size(),
               "highlightAtomColors");
  check_length(data.has_bond_colours, data.bond_colours.size(),
               "highlightBondColors");
  check_length(data.has_radii, data.radii.size(), "highlightAtomRadii");
  check_length(data.has_conf_ids, data.conf_ids.size(), "confIds");
  check_length(data.has_legends, data.legends.size(), "legends");
  auto &draw = *d.drawer;
  if (!(draw.width() > 0 && draw.height() > 0)) {
    throw ValueErrorException("DrawMolecules needs a canvas of fixed size");
  }
  if (d.num_mols || d.global_transform) {
    throw ValueErrorException(
        "DrawMolecules lays out again every molecule the drawer holds; draw "
        "the grid on a new drawer");
  }
  const int pw = draw.panelWidth();
  const int ph = draw.panelHeight();
  if (pw > draw.width() || ph > draw.height()) {
    throw ValueErrorException("the panels must fit on the canvas");
  }
  const int n_cols = draw.width() / pw;
  const int n_rows = draw.height() / ph;
  check_pixels(double(n_rows > 1 ? (n - 1) / n_cols : 0) * ph);
  const auto &opts = draw.drawOptions();
  bool depicts = false;
  for (const auto &m : list) {
    depicts = depicts || (m && !m->getNumConformers());
  }
  const DrawingLocks locks(depicts, comic(d));

  struct Entry {
    std::unique_ptr<DrawnMol> drawn;
    Highlights h;
    std::vector<int> bonds;
    std::string legend;
    int conf_id;
  };
  std::vector<Entry> entries(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (!list[i]) {
      continue;
    }
    const auto &mol = static_cast<const RWMol &>(*list[i]);
    auto &e = entries[i];
    e.h.has_atoms = data.has_highlight_atoms;
    if (e.h.has_atoms) {
      const auto &v = data.highlight_atoms[i].indices;
      e.h.atoms.assign(v.begin(), v.end());
    }
    e.h.has_bonds = data.has_highlight_bonds;
    if (e.h.has_bonds) {
      const auto &v = data.highlight_bonds[i].indices;
      e.h.bonds.assign(v.begin(), v.end());
    }
    check_highlight_atoms(mol, e.h.atoms);
    check_highlight_bonds(mol, e.h.bonds);
    auto read_colours = [](const IndexColourMap &m,
                           std::map<int, RDKit::DrawColour> &res) {
      for (std::size_t k = 0; k < m.keys.size(); ++k) {
        res[static_cast<int>(m.keys[k])] = colour_from_data(m.values[k]);
      }
      check_highlight_colours(res);
    };
    if (data.has_atom_colours) {
      read_colours(data.atom_colours[i], e.h.atom_colours);
    }
    if (data.has_bond_colours) {
      read_colours(data.bond_colours[i], e.h.bond_colours);
    }
    if (data.has_radii) {
      const auto &m = data.radii[i];
      for (std::size_t k = 0; k < m.keys.size(); ++k) {
        e.h.radii[static_cast<int>(m.keys[k])] = m.values[k];
      }
      check_highlight_radii(e.h.radii);
    }
    e.conf_id = data.has_conf_ids ? data.conf_ids[i] : -1;
    e.legend = data.has_legends ? std::string(data.legends[i]) : "";
    // RDKit draws each molecule as given, without symmetrizing its rings.
    const bool symmetrize = false;
    e.drawn = std::make_unique<DrawnMol>(mol, symmetrize);
    check_drawing(d, e.drawn->mol(), e.conf_id, e.legend,
                  draw.supportsAnnotations(), symmetrize);
    check_metadata(d, e.drawn->mol());
    if (e.h.has_bonds) {
      e.bonds = e.h.bonds;
    } else if (opts.continuousHighlight && e.h.has_atoms) {
      detail::getBondHighlightsForAtoms(e.drawn->mol(), e.h.atoms, e.bonds);
    }
  }
  set_up_text_drawer(d);
  std::vector<std::unique_ptr<detail::DrawMol>> probes;
  std::size_t min_scale = 0;
  std::size_t min_font_scale = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto &e = entries[i];
    if (!e.drawn) {
      continue;
    }
    const bool bonds =
        e.h.has_bonds || (opts.continuousHighlight && e.h.has_atoms);
    probes.push_back(std::make_unique<detail::DrawMol>(
        e.drawn->mol(), e.legend, pw, ph, opts, d.text->text_drawer(),
        e.h.has_atoms ? &e.h.atoms : nullptr, bonds ? &e.bonds : nullptr,
        data.has_atom_colours ? &e.h.atom_colours : nullptr,
        data.has_bond_colours ? &e.h.bond_colours : nullptr, nullptr,
        data.has_radii ? &e.h.radii : nullptr, draw.supportsAnnotations(),
        e.conf_id));
    const int row = n_rows > 1 ? static_cast<int>(i) / n_cols : 0;
    const int col = n_cols > 1 ? static_cast<int>(i) % n_cols : 0;
    probes.back()->setOffsets(col * pw, row * ph);
    checked_layout(*probes.back());
    // Indices into the probes, which skip missing molecules as RDKit's
    // drawMols_ do.
    if (probes.back()->getScale() < probes[min_scale]->getScale()) {
      min_scale = probes.size() - 1;
    }
    if (probes.back()->getFontScale() <
        probes[min_font_scale]->getFontScale()) {
      min_font_scale = probes.size() - 1;
    }
  }
  const int drawn = static_cast<int>(probes.size());
  if (opts.drawMolsSameScale && drawn) {
    const double scale = probes[min_scale]->getScale();
    const double font_scale = probes[min_font_scale]->getFontScale();
    for (auto &probe : probes) {
      checked_set_scale(*probe, scale, font_scale, true);
    }
  }
  probes.clear();

  std::vector<RDKit::ROMol *> ptrs;
  std::vector<std::vector<int>> atoms, bonds;
  std::vector<std::map<int, RDKit::DrawColour>> atom_colours, bond_colours;
  std::vector<std::map<int, double>> radii;
  std::vector<int> conf_ids;
  std::vector<std::string> legends;
  for (auto &e : entries) {
    // drawMolecules takes non-const pointers and reads through them.
    ptrs.push_back(e.drawn ? const_cast<RWMol *>(&e.drawn->mol()) : nullptr);
    atoms.push_back(e.h.atoms);
    bonds.push_back(e.h.bonds);
    atom_colours.push_back(e.h.atom_colours);
    bond_colours.push_back(e.h.bond_colours);
    radii.push_back(e.h.radii);
    conf_ids.push_back(e.conf_id);
    legends.push_back(e.legend);
  }
  drawing_call(d, [&] {
    draw.drawMolecules(ptrs, data.has_legends ? &legends : nullptr,
                       data.has_highlight_atoms ? &atoms : nullptr,
                       data.has_highlight_bonds ? &bonds : nullptr,
                       data.has_atom_colours ? &atom_colours : nullptr,
                       data.has_bond_colours ? &bond_colours : nullptr,
                       data.has_radii ? &radii : nullptr,
                       data.has_conf_ids ? &conf_ids : nullptr);
  });
  d.num_mols = drawn;
  d.active_mol = drawn - 1;
}

void borrowed_mols_push_none(BorrowedMols &list) { list.push_none(); }

void drawer_draw_reaction(Drawer &d, const ChemicalReaction &rxn,
                          bool highlight_by_reactant,
                          const rust::Vec<DrawColour> &colours,
                          bool has_colours,
                          const rust::Vec<std::int32_t> &conf_ids,
                          bool has_conf_ids) {
  check_usable(d);
  if (d.num_mols || d.global_transform) {
    throw ValueErrorException(
        "DrawReaction lays out the reaction against every molecule the "
        "drawer holds; draw it on a new drawer");
  }
  std::vector<RDKit::DrawColour> cols;
  for (const auto &c : colours) {
    cols.push_back(colour_from_data(c));
    check_colour(cols.back(), "highlightColorsReactants");
  }
  const std::vector<int> ids(conf_ids.begin(), conf_ids.end());
  // makeReactionComponents indexes confIds from 0 for the reactants, the
  // products, and the agents alike.
  const std::size_t groups =
      std::max({rxn.getReactants().size(), rxn.getProducts().size(),
                rxn.getAgents().size()});
  if (has_conf_ids && ids.size() < groups) {
    throw ValueErrorException(
        "confIds needs an entry for each reactant, product, and agent "
        "position");
  }
  // RDKit depicts the templates without coordinates (makeReactionDrawMol);
  // the wrappers' layout does so first, on a copy.
  bool depicts = false;
  for (const auto *group :
       {&rxn.getReactants(), &rxn.getProducts(), &rxn.getAgents()}) {
    for (const auto &t : *group) {
      depicts = depicts || !t->getNumConformers();
    }
  }
  const DrawingLocks locks(depicts, comic(d));
  check_draw_area(d.drawer->panelWidth(), d.drawer->panelHeight(),
                  d.drawer->drawOptions());
  check_font(d);
  if (!d.svg && d.drawer->drawOptions().includeMetadata) {
    // As check_metadata, for the reaction's pickle.
    check_picklable_templates(rxn);
  }
  set_up_text_drawer(d);
  const ReactionLayout layout(d, rxn, highlight_by_reactant,
                              has_colours ? &cols : nullptr,
                              has_conf_ids ? &ids : nullptr);
  check_canvas_fits(d, layout.canvas_size());
  if (comic(d)) {
    layout.check_comic(*d.drawer);
  }
  const auto components = static_cast<int>(layout.num_components());
  drawing_call(d, [&] {
    d.drawer->drawReaction(rxn, highlight_by_reactant,
                           has_colours ? &cols : nullptr,
                           has_conf_ids ? &ids : nullptr);
  });
  d.num_mols = components;
  d.active_mol = components - 1;
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
  // The copy has coordinates.
  const DrawingLocks locks(false, comic(d));
  try {
    draw_molecule(d, *copy, data, false);
  } catch (...) {
    opts.prepareMolsBeforeDrawing = prepare;
    throw;
  }
  opts.prepareMolsBeforeDrawing = prepare;
}

// Free drawing: lines, shapes, and text the caller places, in pixels
// (rawCoords) or in the molecule's frame, which MolDraw2D::getDrawCoords
// takes to pixels.
namespace {

// The largest coordinate, in pixels or in the molecule's frame, free
// drawing takes: finite, and within float's range, in which arcPoints
// (MolDraw2DDetails.cpp) computes ellipses (a double outside it converts
// to float undefined).
constexpr double kMaxFreeCoordinate = 1e9;

// Point2D::normalize's zero_tolerance: it throws for a shorter vector.
constexpr double kZeroTolerance = 1e-16;

// The most segments a wavy line takes. getWavyLineSegments rounds an odd
// count up, which makes UINT_MAX no segments, and the drawers then read the
// first.
constexpr std::uint32_t kMaxWavySegments = 100000;

// The largest font size free text is drawn at: DrawTextSVG::drawChar
// converts the size to unsigned.
constexpr double kMaxFreeFontSize = 1e6;

// The largest arc angle, in degrees: drawArc steps through the arc in an
// int count of 5 degree steps.
constexpr double kMaxArcAngle = 1e4;

RDGeom::Point2D to_rdkit(const Point2D &p) { return RDGeom::Point2D(p.x, p.y); }

Point2D from_rdkit(const RDGeom::Point2D &p) { return Point2D{p.x, p.y}; }

// What every free-drawing call needs: a drawer that can draw, on a canvas
// it has set up. Before the first drawing on a flexible canvas sizes it,
// MolDraw2DSVG would write ahead of its header, and MolDraw2DCairo has no
// Cairo context (and FreeType's text drawer draws to it without checking).
void check_canvas(const Drawer &d) {
  check_usable(d);
  if (d.drawer->width() < 0) {
    throw ValueErrorException(
        "a flexible canvas has nothing to draw on until its first drawing "
        "sizes it");
  }
}

// Whether getDrawCoords has a transformation to take molecule coordinates
// by: SetScale's, or the last molecule drawn's (drawMols_[activeMolIdx_],
// which it reads without checking the index).
bool has_transform(const Drawer &d) {
  return d.global_transform || (d.active_mol >= 0 && d.active_mol < d.num_mols);
}

void check_free_coordinate(double v) {
  if (!(std::abs(v) <= kMaxFreeCoordinate)) {
    throw ValueErrorException(
        "free drawing coordinates must be finite and within 1e9 (in the "
        "units given, and in pixels)");
  }
}

// A point as the drawer will draw it, in pixels, checked.
RDGeom::Point2D drawn_point(const Drawer &d, const RDGeom::Point2D &p,
                            bool raw) {
  check_free_coordinate(p.x);
  check_free_coordinate(p.y);
  if (raw) {
    return p;
  }
  if (!has_transform(d)) {
    // RDKit's PRECONDITION in getDrawCoords.
    throw std::runtime_error("no scaling info");
  }
  const auto q = d.drawer->getDrawCoords(p);
  check_free_coordinate(q.x);
  check_free_coordinate(q.y);
  return q;
}

// The box from a to b as the drawer will draw it: getDrawCoords scales and
// moves, so the drawn box has the drawn corners.
void check_drawn_box(const Drawer &d, const RDGeom::Point2D &a,
                     const RDGeom::Point2D &b, bool raw) {
  drawn_point(d, a, raw);
  drawn_point(d, b, raw);
}

// The box around the line from a to b, widened by w on each side.
void check_drawn_line_box(const Drawer &d, const RDGeom::Point2D &a,
                          const RDGeom::Point2D &b, double w, bool raw) {
  check_free_coordinate(w);
  const double lo_x = std::min(a.x, b.x) - std::abs(w);
  const double lo_y = std::min(a.y, b.y) - std::abs(w);
  const double hi_x = std::max(a.x, b.x) + std::abs(w);
  const double hi_y = std::max(a.y, b.y) + std::abs(w);
  check_drawn_box(d, RDGeom::Point2D(lo_x, lo_y), RDGeom::Point2D(hi_x, hi_y),
                  raw);
}

// getDrawLineWidth, which each shape's line takes, scales the width by the
// last molecule drawn (scale(), reading drawMols_[activeMolIdx_]) when
// scaleBondWidth is set.
void check_line_width(const Drawer &d) {
  if (d.drawer->drawOptions().scaleBondWidth &&
      !(d.active_mol >= 0 && d.active_mol < d.num_mols)) {
    // RDKit's PRECONDITION in scale().
    throw std::runtime_error("bad active mol index");
  }
}

void check_comic_line(const Drawer &d, const RDGeom::Point2D &a,
                      const RDGeom::Point2D &b) {
  if (d.drawer->drawOptions().comicMode &&
      !((b - a).length() >= kMinComicLineLength)) {
    throw ValueErrorException(
        "comicMode cannot draw a line shorter than 0.06 (in the units it is "
        "given in) by hand: RDKit's hand-drawn line can never end");
  }
}

// Comic mode draws its lines with std::rand (DrawnMol).
std::unique_lock<std::mutex> comic_lock(const Drawer &d) {
  if (d.drawer->drawOptions().comicMode) {
    return std::unique_lock<std::mutex>(random_mutex());
  }
  return {};
}

RDKit::DrawColour checked_colour(const DrawColour &colour, const char *what) {
  const auto c = colour_from_data(colour);
  check_colour(c, what);
  return c;
}

void check_wavy_line(const Drawer &d, const RDGeom::Point2D &a,
                     const RDGeom::Point2D &b, std::uint32_t n_segments,
                     double vert_offset, bool raw) {
  // RDKit's PRECONDITION, then the wrappers' bound.
  if (n_segments < 2) {
    throw std::runtime_error("too few segments");
  }
  if (n_segments > kMaxWavySegments) {
    throw ValueErrorException(
        "a wavy line takes at most 100000 segments (nSegments)");
  }
  check_free_coordinate(a.x);
  check_free_coordinate(a.y);
  check_free_coordinate(b.x);
  check_free_coordinate(b.y);
  // getWavyLineSegments normalizes the line's perpendicular.
  if (!((b - a).length() >= kZeroTolerance)) {
    throw ValueErrorException("a wavy line needs length");
  }
  check_drawn_line_box(d, a, b, vert_offset, raw);
}

// MolDraw2D::setScale, as RDKit computes it: the transformation it keeps
// (globalDrawTrans_), with the molecule laid out first when given.
void set_scale(Drawer &d, int width, int height, const Point2D &minv_data,
               const Point2D &maxv_data, const RWMol *mol) {
  check_usable(d);
  auto &draw = *d.drawer;
  // RDKit's PRECONDITIONs ask for positive sizes, from which it computes
  // the scale; the wrappers bound them as canvases are.
  if (width < 1) {
    throw std::runtime_error("bad width");
  }
  if (height < 1) {
    throw std::runtime_error("bad height");
  }
  if (width > kMaxCanvasSize || height > kMaxCanvasSize) {
    throw ValueErrorException(
        "SetScale's width and height must be at most 32767");
  }
  // It lays the transformation out at the panel size, which a flexible
  // canvas has none of (and in flexi mode keeps none of).
  if (draw.panelWidth() < 0 || draw.panelHeight() < 0) {
    throw ValueErrorException(
        "SetScale needs panels of fixed size: draw on the flexible canvas "
        "first, outside flexi mode");
  }
  const auto minv = to_rdkit(minv_data);
  const auto maxv = to_rdkit(maxv_data);
  for (double v : {minv.x, minv.y, maxv.x, maxv.y}) {
    if (!(std::abs(v) <= kMaxDepictionCoordinate)) {
      throw ValueErrorException(
          "SetScale's minv and maxv must be finite and within 1e5");
    }
  }
  const auto &opts = draw.drawOptions();
  auto &text = d.text->text_drawer();
  double x_min, x_max, x_range, y_min, y_max, y_range;
  double font_scale = 0;
  bool set_font_scale = false;
  std::unique_ptr<DrawingLocks> locks;
  std::unique_ptr<DrawnMol> drawn;
  if (mol) {
    // RDKit lays the molecule out in a DrawMol of the panel's size, as
    // drawMoleculeWithHighlights does (rings as given, the first
    // conformer), but draws nothing, so no std::rand.
    locks = std::make_unique<DrawingLocks>(lacks_conformers(*mol), false);
    const bool symmetrize = false;
    drawn = std::make_unique<DrawnMol>(*mol, symmetrize);
    // DrawMol's defaults: annotations, and the first conformer.
    const bool annotations = true;
    check_drawing(d, drawn->mol(), -1, "", annotations, symmetrize);
    set_up_text_drawer(d);
    detail::DrawMol probe(drawn->mol(), "", draw.panelWidth(),
                          draw.panelHeight(), opts, text);
    checked_layout(probe);
    x_min = std::min(minv.x, probe.xMin_);
    y_min = std::min(minv.y, -probe.yMax_);
    x_max = std::max(maxv.x, probe.xMax_);
    y_max = std::max(maxv.y, -probe.yMin_);
    font_scale = probe.getFontScale();
    set_font_scale = true;
  } else {
    x_min = minv.x;
    y_min = minv.y;
    x_max = maxv.x;
    y_max = maxv.y;
  }
  x_range = x_max - x_min;
  y_range = y_max - y_min;
  if (x_range < 1.0e-4) {
    x_range = 1.0;
    x_min = -0.5;
  }
  if (y_range < 1.0e-4) {
    y_range = 1.0;
    y_min = -0.5;
  }
  x_min -= opts.padding * x_range;
  x_range *= 1 + 2 * opts.padding;
  x_max = x_min + x_range;
  y_min -= opts.padding * y_range;
  y_range *= 1 + 2 * opts.padding;
  y_max = y_min + y_range;
  const double draw_width = width * (1 - 2 * opts.padding);
  const double draw_height = height * (1 - 2 * opts.padding);
  const double scale =
      std::min(double(draw_width) / x_range, double(draw_height) / y_range);
  if (!set_font_scale) {
    const double ofs = text.fontScale();
    text.setFontScale(scale);
    font_scale = text.fontScale();
    text.setFontScale(ofs, true);
  } else {
    font_scale *= scale / font_scale;
  }
  draw.setScale(width, height, minv, maxv, drawn ? &drawn->mol() : nullptr);
  auto global = std::make_unique<detail::DrawMol>(
      draw.panelWidth(), draw.panelHeight(), opts, text, x_min, x_max, y_min,
      y_max, scale, font_scale);
  global->xOffset_ -= width * opts.padding;
  global->yOffset_ -= height * opts.padding;
  d.global = std::move(global);
  d.global_transform = true;
}

}  // namespace

void drawer_draw_line(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                      bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto a = to_rdkit(cds1);
  const auto b = to_rdkit(cds2);
  check_drawn_box(d, a, b, raw_coords);
  d.drawer->drawLine(a, b, raw_coords);
}

void drawer_draw_arrow(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                       bool as_polygon, double frac, double angle,
                       const DrawColour &colour, bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto col = checked_colour(colour, "color");
  // RDKit's PRECONDITION.
  if (!(angle > 1.0e-6)) {
    throw std::runtime_error(
        "Arrow angle must be positive and greater than 0.0.");
  }
  if (!std::isfinite(frac) || !std::isfinite(angle)) {
    throw ValueErrorException("an arrow's frac and angle must be finite");
  }
  const auto begin = to_rdkit(cds1);
  const auto end = to_rdkit(cds2);
  drawn_point(d, begin, raw_coords);
  drawn_point(d, end, raw_coords);
  // The head, as drawArrow computes it: the line ends short of the tip by
  // the miter, and the head's sides divide by the angle's sine and cosine.
  auto tip = end;
  RDGeom::Point2D p1, p2;
  detail::calcArrowHead(tip, p1, p2, begin, frac, d.drawer->getDrawLineWidth(),
                        angle);
  for (const auto &p : {tip, p1, p2}) {
    drawn_point(d, p, raw_coords);
  }
  check_comic_line(d, begin, tip);
  const auto random = comic_lock(d);
  d.drawer->drawArrow(begin, end, as_polygon, frac, angle, col, raw_coords);
}

void drawer_draw_triangle(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                          const Point2D &cds3, bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto a = to_rdkit(cds1);
  const auto b = to_rdkit(cds2);
  const auto c = to_rdkit(cds3);
  for (const auto &p : {a, b, c}) {
    drawn_point(d, p, raw_coords);
  }
  check_comic_line(d, a, b);
  check_comic_line(d, b, c);
  check_comic_line(d, c, a);
  const auto random = comic_lock(d);
  d.drawer->drawTriangle(a, b, c, raw_coords);
}

void drawer_draw_polygon(Drawer &d, rust::Slice<const Point2D> cds,
                         bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  // Python's helper refuses no points, and the drawers' PRECONDITION fewer
  // than three.
  if (cds.empty()) {
    throw ValueErrorException("cds argument must be non-empty");
  }
  if (cds.size() < 3) {
    throw std::runtime_error("must have at least three points");
  }
  std::vector<RDGeom::Point2D> pts;
  pts.reserve(cds.size());
  for (const auto &p : cds) {
    pts.push_back(to_rdkit(p));
    drawn_point(d, pts.back(), raw_coords);
  }
  d.drawer->drawPolygon(pts, raw_coords);
}

void drawer_draw_ellipse(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                         bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto a = to_rdkit(cds1);
  const auto b = to_rdkit(cds2);
  check_drawn_box(d, a, b, raw_coords);
  d.drawer->drawEllipse(a, b, raw_coords);
}

void drawer_draw_rect(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                      bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto a = to_rdkit(cds1);
  const auto b = to_rdkit(cds2);
  check_drawn_box(d, a, b, raw_coords);
  d.drawer->drawRect(a, b, raw_coords);
}

void drawer_draw_arc(Drawer &d, const Point2D &centre, double radius,
                     double angle1, double angle2, bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  // drawArc takes 1 + int((angle2 - angle1) / 5) steps: below one, its
  // step divides by zero, or it has too few points for its polygon.
  if (!(std::abs(angle1) <= kMaxArcAngle && std::abs(angle2) <= kMaxArcAngle &&
        angle2 >= angle1)) {
    throw ValueErrorException(
        "an arc's angles must be within 10000 degrees, with angle2 no less "
        "than angle1");
  }
  const auto c = to_rdkit(centre);
  check_drawn_line_box(d, c, c, radius, raw_coords);
  d.drawer->drawArc(c, radius, angle1, angle2, raw_coords);
}

void drawer_draw_attachment_line(Drawer &d, const Point2D &cds1,
                                 const Point2D &cds2, const DrawColour &colour,
                                 double len, std::uint32_t n_segments,
                                 bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto col = checked_colour(colour, "color");
  const auto a = to_rdkit(cds1);
  const auto b = to_rdkit(cds2);
  check_free_coordinate(a.x);
  check_free_coordinate(a.y);
  check_free_coordinate(b.x);
  check_free_coordinate(b.y);
  check_free_coordinate(len);
  // drawAttachmentLine draws a wavy line of length len across the line's
  // end, along its perpendicular (calcPerpendicular, NaN for a line of no
  // length), and passes rawCoords where drawWavyLine takes its vertical
  // offset: the wavy line is in the molecule's frame, waving by 1 with
  // rawCoords and not at all without.
  const double bx = a.x - b.x;
  const double by = a.y - b.y;
  const double perp_len = std::sqrt(by * by + bx * bx);
  const RDGeom::Point2D perp(-by / perp_len, bx / perp_len);
  if (!(std::isfinite(perp.x) && std::isfinite(perp.y))) {
    throw ValueErrorException("an attachment line needs a line with length");
  }
  const RDGeom::Point2D p1(b.x - perp.x * len / 2, b.y - perp.y * len / 2);
  const RDGeom::Point2D p2(b.x + perp.x * len / 2, b.y + perp.y * len / 2);
  const double vert_offset = raw_coords ? 1.0 : 0.0;
  const bool wavy_raw = false;
  check_wavy_line(d, p1, p2, n_segments, vert_offset, wavy_raw);
  d.drawer->drawAttachmentLine(a, b, col, len, n_segments, raw_coords);
}

void drawer_draw_wavy_line(Drawer &d, const Point2D &cds1, const Point2D &cds2,
                           const DrawColour &colour1, const DrawColour &colour2,
                           std::uint32_t n_segments, double vert_offset,
                           bool raw_coords) {
  check_canvas(d);
  check_line_width(d);
  const auto col1 = checked_colour(colour1, "color1");
  const auto col2 = checked_colour(colour2, "color2");
  const auto a = to_rdkit(cds1);
  const auto b = to_rdkit(cds2);
  check_wavy_line(d, a, b, n_segments, vert_offset, raw_coords);
  d.drawer->drawWavyLine(a, b, col1, col2, n_segments, vert_offset, raw_coords);
}

void drawer_draw_string(Drawer &d, rust::Str text, const Point2D &pos,
                        std::int32_t align, bool has_align, bool raw_coords) {
  check_canvas(d);
  const std::string str(text);
  // Aligned text is laid out as a run (as notes are), and the rest as an
  // atom label, split into pieces (atomLabelToPieces).
  if (has_align) {
    check_annotation(str, d.freetype, "string");
  } else {
    check_label(str, d.freetype, "string");
  }
  auto talign = detail::TextAlignType::MIDDLE;
  if (has_align) {
    switch (align) {
      case 0:
        talign = detail::TextAlignType::MIDDLE;
        break;
      case 1:
        talign = detail::TextAlignType::START;
        break;
      case 2:
        talign = detail::TextAlignType::END;
        break;
      default:
        throw ValueErrorException("align must be 0, 1, or 2");
    }
  }
  const double size = d.text->text_drawer().fontSize();
  if (!(size > 0 && size <= kMaxFreeFontSize)) {
    throw ValueErrorException(
        "the text drawer's font size is out of range for free text");
  }
  const auto p = to_rdkit(pos);
  drawn_point(d, p, raw_coords);
  if (has_align) {
    d.drawer->drawString(str, p, talign, raw_coords);
  } else {
    d.drawer->drawString(str, p, raw_coords);
  }
}

Point2D drawer_get_draw_coords(const Drawer &d, const Point2D &point) {
  check_usable(d);
  if (!has_transform(d)) {
    throw std::runtime_error("no scaling info");
  }
  return from_rdkit(d.drawer->getDrawCoords(to_rdkit(point)));
}

Point2D drawer_get_atom_draw_coords(const Drawer &d, std::uint32_t atom_index) {
  check_usable(d);
  // RDKit's PRECONDITIONs, for the active molecule; it checks the atom.
  if (!(d.active_mol >= 0 && d.active_mol < d.num_mols)) {
    throw std::runtime_error("bad active mol index");
  }
  if (atom_index > static_cast<std::uint32_t>(INT_MAX)) {
    throw std::runtime_error("bad atom number");
  }
  return from_rdkit(d.drawer->getDrawCoords(static_cast<int>(atom_index)));
}

void drawer_set_scale(Drawer &d, int width, int height, const Point2D &minv,
                      const Point2D &maxv) {
  set_scale(d, width, height, minv, maxv, nullptr);
}

void drawer_set_scale_with_mol(Drawer &d, int width, int height,
                               const Point2D &minv, const Point2D &maxv,
                               const RWMol &mol) {
  set_scale(d, width, height, minv, maxv, &mol);
}

void drawer_draw_molecule_acs1996(Drawer &d, const RWMol &mol,
                                  const DrawMoleculeData &data) {
  // MolDraw2DUtils::drawMolACS1996, step by step: it depicts a copy of a
  // molecule without conformers (its rings as given), sets the drawer's
  // options for the molecule's mean bond length, unchecked (a tiny one
  // makes the scaling factor huge), and draws the copy with drawMolecule,
  // highlighting only the bonds given.
  check_usable(d);
  auto &draw = *d.drawer;
  const DrawingLocks locks(lacks_conformers(mol), comic(d));
  if (draw.width() != -1 || draw.height() != -1) {
    BOOST_LOG(rdWarningLog)
        << "ACS drawing mode works best with a flexiCanvas i.e. a drawer"
        << " created with width and height of -1.  The scale will be fixed,"
        << " and that may not look great with a pre-determined size."
        << std::endl;
  }
  const bool symmetrize = false;
  const DrawnMol depicted(mol, symmetrize);
  auto opts = draw.drawOptions();
  RDKit::MolDraw2DUtils::setACS1996Options(
      opts,
      RDKit::MolDraw2DUtils::meanBondLength(depicted.mol(), data.conf_id));
  check_draw_options(opts);
  // The font file can come from $RDBASE, which need not be UTF-8.
  check_utf8(opts.fontFile, "fontFile");
  draw.drawOptions() = opts;
  draw_molecule(d, depicted.mol(), data, false);
}

std::unique_ptr<std::string> mol_to_acs1996_svg(const RWMol &mol,
                                                const DrawMoleculeData &data) {
  // As molToACS1996SVG in Wrap/rdMolDraw2D.cpp, on a drawer of the
  // wrappers' own.
  auto d = mol_draw_2d_svg_new(-1, -1, -1, -1, false);
  drawer_draw_molecule_acs1996(*d, mol, data);
  svg_finish_drawing(*d);
  return svg_drawing_text(*d);
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

std::unique_ptr<Drawer> mol_draw_2d_cairo_new(int width, int height,
                                              int panel_width, int panel_height,
                                              bool no_freetype) {
  check_canvas_size(width, height, panel_width, panel_height);
#ifdef RDK_BUILD_CAIRO_SUPPORT
  if (width > 0) {
    check_cairo_size(width, height);
  }
  auto cairo = std::make_unique<CairoDrawer>(width, height, panel_width,
                                             panel_height, no_freetype);
  TextDrawerAccess &text = *cairo;
  const bool freetype = uses_freetype(text.text_drawer());
#ifdef RDK_BUILD_FREETYPE_SUPPORT
  // When FreeType fails to start, MolDraw2DCairo falls back to Cairo's own
  // text drawer, and its initDrawing, which sets a flexible canvas up at
  // the first drawing, then takes that drawer for FreeType's (a null
  // dynamic_cast it dereferences).
  if (!no_freetype && !freetype && width < 0) {
    throw std::runtime_error(
        "FreeType could not start, and RDKit's Cairo drawer then cannot draw "
        "on a flexible canvas (draw with noFreetype, or on a canvas of fixed "
        "size)");
  }
#endif
  return std::make_unique<Drawer>(std::move(cairo), text, freetype, false);
#else
  static_cast<void>(no_freetype);
  throw std::runtime_error("RDKit was built without Cairo support");
#endif
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
namespace {
CairoDrawer &cairo_drawer(const Drawer &d) {
  return static_cast<CairoDrawer &>(*d.drawer);
}
}  // namespace
#endif

void cairo_finish_drawing(Drawer &d) {
#ifdef RDK_BUILD_CAIRO_SUPPORT
  cairo_drawer(d).finishDrawing();
#else
  static_cast<void>(d);
  throw std::runtime_error("RDKit was built without Cairo support");
#endif
}

rust::Vec<std::uint8_t> cairo_drawing_text(const Drawer &d) {
#ifdef RDK_BUILD_CAIRO_SUPPORT
  // The drawer's state is RDKit's, behind the pointer; png() sets its
  // metadata aside and back, and the drawer is not Sync.
  return to_bytes(cairo_drawer(d).png());
#else
  static_cast<void>(d);
  throw std::runtime_error("RDKit was built without Cairo support");
#endif
}

void cairo_write_drawing_text(const Drawer &d,
                              rust::Slice<const std::uint8_t> path) {
#ifdef RDK_BUILD_CAIRO_SUPPORT
  // As MolDraw2DCairo::writeDrawingText, with the checked PNG, and an
  // error (where RDKit logs one) for a file it cannot write.
  const auto name = to_path(path);
  const auto png = cairo_drawer(d).png();
  std::ofstream out(name, std::ios_base::binary | std::ios_base::out);
  if (out) {
    out.write(png.data(), static_cast<std::streamsize>(png.size()));
    out.close();
  }
  if (!out) {
    throw RDKit::BadFileException("Failed to write PNG file " + name);
  }
#else
  static_cast<void>(d);
  static_cast<void>(path);
  throw std::runtime_error("RDKit was built without Cairo support");
#endif
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
  // atoms as isAtomCandForChiralH in MolDraw2DUtils.cpp picks them. It
  // blocks the logs to kekulize (DrawingLocks).
  const bool depicts = force_coords || lacks_conformers(mol);
  const DrawingLocks locks(depicts, false);
  if (depicts) {
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
  // As molToSVG in Wrap/rdMolDraw2D.cpp, on a drawer of the wrappers'
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
  {
    const DrawingLocks locks(lacks_conformers(mol), comic(*d));
    draw_molecule(*d, mol, data, true);
  }
  svg_finish_drawing(*d);
  return svg_drawing_text(*d);
}

// Contours: MolDraw2DUtils::contourAndDrawGrid draws a grid's values as
// filled boxes (fillGrid) and contour lines, in the molecule's frame, after
// giving the drawer the grid's scale (setScale), which the wrappers do with
// set_scale, so that they keep its transformation.
namespace {

namespace utils = RDKit::MolDraw2DUtils;

// The largest grid value, and level: within it, the differences
// contourAndDrawGrid takes stay finite (an infinite one makes the fill
// colors NaN, which the SVG drawer converts to int).
constexpr double kMaxContourValue = 1e300;

// The most contour levels RDKit computes for itself (it makes nContours of
// them when given none).
constexpr std::uint32_t kMaxContours = 100000;

// The most grid points ContourAndDrawGaussians computes: it converts the
// count to size_t from double, and allocates and fills them all.
constexpr double kMaxGridPoints = 1e8;

utils::ContourParams contour_params_from_data(const ContourParamsData &data) {
  utils::ContourParams ps;
  ps.setScale = data.set_scale;
  ps.dashNegative = data.dash_negative;
  ps.fillGrid = data.fill_grid;
  ps.gridResolution = data.grid_resolution;
  ps.contourWidth = data.contour_width;
  ps.extraGridPadding = data.extra_grid_padding;
  ps.contourColour = colour_from_data(data.contour_colour);
  ps.colourMap.clear();
  for (const auto &c : data.colour_map) {
    ps.colourMap.push_back(colour_from_data(c));
  }
  ps.drawAsLines = data.draw_as_lines;
  ps.coordScaleForQuantization = data.coord_scale_for_quantization;
  ps.isovalScaleForQuantization = data.isoval_scale_for_quantization;
  ps.useFillThreshold = data.use_fill_threshold;
  ps.fillThreshold = data.fill_threshold;
  ps.fillThresholdIsFraction = data.fill_threshold_is_fraction;
  return ps;
}

const RWMol *optional_mol(BorrowedMols &mols) {
  return mols.mols().empty() ? nullptr
                             : static_cast<const RWMol *>(mols.mols()[0].get());
}

void check_contour_value(double v, const char *what) {
  if (!(std::abs(v) <= kMaxContourValue)) {
    throw ValueErrorException(std::string(what) +
                              " must be finite and within 1e300");
  }
}

// What contourAndDrawGrid trusts of its parameters: its colors (the fill
// interpolates between colourMap's), and the contour width, which it sets
// as the drawer's line width (and leaves there if nothing is contoured).
void check_contour_params(const Drawer &d, const utils::ContourParams &ps) {
  // RDKit's PRECONDITION.
  if (ps.colourMap.size() < 2) {
    throw std::runtime_error("colourMap must have at least two entries");
  }
  for (const auto &c : ps.colourMap) {
    check_colour(c, "colourMap");
  }
  check_colour(ps.contourColour, "contourColour");
  auto opts = d.drawer->drawOptions();
  opts.bondLineWidth = ps.contourWidth;
  check_draw_options(opts);
  for (double v : {ps.coordScaleForQuantization, ps.isovalScaleForQuantization,
                   ps.fillThreshold}) {
    if (!std::isfinite(v)) {
      throw ValueErrorException(
          "the contour parameters must be finite (coordScaleForQuantization, "
          "isovalScaleForQuantization, fillThreshold)");
    }
  }
}

// contourAndDrawGrid, checked: the grid as given, then the drawing as RDKit
// will make it (scaled first when setScale is set).
void draw_contours(Drawer &d, const double *grid, std::size_t n_grid,
                   const std::vector<double> &xs, const std::vector<double> &ys,
                   std::uint32_t n_contours, std::vector<double> levels,
                   utils::ContourParams ps, const RWMol *mol) {
  check_canvas(d);
  // Python's checks, then RDKit's.
  if (xs.empty()) {
    throw ValueErrorException("xcoords argument must be non-empty");
  }
  if (ys.empty()) {
    throw ValueErrorException("ycoords argument must be non-empty");
  }
  if (n_grid != xs.size() * ys.size()) {
    throw ValueErrorException(
        "data must have a value for each x and y coordinate (x-major)");
  }
  check_contour_params(d, ps);
  for (double v : xs) {
    check_free_coordinate(v);
  }
  for (double v : ys) {
    check_free_coordinate(v);
  }
  for (std::size_t i = 0; i < n_grid; ++i) {
    check_contour_value(grid[i], "the grid's values");
  }
  for (double v : levels) {
    check_contour_value(v, "levels");
  }
  if (n_contours > kMaxContours) {
    throw ValueErrorException("nContours must be at most 100000");
  }
  if (ps.setScale) {
    set_scale(d, d.drawer->width(), d.drawer->height(),
              Point2D{xs.front(), ys.front()}, Point2D{xs.back(), ys.back()},
              mol);
    ps.setScale = false;
  }
  // RDKit draws nothing unless it finds the grid's range (with no levels,
  // or to fill), and the range is not empty.
  double min_v = std::numeric_limits<double>::max();
  double max_v = std::numeric_limits<double>::lowest();
  if (levels.empty() || ps.fillGrid) {
    for (std::size_t i = 0; i < n_grid; ++i) {
      min_v = std::min(min_v, grid[i]);
      max_v = std::max(max_v, grid[i]);
    }
  }
  if (max_v > min_v) {
    // What it then draws: boxes and lines in the molecule's frame, within
    // the grid's.
    check_line_width(d);
    const auto [x_lo, x_hi] = std::minmax_element(xs.begin(), xs.end());
    const auto [y_lo, y_hi] = std::minmax_element(ys.begin(), ys.end());
    const bool raw = false;
    check_drawn_box(d, RDGeom::Point2D(*x_lo, *y_lo),
                    RDGeom::Point2D(*x_hi, *y_hi), raw);
    if (n_contours) {
      // It throws for these after filling the grid, leaving the drawer's
      // fill and line width as the fill set them.
      if (!levels.empty() && n_contours > levels.size()) {
        throw ValueErrorException(
            "nContours larger than the size of the level list");
      }
      // Contour's PRECONDITIONs.
      if (xs.size() < 2 || ys.size() < 2) {
        throw std::runtime_error("bad bounds");
      }
    }
  }
  utils::contourAndDrawGrid(*d.drawer, grid, xs, ys, n_contours, levels, ps);
}

}  // namespace

ContourParamsData contour_params_default() {
  const utils::ContourParams ps;
  ContourParamsData data;
  data.set_scale = ps.setScale;
  data.dash_negative = ps.dashNegative;
  data.fill_grid = ps.fillGrid;
  data.grid_resolution = ps.gridResolution;
  data.contour_width = ps.contourWidth;
  data.extra_grid_padding = ps.extraGridPadding;
  data.contour_colour = colour_to_data(ps.contourColour);
  for (const auto &c : ps.colourMap) {
    data.colour_map.push_back(colour_to_data(c));
  }
  data.draw_as_lines = ps.drawAsLines;
  data.coord_scale_for_quantization = ps.coordScaleForQuantization;
  data.isoval_scale_for_quantization = ps.isovalScaleForQuantization;
  data.use_fill_threshold = ps.useFillThreshold;
  data.fill_threshold = ps.fillThreshold;
  data.fill_threshold_is_fraction = ps.fillThresholdIsFraction;
  return data;
}

void contour_and_draw_grid(Drawer &d, rust::Slice<const double> data,
                           rust::Slice<const double> xcoords,
                           rust::Slice<const double> ycoords,
                           std::uint32_t n_contours,
                           rust::Slice<const double> levels,
                           const ContourParamsData &params,
                           BorrowedMols &mols) {
  draw_contours(d, data.data(), data.size(),
                std::vector<double>(xcoords.begin(), xcoords.end()),
                std::vector<double>(ycoords.begin(), ycoords.end()), n_contours,
                std::vector<double>(levels.begin(), levels.end()),
                contour_params_from_data(params), optional_mol(mols));
}

void contour_and_draw_gaussians(Drawer &d, rust::Slice<const Point2D> locs,
                                rust::Slice<const double> heights,
                                rust::Slice<const double> widths,
                                std::uint32_t n_contours,
                                rust::Slice<const double> levels,
                                const ContourParamsData &params,
                                BorrowedMols &mols) {
  // MolDraw2DUtils::contourAndDrawGaussians, step by step: it sets the
  // scale from the gaussians, computes their sum on a grid of the drawing's
  // range, and contours that grid. Each step is checked before RDKit's
  // arithmetic: the grid's size, which it converts from double to size_t
  // and allocates, and the gaussians, whose widths it divides by.
  check_canvas(d);
  auto ps = contour_params_from_data(params);
  // Python's checks, then RDKit's PRECONDITIONs.
  if (locs.empty()) {
    throw ValueErrorException("locs argument must be non-empty");
  }
  if (heights.empty()) {
    throw ValueErrorException("heights argument must be non-empty");
  }
  if (widths.empty()) {
    throw ValueErrorException("widths argument must be non-empty");
  }
  if (locs.size() != heights.size() || locs.size() != widths.size()) {
    throw std::runtime_error("size mismatch");
  }
  std::vector<RDGeom::Point2D> points;
  for (const auto &p : locs) {
    points.push_back(to_rdkit(p));
    if (!(std::abs(p.x) <= kMaxDepictionCoordinate &&
          std::abs(p.y) <= kMaxDepictionCoordinate)) {
      throw ValueErrorException(
          "the gaussians' locations must be finite and within 1e5");
    }
  }
  for (double h : heights) {
    if (!(std::abs(h) <= 1e100)) {
      throw ValueErrorException(
          "the gaussians' heights must be finite and within 1e100");
    }
  }
  for (double w : widths) {
    // Each value divides by the width, and by its square.
    if (!(w >= 1e-100 && w <= 1e100)) {
      throw ValueErrorException(
          "the gaussians' widths must be from 1e-100 to 1e100");
    }
  }
  check_contour_params(d, ps);
  if (!(ps.gridResolution > 0 && std::isfinite(ps.gridResolution))) {
    throw ValueErrorException("gridResolution must be positive and finite");
  }
  auto &draw = *d.drawer;
  if (ps.setScale) {
    if (!(std::abs(ps.extraGridPadding) <= kMaxDepictionCoordinate)) {
      throw ValueErrorException(
          "extraGridPadding must be finite and within 1e5");
    }
    RDGeom::Point2D min_p, max_p;
    min_p.x = min_p.y = std::numeric_limits<double>::max();
    max_p.x = max_p.y = std::numeric_limits<double>::lowest();
    for (const auto &loc : points) {
      min_p.x = std::min(loc.x, min_p.x);
      min_p.y = std::min(loc.y, min_p.y);
      max_p.x = std::max(loc.x, max_p.x);
      max_p.y = std::max(loc.y, max_p.y);
    }
    const auto dims = max_p - min_p;
    const double padding = draw.drawOptions().padding;
    min_p.x -= padding * dims.x;
    min_p.y -= padding * dims.y;
    max_p.x += padding * dims.x;
    max_p.y += padding * dims.y;
    if (ps.extraGridPadding > 0) {
      min_p.x -= ps.extraGridPadding;
      min_p.y -= ps.extraGridPadding;
      max_p.x += ps.extraGridPadding;
      max_p.y += ps.extraGridPadding;
    }
    set_scale(d, draw.width(), draw.height(), from_rdkit(min_p),
              from_rdkit(max_p), optional_mol(mols));
  }
  // range() and minPt() read the transformation, or the last molecule
  // drawn (drawMols_[activeMolIdx_]).
  if (!has_transform(d)) {
    throw std::runtime_error("bad active mol");
  }
  const auto range = draw.range();
  const auto min_pt = draw.minPt();
  const double steps_x = std::ceil(range.x / ps.gridResolution);
  const double steps_y = std::ceil(range.y / ps.gridResolution);
  if (!(steps_x >= 0 && steps_y >= 0 &&
        (steps_x + 1) * (steps_y + 1) <= kMaxGridPoints)) {
    throw ValueErrorException(
        "the gaussians' grid would have more than 1e8 points: increase "
        "gridResolution");
  }
  const std::size_t nx = static_cast<std::size_t>(steps_x) + 1;
  const std::size_t ny = static_cast<std::size_t>(steps_y) + 1;
  std::vector<double> xcoords(nx);
  for (std::size_t i = 0; i < nx; ++i) {
    xcoords[i] = min_pt.x + i * ps.gridResolution;
  }
  std::vector<double> ycoords(ny);
  for (std::size_t i = 0; i < ny; ++i) {
    ycoords[i] = min_pt.y + i * ps.gridResolution;
  }
  std::vector<double> grid(nx * ny);
  for (std::size_t ix = 0; ix < nx; ++ix) {
    const auto px = min_pt.x + ix * ps.gridResolution;
    for (std::size_t iy = 0; iy < ny; ++iy) {
      const auto py = min_pt.y + iy * ps.gridResolution;
      const RDGeom::Point2D pt(px, py);
      double accum = 0.0;
      for (std::size_t ig = 0; ig < points.size(); ++ig) {
        const auto d2 = (pt - points[ig]).lengthSq();
        const auto contrib = heights[ig] / widths[ig] *
                             std::exp(-0.5 * d2 / (widths[ig] * widths[ig]));
        accum += contrib;
      }
      grid[ix * ny + iy] = accum / (2 * M_PI);
    }
  }
  ps.setScale = false;
  draw_contours(d, grid.data(), grid.size(), xcoords, ycoords, n_contours,
                std::vector<double>(levels.begin(), levels.end()), ps, nullptr);
}

double mean_bond_length(const RWMol &mol, int conf_id) {
  // Coordinates are summed as doubles; any conformer is fine.
  return RDKit::MolDraw2DUtils::meanBondLength(mol, conf_id);
}

}  // namespace rdkit_cxx
