// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The drawers: MolDraw2DSVG (and MolDraw2DCairo), and the functions that
// draw with them.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
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
#include <RDGeneral/Exceptions.h>
// The text drawer's and DrawMol's interfaces, which RDKit does not install;
// wrapper/vendor holds RDKit's own copies, from the release the crate
// supports.
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawMolMCHCircleAndLine.h>
#include <GraphMol/MolDraw2D/DrawMolMCHLasso.h>
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
// rings (drawMoleculeWithHighlights does not: symmetrize); this does the
// same on a copy, under the depiction lock, so the
// coordinates can be checked before RDKit draws, and DrawMol depicts no
// more. Depiction matches ring templates to the molecule (control_c_lock.h).
// Comic mode draws its lines with std::rand (handdrawnLine in
// MolDraw2DDetails.cpp), a global with no lock on some platforms, so the
// drawing holds random_mutex() too.
class DrawnMol {
 public:
  DrawnMol(const RWMol &mol, const RDKit::MolDrawOptions &opts,
           bool symmetrize = true, bool lock_random = true)
      : d_mol(&mol) {
    if (!mol.getNumConformers()) {
      auto depict = lock_depiction();
      std::lock_guard<std::mutex> control_c(control_c_mutex());
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
    if (opts.comicMode && lock_random) {
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

// As check_layout, for drawMoleculeWithHighlights's DrawMolMCH, which RDKit
// builds from the molecule as given (it does not symmetrize its rings).
void check_mch_layout(Drawer &d, const RWMol &mol, const std::string &legend,
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
      // It takes each reactant's colour modulo their number.
      if (cols->empty()) {
        throw ValueErrorException(
            "highlightByReactant needs at least one colour");
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
    if (d.height() == -1) {
      int height = -1;
      for (const auto *dm : all()) {
        height = std::max(height, dm->height_);
      }
      check_pixels(height + height * 2.0 * opts.padding);
    }
    const std::int64_t arrow_mult = 2;
    std::int64_t plus_width = d_plus_width;
    std::int64_t width_to_use =
        d.panelWidth() == -1 ? d.width() : d.panelWidth();
    if (width_to_use == -1) {
      width_to_use = reaction_width(d_reagents, d_products, d_agents,
                                    opts.padding, arrow_mult, plus_width);
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
  check_legend(legend, opts, d.freetype);
  check_drawn_mol(mol, conf_id, opts, d.freetype, annotations, symmetrized);
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
  check_drawing(d, drawn.mol(), data.conf_id, legend,
                d.drawer->supportsAnnotations(), true);
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

// DrawMoleculeWithHighlights's maps as RDKit takes them, checked: indices
// (DrawMolMCHLasso sets bits of an atom's index without checking it),
// colours, radii, and line width multipliers, which getHighlightBondWidth
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
      check_colour(res.back(), "highlight colour");
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
          "multi-colour highlight radii must be above 0 and at most 16");
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
        "multi-colour highlights need a positive highlightRadius");
  }
  return h;
}

// DrawMolMCHCircleAndLine draws a highlighted bond's colours along the
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
          "multi-colour highlight of it never ends");
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
  check_drawing(d, drawn.mol(), data.conf_id, legend,
                d.drawer->supportsAnnotations(), true);
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

void drawer_draw_molecule_with_highlights(
    Drawer &d, const RWMol &mol, const DrawMoleculeWithHighlightsData &data) {
  check_usable(d);
  const auto &opts = d.drawer->drawOptions();
  auto h = read_multi_highlights(mol, data, opts);
  const std::string legend(data.legend);
  const bool symmetrize = false;
  DrawnMol drawn(mol, opts, symmetrize);
  check_mch_bonds(drawn.mol(), h, opts);
  // DrawMolMCH passes confId where DrawMol takes includeAnnotations, so
  // RDKit draws the first conformer, with annotations unless confId is 0.
  check_drawing(d, drawn.mol(), -1, legend, data.conf_id != 0, symmetrize);
  check_mch_layout(d, drawn.mol(), legend, h, data.conf_id);
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
// panel wider than the canvas.
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

  struct Entry {
    std::unique_ptr<DrawnMol> drawn;
    Highlights h;
    std::vector<int> bonds;
    std::string legend;
    int conf_id;
  };
  std::vector<Entry> entries(n);
  for (std::size_t i = 0; i < n; ++i) {
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
    // RDKit draws each molecule as given, without symmetrizing its rings;
    // the random lock waits for the depictions to finish (U4's order).
    const bool symmetrize = false;
    const bool lock_random = false;
    e.drawn = std::make_unique<DrawnMol>(mol, opts, symmetrize, lock_random);
    check_drawing(d, e.drawn->mol(), e.conf_id, e.legend,
                  draw.supportsAnnotations(), symmetrize);
    if (e.h.has_bonds) {
      e.bonds = e.h.bonds;
    } else if (opts.continuousHighlight && e.h.has_atoms) {
      detail::getBondHighlightsForAtoms(e.drawn->mol(), e.h.atoms, e.bonds);
    }
  }
  std::unique_lock<std::mutex> random;
  if (opts.comicMode) {
    random = std::unique_lock<std::mutex>(random_mutex());
  }
  set_up_text_drawer(d);
  std::vector<std::unique_ptr<detail::DrawMol>> probes;
  std::size_t min_scale = 0;
  std::size_t min_font_scale = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto &e = entries[i];
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
    if (probes.back()->getScale() < probes[min_scale]->getScale()) {
      min_scale = i;
    }
    if (probes.back()->getFontScale() <
        probes[min_font_scale]->getFontScale()) {
      min_font_scale = i;
    }
  }
  if (opts.drawMolsSameScale) {
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
    ptrs.push_back(const_cast<RWMol *>(&e.drawn->mol()));
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
  d.num_mols = static_cast<int>(n);
  d.active_mol = static_cast<int>(n) - 1;
}

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
  // RDKit depicts the templates without coordinates (makeReactionDrawMol),
  // matching ring templates to them (control_c_lock.h); the wrappers' layout
  // does so first, on a copy.
  bool depicts = false;
  for (const auto *group :
       {&rxn.getReactants(), &rxn.getProducts(), &rxn.getAgents()}) {
    for (const auto &t : *group) {
      depicts = depicts || !t->getNumConformers();
    }
  }
  std::unique_lock<std::mutex> depict, control_c, random;
  if (depicts) {
    depict = lock_depiction();
    control_c = std::unique_lock<std::mutex>(control_c_mutex());
  }
  if (d.drawer->drawOptions().comicMode) {
    random = std::unique_lock<std::mutex>(random_mutex());
  }
  check_draw_area(d.drawer->panelWidth(), d.drawer->panelHeight(),
                  d.drawer->drawOptions());
  check_font(d);
  set_up_text_drawer(d);
  const ReactionLayout layout(d, rxn, highlight_by_reactant,
                              has_colours ? &cols : nullptr,
                              has_conf_ids ? &ids : nullptr);
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
