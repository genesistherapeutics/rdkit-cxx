// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/draw_checks.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <Geometry/Transform2D.h>
#include <Geometry/point.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolEnumerator/LinkNode.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

namespace {

[[noreturn]] void refuse(const std::string &message) {
  throw ValueErrorException(message);
}

void check_range(double v, double lo, double hi, const char *what) {
  if (!(v >= lo && v <= hi)) {
    std::ostringstream msg;
    msg << what << " must be from " << lo << " to " << hi << ", not " << v;
    refuse(msg.str());
  }
}

// A value RDKit reads as "off" when it is not positive, else within bounds.
void check_optional_range(double v, double lo, double hi, const char *what) {
  if (std::isfinite(v) && v <= 0) {
    return;
  }
  if (!(v >= lo && v <= hi)) {
    std::ostringstream msg;
    msg << what << " must be at most 0 (off) or from " << lo << " to " << hi
        << ", not " << v;
    refuse(msg.str());
  }
}

// A coordinate DrawMol reads: finite and within the depiction bound, under
// which the drawing's extents, times the largest scale check_draw_options
// allows, stay in int's range.
void check_coordinate(double v, const char *what) {
  if (!(std::abs(v) <= kMaxDepictionCoordinate)) {
    refuse(std::string(what) + " must be finite and within 1e5 to be drawn");
  }
}

enum class TextMode { Normal, Subscript, Superscript };

// RDKit's setStringDrawMode (DrawText.cpp): reads the markup at text[i], if
// any, and moves i to its last character.
bool read_markup(const std::string &text, std::size_t &i, TextMode &mode) {
  if (text.compare(i, 5, "<sub>") == 0) {
    mode = TextMode::Subscript;
    i += 4;
    return true;
  }
  if (text.compare(i, 5, "<sup>") == 0) {
    mode = TextMode::Superscript;
    i += 4;
    return true;
  }
  if (text.compare(i, 6, "</sub>") == 0 || text.compare(i, 6, "</sup>") == 0) {
    mode = TextMode::Normal;
    i += 5;
    return true;
  }
  return false;
}

// One run of text, as the text drawers' getStringRects lay it out (see
// draw_checks.h).
void check_run(const std::string &text, bool freetype, const char *what) {
  auto mode = TextMode::Normal;
  bool normal = false;
  bool printable = false;
  for (std::size_t i = 0; i < text.size(); ++i) {
    if (text[i] == '<' && read_markup(text, i, mode)) {
      continue;
    }
    const auto c = static_cast<unsigned char>(text[i]);
    if (!freetype) {
      // char_widths[(int)c]: a negative index for a byte above 0x7f.
      if (c > 0x7f) {
        refuse(std::string(what) +
               " must be ASCII for a drawer that does not use FreeType");
      }
      // Control characters have zero width, which RDKit divides by.
      if (c >= 0x20 && c < 0x7f) {
        printable = true;
      }
    }
    if (mode == TextMode::Normal) {
      normal = true;
    } else if (mode == TextMode::Subscript && !normal) {
      // adjustStringRectsForSuperSubScript reads rects[-1].
      refuse(std::string(what) +
             " has a subscript before its first normal character");
    }
  }
  if (!normal) {
    // alignString reads rects[0] of an empty run, and a superscript with no
    // normal character reads rects[-1].
    refuse(std::string(what) +
           " has no normal (not super- or subscript) "
           "character to draw");
  }
  if (!freetype && !printable) {
    refuse(std::string(what) + " has no printable character to draw");
  }
}

void check_text_size(const std::string &text, const char *what) {
  if (text.size() > kMaxDrawTextSize) {
    std::ostringstream msg;
    msg << what << " is longer than " << kMaxDrawTextSize
        << " bytes, which RDKit's layout does not keep in range";
    refuse(msg.str());
  }
}

bool is_upper(char c) {
  return std::isupper(static_cast<unsigned char>(c)) != 0;
}

// RDKit's atomLabelToPieces (DrawText.cpp), for orientation N, E, S, W, or C;
// it calls isupper on the label's bytes, so a label it splits must be ASCII.
std::vector<std::string> label_pieces(const std::string &label, char orient,
                                      const char *what) {
  std::vector<std::string> pieces;
  if (label.empty()) {
    return pieces;
  }
  if (label.compare(0, 5, "<lit>") == 0) {
    std::string lit = label.substr(5);
    const auto idx = lit.find("</lit>");
    if (idx != std::string::npos) {
      lit = lit.substr(0, idx);
    }
    pieces.push_back(lit);
    return pieces;
  }
  for (char c : label) {
    if (static_cast<unsigned char>(c) > 0x7f) {
      refuse(std::string(what) +
             " must be ASCII: RDKit passes its bytes to isupper");
    }
  }
  std::string next;
  for (std::size_t i = 0; i < label.size(); ++i) {
    if (label.compare(i, 2, "<s") == 0 || label[i] == ':' ||
        is_upper(label[i])) {
      if (!next.empty()) {
        pieces.push_back(next);
        next.clear();
      }
    }
    next += label[i];
  }
  if (!next.empty()) {
    pieces.push_back(next);
  }
  if (pieces.size() < 2) {
    return pieces;
  }
  if (orient == 'E' || orient == 'S') {
    for (auto &piece : pieces) {
      if (piece == "<sup>+</sup>" || piece == "<sup>-</sup>") {
        const std::string charge = piece;
        piece.clear();
        pieces.push_back(charge);
        break;
      }
    }
  }
  std::vector<std::string> res;
  std::string current;
  bool had_symbol = false;
  for (const auto &piece : pieces) {
    if (piece.empty()) {
      continue;
    }
    if (!is_upper(piece[0])) {
      current += piece;
    } else if (had_symbol) {
      res.push_back(current);
      current = piece;
    } else {
      current += piece;
      had_symbol = true;
    }
  }
  if (!current.empty()) {
    res.push_back(current);
  }
  return res;
}

// RDKit 2024.09's AtomSymbol::adjustColons strips each <...> from a label
// with a colon, in a loop that never ends once the next '<' has no '>' after
// it or a '>' before it. It then reads the rects of the colon's position in
// what is left, of the character before it, and, unless the colon is the
// label's last byte (counting its markup), of the character after it,
// having checked only that the colon is at most one past the last rect
// (RDKit 2026.03 checks both, rdkit/rdkit#9450). There is a rect for each
// character atomLabelToPieces keeps but the <sub>, <sup>, </sub>, and </sup>
// markup getStringRects skips (setStringDrawMode in DrawText.cpp).
void check_colon_markup(const std::string &label, const char *what) {
  if (label.find(':') == std::string::npos) {
    return;
  }
  std::string tmp = label;
  while (true) {
    const auto lt = tmp.find('<');
    if (lt == std::string::npos) {
      break;
    }
    const auto gt = tmp.find('>');
    if (gt == std::string::npos || gt < lt) {
      refuse(std::string(what) +
             " has a colon and a '<' without a '>' after it, which RDKit "
             "2024.09's label layout loops on forever");
    }
    tmp = tmp.substr(0, lt) + tmp.substr(gt + 1);
  }
  const auto colon = tmp.find(':');
  if (colon == std::string::npos) {
    return;
  }
  std::size_t rects = 0;
  for (const auto &piece : label_pieces(label, 'N', what)) {
    for (std::size_t i = 0; i < piece.size(); ++i) {
      if (piece.compare(i, 5, "<sub>") == 0 ||
          piece.compare(i, 5, "<sup>") == 0) {
        i += 4;
      } else if (piece.compare(i, 6, "</sub>") == 0 ||
                 piece.compare(i, 6, "</sup>") == 0) {
        i += 5;
      } else {
        ++rects;
      }
    }
  }
  if (colon <= rects &&
      (colon == rects || (colon + 1 >= rects && colon + 1 < label.size()))) {
    refuse(std::string(what) +
           " has a colon where RDKit 2024.09's label layout reads a character "
           "past the label's last");
  }
}

std::string literal(const std::string &label) {
  return "<lit>" + label + "</lit>";
}

// The text of a property as getPropIfPresent<std::string> reads it, or
// nothing.
template <class T, class K>
bool string_prop(const T &obj, const K &key, std::string &value) {
  return obj.getPropIfPresent(key, value);
}

// DrawMol::getAtomSymbol's label for the atom, checked as it may come out.
void check_atom_symbol(const RDKit::Atom &atom,
                       const RDKit::MolDrawOptions &opts, bool freetype) {
  const auto idx = static_cast<int>(atom.getIdx());
  const auto label = opts.atomLabels.find(idx);
  if (label != opts.atomLabels.end()) {
    if (!label->second.empty()) {
      check_label(literal(label->second), freetype, "atom label (atomLabels)");
    }
    return;
  }
  using RDKit::common_properties::_displayLabel;
  using RDKit::common_properties::_displayLabelW;
  if (atom.hasProp(_displayLabel) || atom.hasProp(_displayLabelW)) {
    std::string lbl;
    std::string lblw;
    string_prop(atom, _displayLabel, lbl);
    string_prop(atom, _displayLabelW, lblw);
    if (lbl.empty()) {
      lbl = lblw;
    }
    // Which one RDKit draws depends on the label's orientation.
    for (const auto *text : {&lbl, &lblw}) {
      if (!text->empty()) {
        check_label(literal(*text), freetype, "atom label (_displayLabel)");
      }
    }
    return;
  }
  if (atom.hasProp(RDKit::common_properties::atomLabel)) {
    const auto text =
        atom.getProp<std::string>(RDKit::common_properties::atomLabel);
    if (!text.empty()) {
      check_label(literal(text), freetype, "atom label (atomLabel)");
    }
    return;
  }
  if (opts.dummiesAreAttachments && atom.getAtomicNum() == 0 &&
      atom.getDegree() == 1) {
    return;
  }
  // The query symbols (A, Q, X, M, with an H for those that allow hydrogen)
  // and the atom lists RDKit writes from element symbols are ASCII letters
  // and punctuation; a complex query shows "?" and its map number.
  static const std::set<std::string> symbol_query_types{"A",  "Q", "QH", "X",
                                                        "XH", "M", "MH"};
  if (opts.useComplexQueryAtomSymbols &&
      symbol_query_types.count(atom.getQueryType())) {
    return;
  }
  if (RDKit::isAtomListQuery(&atom)) {
    return;
  }
  std::string map_num;
  const bool has_map_num = string_prop(atom, "molAtomMapNumber", map_num);
  if (RDKit::isComplexQuery(&atom)) {
    std::string text = "?";
    if (has_map_num) {
      text += ":" + map_num;
    }
    check_label(literal(text), freetype, "atom label (molAtomMapNumber)");
    return;
  }
  // The generated label: isotope, symbol (an element's, or a dummy's
  // dummyLabel), map number, charge, and hydrogens, as superscripts and
  // subscripts. Whether hydrogens are shown depends on what preparation
  // leaves, so both forms are checked.
  std::string pre;
  const unsigned iso = atom.getIsotope();
  if (iso != 0 && ((opts.isotopeLabels && atom.getAtomicNum() != 0) ||
                   (opts.dummyIsotopeLabels && atom.getAtomicNum() == 0))) {
    pre = "<sup>" + std::to_string(iso) + "</sup>";
  }
  std::string post;
  if (has_map_num) {
    post += ":" + map_num;
  }
  if (atom.getFormalCharge() != 0) {
    const int charge = atom.getFormalCharge();
    std::string sign = charge > 0 ? "+" : "-";
    if (std::abs(charge) > 1) {
      sign = std::to_string(std::abs(charge)) + sign;
    }
    post += "<sup>" + sign + "</sup>";
  }
  const auto base = pre + atom.getSymbol() + post;
  check_label(base, freetype, "atom label");
  check_label(base + "H<sub>2</sub>", freetype, "atom label");
}

void check_note(const std::string &note, bool freetype, const char *what) {
  if (!note.empty()) {
    check_annotation(note, freetype, what);
  }
}

// MolDraw2D_detail::getSGroupDataLabels (MolDraw2DDetails.cpp), for a data
// substance group: it joins DATAFIELDS with '|' and drops the last
// character, which an empty list makes pop_back on an empty string; it
// reads FIELDDISP[25], past the end of a shorter string, and positions the
// text by FIELDDISP's coordinates and the molecule's _centroidx and
// _centroidy.
void check_sgroup_data(const RDKit::ROMol &mol, const RDKit::SubstanceGroup &sg,
                       bool freetype) {
  std::string text;
  if (sg.hasProp("DATAFIELDS")) {
    const auto fields = sg.getProp<RDKit::STR_VECT>("DATAFIELDS");
    if (fields.empty()) {
      refuse(
          "a data substance group's DATAFIELDS must not be empty to be "
          "drawn");
    }
    for (const auto &field : fields) {
      text += field + "|";
    }
    text.pop_back();
  }
  if (text.empty()) {
    return;
  }
  check_annotation(text, freetype, "substance group data (DATAFIELDS)");
  std::string field_disp;
  if (!sg.getPropIfPresent("FIELDDISP", field_disp)) {
    return;
  }
  if (field_disp.size() < 25) {
    refuse(
        "a data substance group's FIELDDISP must have at least 25 "
        "characters to be drawn");
  }
  const auto x = RDKit::FileParserUtils::stripSpacesAndCast<double>(
      field_disp.substr(0, 10));
  const auto y = RDKit::FileParserUtils::stripSpacesAndCast<double>(
      field_disp.substr(10, 10));
  check_coordinate(x, "a FIELDDISP coordinate");
  check_coordinate(y, "a FIELDDISP coordinate");
  if (field_disp[25] != 'R' && mol.hasProp("_centroidx")) {
    const auto cx = mol.getProp<double>("_centroidx");
    const auto cy = mol.getProp<double>("_centroidy");
    check_coordinate(cx, "_centroidx");
    check_coordinate(cy, "_centroidy");
  }
}

// DrawMol::extractBrackets: bracket coordinates, and the CONNECT, LABEL, or
// TYPE text it writes by the last bracket (lowercasing TYPE with ::tolower,
// undefined for a byte above 0x7f).
void check_sgroup_brackets(const RDKit::SubstanceGroup &sg, bool freetype,
                           bool annotations) {
  for (const auto &bracket : sg.getBrackets()) {
    for (const auto &p : bracket) {
      check_coordinate(p.x, "a substance group bracket coordinate");
      check_coordinate(p.y, "a substance group bracket coordinate");
    }
  }
  if (!annotations || sg.getBrackets().empty()) {
    return;
  }
  std::string text;
  // DrawAnnotation refuses empty text with an Invariant; refusing it here
  // keeps the drawer usable.
  if (sg.getPropIfPresent("CONNECT", text)) {
    check_annotation(text, freetype, "substance group CONNECT");
  }
  if (sg.getPropIfPresent("LABEL", text)) {
    check_annotation(text, freetype, "substance group LABEL");
  } else if (sg.getPropIfPresent("TYPE", text) && text != "GEN") {
    for (auto &c : text) {
      if (static_cast<unsigned char>(c) > 0x7f) {
        refuse(
            "substance group TYPE must be ASCII: RDKit lowercases it with "
            "tolower");
      }
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    check_annotation(text, freetype, "substance group TYPE");
  }
}

// The atom coordinates DrawMol draws (extractAtomCoords): the conformer,
// centered when centreMoleculesBeforeDrawing (centerMolForDrawing), with y
// inverted and rotated by the rotate option.
std::vector<RDGeom::Point2D> drawn_coords(const RDKit::Conformer &conf,
                                          const RDKit::MolDrawOptions &opts) {
  RDGeom::Point3D centroid;
  if (opts.centreMoleculesBeforeDrawing) {
    centroid = MolTransforms::computeCentroid(conf);
  }
  const double rot = opts.rotate * M_PI / 180.0;
  RDGeom::Transform2D trans;
  trans.SetTransform(RDGeom::Point2D(0.0, 0.0), rot);
  std::vector<RDGeom::Point2D> res;
  for (auto p : conf.getPositions()) {
    p -= centroid;
    RDGeom::Point2D pt{p.x, -p.y};
    if (rot != 0.0) {
      trans.TransformPoint(pt);
    }
    res.push_back(pt);
  }
  return res;
}

// DrawMol::extractLinkNodes places each link node's label off the bond it
// picks, along the bond's direction divided by its length, which a bond of
// zero length (or none picked) makes NaN.
void check_link_nodes(const RDKit::ROMol &mol,
                      const std::vector<RDGeom::Point2D> &coords) {
  if (!mol.hasProp(RDKit::common_properties::molFileLinkNodes)) {
    return;
  }
  const bool strict = false;
  for (const auto &node :
       RDKit::MolEnumerator::utils::getMolLinkNodes(mol, strict)) {
    RDGeom::Point2D label_pt{-1000, -1000};
    RDGeom::Point2D label_perp{0, 0};
    for (const auto &[first, second] : node.bondAtoms) {
      if (first >= coords.size() || second >= coords.size()) {
        throw IndexError(std::max(first, second),
                         static_cast<unsigned>(coords.size()));
      }
      const auto start = coords[first];
      const auto vect = coords[second] - start;
      const auto crossing = start + vect * 0.333;
      RDGeom::Point2D perp{vect.y, -vect.x};
      perp *= 0.333;
      const auto p1 = crossing + perp / 2.;
      const auto p2 = crossing - perp / 2.;
      if (p1.x > label_pt.x) {
        label_pt = p1;
        label_perp = crossing - start;
      }
      if (p2.x > label_pt.x) {
        label_pt = p2;
        label_perp = crossing - start;
      }
    }
    if (!(label_perp.length() > 0)) {
      refuse(
          "a link node's bond has no length (or lies left of x = -1000) "
          "where RDKit draws it, so it cannot place the label");
    }
  }
}

// DrawMol::bondInsideRing walks each ring a double bond is in, from both of
// the bond's atoms, to the ring's next atom, and reads the coordinates of
// atom -1 when there is none. MolDraw2D::drawMolecule recomputes ring
// information that is not symmetrized SSSR before drawing, so only rings it
// keeps are checked (all of them when it draws them as given): every ring
// bond's atoms must have another bond in the ring.
void check_drawn_rings(const RDKit::ROMol &mol, bool symmetrized) {
  const auto *info = mol.getRingInfo();
  if (!info->isInitialized() || (symmetrized && !info->isSymmSssr())) {
    return;
  }
  const auto num_bonds = mol.getNumBonds();
  for (const auto &ring : info->bondRings()) {
    for (int b : ring) {
      if (b < 0 || static_cast<unsigned>(b) >= num_bonds) {
        refuse("the molecule's ring information names a bond it lacks");
      }
    }
    for (int b : ring) {
      const auto *bond = mol.getBondWithIdx(b);
      for (const auto *atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
        bool found = false;
        for (const auto *other : mol.atomBonds(atom)) {
          if (other->getIdx() != bond->getIdx() &&
              std::find(ring.begin(), ring.end(),
                        static_cast<int>(other->getIdx())) != ring.end()) {
            found = true;
            break;
          }
        }
        if (!found) {
          refuse(
              "the molecule's ring information has a ring that is not a "
              "cycle of its bonds, which RDKit's drawing cannot walk");
        }
      }
    }
  }
}

// DrawMol::extractAttachments draws a wavy line across each attachment point
// (a dummy atom with one bond and no label, with dummiesAreAttachments): a
// line one unit of molecule coordinates long whose waves are
// multipleBondOffset times the mean bond length wide, and
// DrawShapeWavyLine::myDraw converts the length over the waves' width to
// int, which a mean bond length of zero makes NaN and a tiny one overflows.
constexpr double kMinAttachmentBondLength = 1e-3;

void check_attachments(const RDKit::ROMol &mol,
                       const RDKit::MolDrawOptions &opts) {
  if (!opts.dummiesAreAttachments || !mol.getNumConformers()) {
    return;
  }
  bool has_attachment = false;
  for (const auto atom : mol.atoms()) {
    if (atom->getAtomicNum() == 0 && atom->getDegree() == 1 &&
        !atom->hasProp(RDKit::common_properties::atomLabel) &&
        !opts.atomLabels.count(static_cast<int>(atom->getIdx()))) {
      has_attachment = true;
      break;
    }
  }
  if (has_attachment && !(RDKit::MolDraw2DUtils::meanBondLength(mol) >=
                          kMinAttachmentBondLength)) {
    refuse(
        "a molecule with attachment points (dummiesAreAttachments) needs "
        "a mean bond length of at least 0.001 to be drawn");
  }
}

// Comic mode draws each line by hand (handdrawnLine in MolDraw2DDetails.cpp),
// shrinking the wiggles until they are small next to the line, which never
// ends for a line of no length: two bonded atoms at one position.
void check_comic_bonds(const RDKit::ROMol &mol,
                       const std::vector<RDGeom::Point2D> &coords) {
  for (const auto bond : mol.bonds()) {
    const auto &a = coords[bond->getBeginAtomIdx()];
    const auto &b = coords[bond->getEndAtomIdx()];
    if (a.x == b.x && a.y == b.y) {
      refuse(
          "comicMode cannot draw a bond whose atoms are at one position: "
          "RDKit's hand-drawn line never ends");
    }
  }
}

}  // namespace

void check_colour(const RDKit::DrawColour &colour, const char *what) {
  for (double v : {colour.r, colour.g, colour.b, colour.a}) {
    if (!(v >= 0.0 && v <= 1.0)) {
      refuse(std::string(what) +
             ": RGBA color values need to be between 0 and 1");
    }
  }
}

void check_draw_options(const RDKit::MolDrawOptions &opts) {
  // Fractions of the canvas. DrawMol shrinks each dimension by 1 - 2 *
  // padding, and the legend takes legendFraction of what is left; at 0.5
  // and 1 nothing is left to draw in, and the scale it computes from the
  // rest, and everything derived from it, divides by zero.
  check_range(opts.padding, 0.0, 0.45, "padding");
  check_range(opts.componentPadding, 0.0, 0.45, "reagentPadding");
  check_range(opts.legendFraction, 0.0, 0.9, "legendFraction");
  // Font sizes in pixels. DrawText converts its font size to unsigned
  // (DrawTextSVG::drawChar), so it must be positive; the bounds keep the
  // sizes DrawMol derives from them, in int, in range.
  check_range(opts.legendFontSize, 1, 512, "legendFontSize");
  if (opts.fixedFontSize != -1) {
    check_range(opts.fixedFontSize, 1, 512, "fixedFontSize");
  }
  check_range(opts.maxFontSize, -1, 512, "maxFontSize");
  check_range(opts.minFontSize, -1, 512, "minFontSize");
  check_optional_range(opts.baseFontSize, 0.05, 4.0, "baseFontSize");
  check_range(opts.annotationFontScale, 0.05, 4.0, "annotationFontScale");
  // Lengths in molecule coordinates, relative to the mean bond length. The
  // wavy lines and dashes RDKit draws divide by multipleBondOffset times it
  // and convert the result to int.
  check_range(opts.multipleBondOffset, 0.01, 1.0, "multipleBondOffset");
  check_range(opts.highlightRadius, 0.0, kMaxHighlightRadius,
              "highlightRadius");
  check_range(opts.variableAtomRadius, 0.0, kMaxHighlightRadius,
              "variableAtomRadius");
  check_range(opts.additionalAtomLabelPadding, 0.0, 1.0,
              "additionalAtomLabelPadding");
  // Line widths in pixels (and int multipliers of them, which
  // getHighlightBondWidth multiplies and DrawMolMCH converts to int).
  check_range(opts.bondLineWidth, 0.0, 64.0, "bondLineWidth");
  check_range(opts.highlightBondWidthMultiplier, 1, 64,
              "highlightBondWidthMultiplier");
  check_range(opts.variableBondWidthMultiplier, 1, 64,
              "variableBondWidthMultiplier");
  // Scales. DrawMol multiplies the drawing's extent by these and converts
  // the product to int.
  check_range(opts.scalingFactor, 0.1, 1000.0, "scalingFactor");
  check_optional_range(opts.fixedScale, 1e-3, 1.0, "fixedScale");
  check_optional_range(opts.fixedBondLength, 0.1, 1e4, "fixedBondLength");
  if (!(std::abs(opts.rotate) <= 1e6)) {
    refuse("rotate must be finite and within 1e6 degrees");
  }
  // extractCloseContacts squares it in int.
  check_range(opts.flagCloseContactsDist, -1, 1000, "flagCloseContactsDist");
  check_colour(opts.highlightColour, "highlightColour");
  check_colour(opts.backgroundColour, "backgroundColour");
  check_colour(opts.queryColour, "queryColour");
  check_colour(opts.legendColour, "legendColour");
  check_colour(opts.symbolColour, "symbolColour");
  check_colour(opts.annotationColour, "annotationColour");
  check_colour(opts.atomNoteColour, "atomNoteColour");
  check_colour(opts.bondNoteColour, "bondNoteColour");
  check_colour(opts.variableAttachmentColour, "variableAttachmentColour");
  for (const auto &[atomic_num, colour] : opts.atomColourPalette) {
    check_colour(colour, "atom palette");
  }
  for (const auto &colour : opts.highlightColourPalette) {
    check_colour(colour, "highlight palette");
  }
  // FreeType opens the font file by its C string, which a NUL would cut
  // short to another path. The file itself is the user's choice, read by
  // FreeType, as a file a user names to any reader is (CONVENTIONS.md U2);
  // what RDKit makes of its glyphs is checked when a drawer loads it.
  if (opts.fontFile.find('\0') != std::string::npos) {
    refuse("fontFile must not contain a NUL character");
  }
  for (const auto &[idx, label] : opts.atomLabels) {
    check_text_size(label, "atomLabels");
  }
}

void check_json_options(const std::string &json) {
  if (json.empty()) {
    return;
  }
  boost::property_tree::ptree pt;
  std::istringstream ss(json);
  boost::property_tree::read_json(ss, pt);
  static const char *const kColourKeys[] = {
      "highlightColour", "backgroundColour", "queryColour",
      "legendColour",    "symbolColour",     "annotationColour",
      "atomNoteColour",  "bondNoteColour",   "variableAttachmentColour"};
  auto check_rgba = [](const boost::property_tree::ptree &node,
                       const std::string &what) {
    if (node.size() < 3) {
      refuse(what + " needs at least three color components");
    }
  };
  // RDKit reads whichever child find() returns for a repeated key, so every
  // one with the name is checked.
  for (const auto &[key, child] : pt) {
    for (const char *name : kColourKeys) {
      if (key == name) {
        check_rgba(child, key);
      }
    }
    if (key == "atomColourPalette") {
      for (const auto &[atomic_num, colour] : child) {
        check_rgba(colour, "atomColourPalette entry " + atomic_num);
      }
    }
  }
}

void check_utf8(const std::string &text, const char *what) {
  const auto *p = reinterpret_cast<const unsigned char *>(text.data());
  const std::size_t n = text.size();
  std::size_t i = 0;
  while (i < n) {
    const unsigned char c = p[i];
    std::size_t len = 0;
    std::uint32_t cp = 0;
    if (c < 0x80) {
      ++i;
      continue;
    } else if ((c & 0xE0) == 0xC0) {
      len = 2;
      cp = c & 0x1F;
    } else if ((c & 0xF0) == 0xE0) {
      len = 3;
      cp = c & 0x0F;
    } else if ((c & 0xF8) == 0xF0) {
      len = 4;
      cp = c & 0x07;
    } else {
      refuse(std::string("invalid UTF-8 in ") + what);
    }
    if (i + len > n) {
      refuse(std::string("invalid UTF-8 in ") + what);
    }
    for (std::size_t k = 1; k < len; ++k) {
      if ((p[i + k] & 0xC0) != 0x80) {
        refuse(std::string("invalid UTF-8 in ") + what);
      }
      cp = (cp << 6) | (p[i + k] & 0x3F);
    }
    // Overlong forms, surrogates, and code points past U+10FFFF.
    const std::uint32_t min_cp[] = {0, 0, 0x80, 0x800, 0x10000};
    if (cp < min_cp[len] || (cp >= 0xD800 && cp <= 0xDFFF) || cp > 0x10FFFF) {
      refuse(std::string("invalid UTF-8 in ") + what);
    }
    i += len;
  }
}

void check_canvas_size(int width, int height, int panel_width,
                       int panel_height) {
  auto check_dim = [](int v, const char *what) {
    if (v != -1 && !(v >= 1 && v <= kMaxCanvasSize)) {
      std::ostringstream msg;
      msg << what << " must be -1 (flexible) or from 1 to " << kMaxCanvasSize
          << ", not " << v;
      refuse(msg.str());
    }
  };
  check_dim(width, "width");
  check_dim(height, "height");
  check_dim(panel_width, "panelWidth");
  check_dim(panel_height, "panelHeight");
  if ((width == -1) != (height == -1)) {
    refuse(
        "width and height must both be -1 (a flexible canvas) or both "
        "fixed");
  }
  if (width == -1 && (panel_width != -1 || panel_height != -1)) {
    refuse("a flexible canvas takes flexible panels (-1)");
  }
}

void check_draw_area(int width, int height, const RDKit::MolDrawOptions &opts) {
  for (int v : {width, height}) {
    if (v > 0 && !(v * (1 - 2 * opts.padding) >= kMinDrawSize)) {
      std::ostringstream msg;
      msg << "the canvas or panel (" << width << " by " << height
          << ") leaves less than a pixel to draw in once padding is taken";
      refuse(msg.str());
    }
  }
}

void check_annotation(const std::string &text, bool freetype,
                      const char *what) {
  check_text_size(text, what);
  check_run(text, freetype, what);
}

void check_label(const std::string &label, bool freetype, const char *what) {
  check_text_size(label, what);
  check_colon_markup(label, what);
  // N, S, and C draw each piece as a run; E draws the pieces joined, and W
  // joined in reverse. E and S move a charge to the end.
  for (char orient : {'N', 'S'}) {
    for (const auto &piece : label_pieces(label, orient, what)) {
      check_run(piece, freetype, what);
    }
  }
  std::string joined;
  for (const auto &piece : label_pieces(label, 'E', what)) {
    joined += piece;
  }
  check_run(joined, freetype, what);
  const auto pieces = label_pieces(label, 'W', what);
  joined.clear();
  for (auto it = pieces.rbegin(); it != pieces.rend(); ++it) {
    joined += *it;
  }
  check_run(joined, freetype, what);
}

void check_legend(const std::string &legend, const RDKit::MolDrawOptions &opts,
                  bool freetype) {
  if (legend.empty()) {
    return;
  }
  // RDKit 2024.09's DrawMol::extractLegend splits the legend at newlines,
  // dropping empty pieces, and then places the last piece, read from an
  // empty list for a legend of only newlines.
  if (legend.find_first_not_of('\n') == std::string::npos) {
    throw ValueErrorException(
        "a legend of only newlines, which RDKit 2024.09 reads past its pieces");
  }
  check_text_size(legend, "legend");
  std::vector<std::string> lines;
  boost::split(lines, legend, boost::is_any_of("\n"));
  for (const auto &line : lines) {
    if (!line.empty()) {
      check_run(line, freetype, "legend");
    }
  }
}

void check_drawn_mol(const RDKit::ROMol &mol, int conf_id,
                     const RDKit::MolDrawOptions &opts, bool freetype,
                     bool annotations, bool symmetrized) {
  // Coordinates: DrawMol draws conf_id, and positions data substance groups
  // from the first conformer (getSGroupDataLabels).
  std::vector<RDGeom::Point2D> coords;
  if (mol.getNumConformers()) {
    const auto *conf = find_conformer(mol, conf_id);
    if (!conf) {
      throw RDKit::ConformerException("Bad Conformer Id");
    }
    for (const auto *c : {conf, &mol.getConformer()}) {
      for (const auto &p : c->getPositions()) {
        check_coordinate(p.x, "an atom coordinate");
        check_coordinate(p.y, "an atom coordinate");
        check_coordinate(p.z, "an atom coordinate");
      }
    }
    coords = drawn_coords(*conf, opts);
  }
  // extractRegions reads the coordinates of each atom of a region with more
  // than one.
  for (const auto &region : opts.atomRegions) {
    if (region.size() > 1) {
      for (int idx : region) {
        if (idx < 0 || static_cast<unsigned>(idx) >= mol.getNumAtoms()) {
          throw IndexError(idx, mol.getNumAtoms());
        }
      }
    }
  }
  check_drawn_rings(mol, symmetrized);
  check_attachments(mol, opts);
  if (!coords.empty() && annotations) {
    check_link_nodes(mol, coords);
  }
  if (!coords.empty() && opts.comicMode) {
    check_comic_bonds(mol, coords);
  }
  // An atom of two bonds with a neighbor at its position is left to RDKit:
  // 2024.09's isLinearAtom normalizes the direction to each neighbor and
  // throws for it (RDKit 2026.03 skips it and reads a bond type it never set).

  // Text.
  using RDKit::common_properties::atomNote;
  using RDKit::common_properties::bondNote;
  const bool atom_indices =
      opts.addAtomIndices && !mol.hasProp("_atomIndicesAdded");
  const bool bond_indices =
      opts.addBondIndices && !mol.hasProp("_bondIndicesAdded");
  for (const auto atom : mol.atoms()) {
    if (!opts.noAtomLabels) {
      check_atom_symbol(*atom, opts, freetype);
    }
    std::string note;
    const bool has_note = string_prop(*atom, atomNote, note);
    if (atom_indices) {
      note = std::to_string(atom->getIdx()) + (has_note ? "," + note : "");
    }
    check_note(note, freetype, "atom note");
    std::string cip;
    if (opts.addStereoAnnotation &&
        string_prop(*atom, RDKit::common_properties::_CIPCode, cip)) {
      check_annotation("(" + cip + ")", freetype, "atom CIP label");
    }
  }
  for (const auto bond : mol.bonds()) {
    std::string note;
    const bool has_note = string_prop(*bond, bondNote, note);
    if (bond_indices) {
      note = std::to_string(bond->getIdx()) + (has_note ? "," + note : "");
    }
    check_note(note, freetype, "bond note");
    std::string cip;
    if (opts.addStereoAnnotation &&
        string_prop(*bond, RDKit::common_properties::_CIPCode, cip)) {
      check_annotation("(" + cip + ")", freetype, "bond CIP label");
    }
  }
  std::string mol_note;
  if (string_prop(mol, RDKit::common_properties::molNote, mol_note)) {
    check_note(mol_note, freetype, "molecule note (molNote)");
  }
  if (!opts.addStereoAnnotation) {
    // DrawMol labels each atom of a stereo group with the group's type, and
    // numbers the OR and AND groups (extractStereoGroups).
    for (const auto &group : mol.getStereoGroups()) {
      switch (group.getGroupType()) {
        case RDKit::StereoGroupType::STEREO_ABSOLUTE:
          check_annotation("abs", freetype, "stereo group label");
          break;
        case RDKit::StereoGroupType::STEREO_OR:
          check_annotation("or1", freetype, "stereo group label");
          break;
        case RDKit::StereoGroupType::STEREO_AND:
          check_annotation("and1", freetype, "stereo group label");
          break;
      }
    }
  }

  for (const auto &sg : RDKit::getSubstanceGroups(mol)) {
    check_sgroup_brackets(sg, freetype, annotations);
    std::string type;
    if (annotations && mol.getNumConformers() &&
        sg.getPropIfPresent("TYPE", type) && type == "DAT") {
      check_sgroup_data(mol, sg, freetype);
    }
  }
}

void check_highlight_atoms(const RDKit::ROMol &mol,
                           const std::vector<int> &atoms) {
  for (int idx : atoms) {
    check_index(static_cast<unsigned>(idx), mol.getNumAtoms());
  }
}

void check_highlight_bonds(const RDKit::ROMol &mol,
                           const std::vector<int> &bonds) {
  for (int idx : bonds) {
    check_index(static_cast<unsigned>(idx), mol.getNumBonds());
  }
}

void check_highlight_colours(const std::map<int, RDKit::DrawColour> &colours) {
  for (const auto &[idx, colour] : colours) {
    check_colour(colour, "highlight color");
  }
}

void check_highlight_radii(const std::map<int, double> &radii) {
  for (const auto &[idx, radius] : radii) {
    check_range(radius, 0.0, kMaxHighlightRadius, "highlight radius");
  }
}

}  // namespace rdkit_cxx
