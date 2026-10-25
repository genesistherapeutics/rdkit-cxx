// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/draw_checks.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>

#include <RDGeneral/Exceptions.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

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
  check_range(opts.highlightRadius, 0.0, 16.0, "highlightRadius");
  check_range(opts.variableAtomRadius, 0.0, 16.0, "variableAtomRadius");
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
      refuse(what + " needs at least three colour components");
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

}  // namespace rdkit_cxx
