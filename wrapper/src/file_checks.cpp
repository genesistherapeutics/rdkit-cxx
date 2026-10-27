// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/file_checks.h"

#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace rdkit_cxx {

std::string read_file(const std::string &name) {
  std::ifstream in(name, std::ios_base::binary);
  if (!in || in.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  std::ostringstream bytes;
  bytes << in.rdbuf();
  if (in.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  return std::move(bytes).str();
}

std::string read_supplier_file(const std::string &name) {
  auto bytes = read_file(name);
  if (bytes.empty()) {
    throw RDKit::BadFileException("Invalid input file " + name);
  }
  return bytes;
}

void check_precision(unsigned precision) {
  if (precision > kMaxPrecision) {
    throw ValueErrorException("precision must be at most " +
                              std::to_string(kMaxPrecision));
  }
}

void check_debug_parse(int debug_parse) {
  if (debug_parse != 0) {
    throw ValueErrorException(
        "debugParse must be 0: RDKit stores it in a global that every SMILES "
        "and SMARTS parse reads without a lock");
  }
}

namespace {

void check_bond_refs(const boost::property_tree::ptree &tree) {
  for (const auto &[key, child] : tree) {
    if (key == "bondArray") {
      for (const auto &[bkey, bond] : child) {
        if (bkey != "bond") {
          continue;
        }
        std::vector<std::string> refs;
        boost::algorithm::split(
            refs, bond.get<std::string>("<xmlattr>.atomRefs2", ""),
            boost::algorithm::is_space());
        if (refs.size() < 2) {
          throw RDKit::FileParseException(
              "atomRefs2 must contain two atom refs that must appear in the "
              "atoms array in MRV file");
        }
      }
    }
    check_bond_refs(child);
  }
}

}  // namespace

namespace {

// An iterator RDKit's CX parser moves, as an index into the text: at() is
// what `*first` reads, the terminator at end(); past() is set once RDKit
// would form an iterator beyond end() + 1 or dereference end() + 1.
class CxCursor {
 public:
  // The index `first` stands at.
  std::size_t i;

  CxCursor(std::string_view text, std::size_t pos) : i(pos), d_text(text) {}

  std::size_t n() const { return d_text.size(); }
  char at() {
    if (i < n()) {
      return d_text[i];
    }
    if (i > n()) {
      d_past = true;
    }
    return '\0';
  }
  char at(std::size_t j) const { return j < n() ? d_text[j] : '\0'; }
  bool digit() { return i <= n() && at() >= '0' && at() <= '9'; }
  void step() {
    ++i;
    if (i > n() + 1) {
      d_past = true;
    }
  }
  bool past() const { return d_past; }

  // As read_text_to: true if it would throw (an escape "&#" not closed by
  // ";"), which ends the parse. `produced` is whether the text read is
  // nonempty.
  bool read_text_to(std::string_view delims, bool &produced) {
    produced = false;
    std::size_t start = i;
    while (i <= n() && delims.find(at()) == std::string_view::npos) {
      if (at() == '&' && n() - i > 2 && at(i + 1) == '#') {
        if (start != i) {
          produced = true;
        }
        std::size_t next = i + 2;
        while (next != n() && at(next) >= '0' && at(next) <= '9') {
          ++next;
        }
        if (next == n() || at(next) != ';') {
          return true;
        }
        if (next > i + 2) {
          produced = true;
        }
        i = next + 1;
        start = i;
      } else {
        step();
      }
    }
    if (start != i) {
      produced = true;
    }
    return false;
  }

 private:
  std::string_view d_text;
  bool d_past = false;
};

// parse_data_sgroup after "SgD:".
bool data_sgroup_past_end(CxCursor c) {
  // read_int_list
  while (true) {
    while (c.digit()) {
      c.step();
    }
    if (c.i >= c.n() || c.at() != ',') {
      break;
    }
    c.step();
  }
  c.step();
  // parse_data_sgroup_attr, for each of the five fields
  for (int field = 0; field < 5 && !c.past(); ++field) {
    if (c.i != c.n() && c.at() != '|' && !c.past()) {
      bool produced;
      if (c.read_text_to(":", produced)) {
        return false;
      }
      c.step();
    }
  }
  if (!c.past() && c.i < c.n() && c.at() == '(') {
    bool produced;
    if (c.read_text_to(")", produced)) {
      return false;
    }
    c.step();
  }
  return c.past();
}

// parse_atom_props after "atomProp:".
bool atom_props_past_end(CxCursor c) {
  if (c.i >= c.n()) {
    return false;
  }
  auto separator = [&c] { return c.at() == '|' || c.at() == ','; };
  while (c.i <= c.n() && !separator()) {
    std::size_t start = c.i;
    while (c.digit()) {
      c.step();
    }
    if (c.i > start) {
      if (c.i >= c.n() || c.at() != '.') {
        return false;
      }
      c.step();
      bool produced;
      if (c.read_text_to(".", produced)) {
        return false;
      }
      if (produced) {
        if (c.i >= c.n() || c.at() != '.') {
          return false;
        }
        c.step();
        if (c.read_text_to(":|,", produced)) {
          return false;
        }
      }
    }
    if (c.i <= c.n() && !separator()) {
      c.step();
    }
  }
  if (c.i <= c.n()) {
    return false;  // RDKit's own checks return first
  }
  return true;  // *first at end() + 1
}

// parse_polymer_sgroup after "Sg:": its type and the ':' after it.
bool polymer_sgroup_past_end(CxCursor c) {
  bool produced;
  if (c.read_text_to(":", produced)) {
    return false;
  }
  c.step();
  return c.past();
}

}  // namespace

std::string cx_extension_problem(std::string_view cx) {
  const std::size_t n = cx.size();
  for (std::size_t p = 0; p < n; ++p) {
    if (cx[p] == 'S' && p + 3 < n && cx.compare(p, 4, "SgD:") == 0) {
      if (data_sgroup_past_end(CxCursor(cx, p + 4))) {
        return "a CXSMILES data S-group (SgD) has fewer than five ':'-ended "
               "fields; RDKit 2024.03 reads past the end of the extension "
               "text";
      }
    } else if (cx[p] == 'S' && p + 2 < n && cx.compare(p, 3, "Sg:") == 0) {
      if (polymer_sgroup_past_end(CxCursor(cx, p + 3))) {
        return "a CXSMILES polymer S-group (Sg) has no ':' after its type; "
               "RDKit 2024.03 steps past the end of the extension text";
      }
    } else if (cx[p] == 'a' && n - p > 9 &&
               cx.compare(p, 9, "atomProp:") == 0) {
      if (atom_props_past_end(CxCursor(cx, p + 9))) {
        return "a CXSMILES atom property (atomProp) runs to the end of the "
               "extension text; RDKit 2024.03 reads past it";
      }
    }
  }
  return {};
}

std::string smiles_cx_problem(const std::string &text) {
  auto sidx = text.find_first_of(" \t");
  if (sidx == std::string::npos || sidx == 0) {
    return {};
  }
  auto cx = boost::trim_copy(text.substr(sidx));
  if (cx.empty() || cx[0] != '|') {
    return {};
  }
  return cx_extension_problem(cx);
}

void check_smiles_cx(const std::string &text, bool allow_cxsmiles) {
  if (!allow_cxsmiles) {
    return;
  }
  auto why = smiles_cx_problem(text);
  if (!why.empty()) {
    throw ValueErrorException(why);
  }
}

void check_reaction_cx(const std::string &text, bool allow_cxsmiles) {
  if (!allow_cxsmiles) {
    return;
  }
  auto sidx = text.find_first_of('|');
  if (sidx == std::string::npos || sidx == 0) {
    return;
  }
  auto why = cx_extension_problem(boost::trim_copy(text.substr(sidx)));
  if (!why.empty()) {
    throw ValueErrorException(why);
  }
}

void check_mrv(const std::string &text) {
  std::istringstream in(text);
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(in, tree);
  check_bond_refs(tree);
}

}  // namespace rdkit_cxx
