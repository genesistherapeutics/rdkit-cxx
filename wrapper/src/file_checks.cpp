// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/file_checks.h"

#include <algorithm>
#include <climits>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <GraphMol/ROMol.h>
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

// RDKit 2024.03's Marvin reader puts each atom with an mrvStereoGroup into
// its group by the index of the first atom with its id (getAtomIndex in
// MarvinDefs.cpp), so two atoms with one id list an atom twice, which RDKit
// 2024.03's StereoGroup keeps and a removal leaves dangling
// (check_stereo_group_removals in mol_checks.h): say, when the reader
// removes the hydrogens. Ids in an atom array must therefore differ.
void check_atom_ids(const boost::property_tree::ptree &tree) {
  for (const auto &[key, child] : tree) {
    if (key == "atomArray") {
      std::set<std::string> ids;
      auto add = [&](const std::string &id) {
        if (!ids.insert(id).second) {
          throw RDKit::FileParseException(
              "atom id " + id +
              " appears twice in an atom array of the MRV file");
        }
      };
      std::vector<std::string> compact;
      boost::algorithm::split(
          compact, child.get<std::string>("<xmlattr>.atomID", ""),
          boost::algorithm::is_space(), boost::algorithm::token_compress_on);
      for (const auto &id : compact) {
        if (!id.empty()) {
          add(id);
        }
      }
      for (const auto &[akey, atom] : child) {
        if (akey == "atom") {
          add(atom.get<std::string>("<xmlattr>.id", ""));
        }
      }
    }
    check_atom_ids(child);
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

// RDKit 2024.03's parse_enhanced_stereo (CXSmilesOps.cpp) adds each atom
// index listed after "a:", "o<n>:", or "&<n>:" to the group the type and
// number name, merging sections by 10 * number + type, so a group can list
// an atom twice, which RDKit 2024.03's StereoGroup keeps and a removal
// leaves dangling (check_stereo_group_removals in mol_checks.h): say, when
// the parser removes the hydrogens ("[H][C@](C)(F)Cl |&1:0,0,1|" crashes
// Python). Every place a section may start is read as one, which may refuse
// a few extensions where such text is inside another section.
bool cx_stereo_group_repeats(std::string_view ext) {
  auto digit = [](char c) { return c >= '0' && c <= '9'; };
  std::map<unsigned int, std::set<unsigned long long>> groups;
  for (std::size_t p = 0; p < ext.size(); ++p) {
    unsigned int type;
    if (ext[p] == 'a') {
      type = 0;
    } else if (ext[p] == 'o') {
      type = 1;
    } else if (ext[p] == '&' && p + 1 < ext.size() && ext[p + 1] != '#') {
      type = 2;
    } else {
      continue;
    }
    std::size_t q = p + 1;
    unsigned int id = 0;
    if (type != 0) {
      for (; q < ext.size() && digit(ext[q]); ++q) {
        id = id * 10 + static_cast<unsigned int>(ext[q] - '0');
      }
    }
    if (q >= ext.size() || ext[q] != ':') {
      continue;
    }
    auto &atoms = groups[10 * id + type];
    for (++q; q < ext.size() && digit(ext[q]);) {
      unsigned long long index = 0;
      for (; q < ext.size() && digit(ext[q]); ++q) {
        index = std::min(index * 10 + static_cast<unsigned>(ext[q] - '0'),
                         static_cast<unsigned long long>(UINT_MAX) + 1);
      }
      if (!atoms.insert(index).second) {
        return true;
      }
      if (q < ext.size() && ext[q] == ',') {
        ++q;
      }
    }
  }
  return false;
}

}  // namespace

std::string cx_extension_problem(std::string_view cx, bool strict) {
  if (cx_stereo_group_repeats(cx)) {
    return "a CXSMILES stereo group lists an atom twice, which RDKit 2024.03 "
           "keeps, leaving it to dangle once the atom is removed";
  }
  const std::size_t n = cx.size();
  for (std::size_t p = 0; p < n; ++p) {
    if (cx[p] == 'S' && p + 3 < n && cx.compare(p, 4, "SgD:") == 0) {
      if (data_sgroup_past_end(CxCursor(cx, p + 4))) {
        return "a CXSMILES data S-group (SgD) has fewer than five ':'-ended "
               "fields; RDKit 2024.03 reads past the end of the extension "
               "text";
      }
    } else if (!strict && cx[p] == 'S' && p + 2 < n &&
               cx.compare(p, 3, "Sg:") == 0) {
      if (polymer_sgroup_past_end(CxCursor(cx, p + 3))) {
        return "a CXSMILES polymer S-group (Sg) has no ':' after its type; "
               "with strictCXSMILES off, RDKit 2024.03 copies the text past "
               "the end of the extension";
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

std::string smiles_cx_problem(const std::string &text, bool strict) {
  auto sidx = text.find_first_of(" \t");
  if (sidx == std::string::npos || sidx == 0) {
    return {};
  }
  auto cx = boost::trim_copy(text.substr(sidx));
  if (cx.empty() || cx[0] != '|') {
    return {};
  }
  return cx_extension_problem(cx, strict);
}

void check_smiles_cx(const std::string &text, bool allow_cxsmiles,
                     bool strict) {
  if (!allow_cxsmiles) {
    return;
  }
  auto why = smiles_cx_problem(text, strict);
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
  // parseReaction throws when the extension does not parse.
  auto why = cx_extension_problem(boost::trim_copy(text.substr(sidx)), true);
  if (!why.empty()) {
    throw ValueErrorException(why);
  }
}

void check_text_cx_fields(std::string_view text, const char *what) {
  const auto bar = text.find('|');
  if (bar == std::string_view::npos) {
    return;
  }
  for (std::string_view field : {"SgD:", "atomProp:", "Sg:"}) {
    if (text.find(field, bar) != std::string_view::npos) {
      throw ValueErrorException(
          std::string(what) +
          " holds a CXSMILES field (SgD:, atomProp:, or Sg:) after a '|'; "
          "RDKit 2024.03 reads past the end of one that ends the SMARTS it "
          "parses there, so such fields are refused in text RDKit's readers "
          "parse");
    }
  }
}

std::vector<std::pair<std::string, std::string>> read_transformation_lines(
    const std::string &text) {
  // As readTransformations(std::istream &) and getReaction(const
  // std::string &), down to their reading past the last token: a line of one
  // field gives an empty SMARTS (the tokenizer clears its token at the end).
  std::vector<std::pair<std::string, std::string>> result;
  std::istringstream in(text);
  constexpr int kMaxLineLen = 512;
  char line[kMaxLineLen];
  while (!in.eof() && !in.fail()) {
    in.getline(line, kMaxLineLen, '\n');
    std::string tmp(line);
    if (tmp.empty() || tmp.compare(0, 2, "//") == 0) {
      continue;
    }
    std::vector<std::string> fields;
    boost::split(fields, tmp, boost::is_any_of("\t"), boost::token_compress_on);
    fields.erase(std::remove(fields.begin(), fields.end(), ""), fields.end());
    fields.resize(std::max<std::size_t>(fields.size(), 2));
    boost::erase_all(fields[0], " ");
    boost::erase_all(fields[1], " ");
    result.emplace_back(std::move(fields[0]), std::move(fields[1]));
  }
  return result;
}

void check_transformation_lines(const std::string &text) {
  for (const auto &[name, smarts] : read_transformation_lines(text)) {
    check_reaction_cx(smarts);
  }
}

void check_query_props(const RDKit::ROMol &mol, const std::string &prop) {
  for (const auto *atom : mol.atoms()) {
    std::string value;
    try {
      if (!atom->getPropIfPresent(prop, value)) {
        continue;
      }
    } catch (const std::exception &) {
      continue;  // not a string, which RDKit refuses
    }
    check_smiles_cx(value, true, true);
  }
}

void check_mrv(const std::string &text) {
  std::istringstream in(text);
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(in, tree);
  check_bond_refs(tree);
  check_atom_ids(tree);
}

}  // namespace rdkit_cxx
