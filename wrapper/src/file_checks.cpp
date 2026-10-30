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

// RDKit 2024.09's Marvin reader puts each atom with an mrvStereoGroup into
// its group by the index of the first atom with its id (getAtomIndex), so two
// atoms with one id list an atom twice, which dangles once removed
// (check_stereo_group_members in mol_checks.cpp): say, when the reader
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

// Where read_text_to (CXSmilesOps.cpp) stops reading `cx` from `i`: at the
// first of `delims`, or one past the end, having read the null there. Its
// "&#N;" escapes hold none of the delimiters it is given here.
std::size_t text_to(std::string_view cx, std::size_t i, const char *delims) {
  auto k = cx.find_first_of(delims, i);
  return k == std::string_view::npos ? cx.size() + 1 : k;
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// parse_atom_props from `i`, just after "atomProp:": whether it reads past
// the null at the end. At the end it reads the null, which ends nothing.
bool atom_props_overrun(std::string_view cx, std::size_t i) {
  const std::size_t n = cx.size();
  if (i >= n) {
    return false;
  }
  auto stop = [&](std::size_t j) {
    return j < n && (cx[j] == '|' || cx[j] == ',');
  };
  while (i <= n && !stop(i)) {
    std::size_t j = i;
    while (j < n && is_digit(cx[j])) {
      ++j;
    }
    if (j > i) {
      i = j;
      if (i >= n || cx[i] != '.') {
        return false;
      }
      auto k = text_to(cx, ++i, ".");
      bool has_name = k != i;
      i = k;
      if (has_name) {
        if (i >= n || cx[i] != '.') {
          return false;
        }
        i = text_to(cx, i + 1, ":|,");
      }
    }
    if (i <= n && !stop(i)) {
      ++i;
    }
  }
  // RDKit returns false at the null and reads the byte after it past it.
  return i > n;
}

// parse_data_sgroup from `p`, at "SgD:": whether a field or its coordinates
// run past the end, after which RDKit steps once more.
bool data_sgroup_overrun(std::string_view cx, std::size_t p) {
  const std::size_t n = cx.size();
  if (p + 3 >= n || cx[p + 3] != ':') {
    return false;
  }
  std::size_t i = p + 4;
  while (true) {
    while (i < n && is_digit(cx[i])) {
      ++i;
    }
    if (i >= n || cx[i] != ',') {
      break;
    }
    ++i;
  }
  ++i;
  for (int field = 0; field < 5; ++field) {
    if (i >= n) {
      return false;  // parse_data_sgroup_attr's precondition throws
    }
    if (cx[i] != '|') {
      i = text_to(cx, i, ":");
      if (i > n) {
        return true;
      }
      ++i;
    }
  }
  return i < n && cx[i] == '(' && text_to(cx, i, ")") > n;
}

// RDKit 2024.09's parse_enhanced_stereo adds each atom index listed after
// "a:", "o<n>:", or "&<n>:" to the group the type and number name, merging
// sections by 10 * number + type, so a group can list an atom twice, which
// leaves a dangling atom once the atom is removed (check_stereo_group_members
// in mol_checks.cpp): say, when the parser removes the hydrogens. Every place
// a section may start is read as one, which may refuse a few extensions where
// such text is inside another section.
void check_cx_stereo_groups(std::string_view cx) {
  std::map<unsigned int, std::set<unsigned long long>> groups;
  for (std::size_t p = 0; p < cx.size(); ++p) {
    unsigned int type;
    if (cx[p] == 'a') {
      type = 0;
    } else if (cx[p] == 'o') {
      type = 1;
    } else if (cx[p] == '&' && p + 1 < cx.size() && cx[p + 1] != '#') {
      type = 2;
    } else {
      continue;
    }
    std::size_t q = p + 1;
    unsigned int id = 0;
    if (type != 0) {
      for (; q < cx.size() && is_digit(cx[q]); ++q) {
        id = id * 10 + static_cast<unsigned int>(cx[q] - '0');
      }
    }
    if (q >= cx.size() || cx[q] != ':') {
      continue;
    }
    auto &atoms = groups[10 * id + type];
    for (++q; q < cx.size() && is_digit(cx[q]);) {
      unsigned long long index = 0;
      for (; q < cx.size() && is_digit(cx[q]); ++q) {
        index = std::min(index * 10 + static_cast<unsigned>(cx[q] - '0'),
                         static_cast<unsigned long long>(UINT_MAX) + 1);
      }
      if (!atoms.insert(index).second) {
        throw ValueErrorException(
            "a CXSMILES stereo group lists an atom twice, which RDKit 2024.09 "
            "keeps, leaving it to dangle once the atom is removed");
      }
      if (q < cx.size() && cx[q] == ',') {
        ++q;
      }
    }
  }
}

}  // namespace

void check_cx_extensions(std::string_view cx, bool strict) {
  check_cx_stereo_groups(cx);
  for (auto p = cx.find("atomProp:"); p != std::string_view::npos;
       p = cx.find("atomProp:", p + 1)) {
    // parse_it takes it for an "atomProp:" block only with more after it.
    if (cx.size() - p > 9 && atom_props_overrun(cx, p + 9)) {
      throw ValueErrorException(
          "a CXSMILES atomProp block must end with '|' or ',': RDKit 2024.09 "
          "reads past the end of the text otherwise");
    }
  }
  if (strict) {
    return;
  }
  for (auto p = cx.find("Sg"); p != std::string_view::npos;
       p = cx.find("Sg", p + 1)) {
    bool overrun = false;
    if (p + 2 < cx.size() && cx[p + 2] == 'D') {
      overrun = data_sgroup_overrun(cx, p);
    } else if (p + 2 < cx.size() && cx[p + 2] == ':') {
      overrun = text_to(cx, p + 3, ":") > cx.size();
    }
    if (overrun) {
      throw ValueErrorException(
          "a CXSMILES substance group's fields must each end with ':' (and "
          "its coordinates with ')') when strictCXSMILES is off: RDKit 2024.09 "
          "copies the text past its end otherwise");
    }
  }
}

void check_smiles_cx(const std::string &text, bool allow_cx, bool strict) {
  if (!allow_cx) {
    return;
  }
  auto sidx = text.find_first_of(" \t");
  if (sidx == std::string::npos || sidx == 0) {
    return;
  }
  auto cx = boost::trim_copy(text.substr(sidx));
  if (!cx.empty() && cx[0] == '|') {
    check_cx_extensions(cx, strict);
  }
}

void check_reaction_cx(const std::string &text, bool allow_cx) {
  if (!allow_cx) {
    return;
  }
  auto sidx = text.find_first_of('|');
  if (sidx == std::string::npos || sidx == 0) {
    return;
  }
  check_cx_extensions(boost::trim_copy(text.substr(sidx)), true);
}

void check_transformation_lines(const std::string &text) {
  std::istringstream in(text);
  char line[512];
  while (!in.eof() && !in.fail()) {
    in.getline(line, sizeof(line), '\n');
    std::string tmp(line);
    if (tmp.empty() || tmp.compare(0, 2, "//") == 0) {
      continue;
    }
    std::vector<std::string> fields;
    boost::split(fields, tmp, boost::is_any_of("\t"), boost::token_compress_on);
    fields.erase(std::remove(fields.begin(), fields.end(), ""), fields.end());
    if (fields.size() >= 2) {
      boost::erase_all(fields[1], " ");
      check_reaction_cx(fields[1], true);
    }
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
