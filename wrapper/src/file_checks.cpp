// Copyright (c) 2026 Genesis Molecular AI
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

// RDKit 2025.03's Marvin reader puts each atom with an mrvStereoGroup into
// its group by the index of the first atom with its id (getAtomIndex), so two
// atoms with one id list an atom twice, which dangles once removed
// (check_stereo_group_members in mol_checks.h): say, when the reader removes
// the hydrogens. Ids in an atom array must therefore differ.
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

void check_mrv(const std::string &text) {
  std::istringstream in(text);
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(in, tree);
  check_bond_refs(tree);
  check_atom_ids(tree);
}

namespace {

// Whether a run of digits is one boost::lexical_cast throws for: a number
// past `max`.
bool past(std::string_view digits, unsigned long long max) {
  unsigned long long value = 0;
  for (char c : digits) {
    value = value * 10 + static_cast<unsigned long long>(c - '0');
    if (value > max) {
      return true;
    }
  }
  return false;
}

bool digit(char c) { return c >= '0' && c <= '9'; }

// Replays parse_atom_props on the extension `ext` from `first`, just after
// an "atomProp:": whether it reads past the end of `ext`. RDKit's iterators
// run to and read `last`, a std::string's terminator, throughout; this one
// goes past it, where a list ends without "|" or ",".
bool atom_props_overrun(std::string_view ext, std::size_t first) {
  const std::size_t last = ext.size();
  auto at = [&](std::size_t i) { return i < last ? ext[i] : '\0'; };
  // read_text_to: false where it throws; `*text` whether it read any.
  auto read_to = [&](std::string_view delims, bool *text) {
    std::size_t start = first;
    *text = false;
    while (first <= last && delims.find(at(first)) == std::string_view::npos) {
      if (at(first) == '&' && last - first > 2 && at(first + 1) == '#') {
        if (start != first) {
          *text = true;
        }
        std::size_t next = first + 2;
        while (next != last && digit(at(next))) {
          ++next;
        }
        if (next == last || at(next) != ';') {
          return false;
        }
        if (next > first + 2) {
          if (past(ext.substr(first + 2, next - first - 2), INT_MAX)) {
            return false;
          }
          *text = true;
        }
        first = next + 1;
        start = first;
      } else {
        ++first;
      }
    }
    if (start != first) {
      *text = true;
    }
    return true;
  };
  if (first >= last) {
    return false;
  }
  while (first <= last && at(first) != '|' && at(first) != ',') {
    std::size_t end = first;
    while (end <= last && digit(at(end))) {
      ++end;
    }
    if (end > first) {
      if (past(ext.substr(first, end - first), UINT_MAX)) {
        return false;
      }
      first = end;
      if (first >= last || at(first) != '.') {
        return false;
      }
      ++first;
      bool name;
      if (!read_to(".", &name)) {
        return false;
      }
      if (name) {
        if (first >= last || at(first) != '.') {
          return false;
        }
        ++first;
        bool value;
        if (!read_to(":|,", &value)) {
          return false;
        }
      }
    }
    if (first <= last && at(first) != '|' && at(first) != ',') {
      ++first;
    }
  }
  // RDKit then reads *first.
  return first > last;
}

// RDKit 2025.03's parse_enhanced_stereo adds each atom index listed after
// "a:", "o<n>:", or "&<n>:" to the group the type and number name, merging
// sections by 10 * number + type, so a group can list an atom twice, which
// leaves a dangling atom once the atom is removed (check_stereo_group_members
// in mol_checks.h): say, when the parser removes the hydrogens. Every place a
// section may start is read as one, which may refuse a few extensions where
// such text is inside another section.
void check_cx_stereo_groups(std::string_view ext) {
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
        throw ValueErrorException(
            "a CXSMILES stereo group lists an atom twice, which RDKit 2025.03 "
            "keeps, leaving it to dangle once the atom is removed");
      }
      if (q < ext.size() && ext[q] == ',') {
        ++q;
      }
    }
  }
}

// Where read_text_to (CXSmilesOps.cpp) stops reading `cx` from `i`: at the
// first of `delims`, or one past the end, having read the terminator there.
// Its "&#N;" escapes hold only digits, none of the delimiters given here.
std::size_t text_to(std::string_view cx, std::size_t i, const char *delims) {
  auto k = cx.find_first_of(delims, i);
  return k == std::string_view::npos ? cx.size() + 1 : k;
}

// Replays parse_data_sgroup on `cx` from `p`, at "SgD:": whether a field or
// its coordinates run past the end, after which RDKit steps once more.
bool data_sgroup_overrun(std::string_view cx, std::size_t p) {
  const std::size_t n = cx.size();
  if (p + 3 >= n || cx[p + 3] != ':') {
    return false;
  }
  std::size_t i = p + 4;
  while (true) {
    while (i < n && digit(cx[i])) {
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

// RDKit's trim (boost::trim_copy, in the classic locale).
std::string_view trimmed(std::string_view text) {
  auto space = [](char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' ||
           c == '\r';
  };
  while (!text.empty() && space(text.front())) {
    text.remove_prefix(1);
  }
  while (!text.empty() && space(text.back())) {
    text.remove_suffix(1);
  }
  return text;
}

}  // namespace

// RDKit 2025.03's CXSMILES parser (CXSmilesOps.cpp, unchanged here since
// 2024.09) reads past the end of the extension in two places:
// - an "atomProp:" list (parse_atom_props) that runs to the end without a
//   "|" or "," reads the byte after the text's terminator, whatever
//   strictCXSMILES says (rdkit/rdkit#9361);
// - a data substance group's ("SgD:") field without a ":" after it, its
//   coordinates without a ")", and a polymer group's ("Sg:") type without a
//   ":" step past the end, and the parse fails there; without
//   strictCXSMILES, MolFromSmiles and MolFromSmarts then copy the text up to
//   where it stopped into _CXSMILES_Data, past the terminator.
// Each such block is read from wherever it appears in the extension, as the
// parser would read it there, which covers wherever the parser does. Stereo
// groups that list an atom twice are refused too (check_cx_stereo_groups).
void check_cx_extensions(std::string_view cx, bool strict) {
  check_cx_stereo_groups(cx);
  constexpr std::string_view key = "atomProp:";
  for (auto at = cx.find(key); at != std::string_view::npos;
       at = cx.find(key, at + 1)) {
    if (atom_props_overrun(cx, at + key.size())) {
      throw ValueErrorException(
          "a CXSMILES atomProp list runs to the end of the extension without "
          "a \"|\" or \",\", where RDKit 2025.03 reads past its end");
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
          "a CXSMILES substance group's fields must each end with \":\" (and "
          "its coordinates with \")\") when strictCXSMILES is off: RDKit "
          "2025.03 copies the text past its end otherwise");
    }
  }
}

// As preprocessSmiles in SmilesParse.cpp.
void check_smiles_cx(std::string_view text, bool allow_cx, bool strict) {
  if (!allow_cx) {
    return;
  }
  auto split = text.find_first_of(" \t");
  if (split == std::string_view::npos || split == 0) {
    return;
  }
  auto cx = trimmed(text.substr(split));
  if (!cx.empty() && cx.front() == '|') {
    check_cx_extensions(cx, strict);
  }
}

// As DaylightParserUtils::reactionFromSmarts in DaylightParser.cpp, which
// parses the extension strictly and throws when it fails.
void check_reaction_cx(std::string_view text, bool allow_cx) {
  if (!allow_cx) {
    return;
  }
  auto split = text.find('|');
  if (split != std::string_view::npos && split != 0) {
    check_cx_extensions(trimmed(text.substr(split)), true);
  }
}

// As readTransformations in TransformCatalogUtils.cpp.
std::vector<std::pair<std::string, std::string>> transformation_data(
    const std::string &text) {
  std::vector<std::pair<std::string, std::string>> data;
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
    if (fields.size() < 2) {
      throw ValueErrorException(
          "a transformation line needs a name and a reaction SMARTS, "
          "separated by a tab");
    }
    boost::erase_all(fields[0], " ");
    boost::erase_all(fields[1], " ");
    data.emplace_back(std::move(fields[0]), std::move(fields[1]));
  }
  return data;
}

// As readTransformations in TransformCatalogUtils.cpp.
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

}  // namespace rdkit_cxx
