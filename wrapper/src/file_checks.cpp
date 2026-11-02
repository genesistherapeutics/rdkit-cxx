// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/file_checks.h"

#include <algorithm>
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

// Where read_text_to (CXSmilesOps.cpp) stops reading `cx` from `i`: at the
// first of `delims`, or one past the end, having read the null there. Its
// "&#N;" escapes hold none of the delimiters it is given here.
std::size_t text_to(std::string_view cx, std::size_t i, const char *delims) {
  auto k = cx.find_first_of(delims, i);
  return k == std::string_view::npos ? cx.size() + 1 : k;
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// parse_data_sgroup from `p`, at "SgD:": whether a field or its coordinates
// run past the end, after which RDKit steps once more.
bool data_sgroup_overrun(std::string_view cx, std::size_t p) {
  const std::size_t n = cx.size();
  if (p + 3 >= n || cx[p + 3] != ':') {
    return false;
  }
  std::size_t i = p + 4;
  // read_int_list
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
  // parse_data_sgroup_attr, for each of the five fields
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

}  // namespace

void check_cx_extensions(std::string_view cx, bool strict) {
  if (strict) {
    return;
  }
  for (auto p = cx.find("Sg"); p != std::string_view::npos;
       p = cx.find("Sg", p + 1)) {
    bool overrun = false;
    if (p + 2 < cx.size() && cx[p + 2] == 'D') {
      overrun = data_sgroup_overrun(cx, p);
    } else if (p + 2 < cx.size() && cx[p + 2] == ':') {
      // parse_polymer_sgroup: its type and the ':' after it.
      overrun = text_to(cx, p + 3, ":") > cx.size();
    }
    if (overrun) {
      throw ValueErrorException(
          "a CXSMILES substance group's fields must each end with ':' (and "
          "its coordinates with ')') when strictCXSMILES is off: RDKit copies "
          "the text past its end otherwise");
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

void check_mrv(const std::string &text) {
  std::istringstream in(text);
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(in, tree);
  check_bond_refs(tree);
}

// As readTransformations in TransformCatalogUtils.cpp: lines of at most 511
// characters, up to the first longer one, skipping empty lines and "//"
// comments, each a name and a reaction SMARTS separated by tabs, with their
// spaces erased.
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

}  // namespace rdkit_cxx
