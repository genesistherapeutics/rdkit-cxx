// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/file_checks.h"

#include <climits>
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

void check_mrv(const std::string &text) {
  std::istringstream in(text);
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(in, tree);
  check_bond_refs(tree);
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

void check_cx_extension(std::string_view ext) {
  constexpr std::string_view key = "atomProp:";
  for (auto at = ext.find(key); at != std::string_view::npos;
       at = ext.find(key, at + 1)) {
    if (atom_props_overrun(ext, at + key.size())) {
      throw ValueErrorException(
          "a CXSMILES atomProp list runs to the end of the extension without "
          "a \"|\" or \",\", where RDKit 2025.03 reads past its end");
    }
  }
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

// As preprocessSmiles in SmilesParse.cpp.
void check_cx_smiles(std::string_view smiles) {
  auto split = smiles.find_first_of(" \t");
  if (split != std::string_view::npos && split != 0) {
    check_cx_extension(trimmed(smiles.substr(split)));
  }
}

// As DaylightParserUtils::reactionFromSmarts in DaylightParser.cpp.
void check_cx_reaction(std::string_view text) {
  auto split = text.find('|');
  if (split != std::string_view::npos && split != 0) {
    check_cx_extension(trimmed(text.substr(split)));
  }
}

}  // namespace rdkit_cxx
