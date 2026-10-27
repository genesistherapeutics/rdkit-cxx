// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/file_checks.h"

#include <cstdint>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string_view>
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

// What parse_atom_props reads from `first`, an index into `s`, which RDKit
// reads to and including its end (as the terminating character): true if
// the loop runs past the end. A read RDKit throws for stops it, safely.
bool atom_props_overrun(const std::string &s, std::size_t first) {
  const std::size_t last = s.size();
  auto at = [&](std::size_t i) { return i < last ? s[i] : '\0'; };
  auto digit = [&](std::size_t i) { return at(i) >= '0' && at(i) <= '9'; };
  // read_int: false for no digits; a number lexical_cast<unsigned int>
  // refuses sets `throws`.
  bool throws = false;
  auto read_int = [&](std::size_t &f) {
    const std::size_t start = f;
    while (f <= last && digit(f)) {
      ++f;
    }
    if (f == start) {
      return false;
    }
    std::uint64_t v = 0;
    for (std::size_t i = start; i < f && !throws; ++i) {
      v = v * 10 + static_cast<std::uint64_t>(s[i] - '0');
      throws = v > UINT32_MAX;
    }
    return true;
  };
  // read_text_to; whether it read anything.
  auto read_text_to = [&](std::size_t &f, std::string_view delims) {
    bool read = false;
    while (f <= last && (at(f) == '\0' || delims.find(at(f)) == delims.npos)) {
      if (at(f) == '&' && last - f > 2 && at(f + 1) == '#') {
        std::size_t next = f + 2;
        while (next != last && digit(next)) {
          ++next;
        }
        if (next == last || at(next) != ';') {
          throws = true;
          return read;
        }
        read = read || next > f + 2;
        f = next + 1;
      } else {
        read = true;
        ++f;
      }
    }
    return read;
  };
  if (first >= last) {
    return false;
  }
  while (first <= last && at(first) != '|' && at(first) != ',') {
    if (read_int(first)) {
      if (throws || first >= last || at(first) != '.') {
        return false;
      }
      ++first;
      const bool name = read_text_to(first, ".");
      if (throws) {
        return false;
      }
      if (name) {
        if (first >= last || at(first) != '.') {
          return false;
        }
        ++first;
        read_text_to(first, ":|,");
        if (throws) {
          return false;
        }
      }
    }
    if (first <= last && at(first) != '|' && at(first) != ',') {
      ++first;
    }
  }
  return first > last;
}

}  // namespace

void check_cx_atom_props(const std::string &text) {
  constexpr std::string_view key = "atomProp:";
  if (text.find(key) == std::string::npos) {
    return;
  }
  const auto trimmed = boost::trim_right_copy(text);
  for (auto at = trimmed.find(key); at != std::string::npos;
       at = trimmed.find(key, at + 1)) {
    if (atom_props_overrun(trimmed, at + key.size())) {
      throw ValueErrorException(
          "a CXSMILES atomProp block runs to the end of the text, which RDKit "
          "2025.09 reads past");
    }
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

}  // namespace rdkit_cxx
