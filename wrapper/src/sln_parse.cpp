// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/sln_parse.h"

#include <algorithm>
#include <mutex>
#include <regex>
#include <string>

#include <GraphMol/SLNParse/SLNParse.h>
#include <boost/algorithm/string.hpp>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/relevant_cycles.h"

namespace rdkit_cxx {

namespace {

// SLNToMol and SLNQueryToMol set the parser's global yysln_debug, which the
// parser reads, without a lock (U4).
std::mutex sln_mutex;

std::unique_ptr<RWMol> owned(RDKit::RWMol *mol) {
  return std::unique_ptr<RWMol>(mol);
}

// The text SLN's parser reads (SLNParse.cpp's toMol): trimmed, without its
// macro definitions ("{name:value}"), and with every match of each
// definition's name, as a regex, replaced by its value. A copy of
// replaceSLNMacroAtoms, which RDKit's headers don't declare.
std::string expanded_sln(std::string inp) {
  boost::trim_if(inp, boost::is_any_of(" \t\r\n"));
  const std::regex defn("\\{(.+?):(.+?)\\}");
  std::string res = std::regex_replace(
      inp, defn, "",
      std::regex_constants::match_default | std::regex_constants::format_sed);
  if (res != inp) {
    std::string::const_iterator start = inp.begin();
    const std::string::const_iterator end = inp.end();
    std::match_results<std::string::const_iterator> what;
    auto flags = std::regex_constants::match_default;
    while (std::regex_search(start, end, what, defn, flags)) {
      std::string name(what[1].first, what[1].second);
      std::string value(what[2].first, what[2].second);
      res = std::regex_replace(res, std::regex(name), value.c_str(),
                               std::regex_constants::match_default |
                                   std::regex_constants::format_sed);
      start = what[0].second;
      flags |= std::regex_constants::match_prev_avail;
    }
  }
  return res;
}

// SLN closes rings only with "@" (sln.yy's AT_TOKEN), each one bond, so a
// molecule, or a query and the queries nested in it, has no more
// independent cycles than the text the parser reads has "@"s. A macro can
// be one, so they're counted once the macros are expanded.
bool few_ring_closures(const std::string &sln) {
  const auto text = expanded_sln(sln);
  return static_cast<std::size_t>(std::count(text.begin(), text.end(), '@')) <=
         MAX_INDEPENDENT_CYCLES;
}

}  // namespace

std::unique_ptr<RWMol> mol_from_sln(rust::Str sln, bool sanitize,
                                    bool debug_parser) {
  std::lock_guard<std::mutex> lock(sln_mutex);
  std::string text(sln);
  // Sanitizing, SLNToMol removes the hydrogens, which sanitizes; without,
  // it perceives nothing.
  if (sanitize && !few_ring_closures(text)) {
    std::unique_ptr<RDKit::RWMol> raw(
        RDKit::SLNToMol(text, false, debug_parser));
    if (raw) {
      check_remove_hs_relevant_cycles(*raw, RDKit::MolOps::RemoveHsParameters(),
                                      true);
    }
  }
  return owned(RDKit::SLNToMol(text, sanitize, debug_parser));
}

std::unique_ptr<RWMol> mol_from_query_sln(rust::Str sln, bool merge_hs,
                                          bool debug_parser) {
  std::lock_guard<std::mutex> lock(sln_mutex);
  std::string text(sln);
  // finalizeQueryMol symmetrizes every query it parses, and there is no
  // parsing one without (README, Known limitations).
  if (!few_ring_closures(text) && !RDKit::MolOps::getUseLegacyRingFinding()) {
    throw ValueErrorException(
        "a query SLN with more than " + std::to_string(MAX_INDEPENDENT_CYCLES) +
        " ring closures cannot be checked for its relevant cycles before "
        "RDKit perceives them");
  }
  return owned(RDKit::SLNQueryToMol(text, merge_hs, debug_parser));
}

}  // namespace rdkit_cxx
