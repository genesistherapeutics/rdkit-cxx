// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/sln_parse.h"

#include <mutex>
#include <string>

#include <GraphMol/SLNParse/SLNParse.h>

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

namespace {

// SLNToMol and SLNQueryToMol set the parser's global yysln_debug, which the
// parser reads, without a lock (U4).
std::mutex sln_mutex;

std::unique_ptr<RWMol> owned(RDKit::RWMol *mol) {
  return std::unique_ptr<RWMol>(mol);
}

}  // namespace

std::unique_ptr<RWMol> mol_from_sln(rust::Str sln, bool sanitize,
                                    bool debug_parser) {
  std::lock_guard<std::mutex> lock(sln_mutex);
  return owned(RDKit::SLNToMol(std::string(sln), sanitize, debug_parser));
}

std::unique_ptr<RWMol> mol_from_query_sln(rust::Str sln, bool merge_hs,
                                          bool debug_parser) {
  std::lock_guard<std::mutex> lock(sln_mutex);
  return owned(RDKit::SLNQueryToMol(std::string(sln), merge_hs, debug_parser));
}

}  // namespace rdkit_cxx
