// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#include "rdkit-cxx/rd_base.h"

#include <RDGeneral/RDLog.h>
#include <RDGeneral/versions.h>

namespace rdkit_cxx {

std::unique_ptr<std::string> rdkit_version() {
  return std::make_unique<std::string>(RDKit::rdkitVersion);
}

std::unique_ptr<std::string> boost_version() {
  return std::make_unique<std::string>(RDKit::boostVersion);
}

// As DisableLog and EnableLog in Code/RDBoost/nbWrap/RDBase.cpp.
void disable_log(rust::Str spec) {
  boost::logging::disable_logs(std::string(spec));
}

void enable_log(rust::Str spec) {
  boost::logging::enable_logs(std::string(spec));
}

}  // namespace rdkit_cxx
