// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>

#include "rust/cxx.h"

namespace rdkit_cxx {

std::unique_ptr<std::string> rdkit_version();
std::unique_ptr<std::string> boost_version();
void disable_log(rust::Str spec);
void enable_log(rust::Str spec);

}  // namespace rdkit_cxx
