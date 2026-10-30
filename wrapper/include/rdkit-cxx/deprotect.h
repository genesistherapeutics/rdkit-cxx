/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdDeprotect: removing protecting groups.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/Deprotect/Deprotect.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using DeprotectData = RDKit::Deprotect::DeprotectData;

// Deprotections Rust lends to one call, copied (the copies share each
// reaction, which RDKit only reads once it is initialized).
class DeprotectList {
 public:
  std::vector<DeprotectData> items;
};

std::unique_ptr<DeprotectData> deprotect_data_new(rust::Str deprotection_class,
                                                  rust::Str reaction_smarts,
                                                  rust::Str abbreviation,
                                                  rust::Str full_name);
std::unique_ptr<DeprotectData> deprotect_data_clone(const DeprotectData &data);
std::unique_ptr<std::string> deprotect_data_deprotection_class(
    const DeprotectData &data);
std::unique_ptr<std::string> deprotect_data_full_name(
    const DeprotectData &data);
std::unique_ptr<std::string> deprotect_data_abbreviation(
    const DeprotectData &data);
std::unique_ptr<std::string> deprotect_data_reaction_smarts(
    const DeprotectData &data);
std::unique_ptr<std::string> deprotect_data_example(const DeprotectData &data);
bool deprotect_data_is_valid(const DeprotectData &data);

std::size_t num_default_deprotections();
std::unique_ptr<DeprotectData> default_deprotection(std::size_t i);

std::unique_ptr<DeprotectList> deprotect_list_new();
void deprotect_list_push(DeprotectList &list, const DeprotectData &data);

// The _default forms use RDKit's built-in deprotections.
std::unique_ptr<RWMol> deprotect_with(const RWMol &mol,
                                      const DeprotectList &list);
std::unique_ptr<RWMol> deprotect_default(const RWMol &mol);
bool deprotect_in_place_with(RWMol &mol, const DeprotectList &list);
bool deprotect_in_place_default(RWMol &mol);

}  // namespace rdkit_cxx
