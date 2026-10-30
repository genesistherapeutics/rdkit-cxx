/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's FragmentRemover and LargestFragmentChooser
// (Fragment.h).

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolStandardize/Fragment.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using FragmentRemover = RDKit::MolStandardize::FragmentRemover;

std::unique_ptr<FragmentRemover> fragment_remover_default();
std::unique_ptr<FragmentRemover> fragment_remover_new(
    rust::Slice<const std::uint8_t> fragment_filename, bool leave_last,
    bool skip_if_all_match);
std::unique_ptr<FragmentRemover> fragment_remover_from_data(
    rust::Str fragment_data, bool leave_last, bool skip_if_all_match);
std::unique_ptr<RWMol> fragment_remover_remove(const FragmentRemover &remover,
                                               const RWMol &mol);
void fragment_remover_remove_in_place(const FragmentRemover &remover,
                                      RWMol &mol);

// The chooser's three options, as LargestFragmentChooser(params) reads them.
std::unique_ptr<RWMol> largest_fragment_chooser_choose(
    const RWMol &mol, bool prefer_organic, bool use_atom_count,
    bool count_heavy_atoms_only);
void largest_fragment_chooser_choose_in_place(RWMol &mol, bool prefer_organic,
                                              bool use_atom_count,
                                              bool count_heavy_atoms_only);

}  // namespace rdkit_cxx
