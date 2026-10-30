/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdCoordGen: 2D coordinates from Schrödinger's CoordGen. CoordGen keeps its
// templates and template directory in globals, so every function holds
// lock_depiction() (CONVENTIONS.md U4).

#pragma once

#include <cstdint>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

struct CoordGenParamsData;

CoordGenParamsData coord_gen_params_default();
// template_mol is used only with has_template_mol.
void add_coords(RWMol &mol, const CoordGenParamsData &data,
                const RWMol &template_mol, bool has_template_mol);
// With the parameters whose template directory
// set_default_template_file_dir sets.
void add_coords_default(RWMol &mol);
void set_default_template_file_dir(rust::Slice<const std::uint8_t> dir);

}  // namespace rdkit_cxx
