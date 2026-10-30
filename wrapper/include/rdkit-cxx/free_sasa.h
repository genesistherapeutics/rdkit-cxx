/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for External/FreeSASA (rdkit.Chem.rdFreeSASA): solvent
// accessible surface areas from the freesasa C library, and its atom
// classifiers. SASAOpts' enums have no fixed underlying type, so cxx cannot
// check them; they cross as int32 and are checked here.

#pragma once

#include <cstdint>
#include <memory>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A shared struct, defined by the bridge (src/ffi/free_sasa.rs).
struct SASAOptsData;

SASAOptsData sasa_opts_default();
// The radii, or none when some atom could not be classified.
rust::Vec<double> classify_atoms(RWMol &mol, const SASAOptsData &opts);
// `query` is null or an atom.
double calc_sasa(RWMol &mol, rust::Slice<const double> radii,
                 std::int32_t conf_idx, const Atom *query,
                 const SASAOptsData &opts);
std::unique_ptr<Atom> make_free_sasa_a_polar_atom_query();
std::unique_ptr<Atom> make_free_sasa_polar_atom_query();

}  // namespace rdkit_cxx
