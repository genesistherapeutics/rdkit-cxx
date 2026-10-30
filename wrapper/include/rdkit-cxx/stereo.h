/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Stereo groups and stereo info, which cross as index-based data
// (src/ffi/stereo.rs) because RDKit's hold pointers into their molecule.

#pragma once

#include <cstdint>

#include <GraphMol/Chirality.h>
#include <GraphMol/StereoGroup.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using StereoGroupType = RDKit::StereoGroupType;
using StereoType = RDKit::Chirality::StereoType;
using StereoSpecified = RDKit::Chirality::StereoSpecified;
using StereoDescriptor = RDKit::Chirality::StereoDescriptor;

struct StereoGroupData;
struct StereoInfoData;

rust::Vec<StereoGroupData> mol_get_stereo_groups(const RWMol &mol);
// As ReadWriteMol::SetStereoGroups in Code/GraphMol/Wrap/Mol.cpp, with
// indices checked against the molecule (CONVENTIONS.md E3).
void rwmol_set_stereo_groups(RWMol &mol,
                             rust::Slice<const StereoGroupData> groups);
// As createStereoGroup in Code/GraphMol/Wrap/StereoGroup.cpp, whose index
// checks have their own messages; here the indices are checked as for
// rwmol_set_stereo_groups.
StereoGroupData create_stereo_group(StereoGroupType group_type,
                                    const RWMol &mol,
                                    rust::Slice<const std::uint32_t> atom_ids,
                                    rust::Slice<const std::uint32_t> bond_ids,
                                    unsigned read_id);
void forward_stereo_group_ids(RWMol &mol);
rust::Vec<StereoInfoData> find_potential_stereo(RWMol &mol, bool clean_it,
                                                bool flag_possible);

}  // namespace rdkit_cxx
