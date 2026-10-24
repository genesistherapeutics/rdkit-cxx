// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause
//
// Value types shared across areas. The bridge (src/ffi/geometry.rs) defines
// them; wrappers that use them include "rdkit-cxx/ffi/geometry.rs.h".

#pragma once

namespace rdkit_cxx {

// rdGeometry.Point3D as plain data. RDGeom::Point3D has a vtable, so it
// cannot cross the bridge itself.
struct Point3D;

// A list of atom or bond indices. cxx has no Vec<Vec<T>>, so nested lists
// cross as Vec<IndexList>.
struct IndexList;

}  // namespace rdkit_cxx
