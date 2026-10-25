/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit's MolWriters.h includes maeparser's Writer.hpp for MaeWriter, which
// holds a std::unique_ptr<schrodinger::mae::Writer>, but RDKit installs
// without maeparser's headers. The wrappers never use the class itself, so a
// declaration serves, as MolSupplier.h declares maeparser's Reader and Block.
// This directory comes after RDKit's in the include path.

#pragma once

namespace schrodinger {
namespace mae {
class Writer;
}  // namespace mae
}  // namespace schrodinger
