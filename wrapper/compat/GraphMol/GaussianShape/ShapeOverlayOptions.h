/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit 2026.09 installs GaussianShape's headers (GaussianShape.h,
// ShapeInput.h, ShapeOverlayOptions.h) at the top of include/rdkit, since
// its CMakeLists.txt gives rdkit_headers no DEST, but they include each other
// as GraphMol/GaussianShape/<name>. This forwards that name to the installed
// header. This directory comes after RDKit's in the include path, so an
// install with the headers in GraphMol/GaussianShape finds its own.

#pragma once

#include <ShapeOverlayOptions.h>
