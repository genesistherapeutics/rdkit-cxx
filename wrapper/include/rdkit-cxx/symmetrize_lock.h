/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// MolAlign::details::symmetrizeTerminalAtoms (GraphMol/MolAlign/
// AlignMolecules.cpp), which the alignment functions and the embedder's
// pruning call, assigns a static SmartsParserParams on every call, so two
// calls at once race. Every wrapper whose RDKit call can reach it holds this
// mutex for that call (CONVENTIONS.md U4). A wrapper that also needs
// control_c_mutex() takes that first, and random_mutex() after this one.

#pragma once

#include <mutex>

namespace rdkit_cxx {

// Defined inline so every wrapper source shares one mutex.
inline std::mutex &symmetrize_mutex() {
  static std::mutex mutex;
  return mutex;
}

}  // namespace rdkit_cxx
