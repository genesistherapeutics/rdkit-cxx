// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/control_c_lock.h"

#include <GraphMol/Chirality.h>

namespace rdkit_cxx {

std::unique_lock<std::mutex> lock_control_c_for_cip() {
  if (RDKit::Chirality::getUseLegacyStereoPerception()) {
    return {};
  }
  return std::unique_lock<std::mutex>(control_c_mutex());
}

}  // namespace rdkit_cxx
