// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/global_state.h"

namespace rdkit_cxx {

std::mutex &control_c_mutex() {
  static std::mutex mutex;
  return mutex;
}

std::mutex &random_mutex() {
  static std::mutex mutex;
  return mutex;
}

}  // namespace rdkit_cxx
