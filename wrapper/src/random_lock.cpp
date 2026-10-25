// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/random_lock.h"

namespace rdkit_cxx {

std::unique_lock<std::mutex> lock_random_generator() {
  static std::mutex mutex;
  return std::unique_lock<std::mutex>(mutex);
}

}  // namespace rdkit_cxx
