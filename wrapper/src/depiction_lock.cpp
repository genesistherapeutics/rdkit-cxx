// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/depiction_lock.h"

#include "rdkit-cxx/control_c_lock.h"

namespace rdkit_cxx {

// The depiction lock comes before control_c_mutex (U4), so a closure that
// holds the latter must not wait for it.
std::unique_lock<std::mutex> lock_depiction() {
  check_not_in_control_c_callback();
  static std::mutex mutex;
  return std::unique_lock<std::mutex>(mutex);
}

}  // namespace rdkit_cxx
