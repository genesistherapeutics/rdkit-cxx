/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Conversions the wrappers share.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Bytes (a pickle, a file path) as a std::string.
inline std::string to_std(rust::Slice<const std::uint8_t> bytes) {
  return std::string(reinterpret_cast<const char *>(bytes.data()),
                     bytes.size());
}

inline rust::Vec<std::uint8_t> to_bytes(const std::string &text) {
  rust::Vec<std::uint8_t> bytes;
  bytes.reserve(text.size());
  for (char c : text) {
    bytes.push_back(static_cast<std::uint8_t>(c));
  }
  return bytes;
}

template <class T>
rust::Vec<T> to_vec(const std::vector<T> &v) {
  rust::Vec<T> result;
  result.reserve(v.size());
  for (const auto &x : v) {
    result.push_back(x);
  }
  return result;
}

}  // namespace rdkit_cxx
