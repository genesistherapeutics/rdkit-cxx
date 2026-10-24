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

#include <RDGeneral/Exceptions.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Bytes (a pickle, a file path) as a std::string.
inline std::string to_std(rust::Slice<const std::uint8_t> bytes) {
  return std::string(reinterpret_cast<const char *>(bytes.data()),
                     bytes.size());
}

// A file path's bytes (A8). C++ opens a path at its first NUL, which would
// be another file, so a path holding one is refused, as Python's open does.
inline std::string to_path(rust::Slice<const std::uint8_t> bytes) {
  auto path = to_std(bytes);
  if (path.find('\0') != std::string::npos) {
    throw ValueErrorException("embedded null byte in a file path");
  }
  return path;
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
