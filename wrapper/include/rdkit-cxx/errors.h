/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Error handling shared by every bridge. Every bridge includes this header
// first (CONVENTIONS.md E5), so cxx uses the trycatch below for every
// function, and RDKit's exceptions reach Rust with their type and payload
// encoded in the message (see src/error.rs).

#pragma once

#include <cstdint>
#include <exception>
#include <string>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Thrown by wrappers that check an index argument (E3). The index is signed,
// as Python's are, and Rust's `Error::Index` carries both values whole.
class IndexError : public std::exception {
 public:
  IndexError(std::int64_t index, std::uint64_t len);
  const char *what() const noexcept override { return d_msg.c_str(); }
  std::int64_t index() const { return d_index; }
  std::uint64_t len() const { return d_len; }

 private:
  std::int64_t d_index;
  std::uint64_t d_len;
  std::string d_msg;
};

// Throws IndexError unless index < len.
void check_index(unsigned index, unsigned len);

// Encodes the exception being handled as
// "\x1Frdkit\x1F<Variant>\x1F<payload>\x1F<what>". Call only from a catch
// block.
std::string encode_current_exception();

// Throws one exception of each kind, for tests of the mapping.
void raise_for_tests(rust::Str kind);

}  // namespace rdkit_cxx

namespace rust::behavior {

template <typename Try, typename Fail>
static void trycatch(Try &&func, Fail &&fail) noexcept try {
  func();
} catch (...) {
  fail(rdkit_cxx::encode_current_exception().c_str());
}

}  // namespace rust::behavior
