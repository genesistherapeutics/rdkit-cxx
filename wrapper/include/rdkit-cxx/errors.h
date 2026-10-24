// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause
//
// Error handling shared by every bridge. Every bridge includes this header
// first (CONVENTIONS.md E5), so cxx uses the trycatch below for every
// function, and RDKit's exceptions reach Rust with their type and payload
// encoded in the message (see src/error.rs).

#pragma once

#include <exception>
#include <string>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Thrown by wrappers that check an index argument (E3).
class IndexError : public std::exception {
 public:
  IndexError(unsigned index, unsigned len);
  const char *what() const noexcept override { return d_msg.c_str(); }
  unsigned index() const { return d_index; }
  unsigned len() const { return d_len; }

 private:
  unsigned d_index;
  unsigned d_len;
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
