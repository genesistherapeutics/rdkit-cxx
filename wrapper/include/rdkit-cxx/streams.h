/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// The streams RDKit's readers read from. Every input the wrappers hand to
// RDKit's text readers (mol blocks, SD files, SMILES and TDT files) is checked
// as RDKit reads it (CONVENTIONS.md U9):
//
// - RDKit counts lines in int (ForwardSDMolSupplier::d_line,
//   TDTMolSupplier::d_line, SmilesMolSupplier::d_line, and
//   MultithreadedMolSupplier::d_line), which overflows past INT_MAX lines.
// - RDKit matches each V3000 enhanced-stereo line ("MDLV30/STE...", joined
//   across "-" continuations) with std::regex (parseEnhancedStereo in
//   MolFileParser.cpp), whose matcher recurses once per character: about
//   300 bytes of stack each, so a long line overflows the stack.
//
// A CheckedStreambuf scans the bytes of its source in order, once each, and
// stops delivering them (reading as end of file) before the first byte that
// breaks either limit, recording why. Bytes are scanned ahead of what RDKit
// asks for, so the stream fails only once RDKit asks for a refused byte, and
// stays failed. Callers check failed() after each call into RDKit and throw.

#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <memory>
#include <mutex>
#include <streambuf>
#include <string>
#include <string_view>
#include <vector>

namespace rdkit_cxx {

// The most lines an input may have; RDKit's counters stay below INT_MAX.
constexpr std::uint64_t kMaxInputLines = (std::uint64_t{1} << 31) - (1 << 20);
// The longest V3000 enhanced-stereo line, in characters after "M  V30 ".
constexpr std::uint64_t kMaxStereoGroupLine = 4000;

// Scans text in order, one byte at a time across calls.
class InputScanner {
 public:
  // Scans bytes; returns how many are within the limits (all of them, or
  // the bytes before the first that breaks one, after which error() says
  // why).
  std::size_t feed(const char *bytes, std::size_t n);
  const std::string &error() const { return d_error; }

 private:
  void end_line();
  void confirm(char c);
  bool fail(std::string why);

  std::string d_error;
  std::uint64_t d_lines = 0;
  // The physical line being read: its length so far, its first seven bytes,
  // and its last two (unconfirmed) bytes of V3000 content.
  std::uint64_t d_pos = 0;
  char d_head[7] = {};
  int d_pending = 0;
  char d_pend[2] = {};
  // The logical V3000 line, joined across continuations.
  bool d_continuing = false;
  std::uint64_t d_logical = 0;
  std::string d_prefix;
};

// Throws ValueErrorException if `text` breaks the limits above.
void check_input_text(std::string_view text);

// A streambuf over `source` that delivers only bytes InputScanner accepts.
// It reads the source in aligned chunks. Bytes past the furthest point read
// are scanned before they are delivered (a forward seek past it scans the
// skipped bytes first), and bytes read again after a seek must be the bytes
// scanned: each scanned chunk keeps a keyed hash, and a chunk that reads
// differently, or a source that ends before an end it has shown (a file that
// changed underneath), fails the stream.
class CheckedStreambuf : public std::streambuf {
 public:
  explicit CheckedStreambuf(std::unique_ptr<std::streambuf> source);
  ~CheckedStreambuf() override;
  CheckedStreambuf(const CheckedStreambuf &) = delete;
  CheckedStreambuf &operator=(const CheckedStreambuf &) = delete;

  // Whether RDKit asked for refused input; safe to call from another
  // thread.
  bool failed() const { return d_failed.load(std::memory_order_acquire); }
  std::string error() const;

 protected:
  int_type underflow() override;
  pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                   std::ios_base::openmode which) override;
  pos_type seekpos(pos_type pos, std::ios_base::openmode which) override;

 private:
  bool fill(std::uint64_t chunk);
  void refuse(const std::string &why, std::uint64_t at);
  int_type hit();
  void changed();
  pos_type moved(pos_type pos);

  std::unique_ptr<std::streambuf> d_source;
  std::unique_ptr<char[]> d_buffer;
  // Where the source is (UINT64_MAX if unknown), and the chunk in the
  // buffer, with how many of its bytes are there.
  std::uint64_t d_source_pos = 0;
  std::uint64_t d_loaded = UINT64_MAX;
  std::size_t d_have = 0;
  // The offset of the get area's first byte, the bytes scanned so far, and
  // the first byte that may not be delivered.
  std::uint64_t d_offset = 0;
  std::uint64_t d_scanned = 0;
  std::uint64_t d_limit = UINT64_MAX;
  // The size below which the source may not shrink: the bytes it has shown,
  // and the end a seek reported.
  std::uint64_t d_least_end = 0;
  // The hashes of each fully scanned chunk and of the scanned part of the
  // next, under a random key.
  std::vector<std::uint64_t> d_hashes;
  std::uint64_t d_tail_hash = 0;
  std::uint64_t d_key[2] = {};
  InputScanner d_scanner;
  std::atomic<bool> d_failed{false};
  mutable std::mutex d_error_mutex;
  std::string d_error;
};

// An istream that owns its CheckedStreambuf.
class CheckedIStream : public std::istream {
 public:
  explicit CheckedIStream(std::unique_ptr<std::streambuf> source);
  const CheckedStreambuf &checked() const { return d_buf; }
  // Throws ValueErrorException if RDKit asked for refused input.
  void check() const;

 private:
  CheckedStreambuf d_buf;
};

// Opens a file as RDKit's suppliers do (openAndCheckStream in MolSupplier.h):
// BadFileException("Bad input file <name>") if it cannot be opened, and, when
// `nonempty`, BadFileException("Invalid input file <name>") if it is empty.
std::unique_ptr<CheckedIStream> open_checked_file(const std::string &name,
                                                  bool nonempty);
// A checked stream over a copy of `text`.
std::unique_ptr<CheckedIStream> checked_text(std::string text);

}  // namespace rdkit_cxx
