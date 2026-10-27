/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// The streams RDKit's readers read from. Every input the wrappers hand to
// RDKit's text readers (mol blocks, SD files, SMILES and TDT files) is checked
// as RDKit reads it (CONVENTIONS.md U12):
//
// - RDKit counts lines in int (ForwardSDMolSupplier::d_line,
//   TDTMolSupplier::d_line, SmilesMolSupplier::d_line, and
//   MultithreadedMolSupplier::d_line), which overflows past INT_MAX lines.
// - RDKit matches each V3000 enhanced-stereo line ("MDLV30/STE...", joined
//   across "-" continuations) with std::regex (parseEnhancedStereo in
//   MolFileParser.cpp), whose matcher recurses once per character: about
//   300 bytes of stack each, so a long line overflows the stack.
// - RDKit keeps a V3000 bond's ENDPTS value as read, and RWMol's atom
//   removals (removeHs, which the readers run when they sanitize and remove
//   hydrogens) rewrite it on every bond: they read its first character
//   without checking that it has one and, for a value in parentheses, the
//   first space-separated number without checking that there is one, and
//   convert each number's std::stod to unsigned int, which is undefined for
//   one that does not fit. A value after "ENDPTS=" (in any case) must not be
//   empty, and one that starts with "(" must be numbers below 2^32 between
//   spaces and a ")" that ends the token. Only a bond block reads ENDPTS, and
//   RDKit reads an END BOND line after it, so input that ends (or is cut
//   off) before a value is checked is never parsed as a molecule.
// - RDKit 2024.09 reads a V3000 substance group's fields with operator>>,
//   unchecked, so a field with too few values leaves variables
//   uninitialized and pops the last character of an empty string
//   (check_sgroup_line in streams.cpp). Each line RDKit reads as a group's
//   (the first nSgroups lines after a BEGIN SGROUP line, after a DEFAULT
//   line, with the counts of the CTAB's COUNTS line) is read as RDKit reads
//   it when it ends, and a line that would be refused is refused from its
//   start. RDKit parses whatever line the input ends in, so when a refusal
//   leaves RDKit within such a line, the stream first delivers
//   kSGroupSuffix, which gives every field RDKit may then read its values,
//   and only then reads as end of file.
//
// A CheckedStreambuf scans the bytes of its source in order, once each, and
// stops delivering them (reading as end of file) before the first byte that
// breaks a limit, recording why. Bytes are scanned ahead of what RDKit
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
#include <utility>
#include <vector>

namespace rdkit_cxx {

// The most lines an input may have; RDKit's counters stay below INT_MAX.
constexpr std::uint64_t kMaxInputLines = (std::uint64_t{1} << 31) - (1 << 20);
// The longest V3000 enhanced-stereo line, in characters after "M  V30 ".
constexpr std::uint64_t kMaxStereoGroupLine = 4000;
// What a refused stream delivers within a substance group's line: four
// values, each of them a number, for whatever field RDKit is reading.
constexpr char kSGroupSuffix[] = " 0 0 0 0";

// Scans text in order, one byte at a time across calls.
class InputScanner {
 public:
  // Scans bytes; returns how many are within the limits (all of them, or
  // the bytes before the first that breaks one, after which error() says
  // why).
  std::size_t feed(const char *bytes, std::size_t n);
  // The input ended: checks the line it ends in, as RDKit reads it.
  void finish();
  const std::string &error() const { return d_error; }
  // Where the input is refused from: the first byte of the line that broke
  // a limit, or the byte that did.
  std::uint64_t error_at() const { return d_error_at; }
  // Whether `pos` falls after the start of a line RDKit reads as a
  // substance group's.
  bool within_sgroup_line(std::uint64_t pos) const;

 private:
  void end_line();
  void end_logical_line();
  void confirm(char c);
  void end_points(char c);
  bool fail(std::string why);
  bool fail_line(std::string why);
  bool d_finished = false;

  std::string d_error;
  std::uint64_t d_error_at = 0;
  // The bytes scanned so far.
  std::uint64_t d_offset = 0;
  std::uint64_t d_lines = 0;
  // The logical V3000 line's content as RDKit reads it, and where its first
  // physical line starts.
  std::string d_line;
  std::uint64_t d_line_start = 0;
  // The CTAB counts last read, and the substance-group block: whether the
  // next line may be its DEFAULT line, its defaults, and how many group
  // lines are left (0 outside one).
  unsigned d_atoms = 0;
  unsigned d_bonds = 0;
  unsigned d_sgroups = 0;
  bool d_block_start = false;
  bool d_first_group = false;
  std::string d_defaults;
  std::uint64_t d_groups_left = 0;
  // The lines RDKit reads as groups', from their first byte to their last
  // newline (UINT64_MAX while open).
  std::vector<std::pair<std::uint64_t, std::uint64_t>> d_sgroup_lines;
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
  // In the logical line: how much of "ENDPTS=" has matched, and then the
  // value after it, with the number being read and how many were read.
  enum class EndPoints : std::uint8_t { kKey, kStart, kNumbers, kClosed };
  EndPoints d_end_points = EndPoints::kKey;
  std::size_t d_key_matched = 0;
  bool d_in_number = false;
  std::uint64_t d_number = 0;
  std::uint64_t d_numbers = 0;
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
  // Ends the input where RDKit reads next, as a refusal does; safe to call
  // from another thread.
  void stop() { d_stopped.store(true, std::memory_order_release); }

 protected:
  int_type underflow() override;
  pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                   std::ios_base::openmode which) override;
  pos_type seekpos(pos_type pos, std::ios_base::openmode which) override;

 private:
  bool fill(std::uint64_t chunk);
  void refuse(const std::string &why, std::uint64_t at);
  int_type hit(std::uint64_t pos);
  int_type changed(std::uint64_t pos);
  int_type end_at(std::uint64_t pos);
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
  // A copy of kSGroupSuffix to deliver, once.
  char d_suffix[sizeof(kSGroupSuffix)] = {};
  bool d_suffixed = false;
  std::atomic<bool> d_failed{false};
  std::atomic<bool> d_stopped{false};
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
  void stop() { d_buf.stop(); }

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
