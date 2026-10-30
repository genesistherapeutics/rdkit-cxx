/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
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
// - RDKit 2025.09 reads a V3000 SGroup line with stream extractions it does
//   not check (ParseV3000SGroupsBlock and the label parsers in
//   MolSGroupParsing.cpp): it reads the index, external index, CSTATE count
//   and bond, COMPNO, and the first value of an array that an extraction
//   left unwritten when the line ran out, and pops the last character of an
//   empty SAP identifier. Each line it would read as an SGroup (as many as
//   the COUNTS line gives, after BEGIN SGROUP and a DEFAULT line) is
//   replayed as RDKit consumes it, going on where RDKit would stop at a
//   missing bookmark or a value it refuses, and refused where RDKit would
//   read such a value. When COUNTS gives no SGroups, RDKit throws when it
//   parses strictly, and otherwise reads SGroup lines until one's header
//   does not read, which the replay follows. (rdkit-2025.03's replay
//   refuses that case, where RDKit 2025.03 loops forever.)
// - RDKit reads a V3000 enhanced-stereo collection's atom numbers with
//   unchecked extractions, so the first is unset when the collection has
//   fewer than its count, and keeps an atom listed twice, which dangles once
//   the reader removes the atom; both are refused.
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
#include <map>
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
  // The input ends here: checks the unfinished line as RDKit would read it.
  // Scanning may go on if the input turns out to go on.
  void finish();
  const std::string &error() const { return d_error; }
  // The offset of the first byte that may not be delivered: the start of a
  // refused SGroup line, which RDKit must not read any of, or else the byte
  // that broke a limit. UINT64_MAX while nothing is refused.
  std::uint64_t refused_at() const { return d_refused_at; }
  // Whether a line of a V3000 CTAB is unfinished, so that it must be
  // scanned to its end before any of it is delivered.
  bool mid_ctab_line() const {
    return d_ctab != Ctab::kNone && (d_pos > 0 || d_continuing);
  }

 private:
  void end_line();
  void confirm(char c);
  void end_points(char c);
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
  // In the logical line: how much of "ENDPTS=" has matched, and then the
  // value after it, with the number being read and how many were read.
  enum class EndPoints : std::uint8_t { kKey, kStart, kNumbers, kClosed };
  EndPoints d_end_points = EndPoints::kKey;
  std::size_t d_key_matched = 0;
  bool d_in_number = false;
  std::uint64_t d_number = 0;
  std::uint64_t d_numbers = 0;
  // The logical V3000 line, and where it is in a CTAB: right after BEGIN
  // CTAB, after its COUNTS line (with the atom, bond, and SGroup counts),
  // right after BEGIN SGROUP, or among the SGroup lines, of which
  // `d_sgroups_left` are left, with the DEFAULT line's labels. `d_line_start`
  // is the offset of the logical line's first byte, of `d_fed` bytes fed before
  // this call.
  bool v3000_line(const std::string &line, bool advance);
  std::string d_line;
  std::uint64_t d_fed = 0;
  std::uint64_t d_line_start = 0;
  std::uint64_t d_refused_at = UINT64_MAX;
  enum class Ctab : std::uint8_t {
    kNone,
    kCtabStart,
    kCounted,
    kSGroupStart,
    kSGroups
  };
  Ctab d_ctab = Ctab::kNone;
  unsigned int d_atoms = 0;
  unsigned int d_bonds = 0;
  unsigned int d_sgroups = 0;
  unsigned int d_sgroups_left = 0;
  bool d_sgroups_until_fail = false;
  std::string d_defaults;
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
  bool read_ahead(std::uint64_t chunk);
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
  // Chunks scanned ahead of the one in the buffer, kept until it is their
  // turn, since the source may not seek back.
  std::map<std::uint64_t, std::string> d_ahead;
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
