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
// - RDKit 2024.03 reads several V3000 S-group values with operator>> into
//   variables it leaves uninitialized, or reads the result unchecked
//   (ParseV3000ParseLabel and the functions it calls in
//   MolSGroupParsing.cpp): CSTATE's count and bond (`CSTATE=(` and
//   `CSTATE=(4` read neither, `CSTATE=(4 ` no bond), the first value of
//   ATOMS, PATOMS, CBONDS, XBONDS, XBHEAD, XBCORR, and BRKXYZ arrays whose
//   nonzero count ends the line, PARENT's and COMPNO's value at the end of
//   the line, and SAP's last string, which it pops a character from even
//   when a read failed and left it empty (2026.03 initializes and checks
//   them). A value after " CSTATE=", " SAP=", " PARENT=", " COMPNO=", or
//   one of those arrays' names must give RDKit's reads something to read
//   before the logical line ends (SgroupValue in streams.cpp).
// - RDKit parses SMARTS in mol blocks, with CXSMILES on and strict: the
//   rest of an "M  MRV SMA" line and the data of a SMARTSQ data S-group
//   (V2000 "M  SCD" and "M  SED" lines, V3000 FIELDDATA). An "M  MRV SMA"
//   line is refused where its SMARTS's extension holds a "SgD:" or
//   "atomProp:" field RDKit 2024.03 reads past the end of
//   (cx_extension_problem in file_checks.h); V3000 lines (joined) and
//   V2000 SCD and SED data (joined as RDKit joins them), whichever S-group
//   they belong to, may not hold either field name.
// - RDKit 2024.03 reads a V3000 enhanced-stereo collection's atoms with
//   unchecked extractions and keeps an atom listed twice
//   (stereo_collection_error in streams.cpp), so such a line is refused.
//
// RDKit parses each line as soon as it has read it, and reads a line cut off
// by the end of the input (or by a refusal) as a whole one. So a byte is
// delivered only once a cut after it would be safe: the scanner holds back
// the bytes from the start of a value or field whose check needs bytes not
// yet scanned (hold()), and refuses (cut()) at the start of a held value,
// never inside one, nor before what it already delivered.
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
#include <vector>

namespace rdkit_cxx {

// The most lines an input may have; RDKit's counters stay below INT_MAX.
constexpr std::uint64_t kMaxInputLines = (std::uint64_t{1} << 31) - (1 << 20);
// The longest V3000 enhanced-stereo line, in characters after "M  V30 ".
constexpr std::uint64_t kMaxStereoGroupLine = 4000;

// The most bytes held back at once (see above); input that needs more is
// refused.
constexpr std::uint64_t kMaxHold = std::uint64_t{1} << 20;
// The longest V3000 S-group value window: white space RDKit skips inside a
// value whose check is pending.
constexpr std::uint64_t kMaxValueWindow = 4096;

// Finds "SgD:" and "atomProp:" in a stream of bytes.
class CxTokenMatcher {
 public:
  // Takes the byte at `pos`; true if it ends a field name, whose start
  // start() then gives.
  bool push(char c, std::uint64_t pos);
  std::uint64_t start() const { return d_start; }
  // The start of the longest field-name prefix the bytes end with, or
  // UINT64_MAX.
  std::uint64_t partial() const;
  void reset() { d_size = 0; }

 private:
  static constexpr std::size_t kLongest = 9;
  char d_tail[kLongest] = {};
  std::uint64_t d_pos[kLongest] = {};
  std::size_t d_size = 0;
  std::uint64_t d_start = 0;
};

// A V3000 logical line's S-group values (see above): what RDKit's reads of
// the value after a label would find if the line ended at each byte.
class SgroupValue {
 public:
  // Takes the byte at `pos` of the logical line; false if RDKit would
  // misread the value (why() and start() then say why and where it began).
  bool push(char c, std::uint64_t pos);
  // Whether the line ending here leaves RDKit a value it misreads.
  bool open() const { return d_open != UINT64_MAX; }
  std::uint64_t start() const { return d_open; }
  const char *why() const { return d_why; }
  void reset() { *this = SgroupValue(); }

 private:
  enum class Kind : std::uint8_t { kNone, kCState, kSap, kScalar, kArray };
  bool fail(const char *why);
  void close();
  Kind d_kind = Kind::kNone;
  int d_state = 0;
  // The label since the last space (at most 7 characters, then too long).
  char d_label[7] = {};
  std::size_t d_label_size = 0;
  bool d_label_long = false;
  bool d_after_space = false;
  std::uint64_t d_label_end = 0;
  std::uint64_t d_number = 0;
  std::uint64_t d_open = UINT64_MAX;
  std::uint64_t d_window = 0;
  const char *d_why = "";
};

// The start of a V3000 S-group line, which RDKit 2024.03 reads as
// `lineStream >> sequenceId >> type >> externalId`
// (ParseV3000SGroupsBlock in MolSGroupParsing.cpp) into variables it never
// initializes: an extraction that fails stores 0 but leaves every later one
// unread, and one at the end of the line reads nothing, so a line that does
// not start with a number and two more tokens (or a cut one) leaves the ids
// uninitialized, and RDKit stores and compares them.
// The block's first line may instead be its DEFAULT line.
class SgroupHeader {
 public:
  void start(bool first) {
    *this = SgroupHeader();
    d_active = true;
    d_state = first ? kDefault : kLead;
  }
  // Takes a byte of the logical line; false if RDKit would misread the ids.
  bool push(char c);
  // Whether a line ending here leaves RDKit ids it misreads.
  bool open() const { return d_active && d_state != kDone; }
  // Whether the line was the block's DEFAULT line, not an S-group.
  bool is_default() const { return d_is_default; }
  bool active() const { return d_active; }

 private:
  enum State : std::uint8_t {
    kDefault,
    kLead,
    kLeadSign,
    kFirst,
    kBeforeType,
    kType,
    kBeforeId,
    kDone
  };
  bool d_active = false;
  bool d_is_default = false;
  std::uint8_t d_default_matched = 0;
  State d_state = kLead;
};

// What a V3000 logical line holds that a cut-off line could misread.
struct LogicalChecks {
  SgroupValue value;
  CxTokenMatcher tokens;
  SgroupHeader header;
};

// Scans text in order, one byte at a time across calls.
class InputScanner {
 public:
  // Scans bytes; returns how many are within the limits (all of them, or
  // the bytes before the first that breaks one, after which error() says
  // why).
  std::size_t feed(const char *bytes, std::size_t n);
  // The input ended after the bytes fed. False if that breaks a limit.
  bool finish();
  const std::string &error() const { return d_error; }
  // After a failure, the first byte that may not be delivered.
  std::uint64_t cut() const { return d_cut; }
  // The first byte that may not be delivered until more is scanned, or
  // UINT64_MAX.
  std::uint64_t hold() const;

 private:
  void end_line();
  void confirm(char c, std::uint64_t pos);
  void end_points(char c);
  bool fail(std::string why);
  bool fail_at(std::string why, std::uint64_t start);
  std::uint64_t logical_hold(bool dangers_only, const char **why) const;
  void mol_line_byte(char c, std::uint64_t pos);
  void end_logical_line();
  std::string mrv_problem() const;

  std::string d_error;
  std::uint64_t d_lines = 0;
  // The offset of the next byte, where a failure cuts the input, and the
  // first byte not yet delivered as safe at the end of the last feed.
  std::uint64_t d_fed = 0;
  std::uint64_t d_cur = 0;
  std::uint64_t d_cut = UINT64_MAX;
  std::uint64_t d_safe = 0;
  // The physical line being read: its length so far, its first seven bytes,
  // and its last two (unconfirmed) bytes of V3000 content, with their
  // offsets.
  std::uint64_t d_pos = 0;
  char d_head[7] = {};
  int d_pending = 0;
  char d_pend[2] = {};
  std::uint64_t d_pend_pos[2] = {};
  // The V3000 logical line's values and CXSMILES field names, a V2000 "M  "
  // line's field names, and those of the SCD and SED data joined.
  LogicalChecks d_checks;
  CxTokenMatcher d_mline_tokens;
  CxTokenMatcher d_data_tokens;
  // A "\r" in SCD or SED data, which getLine drops if it ends the line.
  bool d_data_cr = false;
  std::uint64_t d_data_cr_pos = 0;
  // An "M  MRV " line: whether it is "M  MRV SMA", its bytes after the
  // first seven (at most kMaxHold of them, and whether there were more), and
  // the start of its first CXSMILES field name, or UINT64_MAX.
  bool d_mrv = false;
  bool d_mrv_long = false;
  std::string d_mrv_line;
  std::uint64_t d_mrv_danger = UINT64_MAX;
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
  // Where the logical line began, its first bytes (for the block keywords),
  // and the S-group block it is in: the S-groups the CTAB's COUNTS line
  // names, whether a BEGIN SGROUP line opened a block of them, and how many
  // S-group lines have been read.
  std::uint64_t d_logical_start = 0;
  std::string d_line;
  std::uint64_t d_num_sgroups = 0;
  bool d_in_sgroups = false;
  bool d_sgroup_first = false;
  std::uint64_t d_sgroups_read = 0;
  // Whether the logical line is in a collection block, whose stereo groups
  // RDKit reads.
  bool d_in_collection = false;
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
  bool release(std::uint64_t pos);
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
  // The first byte held back (InputScanner::hold), and a copy of held bytes
  // from chunks no longer in the buffer, from d_carry_pos; whether the
  // scanner was told that the source ended.
  std::uint64_t d_hold = UINT64_MAX;
  std::string d_carry;
  std::uint64_t d_carry_pos = 0;
  bool d_finished = false;
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
