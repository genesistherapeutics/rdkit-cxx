// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/streams.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstring>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <utility>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>

namespace rdkit_cxx {

namespace {

constexpr std::size_t kChunk = std::size_t{1} << 16;
constexpr char kV30[] = "M  V30 ";
constexpr char kStereoGroup[] = "MDLV30/STE";
constexpr std::size_t kStereoGroupPrefix = sizeof(kStereoGroup) - 1;
constexpr char kEndPoints[] = "ENDPTS=";
constexpr std::size_t kEndPointsKey = sizeof(kEndPoints) - 1;
const char *const kBadEndPoints =
    "a V3000 ENDPTS value is empty, or starts with \"(\" but is not numbers "
    "below 2^32 between spaces and a \")\" that ends the token; RWMol reads it "
    "unchecked when it removes atoms";

std::uint64_t rotl(std::uint64_t x, int b) {
  return (x << b) | (x >> (64 - b));
}

// SipHash-2-4 (Aumasson and Bernstein, 2012): input cannot be made to collide
// without the key.
std::uint64_t siphash(const std::uint64_t key[2], const char *data,
                      std::size_t n) {
  std::uint64_t v0 = 0x736f6d6570736575ULL ^ key[0];
  std::uint64_t v1 = 0x646f72616e646f6dULL ^ key[1];
  std::uint64_t v2 = 0x6c7967656e657261ULL ^ key[0];
  std::uint64_t v3 = 0x7465646279746573ULL ^ key[1];
  auto round = [&] {
    v0 += v1;
    v1 = rotl(v1, 13) ^ v0;
    v0 = rotl(v0, 32);
    v2 += v3;
    v3 = rotl(v3, 16) ^ v2;
    v0 += v3;
    v3 = rotl(v3, 21) ^ v0;
    v2 += v1;
    v1 = rotl(v1, 17) ^ v2;
    v2 = rotl(v2, 32);
  };
  const auto *p = reinterpret_cast<const unsigned char *>(data);
  const std::size_t whole = n - n % 8;
  for (std::size_t i = 0; i < whole; i += 8) {
    std::uint64_t m = 0;
    for (int j = 7; j >= 0; --j) {
      m = (m << 8) | p[i + j];
    }
    v3 ^= m;
    round();
    round();
    v0 ^= m;
  }
  std::uint64_t last = static_cast<std::uint64_t>(n) << 56;
  for (std::size_t j = 0; j < n % 8; ++j) {
    last |= static_cast<std::uint64_t>(p[whole + j]) << (8 * j);
  }
  v3 ^= last;
  round();
  round();
  v0 ^= last;
  v2 ^= 0xff;
  for (int i = 0; i < 4; ++i) {
    round();
  }
  return v0 ^ v1 ^ v2 ^ v3;
}

// RDKit 2024.09's V3000 substance-group reader (ParseV3000SGroupsBlock and
// ParseV3000ParseLabel in MolSGroupParsing.cpp) reads its fields with
// operator>> and checks few of the reads: a read that finds nothing left
// leaves its variable as it was, uninitialized for the group's sequence and
// external ids, a CSTATE's count and bond, a COMPNO, and an array's first
// value, and a SAP field read that way pops the last character of an empty
// string. RDKit 2025.09 checks them (rdkit/rdkit#9417). These read a group's
// line as RDKit does, field by field, and say whether RDKit would read such
// a value. Where RDKit would stop reading the line (a check that throws, or
// a lookup of an atom or bond that may), they read on, so they check at
// least the fields RDKit reads.
struct SGroupLimits {
  unsigned atoms = 0;
  unsigned bonds = 0;
};

// Whether a read with operator>> would find something: its sentry skips
// whitespace and fails at the end (or on a failed stream), leaving the
// variable unset.
bool sgroup_readable(std::stringstream &s) {
  if (!s.good()) {
    return false;
  }
  s >> std::ws;
  return s.good();
}

// As ParseV3000Array: false if its first value would be left unset.
template <class T>
bool sgroup_array(std::stringstream &s, int max) {
  s.get();
  unsigned count = 0;
  s >> count;
  if (max >= 0 && count > static_cast<unsigned>(max)) {
    return true;  // RDKit returns before it reads values (non-strict)
  }
  for (unsigned i = 0; i < count; ++i) {
    if (!sgroup_readable(s)) {
      // Later reads keep the value an earlier one set.
      if (i == 0) {
        return false;
      }
      break;
    }
    T value;
    s >> value;
    if (s.fail()) {
      break;  // a failed conversion sets 0, which later reads keep
    }
  }
  s.get();
  return true;
}

// As ParseV3000StringPropLabel, for the read alone.
void sgroup_string(std::stringstream &s) {
  auto next = s.peek();
  if (next == ' ') {
    return;
  }
  std::string value;
  if (next == '"') {
    s.get();
    char chr;
    while (s.get(chr)) {
      if (chr == '"') {
        if (s.peek() != '"') {
          break;
        }
        s.get();
      }
    }
  } else if (next == '\'') {
    std::getline(s, value, '\'');
  } else {
    s >> value;
  }
}

// As ParseV3000ParseLabel: an error, or empty.
std::string sgroup_field(const std::string &label, std::stringstream &s,
                         const std::string &type, const SGroupLimits &limits) {
  const std::string bad = "a V3000 substance group's " + label + " field ";
  if (label == "XBHEAD" || label == "XBCORR" || label == "CBONDS" ||
      label == "XBONDS") {
    if (!sgroup_array<unsigned>(s, static_cast<int>(limits.bonds))) {
      return bad + "lists no value after its count";
    }
  } else if (label == "ATOMS" || label == "PATOMS") {
    if (!sgroup_array<unsigned>(s, static_cast<int>(limits.atoms))) {
      return bad + "lists no value after its count";
    }
  } else if (label == "BRKXYZ") {
    if (!sgroup_array<double>(s, 9)) {
      return bad + "lists no value after its count";
    }
  } else if (label == "CSTATE") {
    s.get();
    if (!sgroup_readable(s)) {
      return bad + "has no count";
    }
    unsigned count = 0;
    s >> count;
    bool expected = type == "SUP" ? count == 4 : count == 1;
    if (!sgroup_readable(s) && expected) {
      return bad + "has no bond";
    }
    unsigned bond = 0;
    s >> bond;
    if (type == "SUP") {
      double x, y, z;
      s >> x >> y >> z;
    }
    s.get();
  } else if (label == "SAP") {
    s.get();
    unsigned count = 0;
    unsigned atom = 0;
    std::string leaving;
    std::string id;
    s >> count >> atom >> leaving >> id;
    if (id.empty()) {
      return bad + "has fewer than four values";
    }
  } else if (label == "PARENT") {
    if (!s.eof()) {
      unsigned parent = 0;
      s >> parent;
    }
  } else if (label == "COMPNO") {
    if (!sgroup_readable(s)) {
      return bad + "has no value";
    }
    unsigned compno = 0;
    s >> compno;
  } else {
    sgroup_string(s);
  }
  return {};
}

// As the label loops of ParseV3000SGroupsBlock, non-strict (which reads on
// past a separator other than a space only to stop): the labels it reads go
// into `parsed`, and with `skip`, the defaults' labels in it are skipped as
// RDKit skips them.
std::string sgroup_labels(std::stringstream &s, const std::string &type,
                          const SGroupLimits &limits,
                          std::set<std::string> &parsed,
                          const std::set<std::string> *skip) {
  while (!s.eof() && !s.fail()) {
    char spacer;
    s.get(spacer);
    if (s.gcount() == 0) {
      continue;
    }
    if (spacer != ' ') {
      return {};
    }
    std::string label;
    std::getline(s, label, '=');
    if (skip && skip->count(label)) {
      auto next = s.peek();
      std::string rest;
      if (next == ' ') {
        return {};
      } else if (next == '(') {
        std::getline(s, rest, ')');
        s.get(spacer);
      } else if (next == '"') {
        s.get(spacer);
        std::getline(s, rest, '"');
      } else {
        std::getline(s, rest, ' ');
        s.putback(' ');
      }
      continue;
    }
    auto error = sgroup_field(label, s, type, limits);
    if (!error.empty()) {
      return error;
    }
    parsed.insert(label);
  }
  return {};
}

// A group's line, with its defaults' line (empty for none): an error, or
// empty.
std::string check_sgroup_line(const std::string &line,
                              const std::string &defaults,
                              const SGroupLimits &limits) {
  std::stringstream s(line);
  if (!sgroup_readable(s)) {
    return "a V3000 substance group's line has no sequence id";
  }
  unsigned sequence = 0;
  s >> sequence;
  std::string type;
  s >> type;
  if (!sgroup_readable(s)) {
    return "a V3000 substance group's line has no external id";
  }
  unsigned external = 0;
  s >> external;
  std::set<std::string> parsed;
  auto error = sgroup_labels(s, type, limits, parsed, nullptr);
  if (!error.empty()) {
    return error;
  }
  std::stringstream d(defaults);
  std::set<std::string> ignored;
  return sgroup_labels(d, type, limits, ignored, &parsed);
}

}  // namespace

bool InputScanner::fail(std::string why) {
  if (d_error.empty()) {
    d_error = std::move(why);
    d_error_at = d_offset;
  }
  return false;
}

// Refuses the logical line from its first byte.
bool InputScanner::fail_line(std::string why) {
  if (d_error.empty()) {
    d_error = std::move(why);
    d_error_at = d_line_start;
  }
  return false;
}

bool InputScanner::within_sgroup_line(std::uint64_t pos) const {
  // The lines are in order; find the last starting before pos.
  auto it = std::upper_bound(
      d_sgroup_lines.begin(), d_sgroup_lines.end(), pos,
      [](std::uint64_t p, const std::pair<std::uint64_t, std::uint64_t> &l) {
        return p <= l.first;
      });
  if (it == d_sgroup_lines.begin()) {
    return false;
  }
  --it;
  return pos > it->first && pos <= it->second;
}

// A logical V3000 line ended (or the input did, within one), with its
// content in d_line: the CTAB counts, and the substance-group block, as
// ParseV3000CTAB and ParseV3000SGroupsBlock read them.
void InputScanner::end_logical_line() {
  std::string upper = boost::to_upper_copy(d_line);
  if (d_groups_left > 0) {
    if (d_block_start && d_line.compare(0, 7, "DEFAULT") == 0 &&
        d_line.size() > 8) {
      d_defaults = boost::trim_right_copy(d_line.substr(7));
    } else {
      // The first group's line is not trimmed after no DEFAULT line.
      std::string line = d_first_group && d_block_start
                             ? d_line
                             : boost::trim_right_copy(d_line);
      auto error = check_sgroup_line(line, d_defaults, {d_atoms, d_bonds});
      if (!error.empty()) {
        fail_line(error + "; RDKit 2024.09 reads it unchecked");
      }
      --d_groups_left;
      d_first_group = false;
    }
    d_block_start = false;
  } else if (upper.compare(0, 7, "COUNTS ") == 0 && upper.size() >= 8) {
    // As ParseV3000CTAB, with FileParserUtils::toUnsigned, which throws
    // for a character other than a digit or "+" (stopping RDKit) and
    // converts with std::from_chars, leaving 0 where it fails.
    std::vector<std::string> counts;
    auto trimmed = boost::trim_copy(upper.substr(7));
    boost::split(counts, trimmed, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    auto number = [&](std::size_t i) -> unsigned {
      if (i >= counts.size()) {
        return 0;
      }
      const auto &text = counts[i];
      auto end = text.find('\0');
      if (end == std::string::npos) {
        end = text.size();
      }
      if (text.substr(0, end).find_first_not_of("0123456789+") !=
          std::string::npos) {
        return 0;
      }
      unsigned res = 0;
      std::from_chars(text.data(), text.data() + text.size(), res);
      return res;
    };
    d_atoms = number(0);
    d_bonds = number(1);
    d_sgroups = number(2);
  } else if (d_sgroups > 0 && upper.compare(0, 12, "BEGIN SGROUP") == 0) {
    d_groups_left = d_sgroups;
    d_sgroups = 0;
    d_block_start = true;
    d_first_group = true;
    d_defaults.clear();
  }
  d_line.clear();
}

void InputScanner::finish() {
  if (!d_error.empty() || d_finished) {
    return;
  }
  d_finished = true;
  bool v30 = d_pos >= 7 && std::memcmp(d_head, kV30, 7) == 0;
  if (v30 || d_continuing) {
    end_line();
  }
}

// A byte of the logical V3000 line, known not to be the line's trailing "\r"
// or continuation "-". RDKit upper-cases the joined line before matching it.
void InputScanner::confirm(char c) {
  ++d_logical;
  d_line += c;
  if (d_prefix.size() < kStereoGroupPrefix) {
    d_prefix += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  }
  if (d_logical > kMaxStereoGroupLine && d_prefix == kStereoGroup) {
    fail("a V3000 enhanced-stereo line is longer than " +
         std::to_string(kMaxStereoGroupLine) +
         " characters; RDKit matches it with a regex that recurses once per "
         "character");
  }
  end_points(c);
}

// A byte of the logical V3000 line, for the ENDPTS values in it. Values are
// checked wherever "ENDPTS=" appears (splitAssignToken upper-cases the
// property's name), not only at the start of a bond line's token.
void InputScanner::end_points(char c) {
  switch (d_end_points) {
    case EndPoints::kStart:
      if (c == '(') {
        d_end_points = EndPoints::kNumbers;
        d_in_number = false;
        d_number = 0;
        d_numbers = 0;
        return;
      }
      // A space or tab ends the token, and a quote may close a quoted one
      // (tokenizeV3000Line), so the value is empty; so it is if only "\r"s,
      // which RDKit strips, are left on the line.
      if (c == ' ' || c == '\t' || c == '"') {
        fail(kBadEndPoints);
        return;
      }
      if (c == '\r') {
        return;
      }
      // RWMol skips a value that does not start with "(".
      d_end_points = EndPoints::kKey;
      break;
    case EndPoints::kNumbers:
      if (c >= '0' && c <= '9') {
        d_in_number = true;
        d_number = d_number * 10 + static_cast<std::uint64_t>(c - '0');
        if (d_number > UINT32_MAX) {
          fail(kBadEndPoints);
        }
      } else if (c == ' ' || c == ')') {
        if (d_in_number) {
          ++d_numbers;
          d_in_number = false;
          d_number = 0;
        }
        if (c == ')') {
          if (d_numbers == 0) {
            fail(kBadEndPoints);
          }
          d_end_points = EndPoints::kClosed;
        }
      } else {
        fail(kBadEndPoints);
      }
      return;
    case EndPoints::kClosed:
      if (c != ' ' && c != '\t') {
        fail(kBadEndPoints);
      }
      d_end_points = EndPoints::kKey;
      return;
    case EndPoints::kKey:
      break;
  }
  char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (upper == kEndPoints[d_key_matched]) {
    if (++d_key_matched == kEndPointsKey) {
      d_end_points = EndPoints::kStart;
      d_key_matched = 0;
    }
  } else {
    d_key_matched = upper == kEndPoints[0] ? 1 : 0;
  }
}

// As getV3000Line in MolFileParser.cpp: strip one "\r"; a trailing "-"
// continues the logical line on the next physical one.
void InputScanner::end_line() {
  bool v30 = d_pos >= 7 && std::memcmp(d_head, kV30, 7) == 0;
  bool continues = false;
  if (v30) {
    if (d_pending > 0 && d_pend[d_pending - 1] == '\r') {
      --d_pending;
    }
    if (d_pending > 0 && d_pend[d_pending - 1] == '-') {
      --d_pending;
      continues = true;
    }
    for (int i = 0; i < d_pending; ++i) {
      confirm(d_pend[i]);
    }
  }
  d_pending = 0;
  d_pos = 0;
  d_continuing = continues;
  if (!continues) {
    // A value that ends the line: one that is empty, or ends in ")" or not.
    if (d_end_points == EndPoints::kStart) {
      fail(kBadEndPoints);
    }
    if (!d_sgroup_lines.empty() && d_sgroup_lines.back().second == UINT64_MAX) {
      d_sgroup_lines.back().second = d_offset;
    }
    if (v30) {
      end_logical_line();
    } else {
      // RDKit reads a CTAB's lines as V3000 lines, throwing at any other,
      // so one ends what it was reading.
      d_line.clear();
      d_sgroups = 0;
      d_groups_left = 0;
      d_block_start = false;
    }
    d_logical = 0;
    d_prefix.clear();
    d_end_points = EndPoints::kKey;
    d_key_matched = 0;
  }
}

std::size_t InputScanner::feed(const char *bytes, std::size_t n) {
  if (!d_error.empty()) {
    return 0;
  }
  for (std::size_t i = 0; i < n; ++i) {
    char c = bytes[i];
    if (d_pos == 0 && !d_continuing) {
      d_line_start = d_offset;
      if (d_groups_left > 0) {
        d_sgroup_lines.emplace_back(d_offset, UINT64_MAX);
      }
    }
    if (c == '\n') {
      end_line();
      if (++d_lines > kMaxInputLines) {
        fail("the input has more than " + std::to_string(kMaxInputLines) +
             " lines; RDKit counts lines in int");
      }
      if (!d_error.empty()) {
        return i;
      }
      ++d_offset;
      continue;
    }
    if (d_pos < 7) {
      d_head[d_pos] = c;
    } else if (std::memcmp(d_head, kV30, 7) == 0) {
      // Content of a V3000 line; the last two bytes may yet be stripped.
      if (d_pending == 2) {
        confirm(d_pend[0]);
        d_pend[0] = d_pend[1];
        d_pend[1] = c;
      } else {
        d_pend[d_pending++] = c;
      }
      if (!d_error.empty()) {
        return i;
      }
    }
    ++d_pos;
    ++d_offset;
  }
  return n;
}

void check_input_text(std::string_view text) {
  InputScanner scanner;
  if (scanner.feed(text.data(), text.size()) < text.size()) {
    throw ValueErrorException(scanner.error());
  }
  scanner.finish();
  if (!scanner.error().empty()) {
    throw ValueErrorException(scanner.error());
  }
}

CheckedStreambuf::CheckedStreambuf(std::unique_ptr<std::streambuf> source)
    : d_source(std::move(source)), d_buffer(new char[kChunk]) {
  std::random_device random;
  for (auto &k : d_key) {
    k = (std::uint64_t{random()} << 32) ^ random();
  }
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
}

CheckedStreambuf::~CheckedStreambuf() = default;

std::string CheckedStreambuf::error() const {
  std::lock_guard<std::mutex> lock(d_error_mutex);
  return d_error;
}

// Bytes from `at` on may not be delivered, because of `why`.
void CheckedStreambuf::refuse(const std::string &why, std::uint64_t at) {
  if (at >= d_limit) {
    return;
  }
  d_limit = at;
  std::lock_guard<std::mutex> lock(d_error_mutex);
  d_error = why;
}

// The failed stream ends where RDKit reads, at `pos`: after kSGroupSuffix
// if that is within a substance group's line (InputScanner), which RDKit
// would otherwise parse cut off.
CheckedStreambuf::int_type CheckedStreambuf::end_at(std::uint64_t pos) {
  if (!d_suffixed && d_scanner.within_sgroup_line(pos)) {
    d_suffixed = true;
    std::memcpy(d_suffix, kSGroupSuffix, sizeof(kSGroupSuffix));
    setg(d_suffix, d_suffix, d_suffix + sizeof(kSGroupSuffix) - 1);
    return traits_type::to_int_type(*gptr());
  }
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
  return traits_type::eof();
}

// RDKit asked for a byte it may not have, at `pos`. The stream stays failed.
CheckedStreambuf::int_type CheckedStreambuf::hit(std::uint64_t pos) {
  d_failed.store(true, std::memory_order_release);
  return end_at(pos);
}

CheckedStreambuf::int_type CheckedStreambuf::changed(std::uint64_t pos) {
  refuse("the input changed while RDKit read it", 0);
  d_failed.store(true, std::memory_order_release);
  return end_at(pos);
}

// Reads chunk `chunk` into the buffer, up to the chunk's end or the source's.
// Bytes of it scanned before must read as they did; the rest are scanned now.
// False if the source changed.
bool CheckedStreambuf::fill(std::uint64_t chunk) {
  using off = std::streambuf::off_type;
  const std::uint64_t base = chunk * kChunk;
  if (d_loaded != chunk) {
    d_loaded = chunk;
    d_have = 0;
  }
  if (d_have == kChunk) {
    return true;
  }
  const std::size_t had = d_have;
  const std::uint64_t from = base + had;
  if (d_source_pos != from) {
    if (d_source->pubseekpos(static_cast<off>(from), std::ios_base::in) ==
        pos_type(off(-1))) {
      d_source_pos = UINT64_MAX;
      return true;  // nothing more to read
    }
    d_source_pos = from;
  }
  auto n = d_source->sgetn(d_buffer.get() + had,
                           static_cast<std::streamsize>(kChunk - had));
  if (n > 0) {
    d_source_pos += static_cast<std::uint64_t>(n);
    d_have += static_cast<std::size_t>(n);
  }
  if (d_have < kChunk && base + d_have < d_least_end) {
    return false;  // it ends sooner than it did
  }
  d_least_end = std::max(d_least_end, base + d_have);
  std::size_t known = 0;
  std::uint64_t expected = 0;
  if (chunk < d_hashes.size()) {
    known = kChunk;
    expected = d_hashes[chunk];
  } else if (chunk == d_hashes.size() && d_scanned > base) {
    known = static_cast<std::size_t>(d_scanned - base);
    expected = d_tail_hash;
  }
  if (d_have < known) {
    return false;  // shorter than it was
  }
  if (had < known && siphash(d_key, d_buffer.get(), known) != expected) {
    return false;
  }
  if (d_scanned >= base && d_scanned < base + d_have && d_scanned < d_limit) {
    auto seen = static_cast<std::size_t>(d_scanned - base);
    std::size_t fresh = d_have - seen;
    auto ok = d_scanner.feed(d_buffer.get() + seen, fresh);
    d_scanned += ok;
    if (!d_scanner.error().empty()) {
      refuse(d_scanner.error(), d_scanner.error_at());
    }
    if (d_scanned == base + kChunk) {
      d_hashes.push_back(siphash(d_key, d_buffer.get(), kChunk));
    } else {
      d_tail_hash = siphash(d_key, d_buffer.get(),
                            static_cast<std::size_t>(d_scanned - base));
    }
  }
  if (d_have < kChunk && d_scanned == base + d_have && d_scanned < d_limit) {
    // The source ended: RDKit reads the line it ends in.
    d_scanner.finish();
    if (!d_scanner.error().empty()) {
      refuse(d_scanner.error(), d_scanner.error_at());
    }
  }
  return true;
}

CheckedStreambuf::int_type CheckedStreambuf::underflow() {
  if (gptr() < egptr()) {
    return traits_type::to_int_type(*gptr());
  }
  if (failed()) {
    return traits_type::eof();
  }
  const std::uint64_t pos =
      d_offset + static_cast<std::uint64_t>(egptr() - eback());
  d_offset = pos;
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
  if (d_stopped.load(std::memory_order_acquire)) {
    refuse("the reader was stopped", pos);
    return hit(pos);
  }
  const std::uint64_t chunk = pos / kChunk;
  // First scan the chunks before it that a seek skipped.
  for (std::uint64_t c = d_scanned / kChunk; c < chunk && d_scanned < d_limit;
       ++c) {
    if (!fill(c)) {
      return changed(pos);
    }
    if (d_have < kChunk) {
      return traits_type::eof();  // the source ends before pos
    }
  }
  if (pos >= d_limit) {
    return hit(pos);
  }
  if (!fill(chunk)) {
    return changed(pos);
  }
  const std::uint64_t base = chunk * kChunk;
  const std::uint64_t end = std::min(base + d_have, d_limit);
  if (pos >= end) {
    return pos >= d_limit ? hit(pos) : traits_type::eof();
  }
  d_offset = base;
  setg(d_buffer.get(), d_buffer.get() + (pos - base),
       d_buffer.get() + (end - base));
  return traits_type::to_int_type(*gptr());
}

// The source moved to `pos`. A position in the chunk in the buffer is
// delivered from it.
CheckedStreambuf::pos_type CheckedStreambuf::moved(pos_type pos) {
  if (pos == pos_type(off_type(-1))) {
    d_source_pos = UINT64_MAX;
    return pos;
  }
  const auto to = static_cast<std::uint64_t>(off_type(pos));
  d_source_pos = to;
  if (d_loaded != UINT64_MAX && !failed()) {
    const std::uint64_t base = d_loaded * kChunk;
    const std::uint64_t end = std::min(base + d_have, d_limit);
    if (to >= base && to < end) {
      d_offset = base;
      setg(d_buffer.get(), d_buffer.get() + (to - base),
           d_buffer.get() + (end - base));
      return pos;
    }
  }
  d_offset = to;
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
  return pos;
}

CheckedStreambuf::pos_type CheckedStreambuf::seekoff(
    off_type off, std::ios_base::seekdir dir, std::ios_base::openmode which) {
  if (!(which & std::ios_base::in)) {
    return pos_type(off_type(-1));
  }
  if (dir == std::ios_base::cur) {
    auto here = static_cast<off_type>(d_offset) + (gptr() - eback());
    if (off == 0) {
      return pos_type(here);
    }
    return seekpos(pos_type(here + off), which);
  }
  if (dir == std::ios_base::end) {
    // RDKit reads up to the end it finds, so the source may not have shrunk.
    auto end = d_source->pubseekoff(0, std::ios_base::end, std::ios_base::in);
    if (end != pos_type(off_type(-1))) {
      if (static_cast<std::uint64_t>(off_type(end)) < d_least_end) {
        changed(0);
        return moved(pos_type(off_type(-1)));
      }
      d_least_end = static_cast<std::uint64_t>(off_type(end));
      if (off == 0) {
        return moved(end);
      }
    }
  }
  return moved(d_source->pubseekoff(off, dir, std::ios_base::in));
}

CheckedStreambuf::pos_type CheckedStreambuf::seekpos(
    pos_type pos, std::ios_base::openmode which) {
  if (!(which & std::ios_base::in) || off_type(pos) < 0) {
    return pos_type(off_type(-1));
  }
  return moved(d_source->pubseekpos(pos, std::ios_base::in));
}

CheckedIStream::CheckedIStream(std::unique_ptr<std::streambuf> source)
    : std::istream(nullptr), d_buf(std::move(source)) {
  rdbuf(&d_buf);
}

void CheckedIStream::check() const {
  if (d_buf.failed()) {
    throw ValueErrorException(d_buf.error());
  }
}

std::unique_ptr<CheckedIStream> open_checked_file(const std::string &name,
                                                  bool nonempty) {
  auto file = std::make_unique<std::filebuf>();
  if (!file->open(name, std::ios_base::in | std::ios_base::binary)) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  if (nonempty && file->sgetc() == std::filebuf::traits_type::eof()) {
    throw RDKit::BadFileException("Invalid input file " + name);
  }
  return std::make_unique<CheckedIStream>(std::move(file));
}

std::unique_ptr<CheckedIStream> checked_text(std::string text) {
  return std::make_unique<CheckedIStream>(std::make_unique<std::stringbuf>(
      std::move(text), std::ios_base::in | std::ios_base::binary));
}

}  // namespace rdkit_cxx
