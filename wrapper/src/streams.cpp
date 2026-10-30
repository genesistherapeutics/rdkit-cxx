// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/streams.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <string_view>
#include <utility>

#include <GraphMol/FileParsers/FileParserUtils.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

namespace {

constexpr std::size_t kChunk = std::size_t{1} << 16;
constexpr char kV30[] = "M  V30 ";
constexpr char kStereoGroup[] = "MDLV30/STE";
constexpr std::size_t kStereoGroupPrefix = sizeof(kStereoGroup) - 1;
constexpr char kEndPoints[] = "ENDPTS=";
constexpr std::size_t kEndPointsKey = sizeof(kEndPoints) - 1;
const char *const kBadCompno =
    "a V3000 line ends after a COMPNO label or \"COMPNO=\"; RDKit reads its "
    "value with an extraction that writes nothing there, and compares the "
    "unset value";
const char *const kBadStereoCollection =
    "a V3000 enhanced-stereo collection has fewer atoms than its count, so "
    "RDKit would read an atom number no extraction wrote";
// A key "COMPNO" or its "=", as the last token of a logical line.
constexpr std::string_view kCompno = "COMPNO=";

bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' ||
         c == '\r';
}

char upper(char c) {
  return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

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

// RDKit matches an enhanced-stereo collection line (parseEnhancedStereo in
// MolFileParser.cpp, on the upper-cased line) with
//   MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *
// which this matches by hand, as rdkit-2025.03's wrapper does: std::regex's
// lazily built tables are not ordered between threads where DRD (the races
// job) can see it.
bool match_stereo_collection(const std::string &line, std::string &type,
                             std::string &count, std::string &atoms) {
  static constexpr std::string_view prefix = "MDLV30/STE";
  auto digit = [](char c) { return c >= '0' && c <= '9'; };
  if (line.compare(0, prefix.size(), prefix) != 0 ||
      line.size() < prefix.size() + 3) {
    return false;
  }
  std::size_t pos = prefix.size();
  type = line.substr(pos, 3);
  pos += 3;
  while (pos < line.size() && digit(line[pos])) {
    ++pos;
  }
  const auto spaces = pos;
  while (pos < line.size() && line[pos] == ' ') {
    ++pos;
  }
  static constexpr std::string_view atoms_tag = "ATOMS=(";
  if (pos == spaces || line.compare(pos, atoms_tag.size(), atoms_tag) != 0) {
    return false;
  }
  pos += atoms_tag.size();
  const auto count_start = pos;
  while (pos < line.size() && digit(line[pos])) {
    ++pos;
  }
  if (pos == count_start) {
    return false;
  }
  count = line.substr(count_start, pos - count_start);
  const auto count_end = pos;
  while (pos < line.size() && line[pos] == ' ') {
    ++pos;
  }
  if (pos == count_end) {
    return false;
  }
  auto end = line.find_last_not_of(' ');
  if (end == std::string::npos || end < pos || line[end] != ')') {
    return false;
  }
  atoms = line.substr(pos, end - pos);
  return true;
}

// Whether RDKit 2026.09 would read an enhanced-stereo collection's first
// atom number unset: `ss >> index` writes nothing when only whitespace is
// left, so a count above 0 with no atoms (ATOMS=(1 )) reads an
// uninitialized index. A shorter list that has a first atom repeats it,
// which RDKit 2026.09 refuses as a duplicate.
bool stereo_collection_unset(const std::string &line) {
  std::string type;
  std::string count_text;
  std::string atoms;
  if (!match_stereo_collection(line, type, count_text, atoms) ||
      (type != "ABS" && type != "REL" && type != "RAC")) {
    return false;
  }
  unsigned int count;
  try {
    count = RDKit::FileParserUtils::toUnsigned(count_text, true);
  } catch (const std::exception &) {
    return false;  // RDKit throws
  }
  return count > 0 && std::all_of(atoms.begin(), atoms.end(),
                                  [](char c) { return is_space(c); });
}

}  // namespace

bool InputScanner::fail(std::string why, std::uint64_t at) {
  if (d_error.empty()) {
    d_error = std::move(why);
    d_refuse_at = at;
  }
  return false;
}

std::uint64_t InputScanner::held_from() const {
  std::uint64_t held = UINT64_MAX;
  if (d_logical > 0 && d_prefix == kStereoGroup) {
    held = d_logical_at;
  }
  if (d_token && d_token_spaced && !d_token_long &&
      kCompno.substr(0, d_token_text.size()) == d_token_text) {
    held = std::min(held, d_token_at);
  }
  return held;
}

std::uint64_t InputScanner::safe_end() const {
  return std::min(d_fed, held_from());
}

// A byte of the logical V3000 line, for COMPNO (streams.h): the line's last
// token, the whitespace before it, and whether a space is the last of that
// whitespace (the separator RDKit requires before a label). A COMPNO label
// is "COMPNO" or "COMPNO=" then, followed only by whitespace.
void InputScanner::compno(char c, std::uint64_t at) {
  if (is_space(c)) {
    d_in_space = true;
    d_space = c;
    d_space_at = at;
    return;
  }
  if (d_in_space || !d_token) {
    d_token = true;
    d_token_spaced = d_in_space && d_space == ' ';
    d_token_at = d_token_spaced ? d_space_at : at;
    d_token_long = false;
    d_token_text.clear();
    d_in_space = false;
  }
  if (d_token_text.size() < kCompno.size()) {
    d_token_text += upper(c);
  } else {
    d_token_long = true;
  }
}

// The end of a logical V3000 line: a COMPNO label may not end it, nor may
// an enhanced-stereo collection RDKit reads an unset atom number from.
void InputScanner::end_fields() {
  if (d_token && d_token_spaced && !d_token_long &&
      (d_token_text == kCompno ||
       d_token_text == kCompno.substr(0, kCompno.size() - 1))) {
    fail(kBadCompno, d_token_at);
  }
  if (d_logical > 0 && d_prefix == kStereoGroup &&
      stereo_collection_unset(d_stereo_line)) {
    fail(kBadStereoCollection, d_logical_at);
  }
  d_token = false;
  d_in_space = false;
  d_stereo_line.clear();
}

// A byte of the logical V3000 line, known not to be the line's trailing "\r"
// or continuation "-". RDKit upper-cases the joined line before matching it.
void InputScanner::confirm(char c, std::uint64_t at) {
  if (d_logical++ == 0) {
    d_logical_at = at;
  }
  if (d_prefix.size() < kStereoGroupPrefix) {
    d_prefix += upper(c);
  }
  if (d_logical > kMaxStereoGroupLine && d_prefix == kStereoGroup) {
    fail("a V3000 enhanced-stereo line is longer than " +
             std::to_string(kMaxStereoGroupLine) +
             " characters; RDKit matches it with a regex that recurses once "
             "per character",
         d_logical_at);
  }
  if (std::string_view(kStereoGroup).substr(0, d_prefix.size()) == d_prefix) {
    d_stereo_line += upper(c);
  }
  end_points(c);
  compno(c, at);
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
      confirm(d_pend[i], d_pend_at[i]);
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
    end_fields();
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
    if (c == '\n') {
      end_line();
      if (++d_lines > kMaxInputLines) {
        fail("the input has more than " + std::to_string(kMaxInputLines) +
             " lines; RDKit counts lines in int");
      }
      if (!d_error.empty()) {
        d_fed += i;
        return i;
      }
      continue;
    }
    if (d_pos < 7) {
      d_head[d_pos] = c;
    } else if (std::memcmp(d_head, kV30, 7) == 0) {
      // Content of a V3000 line; the last two bytes may yet be stripped.
      if (d_pending == 2) {
        confirm(d_pend[0], d_pend_at[0]);
        d_pend[0] = d_pend[1];
        d_pend_at[0] = d_pend_at[1];
        d_pend[1] = c;
        d_pend_at[1] = d_fed + i;
      } else {
        d_pend_at[d_pending] = d_fed + i;
        d_pend[d_pending++] = c;
      }
      // A held field is delivered only once it is read past, so it may not
      // grow without bound.
      const auto held = held_from();
      if (held != UINT64_MAX && d_fed + i - held > kMaxHeld) {
        fail("a V3000 COMPNO field or enhanced-stereo line spans more than " +
                 std::to_string(kMaxHeld) +
                 " bytes, which the checked stream holds back whole",
             held);
      }
      if (!d_error.empty()) {
        d_fed += i;
        return i;
      }
    }
    ++d_pos;
  }
  d_fed += n;
  return n;
}

// RDKit reads a last line without its "\n" as it reads any other.
bool InputScanner::finish() {
  if (d_error.empty() && d_pos > 0) {
    end_line();
  }
  return d_error.empty();
}

void check_input_text(std::string_view text) {
  InputScanner scanner;
  if (scanner.feed(text.data(), text.size()) < text.size() ||
      !scanner.finish()) {
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

// RDKit asked for a byte it may not have. The stream stays failed.
CheckedStreambuf::int_type CheckedStreambuf::hit() {
  d_failed.store(true, std::memory_order_release);
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
  return traits_type::eof();
}

void CheckedStreambuf::changed() {
  refuse("the input changed while RDKit read it", 0);
  d_failed.store(true, std::memory_order_release);
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
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
    // Chunks read ahead before this one are read from the source again, if
    // ever, and checked against their hashes.
    d_ahead.erase(d_ahead.begin(), d_ahead.lower_bound(chunk));
  }
  auto ahead = d_ahead.find(chunk);
  if (ahead != d_ahead.end()) {
    // Bytes read_ahead scanned, which a source that cannot seek back would
    // not give again.
    if (d_have == 0) {
      std::memcpy(d_buffer.get(), ahead->second.data(), ahead->second.size());
      d_have = ahead->second.size();
    }
    d_ahead.erase(ahead);
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
    if (ok < fresh) {
      refuse(d_scanner.error(), std::min(d_scanned, d_scanner.refuse_at()));
    }
    if (d_scanned == base + kChunk) {
      d_hashes.push_back(siphash(d_key, d_buffer.get(), kChunk));
    } else {
      d_tail_hash = siphash(d_key, d_buffer.get(),
                            static_cast<std::size_t>(d_scanned - base));
    }
  }
  if (d_have < kChunk && d_scanned == base + d_have) {
    finish();
  }
  return true;
}

// The source ends at d_scanned: RDKit reads the unfinished line whole.
void CheckedStreambuf::finish() {
  if (!d_finished && d_scanned < d_limit) {
    d_finished = true;
    if (!d_scanner.finish()) {
      refuse(d_scanner.error(), std::min(d_scanned, d_scanner.refuse_at()));
    }
  }
}

// Scans chunk `chunk`, the one after the last scanned, into d_ahead, leaving
// the buffer as it is. False if the source changed.
bool CheckedStreambuf::read_ahead(std::uint64_t chunk) {
  using off = std::streambuf::off_type;
  const std::uint64_t base = chunk * kChunk;
  if (d_source_pos != base) {
    if (d_source->pubseekpos(static_cast<off>(base), std::ios_base::in) ==
        pos_type(off(-1))) {
      d_source_pos = UINT64_MAX;
      return true;  // nothing more to read
    }
    d_source_pos = base;
  }
  std::string bytes(kChunk, '\0');
  std::size_t have = 0;
  auto n = d_source->sgetn(bytes.data(), static_cast<std::streamsize>(kChunk));
  if (n > 0) {
    d_source_pos += static_cast<std::uint64_t>(n);
    have = static_cast<std::size_t>(n);
  }
  bytes.resize(have);
  if (have < kChunk && base + have < d_least_end) {
    return false;  // it ends sooner than it did
  }
  d_least_end = std::max(d_least_end, base + have);
  auto ok = d_scanner.feed(bytes.data(), have);
  d_scanned += ok;
  if (ok < have) {
    refuse(d_scanner.error(), std::min(d_scanned, d_scanner.refuse_at()));
  }
  if (d_scanned == base + kChunk) {
    d_hashes.push_back(siphash(d_key, bytes.data(), kChunk));
  } else {
    d_tail_hash = siphash(d_key, bytes.data(),
                          static_cast<std::size_t>(d_scanned - base));
  }
  if (have < kChunk && d_scanned == base + have) {
    finish();
  }
  d_ahead[chunk] = std::move(bytes);
  return true;
}

// The first byte that may not be delivered yet: past the limit, or in a
// field the scanner has not read past (until the input is known to end).
std::uint64_t CheckedStreambuf::deliverable() const {
  return d_finished ? d_limit : std::min(d_limit, d_scanner.safe_end());
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
  const std::uint64_t chunk = pos / kChunk;
  // First scan the chunks before it that a seek skipped.
  for (std::uint64_t c = d_scanned / kChunk; c < chunk && d_scanned < d_limit;
       ++c) {
    if (!fill(c)) {
      changed();
      return traits_type::eof();
    }
    if (d_have < kChunk) {
      return traits_type::eof();  // the source ends before pos
    }
  }
  if (pos >= d_limit) {
    return hit();
  }
  if (!fill(chunk)) {
    changed();
    return traits_type::eof();
  }
  const std::uint64_t base = chunk * kChunk;
  std::uint64_t end = std::min(base + d_have, deliverable());
  const bool held = pos >= end && pos < base + d_have;
  if (held) {
    // Scan ahead until the field the bytes from `pos` are in is read past,
    // the input ends, or it is refused.
    for (std::uint64_t c = d_scanned / kChunk;
         !d_finished && pos >= deliverable() && d_scanned == c * kChunk &&
         d_scanned < d_limit;
         ++c) {
      if (!read_ahead(c)) {
        changed();
        return traits_type::eof();
      }
    }
    end = std::min(base + d_have, deliverable());
  }
  if (pos >= end) {
    // Bytes still held were refused past them.
    return pos >= d_limit || held ? hit() : traits_type::eof();
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
    const std::uint64_t end = std::min(base + d_have, deliverable());
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
        changed();
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
