// Copyright (c) 2026 Steven Kearnes
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
#include <utility>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>

#include "rdkit-cxx/file_checks.h"

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

const char *const kCxField =
    "a mol block's \"M  \" line holds a CXSMILES field (SgD:, atomProp:, or "
    "Sg:) that RDKit 2024.03 would read past the end of if it ended the "
    "SMARTS it parses there";
constexpr std::uint64_t kNone = UINT64_MAX;

bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f';
}
bool is_digit(char c) { return c >= '0' && c <= '9'; }

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

}  // namespace

bool CxTokenMatcher::push(char c, std::uint64_t pos) {
  if (d_size == kLongest) {
    std::memmove(d_tail, d_tail + 1, kLongest - 1);
    std::memmove(d_pos, d_pos + 1, (kLongest - 1) * sizeof(std::uint64_t));
    --d_size;
  }
  d_tail[d_size] = c;
  d_pos[d_size] = pos;
  ++d_size;
  for (std::string_view name : {"SgD:", "atomProp:", "Sg:"}) {
    if (d_size >= name.size() &&
        std::string_view(d_tail + d_size - name.size(), name.size()) == name) {
      d_start = d_pos[d_size - name.size()];
      return true;
    }
  }
  return false;
}

std::uint64_t CxTokenMatcher::partial() const {
  for (std::size_t k = d_size; k > 0; --k) {
    std::string_view tail(d_tail + d_size - k, k);
    for (std::string_view name : {"SgD:", "atomProp:", "Sg:"}) {
      if (k < name.size() && name.substr(0, k) == tail) {
        return d_pos[d_size - k];
      }
    }
  }
  return kNone;
}

// The labels whose values are checked, and how RDKit 2024.03 reads each
// (ParseV3000ParseLabel in MolSGroupParsing.cpp). Each reads with
// operator>>, which leaves its variable as it was when it finds the end of
// the line before anything to read (and sets it to 0 when what it finds is
// not a number):
// - CSTATE (ParseV3000CStateLabel): one character, then the uninitialized
//   count and bond. The count must have something after its spaces, and,
//   when it is a number, the bond too.
// - SAP (ParseV3000SAPLabel): one character, two unsigned numbers, and two
//   strings, the last of which it pops a character from: a failed read,
//   of a number that is not digits or does not fit, leaves it empty.
// - PARENT and COMPNO: an uninitialized unsigned, with no other character.
// - ATOMS, PATOMS, CBONDS, XBONDS, XBHEAD, XBCORR, and BRKXYZ
//   (ParseV3000Array): one character and a count, then as many values
//   into one uninitialized variable: a nonzero count that ends the line
//   leaves the first one unread.
// A cut-off line is a line too, so the value is held from where a cut
// would leave a misread (the label's last character for the first three,
// the count's first nonzero digit for arrays) until it is safe.
bool SgroupValue::push(char c, std::uint64_t pos) {
  if (d_open != kNone && ++d_window > kMaxValueWindow) {
    return fail(
        "a V3000 S-group value has more white space than the scanner holds "
        "while it checks it");
  }
  switch (d_kind) {
    case Kind::kNone:
      break;
    case Kind::kCState:
      switch (d_state) {
        case 0:
          d_state = 1;
          return true;
        case 1:
          if (is_space(c)) {
            return true;
          }
          if (is_digit(c) || c == '+' || c == '-') {
            d_state = 2;
          } else {
            close();
          }
          return true;
        case 2:
          if (is_digit(c)) {
            return true;
          }
          if (is_space(c)) {
            d_state = 3;
          } else {
            close();
          }
          return true;
        default:
          if (!is_space(c)) {
            close();
          }
          return true;
      }
    case Kind::kSap: {
      const char *bad =
          "a V3000 SAP value is not a character, two unsigned numbers, and "
          "two words; RDKit 2024.03 then reads past a string it left empty";
      switch (d_state) {
        case 0:
          d_state = 1;
          return true;
        case 1:
        case 3:
          if (is_space(c)) {
            return true;
          }
          if (!is_digit(c)) {
            return fail(bad);
          }
          d_number = static_cast<std::uint64_t>(c - '0');
          ++d_state;
          return true;
        case 2:
        case 4:
          if (is_digit(c)) {
            d_number = d_number * 10 + static_cast<std::uint64_t>(c - '0');
            if (d_number > UINT32_MAX) {
              return fail(bad);
            }
            return true;
          }
          if (!is_space(c)) {
            return fail(bad);
          }
          ++d_state;
          return true;
        case 5:
          if (!is_space(c)) {
            d_state = 6;
          }
          return true;
        case 6:
          if (is_space(c)) {
            d_state = 7;
          }
          return true;
        default:
          if (!is_space(c)) {
            close();
          }
          return true;
      }
    }
    case Kind::kScalar:
      if (!is_space(c)) {
        close();
      }
      return true;
    case Kind::kArray:
      switch (d_state) {
        case 0:
          d_state = 1;
          return true;
        case 1:
          if (is_space(c)) {
            return true;
          }
          if (is_digit(c) || c == '+' || c == '-') {
            d_state = 2;
            if (c != '0') {
              d_open = pos;
            }
          } else {
            close();
          }
          return true;
        case 2:
          if (is_digit(c)) {
            if (c != '0' && d_open == kNone) {
              d_open = pos;
            }
            return true;
          }
          if (is_space(c)) {
            d_state = 3;
          } else {
            close();
          }
          return true;
        default:
          if (!is_space(c)) {
            close();
          }
          return true;
      }
  }
  // A label: the text after a space up to "=".
  if (c == ' ') {
    d_after_space = true;
    d_label_size = 0;
    d_label_long = false;
    d_open = kNone;
    return true;
  }
  if (!d_after_space) {
    return true;
  }
  std::string_view label(d_label, d_label_size);
  if (c == '=' && !d_label_long) {
    d_after_space = false;
    if (label == "CSTATE") {
      d_kind = Kind::kCState;
    } else if (label == "SAP") {
      d_kind = Kind::kSap;
    } else if (label == "PARENT" || label == "COMPNO") {
      d_kind = Kind::kScalar;
    } else if (label == "ATOMS" || label == "PATOMS" || label == "CBONDS" ||
               label == "XBONDS" || label == "XBHEAD" || label == "XBCORR" ||
               label == "BRKXYZ") {
      d_kind = Kind::kArray;
      d_why =
          "a V3000 S-group array's nonzero count ends the line; RDKit 2024.03 "
          "reads its first value uninitialized";
    }
    d_state = 0;
    d_window = 0;
    return true;
  }
  if (d_label_size < sizeof(d_label)) {
    d_label[d_label_size++] = c;
  } else {
    d_label_long = true;
  }
  label = std::string_view(d_label, d_label_size);
  d_open = kNone;
  d_window = 0;
  if (!d_label_long) {
    if (label == "CSTATE") {
      d_open = pos;
      d_why =
          "a V3000 CSTATE value ends before its count and bond; RDKit 2024.03 "
          "reads them uninitialized";
    } else if (label == "SAP") {
      d_open = pos;
      d_why =
          "a V3000 SAP value ends before two numbers and two words; RDKit "
          "2024.03 then reads past a string it left empty";
    } else if (label == "PARENT" || label == "COMPNO") {
      d_open = pos;
      d_why =
          "a V3000 PARENT or COMPNO value is empty; RDKit 2024.03 reads it "
          "uninitialized";
    }
  }
  return true;
}

bool SgroupValue::fail(const char *why) {
  d_why = why;
  if (d_open == kNone) {
    d_open = 0;
  }
  return false;
}

void SgroupValue::close() {
  d_kind = Kind::kNone;
  d_state = 0;
  d_open = kNone;
  d_window = 0;
  d_after_space = false;
  d_label_size = 0;
  d_label_long = false;
}

bool InputScanner::fail(std::string why) {
  if (d_error.empty()) {
    d_error = std::move(why);
    // Not inside a value or field a cut would leave misread; otherwise
    // where the last feed ended safely.
    d_cut = hold() == kNone ? d_cur : d_safe;
  }
  return false;
}

bool InputScanner::fail_at(std::string why, std::uint64_t start) {
  if (d_error.empty()) {
    d_error = std::move(why);
    d_cut = std::max(d_safe, start);
  }
  return false;
}

// Where the logical V3000 line would be misread if the input ended here,
// with its pending bytes read as RDKit would read them then; with
// !dangers_only, also the start of a CXSMILES field name it may be
// completing. `why` says why for a danger.
std::uint64_t InputScanner::logical_hold(bool dangers_only,
                                         const char **why) const {
  bool v30 = d_pos >= 7 && std::memcmp(d_head, kV30, 7) == 0;
  if (!v30) {
    // Mid-way through a continuation's "M  V30 ", RDKit would throw.
    return kNone;
  }
  int pending = d_pending;
  if (pending > 0 && d_pend[pending - 1] == '\r') {
    --pending;
  }
  if (pending > 0 && d_pend[pending - 1] == '-') {
    return kNone;  // a continuation with no line after it: RDKit throws
  }
  LogicalChecks checks = d_checks;
  for (int i = 0; i < pending; ++i) {
    if (!checks.value.push(d_pend[i], d_pend_pos[i])) {
      *why = checks.value.why();
      return checks.value.start();
    }
    if (checks.tokens.push(d_pend[i], d_pend_pos[i])) {
      *why = kCxField;
      return checks.tokens.start();
    }
  }
  std::uint64_t result = kNone;
  if (checks.value.open()) {
    *why = checks.value.why();
    result = checks.value.start();
  }
  if (!dangers_only) {
    result = std::min(result, checks.tokens.partial());
  }
  return result;
}

std::uint64_t InputScanner::hold() const {
  const char *why = "";
  std::uint64_t result = logical_hold(false, &why);
  result = std::min(result, d_mline_tokens.partial());
  result = std::min(result, d_data_tokens.partial());
  return result;
}

bool InputScanner::finish() {
  if (!d_error.empty()) {
    return false;
  }
  const char *why = "";
  std::uint64_t start = logical_hold(true, &why);
  if (start != kNone) {
    return fail_at(why, start);
  }
  return true;
}

// A byte of a V2000 "M  " line after its first seven. RDKit joins the data
// of SCD and SED lines (columns 12 to 80 of each, after getLine drops a
// trailing "\r") until an SED line (ParseSGroupV2000SCDSEDLine).
void InputScanner::mol_line_byte(char c, std::uint64_t pos) {
  bool data = std::memcmp(d_head, "M  SCD", 6) == 0 ||
              std::memcmp(d_head, "M  SED", 6) == 0;
  if (!data) {
    if (d_mline_tokens.push(c, pos)) {
      fail_at(kCxField, d_mline_tokens.start());
    }
    return;
  }
  if (d_pos < 11 || d_pos > 79) {
    return;
  }
  if (d_data_cr) {
    d_data_cr = false;
    if (d_data_tokens.push('\r', d_data_cr_pos)) {
      fail_at(kCxField, d_data_tokens.start());
      return;
    }
  }
  if (c == '\r') {
    d_data_cr = true;
    d_data_cr_pos = pos;
    return;
  }
  if (d_data_tokens.push(c, pos)) {
    fail_at(kCxField, d_data_tokens.start());
  }
}

// A byte of the logical V3000 line, known not to be the line's trailing "\r"
// or continuation "-". RDKit upper-cases the joined line before matching it.
void InputScanner::confirm(char c, std::uint64_t pos) {
  ++d_logical;
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
  if (!d_checks.value.push(c, pos)) {
    fail_at(d_checks.value.why(), d_checks.value.start());
  } else if (d_checks.tokens.push(c, pos)) {
    fail_at(kCxField, d_checks.tokens.start());
  }
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
      confirm(d_pend[i], d_pend_pos[i]);
    }
  } else if (d_pos >= 7 && std::memcmp(d_head, "M  SED", 6) == 0) {
    d_data_tokens.reset();
  }
  d_mline_tokens.reset();
  d_data_cr = false;
  d_pending = 0;
  d_pos = 0;
  d_continuing = continues;
  if (!continues) {
    // A value that ends the line: one that is empty, or ends in ")" or not.
    if (d_end_points == EndPoints::kStart) {
      fail(kBadEndPoints);
    }
    if (d_checks.value.open()) {
      fail_at(d_checks.value.why(), d_checks.value.start());
    }
    d_checks.value.reset();
    d_checks.tokens.reset();
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
    const std::uint64_t pos = d_fed + i;
    d_cur = pos;
    if (c == '\n') {
      end_line();
      if (++d_lines > kMaxInputLines) {
        fail("the input has more than " + std::to_string(kMaxInputLines) +
             " lines; RDKit counts lines in int");
      }
      if (!d_error.empty()) {
        return i;
      }
      continue;
    }
    if (d_pos < 7) {
      d_head[d_pos] = c;
    } else if (std::memcmp(d_head, kV30, 7) == 0) {
      // Content of a V3000 line; the last two bytes may yet be stripped.
      if (d_pending == 2) {
        const char first = d_pend[0];
        const std::uint64_t first_pos = d_pend_pos[0];
        d_pend[0] = d_pend[1];
        d_pend_pos[0] = d_pend_pos[1];
        d_pend[1] = c;
        d_pend_pos[1] = pos;
        confirm(first, first_pos);
      } else {
        d_pend[d_pending] = c;
        d_pend_pos[d_pending] = pos;
        ++d_pending;
      }
    } else if (std::memcmp(d_head, "M  ", 3) == 0) {
      mol_line_byte(c, pos);
    }
    if (!d_error.empty()) {
      return i;
    }
    ++d_pos;
  }
  d_fed += n;
  const std::uint64_t held = hold();
  d_safe = held == kNone ? d_fed : held;
  return n;
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
    d_finished = false;
    if (ok < fresh) {
      refuse(d_scanner.error(), d_scanner.cut());
    }
    if (d_scanned == base + kChunk) {
      d_hashes.push_back(siphash(d_key, d_buffer.get(), kChunk));
    } else {
      d_tail_hash = siphash(d_key, d_buffer.get(),
                            static_cast<std::size_t>(d_scanned - base));
    }
  }
  // A short read is the end of the source: what the scanner held for bytes
  // to come is checked as RDKit will read it, cut off there.
  if (d_have < kChunk && d_scanned == base + d_have && d_scanned < d_limit &&
      !d_finished) {
    d_finished = true;
    if (!d_scanner.finish()) {
      refuse(d_scanner.error(), d_scanner.cut());
    }
  }
  if (d_scanned < d_limit) {
    // Once the source has ended and passed, nothing is held.
    d_hold = d_finished ? UINT64_MAX : d_scanner.hold();
    if (d_hold != UINT64_MAX && d_scanned - d_hold > kMaxHold) {
      refuse(
          "the input holds a value or CXSMILES field longer than the " +
              std::to_string(kMaxHold) + " bytes the scanner holds back",
          d_hold);
    }
  }
  return true;
}

// RDKit reached held bytes, from `pos`: they are copied out of the buffer
// (with those of earlier chunks the hold covers) and the next chunk is
// scanned, until the hold moves past `pos` or the input is refused. False
// if the source changed.
bool CheckedStreambuf::release(std::uint64_t pos) {
  while (d_hold != UINT64_MAX && pos >= d_hold && d_scanned < d_limit) {
    const std::uint64_t first = d_hold / kChunk;
    const std::uint64_t last = (d_scanned - 1) / kChunk;
    std::string carry;
    for (std::uint64_t c = first; c <= last; ++c) {
      const std::uint64_t base = c * kChunk;
      const std::uint64_t from = std::max(base, d_hold);
      const std::uint64_t to = std::min(base + kChunk, d_scanned);
      if (from >= to) {
        continue;
      }
      if (d_loaded != c && !d_carry.empty() && d_carry_pos <= from &&
          to <= d_carry_pos + d_carry.size()) {
        carry.append(d_carry, static_cast<std::size_t>(from - d_carry_pos),
                     static_cast<std::size_t>(to - from));
        continue;
      }
      if (!fill(c) || base + d_have < to) {
        return false;
      }
      carry.append(d_buffer.get() + (from - base),
                   static_cast<std::size_t>(to - from));
    }
    d_carry = std::move(carry);
    d_carry_pos = d_hold;
    const std::uint64_t scanned = d_scanned;
    if (!fill(d_scanned / kChunk)) {
      return false;
    }
    if (d_scanned == scanned && d_hold != UINT64_MAX && pos >= d_hold &&
        d_scanned < d_limit) {
      // Nothing more to scan, and the end was not seen: refuse.
      refuse("the input ended inside a value the scanner held back", d_hold);
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
  if (!release(pos)) {
    changed();
    return traits_type::eof();
  }
  if (pos >= d_limit) {
    return hit();
  }
  // Held bytes copied out of the buffer.
  if (!d_carry.empty() && pos >= d_carry_pos &&
      pos < d_carry_pos + d_carry.size() && d_loaded != chunk) {
    const std::uint64_t end =
        std::min({d_carry_pos + d_carry.size(), d_limit, d_hold});
    char *begin = d_carry.data();
    d_offset = d_carry_pos;
    setg(begin, begin + (pos - d_carry_pos), begin + (end - d_carry_pos));
    return traits_type::to_int_type(*gptr());
  }
  if (!fill(chunk)) {
    changed();
    return traits_type::eof();
  }
  const std::uint64_t base = chunk * kChunk;
  const std::uint64_t end = std::min({base + d_have, d_limit, d_hold});
  if (pos >= end) {
    return pos >= d_limit ? hit() : traits_type::eof();
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
    const std::uint64_t end = std::min({base + d_have, d_limit, d_hold});
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
