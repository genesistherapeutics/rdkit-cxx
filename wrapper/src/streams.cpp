// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/streams.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

#include <GraphMol/FileParsers/FileParserUtils.h>
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

// Replays how RDKit 2025.09 consumes one V3000 SGroup line, and its DEFAULT
// line's labels (ParseV3000SGroupsBlock, ParseV3000ParseLabel, and the
// parsers they call in MolSGroupParsing.cpp), on the same kind of stream,
// without RDKit's side effects. It follows the lenient parse, which reads
// at least as far as the strict one, and goes on where RDKit would stop at
// a bookmark that is missing or a value it refuses (it does not have the
// molecule), so it may refuse lines RDKit would stop reading in time.
class SGroupReplay {
 public:
  SGroupReplay(std::string type, unsigned int atoms, unsigned int bonds)
      : d_type(std::move(type)), d_atoms(atoms), d_bonds(bonds) {}

  // Whether `s >> value` writes `value`: its sentry fails, and leaves it as
  // it was, when the stream is not good or only whitespace is left. The
  // sentry's whitespace is skipped as it would skip it.
  static bool writes(std::istream &s) {
    if (!s.good()) {
      s.setstate(std::ios_base::failbit);
      return false;
    }
    s >> std::ws;
    if (s.eof()) {
      s.setstate(std::ios_base::failbit);
      return false;
    }
    return true;
  }

  // Replays the labels after the header, then the DEFAULT line's; the
  // error, or an empty string.
  std::string labels(std::stringstream &line, const std::string &defaults) {
    std::set<std::string> parsed;
    while (!line.eof() && !line.fail()) {
      char spacer;
      line.get(spacer);
      if (line.gcount() == 0) {
        continue;
      }
      if (spacer != ' ') {
        return {};  // the SGroup is invalid; its defaults are not read
      }
      std::string label;
      std::getline(line, label, '=');
      if (!label_value(label, line)) {
        return d_error;
      }
      parsed.insert(label);
    }
    std::stringstream stream(defaults);
    while (!stream.eof() && !stream.fail()) {
      char spacer;
      stream.get(spacer);
      if (stream.gcount() == 0) {
        continue;
      }
      if (spacer != ' ') {
        return {};
      }
      std::string label;
      std::getline(stream, label, '=');
      if (parsed.find(label) == parsed.end()) {
        if (!label_value(label, stream)) {
          return d_error;
        }
        continue;
      }
      spacer = static_cast<char>(stream.peek());
      if (spacer == ' ') {
        return {};
      } else if (spacer == '(') {
        std::getline(stream, label, ')');
        stream.get(spacer);
      } else if (spacer == '"') {
        stream.get(spacer);
        std::getline(stream, label, '"');
      } else {
        std::getline(stream, label, ' ');
        stream.putback(' ');
      }
    }
    return {};
  }

 private:
  bool unset(const std::string &what) {
    d_error = "a V3000 SGroup line has no " + what +
              " where RDKit 2025.09 reads one, so it would read a value no "
              "extraction wrote";
    return false;
  }

  // ParseV3000Array: false if it reads its first value unwritten.
  template <class T>
  bool array(std::stringstream &s, int max_v, std::size_t *count_read) {
    s.get();
    unsigned int count = 0;
    s >> count;
    *count_read = 0;
    if (max_v >= 0 && count > static_cast<unsigned int>(max_v)) {
      return true;  // it warns and reads no values
    }
    *count_read = count;
    T value;
    bool written = false;
    for (unsigned int i = 0; i < count; ++i) {
      if (writes(s)) {
        s >> value;
        written = true;
      } else if (!written) {
        return unset("array value");
      } else {
        break;  // the rest repeat the last value, and read nothing
      }
    }
    s.get();
    return true;
  }

  // ParseV3000ParseLabel: false to stop, with d_error set if RDKit would
  // read an unwritten value.
  bool label_value(const std::string &label, std::stringstream &s) {
    std::size_t count = 0;
    if (label == "XBHEAD" || label == "XBCORR" || label == "CBONDS" ||
        label == "XBONDS") {
      return array<unsigned int>(s, static_cast<int>(d_bonds), &count);
    }
    if (label == "ATOMS" || label == "PATOMS") {
      return array<unsigned int>(s, static_cast<int>(d_atoms), &count);
    }
    if (label == "BRKXYZ") {
      // Any count but 9 throws, which stops the line.
      return array<double>(s, 9, &count) && count == 9;
    }
    if (label == "CSTATE") {
      s.get();
      unsigned int number;
      unsigned int bond;
      if (!writes(s)) {
        return unset("CSTATE count");
      }
      s >> number;
      bool has_bond = writes(s);
      if (has_bond) {
        s >> bond;
      }
      if ((d_type != "SUP" && number != 1) ||
          (d_type == "SUP" && number != 4)) {
        return false;
      }
      if (!has_bond) {
        return unset("CSTATE bond");
      }
      if (d_type == "SUP") {
        double x = 0, y = 0, z = 0;
        s >> x >> y >> z;
      }
      s.get();
      return true;
    }
    if (label == "SAP") {
      s.get();
      unsigned int points = 0;
      unsigned int mark = 0;
      std::string leaving;
      std::string id;
      s >> points >> mark >> leaving >> id;
      if (id.empty()) {
        d_error =
            "a V3000 SGroup SAP value has no identifier; RDKit 2025.09 "
            "would remove the last character of an empty string";
        return false;
      }
      return true;
    }
    if (label == "PARENT") {
      if (s.eof()) {
        return false;
      }
      unsigned int parent;
      s >> parent;
      return !s.fail();
    }
    if (label == "COMPNO") {
      if (!writes(s)) {
        return unset("COMPNO value");
      }
      unsigned int compno;
      s >> compno;
      return compno <= 256u;
    }
    // ParseV3000StringPropLabel.
    std::string value;
    auto next = s.peek();
    if (next == ' ') {
    } else if (next == '"') {
      s.get();
      char chr;
      while (s.get(chr)) {
        if (chr == '"') {
          next = s.peek();
          if (next != '"') {
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
    return true;
  }

  std::string d_type;
  unsigned int d_atoms;
  unsigned int d_bonds;
  std::string d_error;
};

// The error RDKit 2025.09's reading of an enhanced-stereo collection line
// (parseEnhancedStereo in MolFileParser.cpp, on the upper-cased line) would
// make, or an empty string. It reads the count's atom numbers with
// unchecked extractions, so the first may be unset and the rest repeat the
// last, and it keeps an atom listed twice, which then dangles once removed
// (check_stereo_group_duplicates in mol_checks.h): say, when the reader removes
// the hydrogens.
// RDKit matches the line with
//   MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *
// which this matches by hand: the wrapper keeps std::regex off this path,
// which every V3000 read takes, since DRD (the races job) cannot see
// libstdc++'s regex compiler's lazily built tables ordered between threads.
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

std::string stereo_collection_error(const std::string &line) {
  std::string type;
  std::string count_text;
  std::string atoms;
  if (!match_stereo_collection(line, type, count_text, atoms) ||
      (type != "ABS" && type != "REL" && type != "RAC")) {
    return {};
  }
  unsigned int count;
  try {
    count = RDKit::FileParserUtils::toUnsigned(count_text, true);
  } catch (const std::exception &) {
    return {};  // RDKit throws
  }
  std::stringstream stream(atoms);
  std::set<unsigned int> seen;
  unsigned int index = 0;
  bool written = false;
  for (unsigned int i = 0; i < count; ++i) {
    if (SGroupReplay::writes(stream)) {
      stream >> index;
      written = true;
    } else if (!written) {
      return "a V3000 enhanced-stereo collection has fewer atoms than its "
             "count, so RDKit 2025.09 would read an atom number no extraction "
             "wrote";
    }
    if (!seen.insert(index).second) {
      return "a V3000 enhanced-stereo collection lists an atom twice (or "
             "fewer atoms than its count, which RDKit 2025.09 reads as the "
             "last repeated); RDKit keeps the duplicate, which dangles once "
             "the atom is removed";
    }
  }
  return {};
}

// Whether RDKit 2025.09's reading of an SGroup header line (its index, type,
// and external index) fails, where, given no SGroup count, it stops reading
// SGroups.
bool sgroup_header_fails(const std::string &line) {
  std::stringstream stream(line);
  unsigned int index = 0;
  std::string type;
  unsigned int external = 0;
  stream >> index;
  stream >> type;
  stream >> external;
  return stream.fail();
}

// The error RDKit 2025.09's reading of an SGroup header line would make, or
// an empty string.
std::string sgroup_line_error(const std::string &line,
                              const std::string &defaults, unsigned int atoms,
                              unsigned int bonds) {
  std::stringstream stream(line);
  if (!SGroupReplay::writes(stream)) {
    return "a V3000 SGroup line has no index where RDKit 2025.09 reads one, "
           "so it would read a value no extraction wrote";
  }
  unsigned int index;
  stream >> index;
  std::string type;
  stream >> type;
  if (!SGroupReplay::writes(stream)) {
    return "a V3000 SGroup line has no type and external index where RDKit "
           "2025.09 reads them, so it would read a value no extraction wrote";
  }
  unsigned int external;
  stream >> external;
  return SGroupReplay(type, atoms, bonds).labels(stream, defaults);
}

}  // namespace

bool InputScanner::fail(std::string why) {
  if (d_error.empty()) {
    d_error = std::move(why);
  }
  return false;
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
    d_logical = 0;
    d_prefix.clear();
    d_end_points = EndPoints::kKey;
    d_key_matched = 0;
    if (!v30) {
      d_ctab = Ctab::kNone;  // RDKit's CTAB reader throws on such a line
    } else if (d_error.empty()) {
      v3000_line(d_line, true);
    }
    d_line.clear();
  }
}

// A logical V3000 line, as getV3000Line returns it. RDKit reads CTABs
// (ParseV3000CTAB) as BEGIN CTAB, then COUNTS with the atom, bond, and
// SGroup counts, then blocks; after BEGIN SGROUP (in any case), it reads
// as many lines as COUNTS gives as SGroups, after a DEFAULT line if one
// comes first. Unless `advance`, the line is only checked: the input ends
// with it.
bool InputScanner::v3000_line(const std::string &line, bool advance) {
  if (d_ctab == Ctab::kSGroupStart && line.compare(0, 7, "DEFAULT") == 0 &&
      line.size() > 8) {
    if (advance) {
      d_defaults = boost::trim_right_copy(line.substr(7));
      d_ctab = Ctab::kSGroups;
    }
    return true;
  }
  if (d_ctab == Ctab::kSGroupStart || d_ctab == Ctab::kSGroups) {
    // Without an SGroup count, RDKit 2025.09 reads SGroup lines until one's
    // header does not read, which it takes for the end of the block.
    if (d_sgroups_until_fail &&
        sgroup_header_fails(d_ctab == Ctab::kSGroupStart
                                ? line
                                : boost::trim_right_copy(line))) {
      if (advance) {
        d_ctab = Ctab::kCounted;
      }
      return true;
    }
    // RDKit trims the right of each SGroup line but a first one.
    auto error = sgroup_line_error(
        d_ctab == Ctab::kSGroupStart ? line : boost::trim_right_copy(line),
        d_defaults, d_atoms, d_bonds);
    if (!error.empty()) {
      d_refused_at = d_line_start;
      return fail(error);
    }
    if (advance) {
      d_ctab = d_sgroups_until_fail || --d_sgroups_left ? Ctab::kSGroups
                                                        : Ctab::kCounted;
    }
    return true;
  }
  auto upper = boost::to_upper_copy(line);
  if (d_ctab == Ctab::kCounted && upper.compare(0, 10, "MDLV30/STE") == 0) {
    auto error = stereo_collection_error(upper);
    if (!error.empty()) {
      d_refused_at = d_line_start;
      return fail(error);
    }
    return true;
  }
  if (upper.compare(0, 10, "BEGIN CTAB") == 0) {
    if (advance) {
      d_ctab = Ctab::kCtabStart;
    }
  } else if (d_ctab == Ctab::kCtabStart) {
    if (!advance) {
      return true;
    }
    // The line after BEGIN CTAB, which RDKit requires to be COUNTS.
    d_ctab = Ctab::kNone;
    if (upper.size() < 8 || upper.compare(0, 7, "COUNTS ") != 0) {
      return true;
    }
    std::vector<std::string> counts;
    auto trimmed = boost::trim_copy(upper.substr(7));
    boost::split(counts, trimmed, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    if (counts.size() < 2) {
      return true;
    }
    try {
      d_atoms = RDKit::FileParserUtils::toUnsigned(counts[0]);
      d_bonds = RDKit::FileParserUtils::toUnsigned(counts[1]);
      d_sgroups =
          counts.size() > 2 ? RDKit::FileParserUtils::toUnsigned(counts[2]) : 0;
    } catch (const std::exception &) {
      return true;  // RDKit throws
    }
    d_ctab = Ctab::kCounted;
  } else if (d_ctab == Ctab::kCounted && upper.size() >= 12 &&
             upper.compare(0, 12, "BEGIN SGROUP") == 0) {
    // With no SGroups on the COUNTS line, RDKit 2025.09 throws when it
    // parses strictly, and otherwise reads SGroup lines until one's header
    // does not read (ParseV3000CTAB and ParseV3000SGroupsBlock).
    if (advance) {
      d_sgroups_until_fail = d_sgroups == 0;
      d_sgroups_left = d_sgroups;
      d_defaults.clear();
      d_ctab = Ctab::kSGroupStart;
    }
  }
  return true;
}

void InputScanner::finish() {
  if (!d_error.empty() || d_ctab == Ctab::kNone ||
      (d_pos == 0 && !d_continuing)) {
    return;
  }
  // RDKit reads the unfinished line as a whole one: getLine strips one "\r",
  // and RDKit throws on a line that is not a V3000 one or ends in "-" (it
  // reads the next).
  if (d_pos < 7 || std::memcmp(d_head, kV30, 7) != 0) {
    return;
  }
  std::string line = d_line;
  line.append(d_pend, static_cast<std::size_t>(d_pending));
  if (!line.empty() && line.back() == '\r') {
    line.pop_back();
  }
  if (!line.empty() && line.back() == '-') {
    return;
  }
  v3000_line(line, false);
}

std::size_t InputScanner::feed(const char *bytes, std::size_t n) {
  if (!d_error.empty()) {
    return 0;
  }
  // A limit broken at byte i refuses from there, unless a refused line's
  // start was recorded.
  auto broken = [&](std::size_t i) {
    if (d_refused_at == UINT64_MAX) {
      d_refused_at = d_fed + i;
    }
    d_fed += i;
    return i;
  };
  for (std::size_t i = 0; i < n; ++i) {
    char c = bytes[i];
    if (c == '\n') {
      end_line();
      if (++d_lines > kMaxInputLines) {
        fail("the input has more than " + std::to_string(kMaxInputLines) +
             " lines; RDKit counts lines in int");
      }
      if (!d_error.empty()) {
        return broken(i);
      }
      if (!d_continuing) {
        d_line_start = d_fed + i + 1;
      }
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
        return broken(i);
      }
    }
    ++d_pos;
  }
  d_fed += n;
  return n;
}

void check_input_text(std::string_view text) {
  InputScanner scanner;
  scanner.feed(text.data(), text.size());
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
      refuse(d_scanner.error(), d_scanner.refused_at());
    }
    if (d_scanned == base + kChunk) {
      d_hashes.push_back(siphash(d_key, d_buffer.get(), kChunk));
    } else {
      d_tail_hash = siphash(d_key, d_buffer.get(),
                            static_cast<std::size_t>(d_scanned - base));
    }
  }
  // The source ends here, for now; RDKit reads an unfinished line whole.
  if (d_have < kChunk && d_scanned == base + d_have && d_scanned < d_limit) {
    d_scanner.finish();
    if (!d_scanner.error().empty()) {
      refuse(d_scanner.error(), d_scanner.refused_at());
    }
  }
  return true;
}

// Reads and scans chunk `chunk`, the next to scan, without the buffer, and
// keeps it for fill. False if the source changed; d_scanned shows how far
// it got.
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
    refuse(d_scanner.error(), d_scanner.refused_at());
  }
  if (d_scanned == base + kChunk) {
    d_hashes.push_back(siphash(d_key, bytes.data(), kChunk));
  } else {
    d_tail_hash = siphash(d_key, bytes.data(),
                          static_cast<std::size_t>(d_scanned - base));
  }
  if (have < kChunk && d_scanned == base + have && d_scanned < d_limit) {
    d_scanner.finish();
    if (!d_scanner.error().empty()) {
      refuse(d_scanner.error(), d_scanner.refused_at());
    }
  }
  d_ahead[chunk] = std::move(bytes);
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
  if (!fill(chunk)) {
    changed();
    return traits_type::eof();
  }
  // A CTAB line RDKit may read as an SGroup is scanned to its end (or the
  // source's) before any of it is delivered, so that a refused one is cut
  // before its first byte (InputScanner::refused_at): RDKit reads the part
  // of a line it is given as a whole line.
  for (std::uint64_t c = chunk + 1;
       d_scanner.mid_ctab_line() && d_scanned == c * kChunk &&
       d_scanned < d_limit;
       ++c) {
    if (!read_ahead(c)) {
      changed();
      return traits_type::eof();
    }
  }
  const std::uint64_t base = chunk * kChunk;
  const std::uint64_t end = std::min(base + d_have, d_limit);
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
