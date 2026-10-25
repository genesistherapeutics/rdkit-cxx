// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/streams.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <utility>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

namespace {

constexpr std::size_t kChunk = std::size_t{1} << 16;
constexpr char kV30[] = "M  V30 ";
constexpr char kStereoGroup[] = "MDLV30/STE";
constexpr std::size_t kStereoGroupPrefix = sizeof(kStereoGroup) - 1;

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
  if (d_prefix.size() < kStereoGroupPrefix) {
    d_prefix += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  }
  if (d_logical > kMaxStereoGroupLine && d_prefix == kStereoGroup) {
    fail("a V3000 enhanced-stereo line is longer than " +
         std::to_string(kMaxStereoGroupLine) +
         " characters; RDKit matches it with a regex that recurses once per "
         "character");
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
    d_logical = 0;
    d_prefix.clear();
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
        return i;
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
        return i;
      }
    }
    ++d_pos;
  }
  return n;
}

void check_input_text(std::string_view text) {
  InputScanner scanner;
  if (scanner.feed(text.data(), text.size()) < text.size()) {
    throw ValueErrorException(scanner.error());
  }
}

CheckedStreambuf::CheckedStreambuf(std::unique_ptr<std::streambuf> source)
    : d_source(std::move(source)), d_buffer(new char[kChunk]) {
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
}

CheckedStreambuf::~CheckedStreambuf() = default;

std::string CheckedStreambuf::error() const {
  std::lock_guard<std::mutex> lock(d_error_mutex);
  return d_error;
}

void CheckedStreambuf::refuse(const std::string &why, std::uint64_t at) {
  {
    std::lock_guard<std::mutex> lock(d_error_mutex);
    if (d_error.empty()) {
      d_error = why;
    }
  }
  d_limit = std::min(d_limit, at);
  d_failed.store(true, std::memory_order_release);
}

// Scans the bytes from d_scanned up to `offset`, which a seek skipped, and
// leaves the source at `offset`.
bool CheckedStreambuf::scan_to(std::uint64_t offset) {
  using off = std::streambuf::off_type;
  if (d_source->pubseekpos(static_cast<off>(d_scanned), std::ios_base::in) ==
      pos_type(off(-1))) {
    return false;
  }
  while (d_scanned < offset) {
    auto want = static_cast<std::streamsize>(
        std::min<std::uint64_t>(kChunk, offset - d_scanned));
    auto n = d_source->sgetn(d_buffer.get(), want);
    if (n <= 0) {
      break;  // the offset is past the end
    }
    auto ok = d_scanner.feed(d_buffer.get(), static_cast<std::size_t>(n));
    d_scanned += ok;
    if (ok < static_cast<std::size_t>(n)) {
      refuse(d_scanner.error(), d_scanned);
      return false;
    }
  }
  return d_source->pubseekpos(static_cast<off>(offset), std::ios_base::in) !=
         pos_type(off(-1));
}

CheckedStreambuf::int_type CheckedStreambuf::underflow() {
  if (gptr() < egptr()) {
    return traits_type::to_int_type(*gptr());
  }
  std::uint64_t start =
      d_offset + static_cast<std::uint64_t>(egptr() - eback());
  d_offset = start;
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
  if (start >= d_limit || (start > d_scanned && !scan_to(start))) {
    return traits_type::eof();
  }
  auto n = d_source->sgetn(d_buffer.get(), kChunk);
  if (n <= 0) {
    return traits_type::eof();
  }
  auto deliver = static_cast<std::uint64_t>(n);
  if (start + deliver > d_scanned) {
    // Bytes not scanned before: scan them in order.
    auto seen = d_scanned - start;
    auto fresh = static_cast<std::size_t>(deliver - seen);
    auto ok = d_scanner.feed(d_buffer.get() + seen, fresh);
    d_scanned += ok;
    if (ok < fresh) {
      refuse(d_scanner.error(), d_scanned);
    }
  }
  deliver = std::min(deliver, d_limit - start);
  setg(d_buffer.get(), d_buffer.get(), d_buffer.get() + deliver);
  if (deliver == 0) {
    return traits_type::eof();
  }
  return traits_type::to_int_type(*gptr());
}

CheckedStreambuf::pos_type CheckedStreambuf::moved(pos_type pos) {
  if (pos != pos_type(off_type(-1))) {
    d_offset = static_cast<std::uint64_t>(off_type(pos));
    setg(d_buffer.get(), d_buffer.get(), d_buffer.get());
  }
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
