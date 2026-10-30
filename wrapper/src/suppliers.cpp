// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/suppliers.h"

#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/suppliers.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/threads.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace FP = RDKit::v2::FileParsers;

// A streambuf that reads a Rust `Read` value (Python's file objects), once,
// without seeking. A failed read reads as the end of the input and is
// reported by SupplierInput::check.
class RustReadBuf : public std::streambuf {
 public:
  explicit RustReadBuf(rust::Box<RustReader> reader)
      : d_reader(std::move(reader)), d_buffer(new char[kSize]) {}
  bool failed() const { return !d_error.empty(); }
  const std::string &error() const { return d_error; }

 protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    if (failed()) {
      return traits_type::eof();
    }
    std::size_t n = 0;
    try {
      n = rust_reader_read(
          *d_reader,
          rust::Slice<std::uint8_t>(
              reinterpret_cast<std::uint8_t *>(d_buffer.get()), kSize));
    } catch (const std::exception &e) {
      d_error = e.what();
      return traits_type::eof();
    }
    if (n == 0) {
      return traits_type::eof();
    }
    // rust_reader_read refuses a count past the buffer; this keeps the
    // stream inside it regardless.
    if (n > kSize) {
      d_error = "the reader returned more bytes than its buffer holds";
      return traits_type::eof();
    }
    setg(d_buffer.get(), d_buffer.get(), d_buffer.get() + n);
    return traits_type::to_int_type(*gptr());
  }

 private:
  static constexpr std::size_t kSize = std::size_t{1} << 16;
  rust::Box<RustReader> d_reader;
  std::unique_ptr<char[]> d_buffer;
  std::string d_error;
};

void SupplierInput::check() const {
  if (reader && reader->failed()) {
    throw RDKit::BadFileException("could not read the input: " +
                                  reader->error());
  }
  if (stream) {
    stream->check();
  }
}

namespace {

// Calls into RDKit, then checks the input: refused input explains whatever
// RDKit did with it, including an exception. Once the input has failed, RDKit
// is not called again: its suppliers trust the positions they recorded, and
// a failed read may have recorded tellg's -1.
template <class F>
auto checked(const SupplierInput &input, F &&f) -> decltype(f()) {
  input.check();
  try {
    if constexpr (std::is_void_v<decltype(f())>) {
      f();
      input.check();
    } else {
      auto result = f();
      input.check();
      return result;
    }
  } catch (...) {
    input.check();
    throw;
  }
}

FP::MolFileParserParams mol_params(bool sanitize, bool remove_hs,
                                   bool strict_parsing) {
  FP::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  return params;
}

// Builds `w.suppl` over `in`, replacing any supplier and stream `w` had.
// RDKit's constructors read, so the input is checked after them. The RDKit
// supplier does not own the stream, which outlives it: an owning supplier
// deletes its stream when destroyed, including by that check or, through its
// base classes' destructors, by its own constructor throwing.
template <class W, class S, class... Args>
void build(W &w, std::unique_ptr<CheckedIStream> in, const RustReadBuf *reader,
           Args &&...args) {
  SupplierInput input{std::move(in), reader};
  auto suppl = checked(input, [&] {
    return std::make_unique<S>(input.stream.get(), false,
                               std::forward<Args>(args)...);
  });
  // The old supplier goes before the stream it reads.
  w.suppl = std::move(suppl);
  w.input = std::move(input);
}

std::unique_ptr<CheckedIStream> text_of(rust::Slice<const std::uint8_t> data) {
  return checked_text(to_std(data));
}

std::unique_ptr<CheckedIStream> text_of(rust::Str data) {
  return checked_text(std::string(data));
}

// The text of record idx, after checking that there is one (E3).
template <class W>
std::unique_ptr<std::string> item_text(W &w, std::uint32_t idx) {
  return checked(w.input, [&] {
    auto len = w.suppl->length();
    if (idx >= len) {
      throw IndexError(idx, len);
    }
    return std::make_unique<std::string>(w.suppl->getItemText(idx));
  });
}

// Record idx, as Python's __getitem__: past the last record is an
// IndexError.
template <class W>
std::unique_ptr<RWMol> get_record(W &w, std::uint32_t idx) {
  return checked(w.input, [&]() -> std::unique_ptr<RWMol> {
    try {
      return (*w.suppl)[idx];
    } catch (const RDKit::FileParseException &) {
      auto len = w.suppl->length();
      if (idx >= len) {
        throw IndexError(idx, len);
      }
      throw;
    }
  });
}

}  // namespace

// SDMolSupplier

std::unique_ptr<SDMolSupplier> sd_mol_supplier_new() {
  auto s = std::make_unique<SDMolSupplier>();
  build<SDMolSupplier, FP::SDMolSupplier>(*s, checked_text(""), nullptr,
                                          FP::MolFileParserParams());
  return s;
}

std::unique_ptr<SDMolSupplier> sd_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing) {
  auto s = std::make_unique<SDMolSupplier>();
  build<SDMolSupplier, FP::SDMolSupplier>(
      *s, open_checked_file(to_path(path), true), nullptr,
      mol_params(sanitize, remove_hs, strict_parsing));
  return s;
}

// As SDMolSupplier::setData, which restarts the supplier on a new stream.
void sd_mol_supplier_set_data(SDMolSupplier &s,
                              rust::Slice<const std::uint8_t> data,
                              bool sanitize, bool remove_hs,
                              bool strict_parsing) {
  build<SDMolSupplier, FP::SDMolSupplier>(
      s, text_of(data), nullptr,
      mol_params(sanitize, remove_hs, strict_parsing));
}

std::unique_ptr<RWMol> sd_mol_supplier_next(SDMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->next(); });
}

// RDKit's SDMolSupplier counts the lines it reads in int (d_line) and only
// reset() restarts the count, so each random access resets first; reset()
// leaves nothing else that operator[] does not set (U12).
std::unique_ptr<RWMol> sd_mol_supplier_get(SDMolSupplier &s,
                                           std::uint32_t idx) {
  s.suppl->reset();
  return get_record(s, idx);
}

std::unique_ptr<std::string> sd_mol_supplier_get_item_text(SDMolSupplier &s,
                                                           std::uint32_t idx) {
  return item_text(s, idx);
}

std::uint32_t sd_mol_supplier_length(SDMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->length(); });
}

// atEnd reads ahead, so a read that fails there is reported here, not taken
// for the end of the input.
bool sd_mol_supplier_at_end(SDMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->atEnd(); });
}

void sd_mol_supplier_reset(SDMolSupplier &s) { s.suppl->reset(); }

bool sd_mol_supplier_get_process_property_lists(const SDMolSupplier &s) {
  return s.suppl->getProcessPropertyLists();
}

void sd_mol_supplier_set_process_property_lists(SDMolSupplier &s, bool val) {
  s.suppl->setProcessPropertyLists(val);
}

// ForwardSDMolSupplier

std::unique_ptr<ForwardSDMolSupplier> forward_sd_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing) {
  auto s = std::make_unique<ForwardSDMolSupplier>();
  build<ForwardSDMolSupplier, FP::ForwardSDMolSupplier>(
      *s, open_checked_file(to_path(path), false), nullptr,
      mol_params(sanitize, remove_hs, strict_parsing));
  return s;
}

std::unique_ptr<ForwardSDMolSupplier> forward_sd_mol_supplier_from_reader(
    rust::Box<RustReader> reader, bool sanitize, bool remove_hs,
    bool strict_parsing) {
  auto buf = std::make_unique<RustReadBuf>(std::move(reader));
  const RustReadBuf *raw = buf.get();
  auto s = std::make_unique<ForwardSDMolSupplier>();
  build<ForwardSDMolSupplier, FP::ForwardSDMolSupplier>(
      *s, std::make_unique<CheckedIStream>(std::move(buf)), raw,
      mol_params(sanitize, remove_hs, strict_parsing));
  return s;
}

std::unique_ptr<RWMol> forward_sd_mol_supplier_next(ForwardSDMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->next(); });
}

bool forward_sd_mol_supplier_at_end(ForwardSDMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->atEnd(); });
}

bool forward_sd_mol_supplier_get_eof_hit_on_read(
    const ForwardSDMolSupplier &s) {
  return s.suppl->getEOFHitOnRead();
}

bool forward_sd_mol_supplier_get_process_property_lists(
    const ForwardSDMolSupplier &s) {
  return s.suppl->getProcessPropertyLists();
}

void forward_sd_mol_supplier_set_process_property_lists(ForwardSDMolSupplier &s,
                                                        bool val) {
  s.suppl->setProcessPropertyLists(val);
}

// SmilesMolSupplier

namespace {

FP::SmilesMolSupplierParams smiles_params(
    const SmilesMolSupplierParamsData &data) {
  // RDKit's SMILES suppliers index each line's fields by the name column
  // whenever it is not -1 and is less than their number
  // (SmilesMolSupplier::processLine and
  // MultithreadedSmilesMolSupplier::processMoleculeRecord), reading before
  // the fields for a column below -1. A smiles column below 0 is refused by
  // RDKit itself.
  if (data.name_column < -1) {
    throw ValueErrorException("nameColumn must be -1 or more");
  }
  FP::SmilesMolSupplierParams params;
  params.delimiter = std::string(data.delimiter);
  params.smilesColumn = data.smiles_column;
  params.nameColumn = data.name_column;
  params.titleLine = data.title_line;
  params.parseParameters.sanitize = data.sanitize;
  return params;
}

}  // namespace

std::unique_ptr<SmilesMolSupplier> smiles_mol_supplier_new() {
  auto s = std::make_unique<SmilesMolSupplier>();
  build<SmilesMolSupplier, FP::SmilesMolSupplier>(
      *s, checked_text(""), nullptr, FP::SmilesMolSupplierParams());
  return s;
}

std::unique_ptr<SmilesMolSupplier> smiles_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path,
    const SmilesMolSupplierParamsData &params) {
  auto p = smiles_params(params);
  auto s = std::make_unique<SmilesMolSupplier>();
  build<SmilesMolSupplier, FP::SmilesMolSupplier>(
      *s, open_checked_file(to_path(path), true), nullptr, p);
  return s;
}

void smiles_mol_supplier_set_data(SmilesMolSupplier &s, rust::Str data,
                                  const SmilesMolSupplierParamsData &params) {
  auto p = smiles_params(params);
  build<SmilesMolSupplier, FP::SmilesMolSupplier>(s, text_of(data), nullptr, p);
}

std::unique_ptr<RWMol> smiles_mol_supplier_next(SmilesMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->next(); });
}

std::unique_ptr<RWMol> smiles_mol_supplier_get(SmilesMolSupplier &s,
                                               std::uint32_t idx) {
  return get_record(s, idx);
}

std::unique_ptr<std::string> smiles_mol_supplier_get_item_text(
    SmilesMolSupplier &s, std::uint32_t idx) {
  return item_text(s, idx);
}

std::uint32_t smiles_mol_supplier_length(SmilesMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->length(); });
}

bool smiles_mol_supplier_at_end(SmilesMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->atEnd(); });
}

void smiles_mol_supplier_reset(SmilesMolSupplier &s) { s.suppl->reset(); }

// TDTMolSupplier
//
// RDKit's TDTMolSupplier trusts the record positions it collects:
// getItemText reads past their end (d_molpos[idx + 1]) once a "|" line ends
// the input without a newline (Python's GetItemText then returns garbage or
// raises MemoryError), length() loops forever on input without records (as
// Python's len does), and it counts lines in int across calls, which
// getItemText advances on every call that reads to the end. So the wrapper
// indexes the records itself, at each line that starts with "$SMI<" (as
// length() counts them), and reads each with a fresh RDKit supplier over the
// text, placed at the record's start: the supplier reads it, and decides
// whether the input ends after it, as RDKit's own does once it has counted the
// records.

namespace {

FP::TDTMolSupplierParams tdt_params(const TDTMolSupplierParamsData &data) {
  FP::TDTMolSupplierParams params;
  params.nameRecord = std::string(data.name_record);
  params.confId2D = data.conf_id_2d;
  params.confId3D = data.conf_id_3d;
  params.parseParameters.sanitize = data.sanitize;
  return params;
}

// A read-only, seekable streambuf over a string, which it does not copy. It
// has no put area and keeps the default pbackfail, so nothing writes the
// characters.
class TextBuf : public std::streambuf {
 public:
  explicit TextBuf(const std::string &text) {
    char *begin = const_cast<char *>(text.data());
    setg(begin, begin, begin + text.size());
  }

 protected:
  pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                   std::ios_base::openmode which) override {
    const off_type size = egptr() - eback();
    off_type base = size;
    if (dir == std::ios_base::beg) {
      base = 0;
    } else if (dir == std::ios_base::cur) {
      base = gptr() - eback();
    }
    if (!(which & std::ios_base::in) || off < -base || off > size - base) {
      return pos_type(off_type(-1));
    }
    setg(eback(), eback() + base + off, egptr());
    return pos_type(base + off);
  }
  pos_type seekpos(pos_type pos, std::ios_base::openmode which) override {
    return seekoff(off_type(pos), std::ios_base::beg, which);
  }
};

void index_tdt(TDTMolSupplier &s, std::string text,
               FP::TDTMolSupplierParams params) {
  // The fresh suppliers count the lines they read in int (U12).
  check_input_text(text);
  std::vector<std::size_t> starts;
  for (std::size_t pos = 0; pos < text.size();) {
    if (text.compare(pos, 5, "$SMI<") == 0) {
      starts.push_back(pos);
    }
    auto newline = text.find('\n', pos);
    pos = newline == std::string::npos ? text.size() : newline + 1;
  }
  s.text = std::move(text);
  s.starts = std::move(starts);
  s.params = std::move(params);
  s.next = 0;
  s.at_end = s.starts.empty();
}

// Reads record idx with a fresh RDKit supplier, then sets where s stands.
std::unique_ptr<RWMol> read_tdt(TDTMolSupplier &s, std::size_t idx) {
  TextBuf buf(s.text);
  std::istream in(&buf);
  in.seekg(static_cast<std::streamoff>(s.starts[idx]));
  FP::TDTMolSupplier parser(&in, false, s.params);
  auto mol = parser.next();
  s.next = idx + 1;
  s.at_end = parser.atEnd();
  return mol;
}

// The bound on idx, as for the other suppliers (E3).
void check_tdt_index(const TDTMolSupplier &s, std::uint32_t idx) {
  if (idx >= s.starts.size()) {
    throw IndexError(idx, static_cast<unsigned>(s.starts.size()));
  }
}

}  // namespace

std::unique_ptr<TDTMolSupplier> tdt_mol_supplier_new() {
  return std::make_unique<TDTMolSupplier>();
}

std::unique_ptr<TDTMolSupplier> tdt_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path,
    const TDTMolSupplierParamsData &params) {
  auto p = tdt_params(params);
  auto s = std::make_unique<TDTMolSupplier>();
  index_tdt(*s, read_supplier_file(to_path(path)), std::move(p));
  return s;
}

void tdt_mol_supplier_set_data(TDTMolSupplier &s, rust::Str data,
                               const TDTMolSupplierParamsData &params) {
  index_tdt(s, std::string(data), tdt_params(params));
}

// As TDTMolSupplier::next. Past the last record (after a reset of input
// without records, or before lines that follow the last record), RDKit reads
// one more record without a molecule and ends.
std::unique_ptr<RWMol> tdt_mol_supplier_next(TDTMolSupplier &s) {
  if (s.at_end) {
    throw RDKit::FileParseException("EOF hit.");
  }
  if (s.next >= s.starts.size()) {
    ++s.next;
    s.at_end = true;
    return nullptr;
  }
  return read_tdt(s, s.next);
}

// As TDTMolSupplier::operator[]: reading continues after the record, which
// an exception leaves next.
std::unique_ptr<RWMol> tdt_mol_supplier_get(TDTMolSupplier &s,
                                            std::uint32_t idx) {
  check_tdt_index(s, idx);
  s.next = idx;
  return read_tdt(s, idx);
}

// As TDTMolSupplier::getItemText once it has counted the records: up to the
// next record's start, or, for the last, through its "|" line.
std::unique_ptr<std::string> tdt_mol_supplier_get_item_text(TDTMolSupplier &s,
                                                            std::uint32_t idx) {
  check_tdt_index(s, idx);
  const std::string &text = s.text;
  const std::size_t start = s.starts[idx];
  std::size_t end = text.size();
  if (idx + std::size_t{1} < s.starts.size()) {
    end = s.starts[idx + 1];
  } else {
    for (auto pos = text.find('\n', start); pos != std::string::npos;
         pos = text.find('\n', pos + 1)) {
      if (text.compare(pos + 1, 1, "|") == 0) {
        auto newline = text.find('\n', pos + 1);
        end = newline == std::string::npos ? text.size() : newline + 1;
        break;
      }
    }
  }
  return std::make_unique<std::string>(text, start, end - start);
}

std::uint32_t tdt_mol_supplier_length(TDTMolSupplier &s) {
  // The input has fewer lines than kMaxInputLines.
  return static_cast<std::uint32_t>(s.starts.size());
}

bool tdt_mol_supplier_at_end(TDTMolSupplier &s) { return s.at_end; }

// As TDTMolSupplier::reset, which always leaves something to read.
void tdt_mol_supplier_reset(TDTMolSupplier &s) {
  s.next = 0;
  s.at_end = false;
}

// The multithreaded suppliers

namespace {

FP::MultithreadedMolSupplier::Parameters mt_params(
    const MultithreadedParamsData &data) {
  // RDKit converts the thread count to int and negates it
  // (getNumThreadsToUse in RDThreads.h), which overflows for 2^31; 0 asks
  // for one per core.
  // A count the system cannot start fails at the first next() (mt_next).
  if (data.num_writer_threads > kMaxSupplierThreads) {
    throw ValueErrorException("numWriterThreads must be at most " +
                              std::to_string(kMaxSupplierThreads));
  }
  // RDKit's ConcurrentQueue blocks forever when its capacity is 0, so the
  // supplier's threads would never finish; kMaxSupplierQueueSize says why
  // large ones are refused.
  for (std::size_t size : {data.size_input_queue, data.size_output_queue}) {
    if (size == 0 || size > kMaxSupplierQueueSize) {
      throw ValueErrorException("queue sizes must be at least 1 and at most " +
                                std::to_string(kMaxSupplierQueueSize));
    }
  }
  FP::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = data.num_writer_threads;
  params.sizeInputQueue = data.size_input_queue;
  params.sizeOutputQueue = data.size_output_queue;
  return params;
}

// next() returns null both at the end and for a record it could not read;
// only a read record changes the last record id.
//
// The first next() starts the reader thread and every writer thread
// (MultithreadedMolSupplier::startThreads). When the system refuses one,
// std::thread throws std::system_error, leaving the threads already started
// running; the output queue is marked done only when every writer the
// parameters name has finished, so a later next() would wait forever. The
// supplier is refused from then on, and its destructor joins the threads
// that did start.
template <class W>
std::unique_ptr<RWMol> mt_next(W &w, std::uint32_t &record) {
  if (w.broken) {
    throw std::runtime_error("the supplier's threads could not all be started");
  }
  auto before = w.suppl->getLastRecordId();
  std::unique_ptr<RWMol> mol;
  try {
    mol = checked(w.input, [&] { return w.suppl->next(); });
  } catch (const std::system_error &) {
    w.broken = true;
    throw;
  }
  auto after = w.suppl->getLastRecordId();
  record = (mol || after != before) ? after : 0;
  return mol;
}

}  // namespace

std::unique_ptr<MultithreadedSDMolSupplier> mt_sd_mol_supplier_new() {
  auto s = std::make_unique<MultithreadedSDMolSupplier>();
  build<MultithreadedSDMolSupplier, FP::MultithreadedSDMolSupplier>(
      *s, checked_text(""), nullptr, FP::MultithreadedMolSupplier::Parameters(),
      FP::MolFileParserParams());
  return s;
}

std::unique_ptr<MultithreadedSDMolSupplier> mt_sd_mol_supplier_from_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing, const MultithreadedParamsData &params) {
  auto p = mt_params(params);
  auto s = std::make_unique<MultithreadedSDMolSupplier>();
  build<MultithreadedSDMolSupplier, FP::MultithreadedSDMolSupplier>(
      *s, open_checked_file(to_path(path), true), nullptr, p,
      mol_params(sanitize, remove_hs, strict_parsing));
  return s;
}

std::unique_ptr<RWMol> mt_sd_mol_supplier_next(MultithreadedSDMolSupplier &s,
                                               std::uint32_t &record) {
  s.started = true;
  return mt_next(s, record);
}

bool mt_sd_mol_supplier_at_end(MultithreadedSDMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->atEnd(); });
}

std::uint32_t mt_sd_mol_supplier_get_last_record_id(
    const MultithreadedSDMolSupplier &s) {
  return s.suppl->getLastRecordId();
}

std::unique_ptr<std::string> mt_sd_mol_supplier_get_last_item_text(
    const MultithreadedSDMolSupplier &s) {
  return std::make_unique<std::string>(s.suppl->getLastItemText());
}

bool mt_sd_mol_supplier_get_process_property_lists(
    const MultithreadedSDMolSupplier &s) {
  return s.suppl->getProcessPropertyLists();
}

// RDKit's parsing threads read the flag without a lock
// (MultithreadedSDMolSupplier::readMolProps), so it may not change once they
// have started, at the first next().
void mt_sd_mol_supplier_set_process_property_lists(
    MultithreadedSDMolSupplier &s, bool val) {
  if (s.started) {
    throw ValueErrorException(
        "the property-list setting cannot change once reading has started");
  }
  s.suppl->setProcessPropertyLists(val);
}

std::unique_ptr<MultithreadedSmilesMolSupplier> mt_smiles_mol_supplier_new() {
  auto s = std::make_unique<MultithreadedSmilesMolSupplier>();
  build<MultithreadedSmilesMolSupplier, FP::MultithreadedSmilesMolSupplier>(
      *s, checked_text(""), nullptr, FP::MultithreadedMolSupplier::Parameters(),
      FP::SmilesMolSupplierParams());
  return s;
}

std::unique_ptr<MultithreadedSmilesMolSupplier>
mt_smiles_mol_supplier_from_file(rust::Slice<const std::uint8_t> path,
                                 const SmilesMolSupplierParamsData &smiles,
                                 const MultithreadedParamsData &params) {
  auto p = mt_params(params);
  auto sp = smiles_params(smiles);
  auto s = std::make_unique<MultithreadedSmilesMolSupplier>();
  build<MultithreadedSmilesMolSupplier, FP::MultithreadedSmilesMolSupplier>(
      *s, open_checked_file(to_path(path), true), nullptr, p, sp);
  return s;
}

std::unique_ptr<RWMol> mt_smiles_mol_supplier_next(
    MultithreadedSmilesMolSupplier &s, std::uint32_t &record) {
  return mt_next(s, record);
}

bool mt_smiles_mol_supplier_at_end(MultithreadedSmilesMolSupplier &s) {
  return checked(s.input, [&] { return s.suppl->atEnd(); });
}

std::uint32_t mt_smiles_mol_supplier_get_last_record_id(
    const MultithreadedSmilesMolSupplier &s) {
  return s.suppl->getLastRecordId();
}

std::unique_ptr<std::string> mt_smiles_mol_supplier_get_last_item_text(
    const MultithreadedSmilesMolSupplier &s) {
  return std::make_unique<std::string>(s.suppl->getLastItemText());
}

}  // namespace rdkit_cxx
