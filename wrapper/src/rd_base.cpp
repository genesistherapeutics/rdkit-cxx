// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/rd_base.h"

#include <cstdlib>

#include <RDGeneral/versions.h>

namespace rdkit_cxx {

std::unique_ptr<std::string> rdkit_version() {
  return std::make_unique<std::string>(RDKit::rdkitVersion);
}

std::unique_ptr<std::string> boost_version() {
  return std::make_unique<std::string>(RDKit::boostVersion);
}

std::unique_ptr<std::string> rdkit_build() {
  return std::make_unique<std::string>(RDKit::rdkitBuild);
}

// As DisableLog, EnableLog, LogStatus, and the Log*Msg functions in
// Code/RDBoost/Wrap/RDBase.cpp.
void disable_log(rust::Str spec) {
  boost::logging::disable_logs(std::string(spec));
}

void enable_log(rust::Str spec) {
  boost::logging::enable_logs(std::string(spec));
}

std::unique_ptr<std::string> log_status() {
  return std::make_unique<std::string>(boost::logging::log_status());
}

void log_to_cpp_streams() { RDLog::InitLogs(); }

void log_debug_msg(rust::Str msg) {
  BOOST_LOG(rdDebugLog) << std::string(msg) << std::endl;
}

void log_info_msg(rust::Str msg) {
  BOOST_LOG(rdInfoLog) << std::string(msg) << std::endl;
}

void log_warning_msg(rust::Str msg) {
  BOOST_LOG(rdWarningLog) << std::string(msg) << std::endl;
}

void log_error_msg(rust::Str msg) {
  BOOST_LOG(rdErrorLog) << std::string(msg) << std::endl;
}

void log_message(rust::Str spec_str, rust::Str msg) {
  std::string spec(spec_str);
  if (spec == "rdApp.error") {
    log_error_msg(msg);
  } else if (spec == "rdApp.warning") {
    log_warning_msg(msg);
  } else if (spec == "rdApp.info") {
    log_info_msg(msg);
  } else if (spec == "rdApp.debug") {
    log_debug_msg(msg);
  }
}

void seed_random_number_generator(unsigned seed) { std::srand(seed); }

std::unique_ptr<BlockLogs> block_logs_new() {
  return std::make_unique<BlockLogs>();
}

std::unique_ptr<CaptureErrorLog> capture_error_log_new() {
  return std::make_unique<CaptureErrorLog>();
}

std::unique_ptr<std::string> capture_error_log_messages(
    const CaptureErrorLog &capture) {
  return std::make_unique<std::string>(capture.messages());
}

}  // namespace rdkit_cxx
