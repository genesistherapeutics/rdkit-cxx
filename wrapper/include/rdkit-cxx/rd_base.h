/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <memory>
#include <string>

#include <RDGeneral/RDLog.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

using BlockLogs = RDLog::LogStateSetter;
using CaptureErrorLog = RDLog::CaptureErrorLog;

std::unique_ptr<std::string> rdkit_version();
std::unique_ptr<std::string> boost_version();
std::unique_ptr<std::string> rdkit_build();
void disable_log(rust::Str spec);
void enable_log(rust::Str spec);
std::unique_ptr<std::string> log_status();
void log_to_cpp_streams();
void log_debug_msg(rust::Str msg);
void log_info_msg(rust::Str msg);
void log_warning_msg(rust::Str msg);
void log_error_msg(rust::Str msg);
void log_message(rust::Str spec, rust::Str msg);
void seed_random_number_generator(unsigned seed);

std::unique_ptr<BlockLogs> block_logs_new();
std::unique_ptr<CaptureErrorLog> capture_error_log_new();
std::unique_ptr<std::string> capture_error_log_messages(
    const CaptureErrorLog &capture);

}  // namespace rdkit_cxx
