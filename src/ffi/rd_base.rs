// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/rd_base.h");

        type BlockLogs;

        fn rdkit_version() -> Result<UniquePtr<CxxString>>;
        fn boost_version() -> Result<UniquePtr<CxxString>>;
        fn rdkit_build() -> Result<UniquePtr<CxxString>>;
        fn disable_log(spec: &str) -> Result<()>;
        fn enable_log(spec: &str) -> Result<()>;
        fn log_status() -> Result<UniquePtr<CxxString>>;
        fn log_to_cpp_streams() -> Result<()>;
        fn log_debug_msg(msg: &str) -> Result<()>;
        fn log_info_msg(msg: &str) -> Result<()>;
        fn log_warning_msg(msg: &str) -> Result<()>;
        fn log_error_msg(msg: &str) -> Result<()>;
        fn log_message(spec: &str, msg: &str) -> Result<()>;
        fn seed_random_number_generator(seed: u32) -> Result<()>;

        fn block_logs_new() -> Result<UniquePtr<BlockLogs>>;
    }
}
