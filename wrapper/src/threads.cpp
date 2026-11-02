// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/threads.h"

#include <cstdint>
#include <exception>
#include <thread>
#include <vector>

namespace rdkit_cxx {

unsigned num_threads_to_use(int target) {
  if (target >= 1) {
    return static_cast<unsigned>(target);
  }
  std::int64_t res =
      static_cast<std::int64_t>(std::thread::hardware_concurrency()) + target;
  return res >= 1 ? static_cast<unsigned>(res) : 1;
}

void run_on_threads(unsigned n, const std::function<void(unsigned)> &work) {
  if (n <= 1) {
    work(0);
    return;
  }
  std::vector<std::exception_ptr> errors(n);
  auto task = [&](unsigned t) {
    try {
      work(t);
    } catch (...) {
      errors[t] = std::current_exception();
    }
  };
  std::vector<std::thread> threads;
  threads.reserve(n - 1);
  try {
    for (unsigned t = 1; t < n; ++t) {
      threads.emplace_back(task, t);
    }
  } catch (...) {
    // A joinable std::thread's destructor calls std::terminate.
    for (auto &thread : threads) {
      thread.join();
    }
    throw;
  }
  task(0);
  for (auto &thread : threads) {
    thread.join();
  }
  for (const auto &error : errors) {
    if (error) {
      std::rethrow_exception(error);
    }
  }
}

}  // namespace rdkit_cxx
