/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// The multithreaded suppliers the wrappers make.

#pragma once

#include <istream>
#include <string>
#include <tuple>
#include <utility>

#include <GraphMol/FileParsers/MultithreadedMolSupplier.h>
#include <GraphMol/RWMol.h>

#include "rdkit-cxx/streams.h"

namespace rdkit_cxx {

// The checked stream under the one a supplier reads.
struct StopOn {
  CheckedIStream *stream;
};

// RDKit 2024.09's MultithreadedMolSupplier joins its threads in its own
// destructor, which runs after the SD or SMILES supplier's: the reader may
// still be in that supplier's extractNextRecord ("pure virtual method
// called"), and the writers may wait to push into a full output queue that
// nothing empties, and the reader into a full input queue (rdkit/rdkit#8556).
// So the suppliers the wrappers make are this subclass, whose destructor,
// run first, ends the input (CheckedIStream::stop) and takes every record
// left from the output queue, which ends once the reader and every writer
// have; their threads then only return. A Rust reader that blocks keeps
// the drop waiting.
template <class S>
class DrainedSupplier : public S {
 public:
  // Reads `in`, which it stops.
  template <class... Args>
  explicit DrainedSupplier(CheckedIStream *in, Args &&...args)
      : S(in, std::forward<Args>(args)...), d_input(in) {}
  // Reads `in`, a stream over `checked`, which it stops.
  template <class... Args>
  DrainedSupplier(StopOn checked, std::istream *in, Args &&...args)
      : S(in, std::forward<Args>(args)...), d_input(checked.stream) {}
  ~DrainedSupplier() override {
    if (!this->df_started) {
      return;
    }
    d_input->stop();
    std::tuple<RDKit::RWMol *, std::string, unsigned int> r;
    while (this->d_outputQueue->pop(r)) {
      delete std::get<0>(r);
    }
  }

 private:
  CheckedIStream *d_input;
};

}  // namespace rdkit_cxx
