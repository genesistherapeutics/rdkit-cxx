/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of the pickles of RDKit's hierarchical catalogs (MolCatalog and
// FragCatalog, RDCatalog::HierarchCatalog in Catalogs/Catalog.h) and of their
// entries, before RDKit reads them (CONVENTIONS.md U2).
//
// RDKit reads a catalog pickle as a header, a fingerprint length, an entry
// count, the catalog's parameters, each entry, and each entry's down edges.
// Each entry is a molecule pickle followed by fields; MolCatalogEntry's and
// FragCatalogEntry's readers allocate a description's int32 length plus one
// bytes and write a NUL at that length, so a negative length writes before
// the buffer, INT_MAX overflows the sum, and a length past the end of the
// pickle leaves bytes unread that the description then holds. The molecule
// pickles are read by MolPickler, which RDKit fuzzes, but an atom's monomer
// info can have a type its class lacks (U10).
//
// The checks read the pickle exactly as RDKit will, with RDKit's own readers
// for the molecules and parameters, so that RDKit then reads the same bytes,
// in memory, the same way.

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <sstream>
#include <string>

#include <RDGeneral/Exceptions.h>
#include <RDGeneral/StreamOps.h>

namespace rdkit_cxx {

// A stream over a pickle's bytes, as the catalogs' initFromString makes.
std::stringstream pickle_stream(const std::string &text);

// Reads an int32 as RDKit's streamRead does.
std::int32_t read_int32(std::istream &ss);

// Reads a molecule pickle as the entries' readers do, and checks its atoms'
// monomer infos (U10).
void check_entry_mol(std::istream &ss);

// Reads an entry's description as RDKit's readers do: an int32 length and
// that many bytes. Throws ValueErrorException for a length below zero, of
// INT_MAX, or past the end of the `size` bytes of the pickle.
void check_entry_text(std::istream &ss, std::size_t size);

// Checks a catalog pickle as HierarchCatalog::initFromStream reads it, with
// `Params` reading the parameters and `check_entry(ss, size)` reading each
// entry. RDKit checks the down edges' entry indices (an invariant that
// leaks the entries already read), so they are checked here first.
template <class Params, class CheckEntry>
void check_catalog_pickle(const std::string &text, CheckEntry check_entry) {
  auto ss = pickle_stream(text);
  for (int i = 0; i < 4; ++i) {
    read_int32(ss);  // the header, which RDKit ignores
  }
  std::uint32_t fp_length;
  RDKit::streamRead(ss, fp_length);
  std::uint32_t num_entries;
  RDKit::streamRead(ss, num_entries);
  Params params;
  params.initFromStream(ss);
  for (std::uint32_t i = 0; i < num_entries; ++i) {
    check_entry(ss, text.size());
  }
  for (std::uint32_t i = 0; i < num_entries; ++i) {
    std::uint32_t num_neighbors;
    RDKit::streamRead(ss, num_neighbors);
    for (std::uint32_t j = 0; j < num_neighbors; ++j) {
      auto id = read_int32(ss);
      if (id < 0 || static_cast<std::uint32_t>(id) >= num_entries) {
        throw ValueErrorException("a catalog pickle's edge names entry " +
                                  std::to_string(id) + " of " +
                                  std::to_string(num_entries));
      }
    }
  }
}

}  // namespace rdkit_cxx
