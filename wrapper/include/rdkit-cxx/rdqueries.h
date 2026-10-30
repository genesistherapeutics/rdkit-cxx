/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdqueries: query atoms and bonds built from RDKit's query factories.

#pragma once

#include <cstdint>
#include <memory>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/bond.h"
#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared enums, defined by the bridge (src/ffi/rdqueries.rs).
enum class AtomValueQuery : std::uint8_t;
enum class QueryComparison : std::uint8_t;
enum class AtomFlagQuery : std::uint8_t;

std::unique_ptr<Atom> value_query_atom(AtomValueQuery kind,
                                       QueryComparison comparison,
                                       std::int32_t val, bool negate);
std::unique_ptr<Atom> flag_query_atom(AtomFlagQuery kind, bool negate);

std::unique_ptr<Atom> has_prop_query_atom(rust::Str propname, bool negate);
std::unique_ptr<Atom> has_int_prop_with_value_query_atom(
    rust::Str propname, std::int32_t val, bool negate, std::int32_t tolerance);
std::unique_ptr<Atom> has_bool_prop_with_value_query_atom(rust::Str propname,
                                                          bool val,
                                                          bool negate);
std::unique_ptr<Atom> has_string_prop_with_value_query_atom(rust::Str propname,
                                                            rust::Str val,
                                                            bool negate);
std::unique_ptr<Atom> has_double_prop_with_value_query_atom(rust::Str propname,
                                                            double val,
                                                            bool negate,
                                                            double tolerance);
std::unique_ptr<Atom> has_bit_vect_prop_with_value_query_atom(
    rust::Str propname, const BitVect &val, bool negate, double tolerance);

std::unique_ptr<Bond> has_prop_query_bond(rust::Str propname, bool negate);
std::unique_ptr<Bond> has_int_prop_with_value_query_bond(
    rust::Str propname, std::int32_t val, bool negate, std::int32_t tolerance);
std::unique_ptr<Bond> has_bool_prop_with_value_query_bond(rust::Str propname,
                                                          bool val,
                                                          bool negate);
std::unique_ptr<Bond> has_string_prop_with_value_query_bond(rust::Str propname,
                                                            rust::Str val,
                                                            bool negate);
std::unique_ptr<Bond> has_double_prop_with_value_query_bond(rust::Str propname,
                                                            double val,
                                                            bool negate,
                                                            double tolerance);

void replace_atom_with_query_atom(RWMol &mol, std::uint32_t idx);

}  // namespace rdkit_cxx
