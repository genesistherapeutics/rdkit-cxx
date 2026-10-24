/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Properties on molecules, atoms, bonds, conformers, and substance groups.
// They all derive from RDKit::RDProps, so one set of wrappers serves every
// class; each class has an upcast (`<class>_as_props`). RDKit's property
// setters are const, so these take const RDProps&; the Rust API still takes
// &mut self to set.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <RDGeneral/RDProps.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/bond.h"
#include "rdkit-cxx/conformer.h"
#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/substance_group.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using RDProps = RDKit::RDProps;

// A shared struct, defined by the bridge (src/ffi/props.rs).
struct PropValueData;

const RDProps &mol_as_props(const RWMol &mol);
const RDProps &atom_as_props(const Atom &atom);
const RDProps &bond_as_props(const Bond &bond);
const RDProps &conformer_as_props(const Conformer &conf);
const RDProps &substance_group_as_props(const SubstanceGroup &sgroup);

bool props_has_prop(const RDProps &obj, rust::Str key);
// As GetProp, GetPyProp, and GetPropsAsDict in
// Code/GraphMol/nbWrap/props.hpp, with nanobind's messages.
std::unique_ptr<std::string> props_get_prop(const RDProps &obj, rust::Str key);
PropValueData props_get_prop_auto_convert(const RDProps &obj, rust::Str key);
int props_get_int_prop(const RDProps &obj, rust::Str key);
unsigned props_get_unsigned_prop(const RDProps &obj, rust::Str key);
double props_get_double_prop(const RDProps &obj, rust::Str key);
bool props_get_bool_prop(const RDProps &obj, rust::Str key);
std::unique_ptr<BitVect> props_get_explicit_bit_vect_prop(const RDProps &obj,
                                                          rust::Str key);
rust::Vec<std::uint32_t> props_get_unsigned_vect_prop(const RDProps &obj,
                                                      rust::Str key);
std::unique_ptr<std::vector<std::string>> props_get_string_vect_prop(
    const RDProps &obj, rust::Str key);
std::unique_ptr<std::vector<std::string>> props_get_prop_names(
    const RDProps &obj, bool include_private, bool include_computed);
rust::Vec<PropValueData> props_get_props_as_dict(const RDProps &obj,
                                                 bool include_private,
                                                 bool include_computed,
                                                 bool auto_convert_strings);

void props_set_prop(const RDProps &obj, rust::Str key, rust::Str val,
                    bool computed);
void props_set_int_prop(const RDProps &obj, rust::Str key, int val,
                        bool computed);
void props_set_unsigned_prop(const RDProps &obj, rust::Str key, unsigned val,
                             bool computed);
void props_set_double_prop(const RDProps &obj, rust::Str key, double val,
                           bool computed);
void props_set_bool_prop(const RDProps &obj, rust::Str key, bool val,
                         bool computed);
void props_set_explicit_bit_vect_prop(const RDProps &obj, rust::Str key,
                                      const BitVect &val);
void props_clear_prop(const RDProps &obj, rust::Str key);
void props_clear_computed_props(const RDProps &obj);
void mol_clear_computed_props(const RWMol &mol, bool include_rings);

}  // namespace rdkit_cxx
