// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/props.h"

#include <algorithm>

#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "rdkit-cxx/ffi/props.rs.h"

namespace rdkit_cxx {

namespace {

// The tags of PropValueData, matching src/ffi/props.rs.
enum Tag : std::uint8_t {
  None = 0,
  Int,
  Unsigned,
  Double,
  Bool,
  String,
  IntVec,
  UnsignedVec,
  DoubleVec,
  StringVec,
};

PropValueData empty(const std::string &key) {
  return PropValueData{std::make_unique<std::string>(key),
                       None,
                       0,
                       0,
                       0.0,
                       false,
                       std::make_unique<std::string>(),
                       {},
                       {},
                       {},
                       std::make_unique<std::vector<std::string>>()};
}

template <class T>
const char *type_name();
template <>
const char *type_name<int>() {
  return "an integer value";
}
template <>
const char *type_name<unsigned>() {
  return "an unsigned integer value";
}
template <>
const char *type_name<double>() {
  return "a double value";
}
template <>
const char *type_name<bool>() {
  return "a True or False value";
}
template <>
const char *type_name<ExplicitBitVect>() {
  return "an ExplicitBitVect";
}
template <>
const char *type_name<std::vector<unsigned>>() {
  return "a list of unsigned integers";
}
template <>
const char *type_name<std::vector<std::string>>() {
  return "a list of strings";
}

// As GetProp<RDOb, T> in props.hpp.
template <class T>
T get_typed(const RDProps &obj, rust::Str key_str) {
  std::string key(key_str);
  T res;
  try {
    obj.getProp(key, res);
  } catch (const KeyErrorException &) {
    throw KeyErrorException(key);
  } catch (const std::exception &e) {
    throw ValueErrorException(std::string("key `") + key +
                              "` exists but does not result in " +
                              type_name<T>() + " reason: " + e.what());
  }
  return res;
}

// Converts one property as GetPropsAsDict does, or returns false for types
// Python skips.
bool convert(const RDKit::Dict::Pair &item, bool auto_convert_strings,
             PropValueData &out) {
  using namespace RDKit;
  const auto &val = item.val;
  switch (val.getTag()) {
    case RDTypeTag::IntTag:
      out.tag = Int;
      out.int_value = from_rdvalue<int>(val);
      return true;
    case RDTypeTag::DoubleTag:
      out.tag = Double;
      out.double_value = from_rdvalue<double>(val);
      return true;
    case RDTypeTag::FloatTag:
      out.tag = Double;
      out.double_value = from_rdvalue<float>(val);
      return true;
    case RDTypeTag::StringTag: {
      auto value = from_rdvalue<std::string>(val);
      if (auto_convert_strings) {
        auto trimmed = boost::trim_copy(value);
        int ivalue;
        if (boost::conversion::try_lexical_convert(trimmed, ivalue)) {
          out.tag = Int;
          out.int_value = ivalue;
          return true;
        }
        double dvalue;
        if (boost::conversion::try_lexical_convert(trimmed, dvalue)) {
          out.tag = Double;
          out.double_value = dvalue;
          return true;
        }
      }
      out.tag = String;
      *out.string_value = value;
      return true;
    }
    case RDTypeTag::BoolTag:
      out.tag = Bool;
      out.bool_value = from_rdvalue<bool>(val);
      return true;
    case RDTypeTag::UnsignedIntTag:
      out.tag = Unsigned;
      out.unsigned_value = from_rdvalue<unsigned>(val);
      return true;
    case RDTypeTag::AnyTag:
      return false;
    case RDTypeTag::VecDoubleTag:
      out.tag = DoubleVec;
      for (auto v : from_rdvalue<std::vector<double>>(val)) {
        out.double_vec.push_back(v);
      }
      return true;
    case RDTypeTag::VecFloatTag:
      out.tag = DoubleVec;
      for (auto v : from_rdvalue<std::vector<float>>(val)) {
        out.double_vec.push_back(v);
      }
      return true;
    case RDTypeTag::VecIntTag:
      out.tag = IntVec;
      for (auto v : from_rdvalue<std::vector<int>>(val)) {
        out.int_vec.push_back(v);
      }
      return true;
    case RDTypeTag::VecUnsignedIntTag:
      out.tag = UnsignedVec;
      for (auto v : from_rdvalue<std::vector<unsigned>>(val)) {
        out.unsigned_vec.push_back(v);
      }
      return true;
    case RDTypeTag::VecStringTag:
      out.tag = StringVec;
      *out.string_vec = from_rdvalue<std::vector<std::string>>(val);
      return true;
    case RDTypeTag::EmptyTag:
      out.tag = None;
      return true;
    default:
      throw std::runtime_error(
          "Unhandled property type encountered for property: " + item.key);
  }
}

}  // namespace

const RDProps &mol_as_props(const RWMol &mol) { return mol; }
const RDProps &atom_as_props(const Atom &atom) { return atom; }
const RDProps &bond_as_props(const Bond &bond) { return bond; }
const RDProps &conformer_as_props(const Conformer &conf) { return conf; }

const RDProps &substance_group_as_props(const SubstanceGroup &sgroup) {
  return sgroup;
}

bool props_has_prop(const RDProps &obj, rust::Str key) {
  return obj.hasProp(std::string(key));
}

std::unique_ptr<std::string> props_get_prop(const RDProps &obj,
                                            rust::Str key_str) {
  std::string key(key_str);
  auto res = std::make_unique<std::string>();
  if (!obj.getPropIfPresent(key, *res)) {
    throw KeyErrorException(key);
  }
  return res;
}

// As GetPyPropImpl with autoConvert in props.hpp: strings become integers or
// doubles when they parse as one.
PropValueData props_get_prop_auto_convert(const RDProps &obj,
                                          rust::Str key_str) {
  std::string key(key_str);
  for (const auto &item : obj.getDict()) {
    if (item.key != key) {
      continue;
    }
    auto out = empty(key);
    if (item.val.getTag() == RDKit::RDTypeTag::StringTag) {
      int ivalue;
      double dvalue;
      try {
        if (obj.getPropIfPresent(key, ivalue)) {
          out.tag = Int;
          out.int_value = ivalue;
          return out;
        }
      } catch (const std::bad_any_cast &) {
      }
      try {
        if (obj.getPropIfPresent(key, dvalue)) {
          out.tag = Double;
          out.double_value = dvalue;
          return out;
        }
      } catch (const std::bad_any_cast &) {
      }
      out.tag = String;
      obj.getProp(key, *out.string_value);
      return out;
    }
    if (!convert(item, false, out)) {
      throw std::runtime_error(
          "Unhandled property type encountered for property: " + key);
    }
    return out;
  }
  throw KeyErrorException(key);
}

int props_get_int_prop(const RDProps &obj, rust::Str key) {
  return get_typed<int>(obj, key);
}

unsigned props_get_unsigned_prop(const RDProps &obj, rust::Str key) {
  return get_typed<unsigned>(obj, key);
}

double props_get_double_prop(const RDProps &obj, rust::Str key) {
  return get_typed<double>(obj, key);
}

bool props_get_bool_prop(const RDProps &obj, rust::Str key) {
  return get_typed<bool>(obj, key);
}

std::unique_ptr<BitVect> props_get_explicit_bit_vect_prop(const RDProps &obj,
                                                          rust::Str key) {
  return std::make_unique<BitVect>(get_typed<ExplicitBitVect>(obj, key));
}

rust::Vec<std::uint32_t> props_get_unsigned_vect_prop(const RDProps &obj,
                                                      rust::Str key) {
  rust::Vec<std::uint32_t> result;
  for (auto val : get_typed<std::vector<unsigned>>(obj, key)) {
    result.push_back(val);
  }
  return result;
}

std::unique_ptr<std::vector<std::string>> props_get_string_vect_prop(
    const RDProps &obj, rust::Str key) {
  return std::make_unique<std::vector<std::string>>(
      get_typed<std::vector<std::string>>(obj, key));
}

std::unique_ptr<std::vector<std::string>> props_get_prop_names(
    const RDProps &obj, bool include_private, bool include_computed) {
  return std::make_unique<std::vector<std::string>>(
      obj.getPropList(include_private, include_computed));
}

rust::Vec<PropValueData> props_get_props_as_dict(const RDProps &obj,
                                                 bool include_private,
                                                 bool include_computed,
                                                 bool auto_convert_strings) {
  auto keys = obj.getPropList(include_private, include_computed);
  rust::Vec<PropValueData> result;
  for (const auto &item : obj.getDict()) {
    if (std::find(keys.begin(), keys.end(), item.key) == keys.end()) {
      continue;
    }
    auto out = empty(item.key);
    if (convert(item, auto_convert_strings, out)) {
      result.push_back(std::move(out));
    }
  }
  return result;
}

void props_set_prop(const RDProps &obj, rust::Str key, rust::Str val,
                    bool computed) {
  obj.setProp(std::string(key), std::string(val), computed);
}

void props_set_int_prop(const RDProps &obj, rust::Str key, int val,
                        bool computed) {
  obj.setProp(std::string(key), val, computed);
}

void props_set_unsigned_prop(const RDProps &obj, rust::Str key, unsigned val,
                             bool computed) {
  obj.setProp(std::string(key), val, computed);
}

void props_set_double_prop(const RDProps &obj, rust::Str key, double val,
                           bool computed) {
  obj.setProp(std::string(key), val, computed);
}

void props_set_bool_prop(const RDProps &obj, rust::Str key, bool val,
                         bool computed) {
  obj.setProp(std::string(key), val, computed);
}

// As MolClearProp in props.hpp: a missing key is not an error.
void props_set_explicit_bit_vect_prop(const RDProps &obj, rust::Str key,
                                      const BitVect &val) {
  obj.setProp(std::string(key), as_explicit(val));
}

void props_clear_prop(const RDProps &obj, rust::Str key_str) {
  std::string key(key_str);
  if (obj.hasProp(key)) {
    obj.clearProp(key);
  }
}

void props_clear_computed_props(const RDProps &obj) {
  obj.clearComputedProps();
}

void mol_clear_computed_props(const RWMol &mol, bool include_rings) {
  mol.clearComputedProps(include_rings);
}

}  // namespace rdkit_cxx
