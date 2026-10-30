/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDescriptors's property registry: Properties, property functors, and
// property range queries. RDKit's registry is a global vector it writes
// without a lock, so every wrapper that reads or writes it holds one (U4).

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/Descriptors/Property.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A Rust closure, defined by the bridge (src/ffi/properties.rs).
struct PropertyCallback;

// A property functor shared with the registry and with Properties objects:
// one of RDKit's, or one that calls a Rust closure.
struct PropertyFunctorHandle {
  boost::shared_ptr<RDKit::Descriptors::PropertyFunctor> ptr;
};

using Properties = RDKit::Descriptors::Properties;
using PropertyRangeQuery = RDKit::Descriptors::PROP_RANGE_QUERY;

std::unique_ptr<PropertyFunctorHandle> property_functor_new(
    rust::Str name, rust::Str version, rust::Box<PropertyCallback> callback);
std::unique_ptr<std::string> property_functor_get_name(
    const PropertyFunctorHandle &functor);
std::unique_ptr<std::string> property_functor_get_version(
    const PropertyFunctorHandle &functor);
double property_functor_call(const PropertyFunctorHandle &functor,
                             const RWMol &mol);

std::unique_ptr<Properties> properties_new();
std::unique_ptr<Properties> properties_from_prop_names(
    rust::Slice<const rust::Str> prop_names);
std::unique_ptr<std::vector<std::string>> properties_get_property_names(
    const Properties &properties);
rust::Vec<double> properties_compute_properties(const Properties &properties,
                                                const RWMol &mol);
rust::Vec<double> properties_compute_properties_annotate(
    const Properties &properties, RWMol &mol);
void properties_annotate_properties(const Properties &properties, RWMol &mol);
std::unique_ptr<std::vector<std::string>> properties_get_available_properties();
std::unique_ptr<PropertyFunctorHandle> properties_get_property(rust::Str name);
std::int32_t properties_register_property(const PropertyFunctorHandle &functor);

std::unique_ptr<PropertyRangeQuery> make_property_range_query(rust::Str name,
                                                              double min,
                                                              double max);
bool property_range_query_match(const PropertyRangeQuery &query,
                                const RWMol &mol);

}  // namespace rdkit_cxx
