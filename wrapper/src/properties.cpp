// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/properties.h"

#include <atomic>
#include <cstdlib>
#include <mutex>
#include <utility>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/properties.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace Desc = RDKit::Descriptors;

// Properties::registry is a global vector that registerProperty writes
// ("XXX Add mutex?" in Property.cpp) while the constructors, getProperty,
// and getAvailableProperties read it (U4). Every wrapper that touches it
// holds this lock; none calls a functor while holding it.
std::mutex &registry_mutex() {
  // A leaf: held for one RDKit call that takes no other lock, so a closure
  // running under control_c_mutex (U4) can take it without deadlock. A
  // wrapper that takes another lock while holding it must first call
  // check_not_in_control_c_callback(), as rgroup_mutex() does.
  static std::mutex mutex;
  return mutex;
}

// RDKit's own functors include descriptors that match substructures
// (control_c_lock.h), so computing properties holds control_c_mutex
// (HoldControlC). A Rust functor among them runs a closure that may call
// wrappers that take it, so it releases the lock while the closure runs; no
// ControlCHandler is alive between two functors.

// Set at exit, before RDKit's registry is destroyed: closures still in it
// are leaked then rather than dropped, since Rust's thread-locals on the
// main thread are already gone.
std::atomic<bool> exiting{false};

// A property functor whose value is a Rust closure. RDKit's operator() takes
// any ROMol; the closure gets its own copy, as a Rust molecule, without
// cached matrices or descriptors (U6, U8). A panic in the closure becomes
// an exception in Rust's half of the bridge.
class RustPropertyFunctor : public Desc::PropertyFunctor {
 public:
  RustPropertyFunctor(std::string name, std::string version,
                      rust::Box<PropertyCallback> callback)
      : Desc::PropertyFunctor(std::move(name), std::move(version)),
        d_callback(std::move(callback)) {
    static std::once_flag registered;
    std::call_once(registered, [] { std::atexit([] { exiting = true; }); });
  }
  ~RustPropertyFunctor() override {
    if (exiting) {
      (void)d_callback.into_raw();
    }
  }

  double operator()(const RDKit::ROMol &mol) const override {
    auto copy = std::make_unique<RWMol>(mol);
    clear_matrix_caches(*copy);
    clear_descriptor_caches(*copy);
    // A computation's lock is released for the closure and taken back for
    // the functors after it.
    ReleaseControlC release;
    return call_property_callback(*d_callback, std::move(copy));
  }

 private:
  rust::Box<PropertyCallback> d_callback;
};

std::unique_ptr<std::vector<std::string>> strings(
    std::vector<std::string> values) {
  return std::make_unique<std::vector<std::string>>(std::move(values));
}

}  // namespace

std::unique_ptr<PropertyFunctorHandle> property_functor_new(
    rust::Str name, rust::Str version, rust::Box<PropertyCallback> callback) {
  auto handle = std::make_unique<PropertyFunctorHandle>();
  handle->ptr = boost::make_shared<RustPropertyFunctor>(
      std::string(name), std::string(version), std::move(callback));
  return handle;
}

std::unique_ptr<std::string> property_functor_get_name(
    const PropertyFunctorHandle &functor) {
  return std::make_unique<std::string>(functor.ptr->getName());
}

std::unique_ptr<std::string> property_functor_get_version(
    const PropertyFunctorHandle &functor) {
  return std::make_unique<std::string>(functor.ptr->getVersion());
}

// RDKit's own functors call their descriptors with the defaults, which read
// the descriptor caches (U8), and may match substructures; a Rust functor
// calls its closure without the lock.
double property_functor_call(const PropertyFunctorHandle &functor,
                             const RWMol &mol) {
  DescriptorCacheGuard guard(mol);
  std::unique_lock<ControlCMutex> control_c;
  if (!dynamic_cast<const RustPropertyFunctor *>(functor.ptr.get())) {
    control_c = std::unique_lock<ControlCMutex>(control_c_mutex());
  }
  return (*functor.ptr)(mol);
}

std::unique_ptr<Properties> properties_new() {
  std::lock_guard<std::mutex> lock(registry_mutex());
  return std::make_unique<Properties>();
}

std::unique_ptr<Properties> properties_from_prop_names(
    rust::Slice<const rust::Str> prop_names) {
  std::vector<std::string> names;
  for (auto name : prop_names) {
    names.emplace_back(name);
  }
  std::lock_guard<std::mutex> lock(registry_mutex());
  return std::make_unique<Properties>(names);
}

std::unique_ptr<std::vector<std::string>> properties_get_property_names(
    const Properties &properties) {
  return strings(properties.getPropertyNames());
}

rust::Vec<double> properties_compute_properties(const Properties &properties,
                                                const RWMol &mol) {
  DescriptorCacheGuard guard(mol);
  HoldControlC control_c;
  return to_vec(properties.computeProperties(mol, false));
}

rust::Vec<double> properties_compute_properties_annotate(
    const Properties &properties, RWMol &mol) {
  DescriptorCacheGuard guard(mol);
  HoldControlC control_c;
  return to_vec(properties.computeProperties(mol, true));
}

void properties_annotate_properties(const Properties &properties, RWMol &mol) {
  DescriptorCacheGuard guard(mol);
  HoldControlC control_c;
  properties.annotateProperties(mol);
}

std::unique_ptr<std::vector<std::string>>
properties_get_available_properties() {
  std::lock_guard<std::mutex> lock(registry_mutex());
  return strings(Properties::getAvailableProperties());
}

std::unique_ptr<PropertyFunctorHandle> properties_get_property(rust::Str name) {
  auto handle = std::make_unique<PropertyFunctorHandle>();
  std::lock_guard<std::mutex> lock(registry_mutex());
  handle->ptr = Properties::getProperty(std::string(name));
  return handle;
}

// registerProperty replaces a functor of the same name, which could drop
// its Rust closure under the lock; a reference to it is held until after.
std::int32_t properties_register_property(
    const PropertyFunctorHandle &functor) {
  boost::shared_ptr<Desc::PropertyFunctor> replaced;
  std::lock_guard<std::mutex> lock(registry_mutex());
  // registerDescriptors fills the registry on first use; registering first
  // would let it replace this functor later.
  Properties::getAvailableProperties();
  for (const auto &registered : Properties::registry) {
    if (registered->getName() == functor.ptr->getName()) {
      replaced = registered;
    }
  }
  auto idx = Properties::registerProperty(functor.ptr);
  return idx;
}

std::unique_ptr<PropertyRangeQuery> make_property_range_query(rust::Str name,
                                                              double min,
                                                              double max) {
  std::lock_guard<std::mutex> lock(registry_mutex());
  return std::unique_ptr<PropertyRangeQuery>(
      Desc::makePropertyRangeQuery(std::string(name), min, max));
}

// The query calls the descriptor function of RDKit's functor, which reads
// the descriptor caches (U8) and may match substructures. A Rust functor has
// none, which RDKit's precondition reports.
bool property_range_query_match(const PropertyRangeQuery &query,
                                const RWMol &mol) {
  DescriptorCacheGuard guard(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return query.Match(mol);
}

}  // namespace rdkit_cxx
