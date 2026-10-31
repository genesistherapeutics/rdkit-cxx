/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/mol_bundle.h"

#include <string>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"

namespace rdkit_cxx {

namespace {

boost::shared_ptr<RDKit::ROMol> borrow(const RWMol &mol) {
  // MolBundle holds non-const pointers but only reads through them here.
  return boost::shared_ptr<RDKit::ROMol>(const_cast<RWMol *>(&mol),
                                         [](RDKit::ROMol *) {});
}

}  // namespace

std::unique_ptr<MolBundle> mol_bundle_new() {
  return std::make_unique<MolBundle>();
}

void mol_bundle_add_borrowed(MolBundle &bundle, const RWMol &mol) {
  bundle.addMol(borrow(mol));
}

void fixed_mol_size_mol_bundle_check(const RWMol &first, const RWMol &nmol) {
  RDKit::FixedMolSizeMolBundle bundle;
  bundle.addMol(borrow(first));
  bundle.addMol(borrow(nmol));
}

std::unique_ptr<MolList> mol_bundle_from_pkl_string(
    rust::Slice<const std::uint8_t> pkl) {
  std::string text(reinterpret_cast<const char *>(pkl.data()), pkl.size());
  check_mol_pickle_versions(text);
  MolBundle bundle(text);
  auto list = std::make_unique<MolList>();
  for (const auto &mol : bundle.getMols()) {
    check_pickled_mol(*mol);
    list->push(std::make_unique<RWMol>(*mol));
  }
  return list;
}

rust::Vec<std::uint8_t> mol_bundle_to_binary(const MolBundle &bundle) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  for (const auto &mol : bundle.getMols()) {
    check_picklable(*mol);
  }
  auto res = bundle.serialize();
  rust::Vec<std::uint8_t> bytes;
  bytes.reserve(res.size());
  for (char c : res) {
    bytes.push_back(static_cast<std::uint8_t>(c));
  }
  return bytes;
#else
  (void)bundle;
  throw ValueErrorException("Pickling of MolBundle instances is not enabled");
#endif
}

bool mol_bundle_can_serialize() { return RDKit::MolBundleCanSerialize(); }

}  // namespace rdkit_cxx
