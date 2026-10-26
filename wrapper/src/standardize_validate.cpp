// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_validate.h"

#include <mutex>
#include <utility>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/standardize_validate.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;

namespace {

// A validation method whose validate is a Rust closure. RDKit calls it from
// the thread that validates; copies (MolVSValidation copies its methods)
// share the closure, which is why it is Send + Sync. It gets its own copy
// of the molecule, as a Rust molecule, without cached matrices or
// descriptors (U6, U8). A panic in it becomes an exception in Rust's half of
// the bridge. It runs without the validation's control_c_mutex, so it can
// call wrappers that take it; RDKit calls it between matches.
class RustValidationMethod : public MS::ValidationMethod {
 public:
  explicit RustValidationMethod(
      std::shared_ptr<rust::Box<ValidationCallback>> callback)
      : d_callback(std::move(callback)) {}

  std::vector<MS::ValidationErrorInfo> validate(
      const RDKit::ROMol &mol, bool reportAllFailures) const override {
    auto copy = std::make_unique<RWMol>(mol);
    clear_matrix_caches(*copy);
    clear_descriptor_caches(*copy);
    ReleaseControlC release;
    auto messages = call_validation_callback(**d_callback, std::move(copy),
                                             reportAllFailures);
    std::vector<MS::ValidationErrorInfo> errors;
    for (const auto &message : messages) {
      errors.emplace_back(message);
    }
    return errors;
  }

  std::shared_ptr<MS::ValidationMethod> copy() const override {
    return std::make_shared<RustValidationMethod>(d_callback);
  }

 private:
  std::shared_ptr<rust::Box<ValidationCallback>> d_callback;
};

std::unique_ptr<ValidationMethodHandle> handle(
    std::shared_ptr<MS::ValidationMethod> ptr, bool sums_charges = false) {
  auto result = std::make_unique<ValidationMethodHandle>();
  result->ptr = std::move(ptr);
  result->sums_charges = sums_charges;
  return result;
}

std::vector<std::shared_ptr<RDKit::Atom>> atoms(const AtomList &list,
                                                const char *what) {
  if (list.items.empty()) {
    throw ValueErrorException(std::string(what) +
                              " argument must be non-empty");
  }
  return list.items;
}

std::unique_ptr<std::vector<std::string>> strings(
    std::vector<std::string> values) {
  return std::make_unique<std::vector<std::string>>(std::move(values));
}

}  // namespace

std::unique_ptr<ValidationMethodHandle> validation_rust(
    rust::Box<ValidationCallback> callback) {
  return handle(std::make_shared<RustValidationMethod>(
      std::make_shared<rust::Box<ValidationCallback>>(std::move(callback))));
}

std::unique_ptr<ValidationMethodHandle> validation_rdkit(
    bool allow_empty_molecules) {
  return handle(std::make_shared<MS::RDKitValidation>(allow_empty_molecules));
}

std::unique_ptr<ValidationMethodHandle> validation_no_atom() {
  return handle(std::make_shared<MS::NoAtomValidation>());
}

std::unique_ptr<ValidationMethodHandle> validation_fragment() {
  return handle(std::make_shared<MS::FragmentValidation>());
}

std::unique_ptr<ValidationMethodHandle> validation_neutral() {
  return handle(std::make_shared<MS::NeutralValidation>(), true);
}

std::unique_ptr<ValidationMethodHandle> validation_isotope(bool strict) {
  return handle(std::make_shared<MS::IsotopeValidation>(strict));
}

std::unique_ptr<ValidationMethodHandle> validation_mol_vs_default() {
  return handle(std::make_shared<MS::MolVSValidation>(), true);
}

std::unique_ptr<ValidationMethodHandle> validation_mol_vs(
    const ValidationList &list) {
  if (list.items.empty()) {
    throw ValueErrorException("validations argument must be non-empty");
  }
  return handle(std::make_shared<MS::MolVSValidation>(list.items),
                list.sums_charges);
}

std::unique_ptr<ValidationMethodHandle> validation_allowed_atoms(
    const AtomList &list) {
  return handle(std::make_shared<MS::AllowedAtomsValidation>(
      atoms(list, "allowedAtoms")));
}

std::unique_ptr<ValidationMethodHandle> validation_disallowed_atoms(
    const AtomList &list) {
  return handle(std::make_shared<MS::DisallowedAtomsValidation>(
      atoms(list, "disallowedAtoms")));
}

std::unique_ptr<ValidationMethodHandle> validation_disallowed_radical() {
  return handle(std::make_shared<MS::DisallowedRadicalValidation>());
}

std::unique_ptr<ValidationMethodHandle> validation_features(
    const FeaturesValidationData &data) {
  return handle(std::make_shared<MS::FeaturesValidation>(
      data.allow_enhanced_stereo, data.allow_aromatic_bond_type,
      data.allow_dative_bond_type, data.allow_queries, data.allow_dummies,
      data.allow_atom_aliases));
}

std::unique_ptr<ValidationMethodHandle> validation_is_2d(double threshold) {
  return handle(std::make_shared<MS::Is2DValidation>(threshold));
}

std::unique_ptr<ValidationMethodHandle> validation_layout_2d(
    const Layout2DValidationData &data) {
  return handle(std::make_shared<MS::Layout2DValidation>(
      data.clash_limit, data.bond_length_limit, data.allow_long_bonds_in_rings,
      data.allow_atom_bond_clash_exemption, data.min_median_bond_length));
}

std::unique_ptr<ValidationMethodHandle> validation_stereo() {
  return handle(std::make_shared<MS::StereoValidation>());
}

std::unique_ptr<ValidationMethodHandle> validation_clone(
    const ValidationMethodHandle &original) {
  return handle(original.ptr, original.sums_charges);
}

// NeutralValidation sums the formal charges in int (MolOps::getFormalCharge),
// so a molecule whose charges could overflow it is refused first (U9).
// FragmentValidation, alone or among MolVSValidation's methods, matches
// substructures (control_c_lock.h); a Rust method runs without the lock.
std::unique_ptr<std::vector<std::string>> validation_validate(
    const ValidationMethodHandle &handle, const RWMol &mol,
    bool report_all_failures) {
  if (handle.sums_charges) {
    check_formal_charge_sum(mol);
  }
  // FragmentValidation splits the molecule into its fragments
  // (getMolFrags), removing atoms from copies (check_removals).
  check_removals(mol);
  HoldControlC control_c;
  return strings(handle.ptr->validate(mol, report_all_failures));
}

std::unique_ptr<ValidationList> validation_list_new() {
  return std::make_unique<ValidationList>();
}

void validation_list_push(ValidationList &list,
                          const ValidationMethodHandle &handle) {
  list.items.push_back(handle.ptr->copy());
  list.sums_charges = list.sums_charges || handle.sums_charges;
}

std::unique_ptr<AtomList> atom_list_new() {
  return std::make_unique<AtomList>();
}

void atom_list_push(AtomList &list, const Atom &atom) {
  list.items.emplace_back(atom.copy());
}

// As validateSmiles in Validate.cpp, with MolVSValidation's charge check.
std::unique_ptr<std::vector<std::string>> validate_smiles(rust::Str smiles) {
  std::string input(smiles);
  std::unique_ptr<RWMol> mol(RDKit::SmilesToMol(input));
  if (!mol) {
    throw ValueErrorException("SMILES Parse Error: syntax error for input: " +
                              input);
  }
  check_formal_charge_sum(*mol);
  // CXSMILES can give the substance groups PARENT chains
  // (validation_validate).
  check_removals(*mol);
  MS::MolVSValidation validation;
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return strings(validation.validate(*mol, true));
}

}  // namespace rdkit_cxx
