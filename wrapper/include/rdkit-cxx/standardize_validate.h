/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's validation methods (Validate.h), and validation methods
// that call Rust closures.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <GraphMol/MolStandardize/Validate.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Defined by the bridge (src/ffi/standardize_validate.rs).
struct ValidationCallback;
struct FeaturesValidationData;
struct Layout2DValidationData;

// A validation method, shared: copies of a handle share one, which RDKit
// only reads. sums_charges says whether its validation sums the molecule's
// formal charges in int (NeutralValidation, alone or in a MolVSValidation),
// which the wrapper checks first (U9).
struct ValidationMethodHandle {
  std::shared_ptr<RDKit::MolStandardize::ValidationMethod> ptr;
  bool sums_charges = false;
};

// Validation methods Rust lends to MolVSValidation's constructor.
struct ValidationList {
  std::vector<std::shared_ptr<RDKit::MolStandardize::ValidationMethod>> items;
  bool sums_charges = false;
};

// Atoms Rust lends to the atom-list validations, copied.
struct AtomList {
  std::vector<std::shared_ptr<RDKit::Atom>> items;
};

std::unique_ptr<ValidationMethodHandle> validation_rust(
    rust::Box<ValidationCallback> callback);
std::unique_ptr<ValidationMethodHandle> validation_rdkit(
    bool allow_empty_molecules);
std::unique_ptr<ValidationMethodHandle> validation_no_atom();
std::unique_ptr<ValidationMethodHandle> validation_fragment();
std::unique_ptr<ValidationMethodHandle> validation_neutral();
std::unique_ptr<ValidationMethodHandle> validation_isotope(bool strict);
std::unique_ptr<ValidationMethodHandle> validation_mol_vs_default();
// Copies each method, as Python's MolVSValidation(validations) does; an
// empty list is ValueErrorException, as there.
std::unique_ptr<ValidationMethodHandle> validation_mol_vs(
    const ValidationList &list);
std::unique_ptr<ValidationMethodHandle> validation_allowed_atoms(
    const AtomList &atoms);
std::unique_ptr<ValidationMethodHandle> validation_disallowed_atoms(
    const AtomList &atoms);
std::unique_ptr<ValidationMethodHandle> validation_disallowed_radical();
std::unique_ptr<ValidationMethodHandle> validation_features(
    const FeaturesValidationData &data);
std::unique_ptr<ValidationMethodHandle> validation_is_2d(double threshold);
std::unique_ptr<ValidationMethodHandle> validation_layout_2d(
    const Layout2DValidationData &data);
std::unique_ptr<ValidationMethodHandle> validation_stereo();

std::unique_ptr<ValidationMethodHandle> validation_clone(
    const ValidationMethodHandle &handle);
std::unique_ptr<std::vector<std::string>> validation_validate(
    const ValidationMethodHandle &handle, const RWMol &mol,
    bool report_all_failures);

std::unique_ptr<ValidationList> validation_list_new();
void validation_list_push(ValidationList &list,
                          const ValidationMethodHandle &handle);
std::unique_ptr<AtomList> atom_list_new();
void atom_list_push(AtomList &list, const Atom &atom);

std::unique_ptr<std::vector<std::string>> validate_smiles(rust::Str smiles);

}  // namespace rdkit_cxx
