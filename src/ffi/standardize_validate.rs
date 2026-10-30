// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::panic::{AssertUnwindSafe, catch_unwind};

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::ffi::mol::ffi::RWMol;

/// A validation method's closure: the molecule and whether to report every
/// failure, to the messages.
pub(crate) type ValidateFn = dyn Fn(&Mol, bool) -> Vec<String> + Send + Sync;

/// The closure behind a validation method made in Rust. Copies of the
/// method share it, so threads holding them may call it at once.
pub(crate) struct ValidationCallback(pub(crate) Box<ValidateFn>);

/// A panic's message, for the exception C++ sees.
pub(crate) fn panic_message(payload: Box<dyn std::any::Any + Send>, what: &str) -> String {
    let message = payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| format!("{what} panicked"));
    format!("{what} panicked: {message}")
}

/// Called by the C++ method with its own copy of the molecule. A panic is
/// caught here, at the boundary, and becomes the exception C++ sees.
fn call_validation_callback(
    callback: &ValidationCallback,
    mol: UniquePtr<RWMol>,
    report_all_failures: bool,
) -> Result<Vec<String>, String> {
    let mol = Mol::from_ffi(mol);
    catch_unwind(AssertUnwindSafe(|| (callback.0)(&mol, report_all_failures)))
        .map_err(|payload| panic_message(payload, "a validation method"))
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolStandardize::FeaturesValidation`'s fields.
    struct FeaturesValidationData {
        allow_enhanced_stereo: bool,
        allow_aromatic_bond_type: bool,
        allow_dative_bond_type: bool,
        allow_queries: bool,
        allow_dummies: bool,
        allow_atom_aliases: bool,
    }

    /// `RDKit::MolStandardize::Layout2DValidation`'s fields.
    struct Layout2DValidationData {
        clash_limit: f64,
        bond_length_limit: f64,
        allow_long_bonds_in_rings: bool,
        allow_atom_bond_clash_exemption: bool,
        min_median_bond_length: f64,
    }

    extern "Rust" {
        type ValidationCallback;

        fn call_validation_callback(
            callback: &ValidationCallback,
            mol: UniquePtr<RWMol>,
            report_all_failures: bool,
        ) -> Result<Vec<String>>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize_validate.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;
        type ValidationMethodHandle;
        type ValidationList;
        type AtomList;

        fn validation_rust(
            callback: Box<ValidationCallback>,
        ) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_rdkit(
            allow_empty_molecules: bool,
        ) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_no_atom() -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_fragment() -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_neutral() -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_isotope(strict: bool) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_mol_vs_default() -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_mol_vs(list: &ValidationList) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_allowed_atoms(atoms: &AtomList) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_disallowed_atoms(
            atoms: &AtomList,
        ) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_disallowed_radical() -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_features(
            data: &FeaturesValidationData,
        ) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_is_2d(threshold: f64) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_layout_2d(
            data: &Layout2DValidationData,
        ) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_stereo() -> Result<UniquePtr<ValidationMethodHandle>>;

        fn validation_clone(
            handle: &ValidationMethodHandle,
        ) -> Result<UniquePtr<ValidationMethodHandle>>;
        fn validation_validate(
            handle: &ValidationMethodHandle,
            mol: &RWMol,
            report_all_failures: bool,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;

        fn validation_list_new() -> Result<UniquePtr<ValidationList>>;
        fn validation_list_push(
            list: Pin<&mut ValidationList>,
            handle: &ValidationMethodHandle,
        ) -> Result<()>;
        fn atom_list_new() -> Result<UniquePtr<AtomList>>;
        fn atom_list_push(list: Pin<&mut AtomList>, atom: &Atom) -> Result<()>;

        fn validate_smiles(smiles: &str) -> Result<UniquePtr<CxxVector<CxxString>>>;
    }
}

// SAFETY: a handle shares (through a shared_ptr, with atomic counts) a
// validation method that nothing changes after it is made: RDKit's methods
// hold their options and copied atoms, whose const validate only reads them,
// and Rust's hold a Send + Sync closure.
unsafe impl Send for ffi::ValidationMethodHandle {}
unsafe impl Sync for ffi::ValidationMethodHandle {}
