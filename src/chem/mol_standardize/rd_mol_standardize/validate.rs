// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Validating molecules.

use cxx::{CxxString, CxxVector, UniquePtr};

use crate::Result;
use crate::chem::{AsAtom, Mol};
use crate::error::expect_ok;
use crate::ffi::standardize_validate::{ValidationCallback, ffi};

/// The messages, which may hold text from the molecule (an atom's alias),
/// as Rust strings: `Error::Value` for any that is not UTF-8.
fn messages(values: UniquePtr<CxxVector<CxxString>>) -> Result<Vec<String>> {
    values
        .iter()
        .map(|value| {
            value
                .to_str()
                .map(str::to_owned)
                .map_err(|err| crate::Error::Value {
                    message: format!("invalid UTF-8 in a validation message: {err}"),
                })
        })
        .collect()
}

/// A validation method, which lists what it finds wrong with a molecule:
/// one of RDKit's, converted with `From`, or one made from a Rust closure.
/// Mirrors `rdMolStandardize.ValidationMethod`, the class Python subclasses
/// to define one. Clones share the method.
pub struct ValidationMethod {
    /// Never null.
    ptr: UniquePtr<ffi::ValidationMethodHandle>,
}

impl std::fmt::Debug for ValidationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ValidationMethod").finish_non_exhaustive()
    }
}

impl Clone for ValidationMethod {
    fn clone(&self) -> ValidationMethod {
        ValidationMethod::from_ffi(expect_ok(ffi::validation_clone(self.raw())))
    }
}

impl ValidationMethod {
    fn from_ffi(ptr: UniquePtr<ffi::ValidationMethodHandle>) -> ValidationMethod {
        assert!(!ptr.is_null(), "RDKit returned a null validation method");
        ValidationMethod { ptr }
    }

    fn raw(&self) -> &ffi::ValidationMethodHandle {
        self.ptr.as_ref().unwrap()
    }

    /// A validation method that calls `validate` with a copy of the
    /// molecule and `report_all_failures`, and returns its messages. Python
    /// subclasses the class and overrides `validate` (and `copy`); here the
    /// closure is the argument. A `MolVSValidation` made from the method
    /// shares it, so the closure must be `Send + Sync`. A panic in it is
    /// caught, and the call that made it returns `Error::Runtime`. Mirrors
    /// `ValidationMethod()`.
    pub fn new<F>(validate: F) -> ValidationMethod
    where
        F: Fn(&Mol, bool) -> Vec<String> + Send + Sync + 'static,
    {
        let callback = Box::new(ValidationCallback(Box::new(validate)));
        ValidationMethod::from_ffi(expect_ok(ffi::validation_rust(callback)))
    }

    /// What the method finds wrong with the molecule, stopping at the first
    /// problem. `Error::Value` for a message that is not UTF-8. Mirrors
    /// `ValidationMethod.validate`.
    pub fn validate(&self, mol: &Mol) -> Result<Vec<String>> {
        self.validate_with_report_all_failures(mol, false)
    }

    /// With `report_all_failures`, every problem. A method that sums the
    /// formal charges (`NeutralValidation`, alone or in a
    /// `MolVSValidation`) returns `Error::Value` for a molecule whose
    /// charges could overflow RDKit's `int` sum (CONVENTIONS.md U9). Mirrors
    /// `ValidationMethod.validate(mol, reportAllFailures)`.
    pub fn validate_with_report_all_failures(
        &self,
        mol: &Mol,
        report_all_failures: bool,
    ) -> Result<Vec<String>> {
        messages(ffi::validation_validate(
            self.raw(),
            mol.as_ffi(),
            report_all_failures,
        )?)
    }
}

/// `validate` and `validate_with_report_all_failures` for a validation,
/// through its `ValidationMethod`.
macro_rules! validate_methods {
    ($name:ident) => {
        impl $name {
            #[doc = concat!("What the validation finds wrong with the molecule, stopping at the first problem. Mirrors `", stringify!($name), ".validate`.")]
            pub fn validate(&self, mol: &Mol) -> Result<Vec<String>> {
                ValidationMethod::from(self.clone()).validate(mol)
            }

            #[doc = concat!("With `report_all_failures`, every problem. Mirrors `", stringify!($name), ".validate(mol, reportAllFailures)`.")]
            pub fn validate_with_report_all_failures(
                &self,
                mol: &Mol,
                report_all_failures: bool,
            ) -> Result<Vec<String>> {
                ValidationMethod::from(self.clone())
                    .validate_with_report_all_failures(mol, report_all_failures)
            }
        }
    };
}

/// A validation with no options: a unit struct, its `Default` Python's
/// constructor.
macro_rules! unit_validation {
    ($(#[$meta:meta])* $name:ident, $new:ident) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
        pub struct $name;

        impl From<$name> for ValidationMethod {
            fn from(_: $name) -> ValidationMethod {
                ValidationMethod::from_ffi(expect_ok(ffi::$new()))
            }
        }

        validate_methods!($name);
    };
}

unit_validation!(
    /// Checks that each atom has a valence RDKit allows. Mirrors
    /// `rdMolStandardize.RDKitValidation`.
    RDKitValidation,
    validation_rdkit
);

unit_validation!(
    /// Checks that the molecule has atoms. Mirrors
    /// `rdMolStandardize.NoAtomValidation`.
    NoAtomValidation,
    validation_no_atom
);

unit_validation!(
    /// Reports the fragments of RDKit's built-in list (solvents,
    /// counter-ions) the molecule has. Mirrors
    /// `rdMolStandardize.FragmentValidation`.
    FragmentValidation,
    validation_fragment
);

unit_validation!(
    /// Reports a net charge. `Error::Value` for a molecule whose charges
    /// could overflow RDKit's `int` sum (CONVENTIONS.md U9). Mirrors
    /// `rdMolStandardize.NeutralValidation`.
    NeutralValidation,
    validation_neutral
);

unit_validation!(
    /// Reports isotopes. Mirrors `rdMolStandardize.IsotopeValidation`.
    IsotopeValidation,
    validation_isotope
);

/// A validation that holds other validation methods, or RDKit's atoms, and
/// so is a `ValidationMethod` itself.
macro_rules! method_validation {
    ($name:ident) => {
        impl From<$name> for ValidationMethod {
            fn from(validation: $name) -> ValidationMethod {
                validation.method
            }
        }

        impl $name {
            #[doc = concat!("What the validation finds wrong with the molecule, stopping at the first problem. Mirrors `", stringify!($name), ".validate`.")]
            pub fn validate(&self, mol: &Mol) -> Result<Vec<String>> {
                self.method.validate(mol)
            }

            #[doc = concat!("With `report_all_failures`, every problem. Mirrors `", stringify!($name), ".validate(mol, reportAllFailures)`.")]
            pub fn validate_with_report_all_failures(
                &self,
                mol: &Mol,
                report_all_failures: bool,
            ) -> Result<Vec<String>> {
                self.method
                    .validate_with_report_all_failures(mol, report_all_failures)
            }
        }
    };
}

/// MolVS's validations, in order, stopping at the first that fails unless
/// every failure is reported. Mirrors `rdMolStandardize.MolVSValidation`.
#[derive(Clone)]
pub struct MolVSValidation {
    method: ValidationMethod,
}

impl std::fmt::Debug for MolVSValidation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MolVSValidation").finish_non_exhaustive()
    }
}

/// `NoAtomValidation`, `FragmentValidation`, `NeutralValidation`, and
/// `IsotopeValidation`. Mirrors `MolVSValidation()`.
impl Default for MolVSValidation {
    fn default() -> MolVSValidation {
        MolVSValidation {
            method: ValidationMethod::from_ffi(expect_ok(ffi::validation_mol_vs_default())),
        }
    }
}

impl MolVSValidation {
    /// These validations, copied (a Rust one shares its closure).
    /// `Error::Value` for none, as in Python. Mirrors
    /// `MolVSValidation(validations)`.
    pub fn new(validations: &[ValidationMethod]) -> Result<MolVSValidation> {
        let mut list = expect_ok(ffi::validation_list_new());
        for validation in validations {
            expect_ok(ffi::validation_list_push(list.pin_mut(), validation.raw()));
        }
        Ok(MolVSValidation {
            method: ValidationMethod::from_ffi(ffi::validation_mol_vs(&list)?),
        })
    }
}

method_validation!(MolVSValidation);

fn atom_list(atoms: &[&dyn AsAtom]) -> UniquePtr<ffi::AtomList> {
    let mut list = expect_ok(ffi::atom_list_new());
    for atom in atoms {
        expect_ok(ffi::atom_list_push(list.pin_mut(), atom.atom()));
    }
    list
}

/// Reports atoms that match none of its atoms (or query atoms). Mirrors
/// `rdMolStandardize.AllowedAtomsValidation`.
#[derive(Clone)]
pub struct AllowedAtomsValidation {
    method: ValidationMethod,
}

impl std::fmt::Debug for AllowedAtomsValidation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AllowedAtomsValidation")
            .finish_non_exhaustive()
    }
}

impl AllowedAtomsValidation {
    /// The allowed atoms, copied. `Error::Value` for none, as in Python.
    /// Mirrors `AllowedAtomsValidation(atoms)`.
    pub fn new(atoms: &[&dyn AsAtom]) -> Result<AllowedAtomsValidation> {
        Ok(AllowedAtomsValidation {
            method: ValidationMethod::from_ffi(ffi::validation_allowed_atoms(&atom_list(atoms))?),
        })
    }
}

method_validation!(AllowedAtomsValidation);

/// Reports atoms that match any of its atoms (or query atoms). Mirrors
/// `rdMolStandardize.DisallowedAtomsValidation`.
#[derive(Clone)]
pub struct DisallowedAtomsValidation {
    method: ValidationMethod,
}

impl std::fmt::Debug for DisallowedAtomsValidation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DisallowedAtomsValidation")
            .finish_non_exhaustive()
    }
}

impl DisallowedAtomsValidation {
    /// The disallowed atoms, copied. `Error::Value` for none, as in Python.
    /// Mirrors `DisallowedAtomsValidation(atoms)`.
    pub fn new(atoms: &[&dyn AsAtom]) -> Result<DisallowedAtomsValidation> {
        Ok(DisallowedAtomsValidation {
            method: ValidationMethod::from_ffi(ffi::validation_disallowed_atoms(&atom_list(
                atoms,
            ))?),
        })
    }
}

method_validation!(DisallowedAtomsValidation);

/// What `MolVSValidation::default()` finds wrong with a molecule parsed from
/// SMILES, reporting every failure. `Error::Value` for SMILES that do not
/// parse, as in Python, and sanitization's errors (`Error::AtomValence` and
/// the like) for a molecule that does not sanitize. Mirrors
/// `rdMolStandardize.ValidateSmiles`.
pub fn validate_smiles(smiles: &str) -> Result<Vec<String>> {
    messages(ffi::validate_smiles(smiles)?)
}
