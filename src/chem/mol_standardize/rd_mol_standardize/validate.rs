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

/// Checks that each atom has a valence RDKit allows, and that the molecule
/// has atoms (unless empty molecules are allowed). Mirrors
/// `rdMolStandardize.RDKitValidation`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RDKitValidation {
    /// Allows molecules with no atoms. Mirrors
    /// `RDKitValidation.allowEmptyMolecules`.
    pub allow_empty_molecules: bool,
}

impl From<RDKitValidation> for ValidationMethod {
    fn from(validation: RDKitValidation) -> ValidationMethod {
        ValidationMethod::from_ffi(expect_ok(ffi::validation_rdkit(
            validation.allow_empty_molecules,
        )))
    }
}

validate_methods!(RDKitValidation);

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
    /// Reports atoms with unpaired electrons, except the radicals of nitric
    /// oxide and aminoxyls. Mirrors
    /// `rdMolStandardize.DisallowedRadicalValidation`.
    DisallowedRadicalValidation,
    validation_disallowed_radical
);

unit_validation!(
    /// Checks the wedged and dashed bonds of stereocenters against the 2D
    /// layout. Mirrors `rdMolStandardize.StereoValidation`.
    StereoValidation,
    validation_stereo
);

/// Reports isotopes: every one, or with `strict`, those RDKit's periodic
/// table does not know. Mirrors `rdMolStandardize.IsotopeValidation`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IsotopeValidation {
    /// Reports only unknown isotopes. Mirrors `IsotopeValidation.strict`.
    pub strict: bool,
}

impl From<IsotopeValidation> for ValidationMethod {
    fn from(validation: IsotopeValidation) -> ValidationMethod {
        ValidationMethod::from_ffi(expect_ok(ffi::validation_isotope(validation.strict)))
    }
}

validate_methods!(IsotopeValidation);

/// Reports features a molecule may not have: query atoms and bonds, dummy
/// atoms, atom aliases, aromatic or dative bonds, enhanced stereo. Mirrors
/// `rdMolStandardize.FeaturesValidation`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FeaturesValidation {
    /// Mirrors `FeaturesValidation.allowEnhancedStereo`.
    pub allow_enhanced_stereo: bool,
    /// Mirrors `FeaturesValidation.allowAromaticBondType`.
    pub allow_aromatic_bond_type: bool,
    /// Mirrors `FeaturesValidation.allowDativeBondType`.
    pub allow_dative_bond_type: bool,
    /// Mirrors `FeaturesValidation.allowQueries`.
    pub allow_queries: bool,
    /// Mirrors `FeaturesValidation.allowDummies`.
    pub allow_dummies: bool,
    /// Mirrors `FeaturesValidation.allowAtomAliases`.
    pub allow_atom_aliases: bool,
}

impl From<FeaturesValidation> for ValidationMethod {
    fn from(validation: FeaturesValidation) -> ValidationMethod {
        ValidationMethod::from_ffi(expect_ok(ffi::validation_features(
            &ffi::FeaturesValidationData {
                allow_enhanced_stereo: validation.allow_enhanced_stereo,
                allow_aromatic_bond_type: validation.allow_aromatic_bond_type,
                allow_dative_bond_type: validation.allow_dative_bond_type,
                allow_queries: validation.allow_queries,
                allow_dummies: validation.allow_dummies,
                allow_atom_aliases: validation.allow_atom_aliases,
            },
        )))
    }
}

validate_methods!(FeaturesValidation);

/// Checks that the molecule has 2D coordinates: a conformer, without Z
/// coordinates above the threshold, whose atoms are not all at one point.
/// Mirrors `rdMolStandardize.Is2DValidation`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Is2DValidation {
    /// The largest Z coordinate, and the smallest X or Y extent, taken as 0.
    /// Mirrors `Is2DValidation.threshold`.
    pub threshold: f64,
}

/// Mirrors `Is2DValidation()`.
impl Default for Is2DValidation {
    fn default() -> Is2DValidation {
        Is2DValidation { threshold: 1e-3 }
    }
}

impl From<Is2DValidation> for ValidationMethod {
    fn from(validation: Is2DValidation) -> ValidationMethod {
        ValidationMethod::from_ffi(expect_ok(ffi::validation_is_2d(validation.threshold)))
    }
}

validate_methods!(Is2DValidation);

/// Checks a 2D layout: for atoms too close to other atoms or to bonds, and
/// bonds too long, relative to the median bond length. Mirrors
/// `rdMolStandardize.Layout2DValidation`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Layout2DValidation {
    /// The closest atoms may come, as a fraction of the median bond length.
    /// Mirrors `Layout2DValidation.clashLimit`.
    pub clash_limit: f64,
    /// The longest a bond may be, as a multiple of the median bond length.
    /// Mirrors `Layout2DValidation.bondLengthLimit`.
    pub bond_length_limit: f64,
    /// Allows long bonds in rings. Mirrors
    /// `Layout2DValidation.allowLongBondsInRings`.
    pub allow_long_bonds_in_rings: bool,
    /// Exempts long ring bonds from the atom-bond clash check. Mirrors
    /// `Layout2DValidation.allowAtomBondClashExemption`.
    pub allow_atom_bond_clash_exemption: bool,
    /// The shortest the median bond length may be. Mirrors
    /// `Layout2DValidation.minMedianBondLength`.
    pub min_median_bond_length: f64,
}

/// Mirrors `Layout2DValidation()`, whose `minMedianBondLength` is 0 in
/// Python (the C++ default is 0.001).
impl Default for Layout2DValidation {
    fn default() -> Layout2DValidation {
        Layout2DValidation {
            clash_limit: 0.15,
            bond_length_limit: 25.0,
            allow_long_bonds_in_rings: true,
            allow_atom_bond_clash_exemption: true,
            min_median_bond_length: 0.0,
        }
    }
}

impl From<Layout2DValidation> for ValidationMethod {
    fn from(validation: Layout2DValidation) -> ValidationMethod {
        ValidationMethod::from_ffi(expect_ok(ffi::validation_layout_2d(
            &ffi::Layout2DValidationData {
                clash_limit: validation.clash_limit,
                bond_length_limit: validation.bond_length_limit,
                allow_long_bonds_in_rings: validation.allow_long_bonds_in_rings,
                allow_atom_bond_clash_exemption: validation.allow_atom_bond_clash_exemption,
                min_median_bond_length: validation.min_median_bond_length,
            },
        )))
    }
}

validate_methods!(Layout2DValidation);

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
/// the like) for a molecule that does not sanitize, or `Error::Value` for one
/// with too many relevant cycles (`sanitize_mol`). Mirrors
/// `rdMolStandardize.ValidateSmiles`.
pub fn validate_smiles(smiles: &str) -> Result<Vec<String>> {
    messages(ffi::validate_smiles(smiles)?)
}
