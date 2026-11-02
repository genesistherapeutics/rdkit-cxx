// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reionizing and uncharging.

use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::standardize_charge::ffi;
use crate::util::path_bytes;

/// A charge to give atoms that match SMARTS before reionizing. Mirrors
/// `rdMolStandardize.ChargeCorrection`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChargeCorrection {
    /// The correction's name. Mirrors `ChargeCorrection.Name`.
    pub name: String,
    /// The atoms to correct. Mirrors `ChargeCorrection.Smarts`.
    pub smarts: String,
    /// The charge to give them. Mirrors `ChargeCorrection.Charge`.
    pub charge: i32,
}

impl ChargeCorrection {
    /// Mirrors `ChargeCorrection(name, smarts, charge)`.
    pub fn new(name: &str, smarts: &str, charge: i32) -> ChargeCorrection {
        ChargeCorrection {
            name: name.to_string(),
            smarts: smarts.to_string(),
            charge,
        }
    }

    fn to_ffi(&self) -> ffi::ChargeCorrectionData {
        ffi::ChargeCorrectionData {
            name: self.name.as_bytes().to_vec(),
            smarts: self.smarts.as_bytes().to_vec(),
            charge: self.charge,
        }
    }
}

fn to_ffi(ccs: &[ChargeCorrection]) -> Vec<ffi::ChargeCorrectionData> {
    ccs.iter().map(ChargeCorrection::to_ffi).collect()
}

/// RDKit's default charge corrections: alkali metals +1, Mg and Ca +2, and
/// chlorine -1, for atoms with no bonds or charge. Mirrors
/// `rdMolStandardize.CHARGE_CORRECTIONS`.
pub fn charge_corrections() -> Vec<ChargeCorrection> {
    let text = |bytes: Vec<u8>| String::from_utf8(bytes).expect("RDKit's corrections are UTF-8");
    expect_ok(ffi::charge_corrections())
        .into_iter()
        .map(|cc| ChargeCorrection {
            name: text(cc.name),
            smarts: text(cc.smarts),
            charge: cc.charge,
        })
        .collect()
}

/// Corrects charges, then moves protons so the strongest acids are the
/// ones ionized. Mirrors `rdMolStandardize.Reionizer`.
pub struct Reionizer {
    /// Never null.
    ptr: UniquePtr<ffi::Reionizer>,
}

impl std::fmt::Debug for Reionizer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reionizer").finish_non_exhaustive()
    }
}

/// RDKit's built-in acid-base pairs and `charge_corrections()`. Mirrors
/// `Reionizer()`.
impl Default for Reionizer {
    fn default() -> Reionizer {
        Reionizer::from_ffi(expect_ok(ffi::reionizer_default()))
    }
}

impl Reionizer {
    fn from_ffi(ptr: UniquePtr<ffi::Reionizer>) -> Reionizer {
        assert!(!ptr.is_null(), "RDKit returned a null Reionizer");
        Reionizer { ptr }
    }

    fn raw(&self) -> &ffi::Reionizer {
        self.ptr.as_ref().unwrap()
    }

    /// The acid-base pairs of a file (name, acid SMARTS, and base SMARTS,
    /// tab-separated, strongest acid first; empty for RDKit's built-in
    /// ones), with `charge_corrections()`. RDKit reads each file once per
    /// process and keeps it, under its name. `Error::Os` for a file that
    /// cannot be read and `Error::Value` for SMARTS that do not parse.
    /// Mirrors `Reionizer(acidbaseFile)`.
    pub fn new<P: AsRef<Path>>(acidbase_file: P) -> Result<Reionizer> {
        Ok(Reionizer::from_ffi(ffi::reionizer_new(
            path_bytes(acidbase_file.as_ref()),
            &[],
            true,
        )?))
    }

    /// With these charge corrections. `Error::Value` for a correction whose
    /// SMARTS do not parse, which RDKit would dereference as null. Mirrors
    /// `Reionizer(acidbaseFile, ccs)`.
    pub fn new_with_ccs<P: AsRef<Path>>(
        acidbase_file: P,
        ccs: &[ChargeCorrection],
    ) -> Result<Reionizer> {
        Ok(Reionizer::from_ffi(ffi::reionizer_new(
            path_bytes(acidbase_file.as_ref()),
            &to_ffi(ccs),
            false,
        )?))
    }

    /// A reionized copy of the molecule. `Error::Value` for a molecule of
    /// more than `i32::MAX / 255` atoms, whose charge sums could overflow
    /// RDKit's `int`s (CONVENTIONS.md U9). Mirrors `Reionizer.reionize`.
    pub fn reionize(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::reionizer_reionize(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// Reionizes the molecule. Mirrors `Reionizer.reionizeInPlace`.
    pub fn reionize_in_place(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::reionizer_reionize_in_place(self.raw(), mol.pin_mut())?)
    }
}

/// A reionizer from acid-base pairs, one per line: a name, acid SMARTS, and
/// base SMARTS, tab-separated, strongest acid first (lines starting `//`
/// are comments), with no charge corrections, as Python's default `None`
/// gives. `Error::Value` for SMARTS that do not parse. Mirrors
/// `rdMolStandardize.ReionizerFromData`.
pub fn reionizer_from_data(param_data: &str) -> Result<Reionizer> {
    reionizer_from_data_with_charge_corrections(param_data, &[])
}

/// Mirrors `rdMolStandardize.ReionizerFromData(paramData,
/// chargeCorrections)`. `Error::Value` for a correction whose SMARTS do not
/// parse, which RDKit would dereference as null.
pub fn reionizer_from_data_with_charge_corrections(
    param_data: &str,
    charge_corrections: &[ChargeCorrection],
) -> Result<Reionizer> {
    Ok(Reionizer::from_ffi(ffi::reionizer_from_data(
        param_data,
        &to_ffi(charge_corrections),
    )?))
}

/// Neutralizes ionized acids and bases by adding or removing hydrogens.
/// Mirrors `rdMolStandardize.Uncharger`.
pub struct Uncharger {
    /// Never null.
    ptr: UniquePtr<ffi::Uncharger>,
}

impl std::fmt::Debug for Uncharger {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Uncharger").finish_non_exhaustive()
    }
}

/// Parameters for `Uncharger::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UnchargerParams {
    /// Neutralizes sites in a canonical order.
    pub canonical_order: bool,
    /// Neutralizes every site it can, even if the molecule ends charged.
    pub force: bool,
    /// Only changes protonation, leaving charges that need a hydride.
    pub protonation_only: bool,
}

impl Default for UnchargerParams {
    fn default() -> UnchargerParams {
        UnchargerParams {
            canonical_order: true,
            force: false,
            protonation_only: false,
        }
    }
}

/// Mirrors `Uncharger()`.
impl Default for Uncharger {
    fn default() -> Uncharger {
        Uncharger::new_with_params(&UnchargerParams::default())
    }
}

impl Uncharger {
    fn raw(&self) -> &ffi::Uncharger {
        self.ptr.as_ref().unwrap()
    }

    /// Mirrors `Uncharger(canonicalOrder, force, protonationOnly)`.
    pub fn new_with_params(params: &UnchargerParams) -> Uncharger {
        let ptr = expect_ok(ffi::uncharger_new(
            params.canonical_order,
            params.force,
            params.protonation_only,
        ));
        assert!(!ptr.is_null(), "RDKit returned a null Uncharger");
        Uncharger { ptr }
    }

    /// An uncharged copy of the molecule. `Error::Value` for a molecule
    /// whose charges could overflow RDKit's `int` sum (CONVENTIONS.md U9).
    /// Mirrors `Uncharger.uncharge`.
    pub fn uncharge(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::uncharger_uncharge(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// Uncharges the molecule. Mirrors `Uncharger.unchargeInPlace`.
    pub fn uncharge_in_place(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::uncharger_uncharge_in_place(self.raw(), mol.pin_mut())?)
    }
}
