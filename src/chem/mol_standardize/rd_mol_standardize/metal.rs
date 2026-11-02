// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Disconnecting metals.

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::standardize_normalize::ffi;

/// Options for disconnecting metals. Mirrors
/// `rdMolStandardize.MetalDisconnectorOptions`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MetalDisconnectorOptions {
    /// Splits Grignard-type complexes.
    pub split_grignards: bool,
    /// Splits metal-aromatic carbon bonds.
    pub split_aromatic_c: bool,
    /// Adjusts the charges of the atoms that were bonded.
    pub adjust_charges: bool,
    /// Removes the dummy atoms that stand for haptic bonds: those bonded to
    /// a metal by a bond with the `_MolFileBondEndPts` property.
    pub remove_haptic_dummies: bool,
}

/// Mirrors `MetalDisconnectorOptions()`.
impl Default for MetalDisconnectorOptions {
    fn default() -> MetalDisconnectorOptions {
        MetalDisconnectorOptions::from_ffi(expect_ok(ffi::metal_disconnector_options_default()))
    }
}

impl MetalDisconnectorOptions {
    fn from_ffi(data: ffi::MetalDisconnectorOptionsData) -> MetalDisconnectorOptions {
        MetalDisconnectorOptions {
            split_grignards: data.split_grignards,
            split_aromatic_c: data.split_aromatic_c,
            adjust_charges: data.adjust_charges,
            remove_haptic_dummies: data.remove_haptic_dummies,
        }
    }

    fn to_ffi(self) -> ffi::MetalDisconnectorOptionsData {
        ffi::MetalDisconnectorOptionsData {
            split_grignards: self.split_grignards,
            split_aromatic_c: self.split_aromatic_c,
            adjust_charges: self.adjust_charges,
            remove_haptic_dummies: self.remove_haptic_dummies,
        }
    }
}

/// Disconnects metals bonded to nitrogen, oxygen, or fluorine, and
/// transition metals bonded to other non-metals, adjusting the charges.
/// Mirrors `rdMolStandardize.MetalDisconnector`.
pub struct MetalDisconnector {
    /// Never null.
    ptr: UniquePtr<ffi::MetalDisconnector>,
}

impl std::fmt::Debug for MetalDisconnector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MetalDisconnector").finish_non_exhaustive()
    }
}

/// Mirrors `MetalDisconnector()`.
impl Default for MetalDisconnector {
    fn default() -> MetalDisconnector {
        MetalDisconnector::new()
    }
}

impl MetalDisconnector {
    fn raw(&self) -> &ffi::MetalDisconnector {
        self.ptr.as_ref().unwrap()
    }

    /// A disconnector with `MetalDisconnectorOptions`' defaults. Mirrors
    /// `MetalDisconnector()`.
    pub fn new() -> MetalDisconnector {
        MetalDisconnector::new_with_options(&MetalDisconnectorOptions::default())
    }

    /// Mirrors `MetalDisconnector(options)`.
    pub fn new_with_options(options: &MetalDisconnectorOptions) -> MetalDisconnector {
        let ptr = expect_ok(ffi::metal_disconnector_new(&options.to_ffi()));
        assert!(!ptr.is_null(), "RDKit returned a null MetalDisconnector");
        MetalDisconnector { ptr }
    }

    /// The SMARTS of the metals to disconnect from nitrogen, oxygen, and
    /// fluorine, and the atom bonded to each. Mirrors
    /// `MetalDisconnector.MetalNof`.
    pub fn metal_nof(&self) -> Result<String> {
        to_string(ffi::metal_disconnector_get_metal_nof(self.raw())?)
    }

    /// The SMARTS of the metals to disconnect from other non-metals, and the
    /// atom bonded to each. Mirrors `MetalDisconnector.MetalNon`.
    pub fn metal_non(&self) -> Result<String> {
        to_string(ffi::metal_disconnector_get_metal_non(self.raw())?)
    }

    /// Sets the query `metal_nof` describes: the metal, then the atom bonded
    /// to it. `Error::Value` unless it has those two atoms, joined by a
    /// bond, or none: RDKit reads each match's second atom and the bond
    /// between the two without checking either, and Python crashes for a
    /// query of one atom, of two unbonded atoms, or of three. Mirrors
    /// `MetalDisconnector.SetMetalNof`.
    pub fn set_metal_nof(&mut self, mol: &Mol) -> Result<()> {
        Ok(ffi::metal_disconnector_set_metal_nof(
            self.ptr.pin_mut(),
            mol.as_ffi(),
        )?)
    }

    /// Sets the query `metal_non` describes, as `set_metal_nof` does.
    /// Mirrors `MetalDisconnector.SetMetalNon`.
    pub fn set_metal_non(&mut self, mol: &Mol) -> Result<()> {
        Ok(ffi::metal_disconnector_set_metal_non(
            self.ptr.pin_mut(),
            mol.as_ffi(),
        )?)
    }

    /// A disconnected copy of the molecule. Mirrors
    /// `MetalDisconnector.Disconnect`.
    pub fn disconnect(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::metal_disconnector_disconnect(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// Disconnects the molecule. Mirrors
    /// `MetalDisconnector.DisconnectInPlace`.
    pub fn disconnect_in_place(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::metal_disconnector_disconnect_in_place(
            self.raw(),
            mol.pin_mut(),
        )?)
    }
}

/// Disconnects every bond to a metal, by the rules Syngenta prefers:
/// Grignards and metal-aromatic carbon bonds too, without adjusting
/// charges, and removing haptic dummies. These are the C++ function's
/// defaults, not `MetalDisconnectorOptions`'. Mirrors
/// `rdMolStandardize.DisconnectOrganometallics`.
pub fn disconnect_organometallics(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::disconnect_organometallics(
        mol.as_ffi(),
    )?))
}

/// Mirrors `rdMolStandardize.DisconnectOrganometallics(mol, params)`.
pub fn disconnect_organometallics_with_params(
    mol: &Mol,
    params: &MetalDisconnectorOptions,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::disconnect_organometallics_with_options(
        mol.as_ffi(),
        &params.to_ffi(),
    )?))
}

/// As `disconnect_organometallics`, in place. Mirrors
/// `rdMolStandardize.DisconnectOrganometallicsInPlace`.
pub fn disconnect_organometallics_in_place(mol: &mut Mol) -> Result<()> {
    Ok(ffi::disconnect_organometallics_in_place(mol.pin_mut())?)
}

/// Mirrors `rdMolStandardize.DisconnectOrganometallicsInPlace(mol,
/// params)`.
pub fn disconnect_organometallics_in_place_with_params(
    mol: &mut Mol,
    params: &MetalDisconnectorOptions,
) -> Result<()> {
    Ok(ffi::disconnect_organometallics_in_place_with_options(
        mol.pin_mut(),
        &params.to_ffi(),
    )?)
}
