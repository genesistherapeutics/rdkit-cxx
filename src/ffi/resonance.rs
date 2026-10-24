// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::ResonanceMolSupplier::ResonanceFlags`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum ResonanceFlags {
        ALLOW_INCOMPLETE_OCTETS = 0x1,
        ALLOW_CHARGE_SEPARATION = 0x2,
        KEKULE_ALL = 0x4,
        UNCONSTRAINED_CATIONS = 0x8,
        UNCONSTRAINED_ANIONS = 0x10,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/resonance.h");

        type ResonanceMolSupplier;
        type ResonanceFlags;
        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn resonance_mol_supplier_new(
            mol: &RWMol,
            flags: u32,
            max_structs: u32,
        ) -> Result<UniquePtr<ResonanceMolSupplier>>;
        fn resonance_mol_supplier_get_num_conj_grps(suppl: &ResonanceMolSupplier) -> Result<u32>;
        fn resonance_mol_supplier_get_bond_conj_grp_idx(
            suppl: &ResonanceMolSupplier,
            bi: u32,
        ) -> Result<i32>;
        fn resonance_mol_supplier_get_atom_conj_grp_idx(
            suppl: &ResonanceMolSupplier,
            ai: u32,
        ) -> Result<i32>;
        fn resonance_mol_supplier_set_num_threads(
            suppl: Pin<&mut ResonanceMolSupplier>,
            num_threads: i32,
        ) -> Result<()>;
        fn resonance_mol_supplier_was_canceled(suppl: &ResonanceMolSupplier) -> Result<bool>;
        fn resonance_mol_supplier_enumerate(suppl: Pin<&mut ResonanceMolSupplier>) -> Result<()>;
        fn resonance_mol_supplier_get_is_enumerated(suppl: &ResonanceMolSupplier) -> Result<bool>;
        fn resonance_mol_supplier_length(suppl: Pin<&mut ResonanceMolSupplier>) -> Result<u32>;
        fn resonance_mol_supplier_reset(suppl: Pin<&mut ResonanceMolSupplier>) -> Result<()>;
        fn resonance_mol_supplier_at_end(suppl: Pin<&mut ResonanceMolSupplier>) -> Result<bool>;
        fn resonance_mol_supplier_next(
            suppl: Pin<&mut ResonanceMolSupplier>,
        ) -> Result<UniquePtr<RWMol>>;
        fn resonance_mol_supplier_get_item(
            suppl: Pin<&mut ResonanceMolSupplier>,
            idx: u32,
        ) -> Result<UniquePtr<RWMol>>;
    }
}

// SAFETY: a supplier owns a copy of its molecule and its enumeration state,
// with no thread affinity; the crate never sets a progress callback. Like
// molecules, it is not Sync (CONVENTIONS.md S1).
unsafe impl Send for ffi::ResonanceMolSupplier {}
