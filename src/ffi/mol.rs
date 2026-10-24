// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol.h");

        type RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;
        type Bond = crate::ffi::bond::ffi::Bond;
        type BondType = crate::ffi::bond::ffi::BondType;

        fn mol_new() -> Result<UniquePtr<RWMol>>;
        fn mol_clone(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn mol_get_num_atoms(mol: &RWMol) -> Result<u32>;
        fn mol_get_num_atoms_with_only_explicit(mol: &RWMol, only_explicit: bool) -> Result<u32>;
        fn mol_get_num_heavy_atoms(mol: &RWMol) -> Result<u32>;
        fn mol_get_num_bonds(mol: &RWMol, only_heavy: bool) -> Result<u32>;
        fn mol_get_atom_with_idx(mol: &RWMol, idx: u32) -> Result<&Atom>;
        fn mol_get_bond_with_idx(mol: &RWMol, idx: u32) -> Result<&Bond>;
        fn mol_get_bond_between_atoms(mol: &RWMol, idx1: u32, idx2: u32) -> Result<i64>;
        fn mol_get_atom_with_idx_mut(mol: Pin<&mut RWMol>, idx: u32) -> Result<Pin<&mut Atom>>;
        fn mol_get_bond_with_idx_mut(mol: Pin<&mut RWMol>, idx: u32) -> Result<Pin<&mut Bond>>;
        fn mol_update_property_cache(mol: Pin<&mut RWMol>, strict: bool) -> Result<()>;

        fn rwmol_add_atom(mol: Pin<&mut RWMol>, atom: &Atom) -> Result<u32>;
        fn rwmol_add_bond(
            mol: Pin<&mut RWMol>,
            begin_atom_idx: u32,
            end_atom_idx: u32,
            order: BondType,
        ) -> Result<u32>;
        fn rwmol_remove_atom(mol: Pin<&mut RWMol>, idx: u32) -> Result<()>;
        fn rwmol_remove_bond(mol: Pin<&mut RWMol>, idx1: u32, idx2: u32) -> Result<()>;
        fn rwmol_replace_atom(
            mol: Pin<&mut RWMol>,
            index: u32,
            new_atom: &Atom,
            update_label: bool,
            preserve_props: bool,
        ) -> Result<()>;
    }
}

// SAFETY: an RWMol owns all of its data and has no thread affinity, so it can
// move between threads. It is not Sync: RDKit writes computed properties
// through const molecules (CONVENTIONS.md S1).
unsafe impl Send for ffi::RWMol {}
