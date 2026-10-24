// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/ring_info.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type IndexList = crate::ffi::geometry::ffi::IndexList;

        fn ring_info_num_rings(mol: &RWMol) -> Result<u32>;
        fn ring_info_atom_rings(mol: &RWMol) -> Result<Vec<IndexList>>;
        fn ring_info_bond_rings(mol: &RWMol) -> Result<Vec<IndexList>>;
        fn ring_info_num_atom_rings(mol: &RWMol, idx: u32) -> Result<u32>;
        fn ring_info_num_bond_rings(mol: &RWMol, idx: u32) -> Result<u32>;
        fn ring_info_is_atom_in_ring_of_size(mol: &RWMol, idx: u32, size: u32) -> Result<bool>;
        fn ring_info_is_bond_in_ring_of_size(mol: &RWMol, idx: u32, size: u32) -> Result<bool>;
        fn ring_info_atom_members(mol: &RWMol, idx: u32) -> Result<Vec<i32>>;
        fn ring_info_bond_members(mol: &RWMol, idx: u32) -> Result<Vec<i32>>;
        fn ring_info_atom_ring_sizes(mol: &RWMol, idx: u32) -> Result<Vec<i32>>;
        fn ring_info_bond_ring_sizes(mol: &RWMol, idx: u32) -> Result<Vec<i32>>;
        fn ring_info_min_atom_ring_size(mol: &RWMol, idx: u32) -> Result<u32>;
        fn ring_info_min_bond_ring_size(mol: &RWMol, idx: u32) -> Result<u32>;
        fn ring_info_are_atoms_in_same_ring(mol: &RWMol, idx1: u32, idx2: u32) -> Result<bool>;
        fn ring_info_are_atoms_in_same_ring_of_size(
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            size: u32,
        ) -> Result<bool>;
        fn ring_info_are_bonds_in_same_ring(mol: &RWMol, idx1: u32, idx2: u32) -> Result<bool>;
        fn ring_info_are_bonds_in_same_ring_of_size(
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            size: u32,
        ) -> Result<bool>;
        fn ring_info_is_ring_fused(mol: &RWMol, ring_idx: u32) -> Result<bool>;
        fn ring_info_are_rings_fused(mol: &RWMol, ring1_idx: u32, ring2_idx: u32) -> Result<bool>;
        fn ring_info_num_fused_bonds(mol: &RWMol, ring_idx: u32) -> Result<u32>;
        fn ring_info_are_ring_families_initialized(mol: &RWMol) -> Result<bool>;
        fn ring_info_num_ring_families(mol: &RWMol) -> Result<u32>;
        fn ring_info_num_relevant_cycles(mol: &RWMol) -> Result<u32>;
        fn ring_info_atom_ring_families(mol: &RWMol) -> Result<Vec<IndexList>>;
        fn ring_info_bond_ring_families(mol: &RWMol) -> Result<Vec<IndexList>>;
        fn ring_info_add_ring(
            mol: Pin<&mut RWMol>,
            atom_ids: &[u32],
            bond_ids: &[u32],
        ) -> Result<()>;
    }
}
