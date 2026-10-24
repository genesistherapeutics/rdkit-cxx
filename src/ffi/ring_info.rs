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
    }
}
