// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/monomer_info.h");

        type AtomMonomerInfo;
        type AtomPDBResidueInfo;
        type Atom = crate::ffi::atom::ffi::Atom;

        fn monomer_info_new() -> Result<UniquePtr<AtomMonomerInfo>>;
        fn monomer_info_with(monomer_type: i32, name: &str) -> Result<UniquePtr<AtomMonomerInfo>>;
        fn monomer_info_clone(info: &AtomMonomerInfo) -> Result<UniquePtr<AtomMonomerInfo>>;
        fn monomer_info_get_name(info: &AtomMonomerInfo) -> Result<UniquePtr<CxxString>>;
        fn monomer_info_get_monomer_type(info: &AtomMonomerInfo) -> Result<i32>;
        fn monomer_info_set_name(info: Pin<&mut AtomMonomerInfo>, nm: &str) -> Result<()>;
        fn monomer_info_set_monomer_type(info: Pin<&mut AtomMonomerInfo>, typ: i32) -> Result<()>;

        fn pdb_residue_info_new() -> Result<UniquePtr<AtomPDBResidueInfo>>;
        #[allow(clippy::too_many_arguments)]
        fn pdb_residue_info_with(
            atom_name: &str,
            serial_number: i32,
            alt_loc: &str,
            residue_name: &str,
            residue_number: i32,
            chain_id: &str,
            insertion_code: &str,
            occupancy: f64,
            temp_factor: f64,
            is_hetero_atom: bool,
            secondary_structure: u32,
            segment_number: u32,
        ) -> Result<UniquePtr<AtomPDBResidueInfo>>;
        fn pdb_residue_info_clone(
            info: &AtomPDBResidueInfo,
        ) -> Result<UniquePtr<AtomPDBResidueInfo>>;
        fn pdb_residue_info_as_monomer_info(info: &AtomPDBResidueInfo) -> Result<&AtomMonomerInfo>;
        fn pdb_residue_info_as_monomer_info_mut(
            info: Pin<&mut AtomPDBResidueInfo>,
        ) -> Result<Pin<&mut AtomMonomerInfo>>;
        fn pdb_residue_info_get_serial_number(info: &AtomPDBResidueInfo) -> Result<i32>;
        fn pdb_residue_info_get_residue_name(
            info: &AtomPDBResidueInfo,
        ) -> Result<UniquePtr<CxxString>>;
        fn pdb_residue_info_get_residue_number(info: &AtomPDBResidueInfo) -> Result<i32>;
        fn pdb_residue_info_get_chain_id(info: &AtomPDBResidueInfo) -> Result<UniquePtr<CxxString>>;
        fn pdb_residue_info_set_residue_name(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: &str,
        ) -> Result<()>;
        fn pdb_residue_info_set_residue_number(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: i32,
        ) -> Result<()>;
        fn pdb_residue_info_set_chain_id(info: Pin<&mut AtomPDBResidueInfo>, val: &str)
        -> Result<()>;
        fn pdb_residue_info_get_alt_loc(info: &AtomPDBResidueInfo) -> Result<UniquePtr<CxxString>>;
        fn pdb_residue_info_get_insertion_code(
            info: &AtomPDBResidueInfo,
        ) -> Result<UniquePtr<CxxString>>;
        fn pdb_residue_info_get_occupancy(info: &AtomPDBResidueInfo) -> Result<f64>;
        fn pdb_residue_info_get_temp_factor(info: &AtomPDBResidueInfo) -> Result<f64>;
        fn pdb_residue_info_get_is_hetero_atom(info: &AtomPDBResidueInfo) -> Result<bool>;
        fn pdb_residue_info_get_secondary_structure(info: &AtomPDBResidueInfo) -> Result<u32>;
        fn pdb_residue_info_get_segment_number(info: &AtomPDBResidueInfo) -> Result<u32>;
        fn pdb_residue_info_set_serial_number(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: i32,
        ) -> Result<()>;
        fn pdb_residue_info_set_alt_loc(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: &str,
        ) -> Result<()>;
        fn pdb_residue_info_set_insertion_code(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: &str,
        ) -> Result<()>;
        fn pdb_residue_info_set_occupancy(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: f64,
        ) -> Result<()>;
        fn pdb_residue_info_set_temp_factor(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: f64,
        ) -> Result<()>;
        fn pdb_residue_info_set_is_hetero_atom(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: bool,
        ) -> Result<()>;
        fn pdb_residue_info_set_secondary_structure(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: u32,
        ) -> Result<()>;
        fn pdb_residue_info_set_segment_number(
            info: Pin<&mut AtomPDBResidueInfo>,
            val: u32,
        ) -> Result<()>;

        fn atom_get_monomer_info(atom: &Atom) -> Result<UniquePtr<AtomMonomerInfo>>;
        fn atom_get_pdb_residue_info(atom: &Atom) -> Result<UniquePtr<AtomPDBResidueInfo>>;
        fn atom_set_monomer_info(atom: Pin<&mut Atom>, info: &AtomMonomerInfo) -> Result<()>;
    }
}

// SAFETY: monomer info is plain data with no thread affinity.
unsafe impl Send for ffi::AtomMonomerInfo {}
unsafe impl Send for ffi::AtomPDBResidueInfo {}
