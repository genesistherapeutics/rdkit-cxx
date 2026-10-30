// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolOps::AromaticityModel`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum AromaticityModel {
        AROMATICITY_DEFAULT = 0x0,
        AROMATICITY_RDKIT = 0x1,
        AROMATICITY_SIMPLE = 0x2,
        AROMATICITY_MDL = 0x4,
        AROMATICITY_CUSTOM = 0xFFFFFFF,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/perception.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type AromaticityModel;

        fn get_sssr(mol: Pin<&mut RWMol>, include_dative_bonds: bool) -> Result<Vec<IndexList>>;
        fn get_symm_sssr(
            mol: Pin<&mut RWMol>,
            include_dative_bonds: bool,
        ) -> Result<Vec<IndexList>>;
        fn fast_find_rings(mol: Pin<&mut RWMol>) -> Result<()>;
        fn find_ring_families(mol: Pin<&mut RWMol>) -> Result<()>;

        fn cleanup(mol: Pin<&mut RWMol>) -> Result<()>;
        fn cleanup_chirality(mol: Pin<&mut RWMol>) -> Result<()>;
        fn cleanup_atropisomers(mol: Pin<&mut RWMol>) -> Result<()>;
        fn cleanup_organometallics(mol: Pin<&mut RWMol>) -> Result<()>;
        fn set_aromaticity(mol: Pin<&mut RWMol>, model: AromaticityModel) -> Result<()>;
        fn set_conjugation(mol: Pin<&mut RWMol>) -> Result<()>;
        fn set_hybridization(mol: Pin<&mut RWMol>) -> Result<()>;
        fn assign_radicals(mol: Pin<&mut RWMol>) -> Result<()>;
        fn kekulize_if_possible(mol: Pin<&mut RWMol>, clear_aromatic_flags: bool) -> Result<()>;
        fn haptic_bonds_to_dative(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn dative_bonds_to_haptic(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn needs_hs(mol: &RWMol) -> Result<bool>;
        fn count_atom_elec(atom: &Atom) -> Result<i32>;
        fn atom_has_conjugated_bond(atom: &Atom) -> Result<bool>;
        fn get_formal_charge(mol: &RWMol) -> Result<i32>;
        fn get_shortest_path(mol: &RWMol, aid1: u32, aid2: u32) -> Result<Vec<u32>>;
        fn set_terminal_atom_coords(mol: Pin<&mut RWMol>, idx: u32, other_idx: u32) -> Result<()>;
    }
}
