// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::Atom::HybridizationType`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum HybridizationType {
        UNSPECIFIED = 0,
        S,
        SP,
        SP2,
        SP3,
        SP2D,
        SP3D,
        SP3D2,
        OTHER,
    }

    /// `RDKit::Atom::ChiralType`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum ChiralType {
        CHI_UNSPECIFIED = 0,
        CHI_TETRAHEDRAL_CW,
        CHI_TETRAHEDRAL_CCW,
        CHI_OTHER,
        CHI_TETRAHEDRAL,
        CHI_ALLENE,
        CHI_SQUAREPLANAR,
        CHI_TRIGONALBIPYRAMIDAL,
        CHI_OCTAHEDRAL,
    }

    /// `RDKit::Atom::ValenceType`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum ValenceType {
        IMPLICIT = 0,
        EXPLICIT,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/atom.h");

        type Atom;
        type HybridizationType;
        type ChiralType;
        type ValenceType;

        fn atom_noatom() -> Result<u32>;

        fn atom_new(what: &str) -> Result<UniquePtr<Atom>>;
        fn atom_from_num(num: u32) -> Result<UniquePtr<Atom>>;
        fn atom_clone(atom: &Atom) -> Result<UniquePtr<Atom>>;

        fn atom_get_idx(atom: &Atom) -> Result<u32>;
        fn atom_get_atomic_num(atom: &Atom) -> Result<i32>;
        fn atom_get_symbol(atom: &Atom) -> Result<UniquePtr<CxxString>>;
        fn atom_get_formal_charge(atom: &Atom) -> Result<i32>;
        fn atom_get_is_aromatic(atom: &Atom) -> Result<bool>;
        fn atom_get_isotope(atom: &Atom) -> Result<u32>;
        fn atom_get_num_explicit_hs(atom: &Atom) -> Result<u32>;
        fn atom_get_hybridization(atom: &Atom) -> Result<HybridizationType>;
        fn atom_get_chiral_tag(atom: &Atom) -> Result<ChiralType>;
        fn atom_get_atom_map_num(atom: &Atom) -> Result<i32>;
        fn atom_get_mass(atom: &Atom) -> Result<f64>;
        fn atom_get_no_implicit(atom: &Atom) -> Result<bool>;
        fn atom_get_num_radical_electrons(atom: &Atom) -> Result<u32>;
        fn atom_has_owning_mol(atom: &Atom) -> Result<bool>;
        fn atom_has_query(atom: &Atom) -> Result<bool>;
        fn atom_needs_update_property_cache(atom: &Atom) -> Result<bool>;
        fn atom_get_valence(atom: &Atom, which: ValenceType) -> Result<u32>;
        fn atom_get_smarts(
            atom: &Atom,
            do_kekule: bool,
            all_hs_explicit: bool,
            isomeric_smiles: bool,
        ) -> Result<UniquePtr<CxxString>>;

        fn atom_set_atomic_num(atom: Pin<&mut Atom>, new_num: i32) -> Result<()>;
        fn atom_set_formal_charge(atom: Pin<&mut Atom>, what: i32) -> Result<()>;
        fn atom_set_is_aromatic(atom: Pin<&mut Atom>, what: bool) -> Result<()>;
        fn atom_set_isotope(atom: Pin<&mut Atom>, what: u32) -> Result<()>;
        fn atom_set_num_explicit_hs(atom: Pin<&mut Atom>, what: u32) -> Result<()>;
        fn atom_set_no_implicit(atom: Pin<&mut Atom>, what: bool) -> Result<()>;
        fn atom_set_hybridization(atom: Pin<&mut Atom>, what: HybridizationType) -> Result<()>;
        fn atom_set_chiral_tag(atom: Pin<&mut Atom>, what: ChiralType) -> Result<()>;
        fn atom_set_atom_map_num(atom: Pin<&mut Atom>, mapno: i32, strict: bool) -> Result<()>;
        fn atom_set_num_radical_electrons(atom: Pin<&mut Atom>, num: u32) -> Result<()>;
        fn atom_clear_property_cache(atom: Pin<&mut Atom>) -> Result<()>;
        fn atom_invert_chirality(atom: Pin<&mut Atom>) -> Result<bool>;

        fn atom_get_degree(atom: &Atom) -> Result<u32>;
        fn atom_get_total_degree(atom: &Atom) -> Result<u32>;
        fn atom_get_total_num_hs(atom: &Atom, include_neighbors: bool) -> Result<u32>;
        fn atom_get_num_implicit_hs(atom: &Atom) -> Result<u32>;
        fn atom_get_total_valence(atom: &Atom) -> Result<u32>;
        fn atom_is_in_ring(atom: &Atom) -> Result<bool>;
        fn atom_is_in_ring_size(atom: &Atom, size: i32) -> Result<bool>;
        fn atom_has_valence_violation(atom: &Atom) -> Result<bool>;
        fn atom_get_alias(atom: &Atom) -> Result<UniquePtr<CxxString>>;
        fn atom_set_alias(atom: Pin<&mut Atom>, rlabel: &str) -> Result<()>;
        fn atom_get_r_label(atom: &Atom) -> Result<i32>;
        fn atom_set_r_label(atom: Pin<&mut Atom>, rlabel: i32) -> Result<()>;
        fn atom_get_value(atom: &Atom) -> Result<UniquePtr<CxxString>>;
        fn atom_set_value(atom: Pin<&mut Atom>, rlabel: &str) -> Result<()>;
        fn atom_get_supplemental_smiles_label(atom: &Atom) -> Result<UniquePtr<CxxString>>;
        fn atom_set_supplemental_smiles_label(atom: Pin<&mut Atom>, label: &str) -> Result<()>;
        fn atom_get_num_pi_electrons(atom: &Atom) -> Result<u32>;
        fn atom_get_neighbors(atom: &Atom) -> Result<Vec<u32>>;
        fn atom_get_bonds(atom: &Atom) -> Result<Vec<u32>>;
        fn atom_update_property_cache(atom: Pin<&mut Atom>, strict: bool) -> Result<()>;
    }
}

// SAFETY: an Atom outside a molecule owns all of its data. Atoms inside a
// molecule are reached only through views that borrow the molecule, which is
// not Sync (CONVENTIONS.md S1).
unsafe impl Send for ffi::Atom {}
