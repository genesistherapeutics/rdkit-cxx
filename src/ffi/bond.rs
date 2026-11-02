// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::Bond::BondType`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum BondType {
        UNSPECIFIED = 0,
        SINGLE,
        DOUBLE,
        TRIPLE,
        QUADRUPLE,
        QUINTUPLE,
        HEXTUPLE,
        ONEANDAHALF,
        TWOANDAHALF,
        THREEANDAHALF,
        FOURANDAHALF,
        FIVEANDAHALF,
        AROMATIC,
        IONIC,
        HYDROGEN,
        THREECENTER,
        DATIVEONE,
        DATIVE,
        DATIVEL,
        DATIVER,
        OTHER,
        ZERO,
    }

    /// `RDKit::Bond::BondDir`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum BondDir {
        NONE = 0,
        BEGINWEDGE,
        BEGINDASH,
        ENDDOWNRIGHT,
        ENDUPRIGHT,
        EITHERDOUBLE,
        UNKNOWN,
    }

    /// `RDKit::Bond::BondStereo`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum BondStereo {
        STEREONONE = 0,
        STEREOANY,
        STEREOZ,
        STEREOE,
        STEREOCIS,
        STEREOTRANS,
        STEREOATROPCW,
        STEREOATROPCCW,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/atom.h");
        include!("rdkit-cxx/bond.h");

        type Bond;
        type BondType;
        type BondDir;
        type BondStereo;
        type Atom = crate::ffi::atom::ffi::Atom;

        fn bond_get_idx(bond: &Bond) -> Result<u32>;
        fn bond_get_bond_type(bond: &Bond) -> Result<BondType>;
        fn bond_get_bond_type_as_double(bond: &Bond) -> Result<f64>;
        fn bond_get_begin_atom_idx(bond: &Bond) -> Result<u32>;
        fn bond_get_end_atom_idx(bond: &Bond) -> Result<u32>;
        fn bond_get_other_atom_idx(bond: &Bond, this_idx: u32) -> Result<u32>;
        fn bond_get_is_aromatic(bond: &Bond) -> Result<bool>;
        fn bond_get_is_conjugated(bond: &Bond) -> Result<bool>;
        fn bond_get_bond_dir(bond: &Bond) -> Result<BondDir>;
        fn bond_get_stereo(bond: &Bond) -> Result<BondStereo>;
        fn bond_get_stereo_atoms(bond: &Bond) -> Result<Vec<i32>>;
        fn bond_has_owning_mol(bond: &Bond) -> Result<bool>;
        fn bond_has_query(bond: &Bond) -> Result<bool>;
        fn bond_is_in_ring(bond: &Bond) -> Result<bool>;
        fn bond_is_in_ring_size(bond: &Bond, size: i32) -> Result<bool>;
        fn bond_get_valence_contrib(bond: &Bond, at: &Atom) -> Result<f64>;
        fn bond_get_smarts(bond: &Bond, all_bonds_explicit: bool) -> Result<UniquePtr<CxxString>>;

        fn bond_set_bond_type(bond: Pin<&mut Bond>, b_t: BondType) -> Result<()>;
        fn bond_set_is_aromatic(bond: Pin<&mut Bond>, what: bool) -> Result<()>;
        fn bond_set_is_conjugated(bond: Pin<&mut Bond>, what: bool) -> Result<()>;
        fn bond_set_bond_dir(bond: Pin<&mut Bond>, what: BondDir) -> Result<()>;
        fn bond_set_stereo(bond: Pin<&mut Bond>, what: BondStereo) -> Result<()>;
        fn bond_set_stereo_atoms(bond: Pin<&mut Bond>, bgn_idx: u32, end_idx: u32) -> Result<()>;
        fn bond_invert_chirality(bond: Pin<&mut Bond>) -> Result<bool>;
        fn bond_clone(bond: &Bond) -> Result<UniquePtr<Bond>>;
    }
}

// SAFETY: as for Atom (atom.rs): a Bond outside a molecule owns all of its
// data (its stereo atoms, macro info, query, and properties; RDKit's only
// static, QueryBond.cpp's table of equality queries, is initialized once and
// only read). Bonds inside a molecule are reached only through views that
// borrow the molecule, which is not Sync (CONVENTIONS.md S1).
unsafe impl Send for ffi::Bond {}
