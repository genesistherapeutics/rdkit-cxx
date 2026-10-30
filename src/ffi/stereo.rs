// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::StereoGroupType`; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum StereoGroupType {
        STEREO_ABSOLUTE = 0,
        STEREO_OR = 1,
        STEREO_AND = 2,
    }

    /// `RDKit::Chirality::StereoType`; cxx checks the values.
    #[repr(i32)]
    enum StereoType {
        Unspecified = 0,
        Atom_Tetrahedral,
        Atom_SquarePlanar,
        Atom_TrigonalBipyramidal,
        Atom_Octahedral,
        Bond_Double,
        Bond_Cumulene_Even,
        Bond_Atropisomer,
    }

    /// `RDKit::Chirality::StereoSpecified`; cxx checks the values.
    #[repr(i32)]
    enum StereoSpecified {
        Unspecified = 0,
        Specified,
        Unknown,
    }

    /// `RDKit::Chirality::StereoDescriptor`; cxx checks the values.
    #[repr(i32)]
    enum StereoDescriptor {
        None = 0,
        Tet_CW,
        Tet_CCW,
        Bond_Cis,
        Bond_Trans,
        Bond_AtropCW,
        Bond_AtropCCW,
    }

    /// A stereo group by atom and bond index.
    struct StereoGroupData {
        group_type: StereoGroupType,
        atoms: Vec<u32>,
        bonds: Vec<u32>,
        read_id: u32,
        write_id: u32,
    }

    /// `RDKit::Chirality::StereoInfo`.
    struct StereoInfoData {
        type_: StereoType,
        specified: StereoSpecified,
        centered_on: u32,
        descriptor: StereoDescriptor,
        permutation: u32,
        controlling_atoms: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/stereo.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type StereoGroupType;
        type StereoType;
        type StereoSpecified;
        type StereoDescriptor;

        fn mol_get_stereo_groups(mol: &RWMol) -> Result<Vec<StereoGroupData>>;
        fn rwmol_set_stereo_groups(mol: Pin<&mut RWMol>, groups: &[StereoGroupData]) -> Result<()>;
        fn create_stereo_group(
            group_type: StereoGroupType,
            mol: &RWMol,
            atom_ids: &[u32],
            read_id: u32,
        ) -> Result<StereoGroupData>;
        fn forward_stereo_group_ids(mol: Pin<&mut RWMol>) -> Result<()>;
        fn find_potential_stereo(
            mol: Pin<&mut RWMol>,
            clean_it: bool,
            flag_possible: bool,
        ) -> Result<Vec<StereoInfoData>>;
    }
}
