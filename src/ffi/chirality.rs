// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::StereoGroupAbsOptions`; cxx checks the values.
    #[repr(i32)]
    enum StereoGroupAbsOptions {
        OnlyIncludeWhenOtherGroupsExist,
        NeverInclude,
        AlwaysInclude,
    }

    /// A pair of atom indices.
    struct IndexPair {
        first: u32,
        second: u32,
    }

    /// `Chirality::BondWedgingParameters`.
    struct BondWedgingParametersData {
        wedge_two_bonds_if_possible: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/chirality.h");
        include!("rdkit-cxx/ffi/stereo.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Bond = crate::ffi::bond::ffi::Bond;
        type Conformer = crate::ffi::conformer::ffi::Conformer;
        type StereoGroupType = crate::ffi::stereo::ffi::StereoGroupType;
        type StereoGroupAbsOptions;

        fn assign_stereochemistry(
            mol: Pin<&mut RWMol>,
            clean_it: bool,
            force: bool,
            flag_possible_stereo_centers: bool,
        ) -> Result<()>;
        fn assign_stereochemistry_from_3d(
            mol: Pin<&mut RWMol>,
            conf_id: i32,
            replace_existing_tags: bool,
        ) -> Result<()>;
        fn assign_chiral_types_from_bond_dirs(
            mol: Pin<&mut RWMol>,
            conf_id: i32,
            replace_existing_tags: bool,
        ) -> Result<()>;
        fn assign_atom_chiral_tags_from_structure(
            mol: Pin<&mut RWMol>,
            conf_id: i32,
            replace_existing_tags: bool,
        ) -> Result<()>;
        fn assign_atom_chiral_tags_from_mol_parity(
            mol: Pin<&mut RWMol>,
            replace_existing_tags: bool,
        ) -> Result<()>;
        fn find_potential_stereo_bonds(mol: Pin<&mut RWMol>, clean_it: bool) -> Result<()>;
        fn remove_stereochemistry(mol: Pin<&mut RWMol>) -> Result<()>;
        fn remove_non_explicit_3d_chirality(mol: Pin<&mut RWMol>) -> Result<()>;
        fn detect_bond_stereochemistry(mol: Pin<&mut RWMol>, conf_id: i32) -> Result<()>;
        fn set_double_bond_neighbor_directions(mol: Pin<&mut RWMol>) -> Result<()>;
        fn set_double_bond_neighbor_directions_with_conf(
            mol: Pin<&mut RWMol>,
            conf: &Conformer,
        ) -> Result<()>;
        fn set_bond_stereo_from_directions(mol: Pin<&mut RWMol>) -> Result<()>;
        fn compute_atom_cip_ranks(mol: Pin<&mut RWMol>) -> Result<Vec<u32>>;
        fn find_meso_centers(
            mol: &RWMol,
            include_isotopes: bool,
            include_atom_maps: bool,
        ) -> Result<Vec<IndexPair>>;

        fn canonicalize_stereo_groups(
            mol: &RWMol,
            output_absolute_groups: StereoGroupAbsOptions,
            max_stereo_groups: u32,
        ) -> Result<UniquePtr<RWMol>>;
        fn cleanup_stereo_groups(mol: Pin<&mut RWMol>) -> Result<()>;
        fn simplify_enhanced_stereo(
            mol: Pin<&mut RWMol>,
            remove_affected_stereo_groups: bool,
        ) -> Result<()>;
        fn translate_chiral_flag_to_stereo_groups(
            mol: Pin<&mut RWMol>,
            zero_flag_group_type: StereoGroupType,
        ) -> Result<()>;
        fn add_stereo_annotations(
            mol: Pin<&mut RWMol>,
            abs_label: &str,
            or_label: &str,
            and_label: &str,
            cip_label: &str,
            bond_label: &str,
        ) -> Result<()>;

        fn add_wavy_bonds_for_stereo_any(
            mol: Pin<&mut RWMol>,
            clear_double_bond_flags: bool,
            add_when_impossible: u32,
        ) -> Result<()>;
        fn wedge_bond(bond: Pin<&mut Bond>, from_atom_idx: u32, conf_id: i32) -> Result<()>;
        fn bond_wedging_parameters_default() -> Result<BondWedgingParametersData>;
        fn wedge_mol_bonds(
            mol: Pin<&mut RWMol>,
            conformer: &Conformer,
            params: &BondWedgingParametersData,
        ) -> Result<()>;
        fn reapply_mol_block_wedging(
            mol: Pin<&mut RWMol>,
            all_bond_types: bool,
            verify: bool,
        ) -> Result<()>;

        fn set_allow_nontetrahedral_chirality(val: bool) -> Result<()>;
        fn get_allow_nontetrahedral_chirality() -> Result<bool>;
        fn set_use_legacy_stereo_perception(val: bool) -> Result<()>;
        fn get_use_legacy_stereo_perception() -> Result<bool>;

        fn expand_attachment_points(
            mol: Pin<&mut RWMol>,
            add_as_queries: bool,
            add_coords: bool,
        ) -> Result<()>;
        fn collapse_attachment_points(mol: Pin<&mut RWMol>, marked_only: bool) -> Result<()>;
        fn detect_chemistry_problems(
            mol: &RWMol,
            sanitize_ops: u32,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
    }
}
