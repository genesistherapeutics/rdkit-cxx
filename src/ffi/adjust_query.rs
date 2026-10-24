// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolOps::AdjustQueryWhichFlags`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum AdjustQueryWhichFlags {
        ADJUST_IGNORENONE = 0x0,
        ADJUST_IGNORECHAINS = 0x1,
        ADJUST_IGNOREDUMMIES = 0x2,
        ADJUST_IGNORERINGS = 0x4,
        ADJUST_IGNORENONDUMMIES = 0x8,
        ADJUST_IGNOREMAPPED = 0x10,
        ADJUST_IGNOREALL = 0xFFFFFFF,
    }

    /// `MolOps::AdjustQueryParameters`.
    struct AdjustQueryParametersData {
        adjust_degree: bool,
        adjust_degree_flags: u32,
        adjust_heavy_degree: bool,
        adjust_heavy_degree_flags: u32,
        adjust_ring_count: bool,
        adjust_ring_count_flags: u32,
        make_dummies_queries: bool,
        aromatize_if_possible: bool,
        make_bonds_generic: bool,
        make_bonds_generic_flags: u32,
        make_atoms_generic: bool,
        make_atoms_generic_flags: u32,
        adjust_ring_chain: bool,
        adjust_ring_chain_flags: u32,
        use_stereo_care_for_bonds: bool,
        adjust_conjugated_five_rings: bool,
        set_mdl_five_ring_aromaticity: bool,
        adjust_single_bonds_to_degree_one_neighbors: bool,
        adjust_single_bonds_between_aromatic_atoms: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/adjust_query.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type AdjustQueryWhichFlags;

        fn adjust_query_parameters_default() -> Result<AdjustQueryParametersData>;
        fn adjust_query_parameters_no_adjustments() -> Result<AdjustQueryParametersData>;
        fn adjust_query_properties(
            mol: &RWMol,
            params: &AdjustQueryParametersData,
        ) -> Result<UniquePtr<RWMol>>;
        fn adjust_query_properties_with_generic_groups(
            mol: &RWMol,
            params: &AdjustQueryParametersData,
        ) -> Result<UniquePtr<RWMol>>;
        fn set_generic_queries_from_properties(
            mol: Pin<&mut RWMol>,
            use_atom_labels: bool,
            use_s_groups: bool,
        ) -> Result<()>;
        fn convert_generic_queries_to_substance_groups(mol: Pin<&mut RWMol>) -> Result<()>;
    }
}
