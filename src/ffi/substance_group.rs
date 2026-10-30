// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::SubstanceGroup::CState`.
    struct SubstanceGroupCStateData {
        bond_idx: u32,
        vector: Point3D,
    }

    /// `RDKit::SubstanceGroup::AttachPoint`.
    struct SubstanceGroupAttachData {
        a_idx: u32,
        lv_idx: i32,
        id: UniquePtr<CxxString>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/substance_group.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type SubstanceGroup;
        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;

        fn mol_get_num_substance_groups(mol: &RWMol) -> Result<u32>;
        fn mol_get_substance_group(mol: &RWMol, idx: u32) -> Result<&SubstanceGroup>;
        fn mol_get_substance_group_mut(
            mol: Pin<&mut RWMol>,
            idx: u32,
        ) -> Result<Pin<&mut SubstanceGroup>>;
        fn mol_clear_substance_groups(mol: Pin<&mut RWMol>) -> Result<()>;
        fn mol_create_substance_group(mol: Pin<&mut RWMol>, type_: &str) -> Result<u32>;
        fn mol_create_data_substance_group(
            mol: Pin<&mut RWMol>,
            field_name: &str,
            value: &str,
        ) -> Result<u32>;
        fn mol_add_substance_group(mol: Pin<&mut RWMol>, sgroup: &SubstanceGroup) -> Result<u32>;

        fn substance_group_get_index_in_mol(sgroup: &SubstanceGroup) -> Result<u32>;
        fn substance_group_get_atoms(sgroup: &SubstanceGroup) -> Result<Vec<u32>>;
        fn substance_group_get_parent_atoms(sgroup: &SubstanceGroup) -> Result<Vec<u32>>;
        fn substance_group_get_bonds(sgroup: &SubstanceGroup) -> Result<Vec<u32>>;
        fn substance_group_get_cstates(
            sgroup: &SubstanceGroup,
        ) -> Result<Vec<SubstanceGroupCStateData>>;
        fn substance_group_get_attach_points(
            sgroup: &SubstanceGroup,
        ) -> Result<Vec<SubstanceGroupAttachData>>;
        fn substance_group_get_brackets(sgroup: &SubstanceGroup) -> Result<Vec<Point3D>>;

        fn substance_group_set_atoms(sgroup: Pin<&mut SubstanceGroup>, atoms: &[u32])
        -> Result<()>;
        fn substance_group_set_parent_atoms(
            sgroup: Pin<&mut SubstanceGroup>,
            patoms: &[u32],
        ) -> Result<()>;
        fn substance_group_set_bonds(sgroup: Pin<&mut SubstanceGroup>, bonds: &[u32])
        -> Result<()>;
        fn substance_group_add_atom_with_idx(
            sgroup: Pin<&mut SubstanceGroup>,
            idx: u32,
        ) -> Result<()>;
        fn substance_group_add_parent_atom_with_idx(
            sgroup: Pin<&mut SubstanceGroup>,
            idx: u32,
        ) -> Result<()>;
        fn substance_group_add_bond_with_idx(
            sgroup: Pin<&mut SubstanceGroup>,
            idx: u32,
        ) -> Result<()>;
        fn substance_group_add_atom_with_bookmark(
            sgroup: Pin<&mut SubstanceGroup>,
            mark: i32,
        ) -> Result<()>;
        fn substance_group_add_parent_atom_with_bookmark(
            sgroup: Pin<&mut SubstanceGroup>,
            mark: i32,
        ) -> Result<()>;
        fn substance_group_add_bond_with_bookmark(
            sgroup: Pin<&mut SubstanceGroup>,
            mark: i32,
        ) -> Result<()>;
        fn substance_group_add_cstate(
            sgroup: Pin<&mut SubstanceGroup>,
            bond_idx: u32,
            vector: &Point3D,
        ) -> Result<()>;
        fn substance_group_add_attach_point(
            sgroup: Pin<&mut SubstanceGroup>,
            a_idx: u32,
            lv_idx: i32,
            id_str: &str,
        ) -> Result<()>;
        fn substance_group_add_bracket(
            sgroup: Pin<&mut SubstanceGroup>,
            pts: &[Point3D],
        ) -> Result<()>;
        fn substance_group_clear_brackets(sgroup: Pin<&mut SubstanceGroup>) -> Result<()>;
        fn substance_group_clear_cstates(sgroup: Pin<&mut SubstanceGroup>) -> Result<()>;
        fn substance_group_clear_attach_points(sgroup: Pin<&mut SubstanceGroup>) -> Result<()>;

        fn substance_group_cstate_default() -> Result<SubstanceGroupCStateData>;
        fn substance_group_attach_default() -> Result<SubstanceGroupAttachData>;
    }
}
