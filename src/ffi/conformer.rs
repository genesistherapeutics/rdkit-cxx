// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/conformer.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type Conformer;
        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;

        fn conformer_new() -> Result<UniquePtr<Conformer>>;
        fn conformer_with_num_atoms(num_atoms: u32) -> Result<UniquePtr<Conformer>>;
        fn conformer_clone(conf: &Conformer) -> Result<UniquePtr<Conformer>>;

        fn conformer_get_id(conf: &Conformer) -> Result<u32>;
        fn conformer_get_num_atoms(conf: &Conformer) -> Result<u32>;
        fn conformer_is_3d(conf: &Conformer) -> Result<bool>;
        fn conformer_has_owning_mol(conf: &Conformer) -> Result<bool>;
        fn conformer_get_atom_position(conf: &Conformer, aid: u32) -> Result<Point3D>;
        fn conformer_get_positions(conf: &Conformer) -> Result<Vec<Point3D>>;
        fn conformer_set_atom_position(
            conf: Pin<&mut Conformer>,
            aid: u32,
            loc: &Point3D,
        ) -> Result<()>;
        fn conformer_set_positions(conf: Pin<&mut Conformer>, positions: &[Point3D]) -> Result<()>;
        fn conformer_set_3d(conf: Pin<&mut Conformer>, v: bool) -> Result<()>;
        fn conformer_set_id(conf: Pin<&mut Conformer>, id: u32) -> Result<()>;

        fn mol_get_num_conformers(mol: &RWMol) -> Result<u32>;
        fn mol_get_conformer(mol: &RWMol, id: i32) -> Result<&Conformer>;
        fn mol_find_conformer(mol: &RWMol, id: u32) -> Result<&Conformer>;
        fn mol_find_conformer_mut(mol: Pin<&mut RWMol>, id: u32) -> Result<Pin<&mut Conformer>>;
        fn mol_set_conformer_id(mol: Pin<&mut RWMol>, current: u32, id: u32) -> Result<()>;
        fn mol_get_conformer_at(mol: &RWMol, position: u32) -> Result<&Conformer>;
        fn mol_add_conformer(
            mol: Pin<&mut RWMol>,
            conf: &Conformer,
            assign_id: bool,
        ) -> Result<u32>;
        fn mol_remove_conformer(mol: Pin<&mut RWMol>, id: u32) -> Result<()>;
        fn mol_remove_all_conformers(mol: Pin<&mut RWMol>) -> Result<()>;
    }
}

// SAFETY: a Conformer outside a molecule owns all of its data. Conformers in a
// molecule are reached only through views that borrow the molecule, which is
// not Sync (CONVENTIONS.md S1).
unsafe impl Send for ffi::Conformer {}
