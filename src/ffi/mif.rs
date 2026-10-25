// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mif.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type UniformRealValueGrid3D = crate::ffi::grid::ffi::UniformRealValueGrid3D;
        type Coulomb;
        type CoulombDielectric;
        type MMFFVdWaals;
        type UFFVdWaals;
        type HBond;
        type Hydrophilic;

        fn construct_grid(
            mol: &RWMol,
            conf_id: i32,
            margin: f64,
            spacing: f64,
        ) -> Result<UniquePtr<UniformRealValueGrid3D>>;

        #[allow(clippy::too_many_arguments)]
        fn coulomb_new(
            mol: &RWMol,
            conf_id: i32,
            probe_charge: f64,
            abs_val: bool,
            charge_key: &str,
            softcore_param: f64,
            cutoff: f64,
        ) -> Result<UniquePtr<Coulomb>>;
        fn coulomb_from_charges(
            charges: &[f64],
            positions: &[Point3D],
            probe_charge: f64,
            abs_val: bool,
            softcore_param: f64,
            cutoff: f64,
        ) -> Result<UniquePtr<Coulomb>>;
        fn coulomb_call(functor: &Coulomb, x: f64, y: f64, z: f64, threshold: f64) -> Result<f64>;

        #[allow(clippy::too_many_arguments)]
        fn coulomb_dielectric_new(
            mol: &RWMol,
            conf_id: i32,
            probe_charge: f64,
            abs_val: bool,
            charge_key: &str,
            softcore_param: f64,
            cutoff: f64,
            epsilon: f64,
            xi: f64,
        ) -> Result<UniquePtr<CoulombDielectric>>;
        #[allow(clippy::too_many_arguments)]
        fn coulomb_dielectric_from_charges(
            charges: &[f64],
            positions: &[Point3D],
            probe_charge: f64,
            abs_val: bool,
            softcore_param: f64,
            cutoff: f64,
            epsilon: f64,
            xi: f64,
        ) -> Result<UniquePtr<CoulombDielectric>>;
        fn coulomb_dielectric_call(
            functor: &CoulombDielectric,
            x: f64,
            y: f64,
            z: f64,
            threshold: f64,
        ) -> Result<f64>;

        fn mmff_vd_waals_new(
            mol: &RWMol,
            conf_id: i32,
            probe_atom_type: u32,
            scaling: bool,
            cutoff: f64,
        ) -> Result<UniquePtr<MMFFVdWaals>>;
        fn mmff_vd_waals_call(
            functor: &MMFFVdWaals,
            x: f64,
            y: f64,
            z: f64,
            threshold: f64,
        ) -> Result<f64>;

        fn uff_vd_waals_new(
            mol: &RWMol,
            conf_id: i32,
            probe_atom_type: &str,
            cutoff: f64,
        ) -> Result<UniquePtr<UFFVdWaals>>;
        fn uff_vd_waals_call(
            functor: &UFFVdWaals,
            x: f64,
            y: f64,
            z: f64,
            threshold: f64,
        ) -> Result<f64>;

        fn h_bond_new(
            mol: &RWMol,
            conf_id: i32,
            probe_atom_type: &str,
            fixed: bool,
            cutoff: f64,
        ) -> Result<UniquePtr<HBond>>;
        fn h_bond_call(functor: &HBond, x: f64, y: f64, z: f64, threshold: f64) -> Result<f64>;

        fn hydrophilic_new(
            mol: &RWMol,
            conf_id: i32,
            fixed: bool,
            cutoff: f64,
        ) -> Result<UniquePtr<Hydrophilic>>;
        fn hydrophilic_call(
            functor: &Hydrophilic,
            x: f64,
            y: f64,
            z: f64,
            threshold: f64,
        ) -> Result<f64>;

        fn calculate_descriptors_coulomb(
            grid: Pin<&mut UniformRealValueGrid3D>,
            descriptor: &Coulomb,
            threshold: f64,
        ) -> Result<()>;
        fn calculate_descriptors_coulomb_dielectric(
            grid: Pin<&mut UniformRealValueGrid3D>,
            descriptor: &CoulombDielectric,
            threshold: f64,
        ) -> Result<()>;
        fn calculate_descriptors_mmff_vd_waals(
            grid: Pin<&mut UniformRealValueGrid3D>,
            descriptor: &MMFFVdWaals,
            threshold: f64,
        ) -> Result<()>;
        fn calculate_descriptors_uff_vd_waals(
            grid: Pin<&mut UniformRealValueGrid3D>,
            descriptor: &UFFVdWaals,
            threshold: f64,
        ) -> Result<()>;
        fn calculate_descriptors_h_bond(
            grid: Pin<&mut UniformRealValueGrid3D>,
            descriptor: &HBond,
            threshold: f64,
        ) -> Result<()>;
        fn calculate_descriptors_hydrophilic(
            grid: Pin<&mut UniformRealValueGrid3D>,
            descriptor: &Hydrophilic,
            threshold: f64,
        ) -> Result<()>;

        /// # Safety
        ///
        /// The molecule is null or valid.
        unsafe fn write_to_cube_file(
            grid: &UniformRealValueGrid3D,
            filename: &[u8],
            mol: *const RWMol,
            conf_id: i32,
        ) -> Result<()>;
        fn read_from_cube_file(
            filename: &[u8],
            grid: Pin<&mut UniformRealValueGrid3D>,
        ) -> Result<UniquePtr<RWMol>>;
    }
}

// SAFETY: each functor owns its data (MMFFVdWaals its own copy of the
// molecule), and the force-field parameters they point to are RDKit's
// immutable globals, so they may move between threads.
unsafe impl Send for ffi::Coulomb {}
unsafe impl Send for ffi::CoulombDielectric {}
unsafe impl Send for ffi::MMFFVdWaals {}
unsafe impl Send for ffi::UFFVdWaals {}
unsafe impl Send for ffi::HBond {}
unsafe impl Send for ffi::Hydrophilic {}

// SAFETY: these functors' const call only reads them. CoulombDielectric,
// HBond, and Hydrophilic are not Sync: their const operator() writes the
// mutable vectors they compute in (d_dists, d_eneContrib,
// d_vectTargetProbe).
unsafe impl Sync for ffi::Coulomb {}
unsafe impl Sync for ffi::MMFFVdWaals {}
unsafe impl Sync for ffi::UFFVdWaals {}
