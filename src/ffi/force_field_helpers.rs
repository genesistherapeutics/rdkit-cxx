// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// One conformer's result: 0 when it converged, 1 when it needs more
    /// iterations, -1 when MMFF could not type the molecule.
    #[derive(Clone, Copy, Default)]
    struct ConfOptimization {
        not_converged: i32,
        energy: f64,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/force_field_helpers.h");
        include!("rdkit-cxx/ffi/force_field.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type ForceField = crate::ffi::force_field::ffi::ForceField;
        type MMFFMolProperties = crate::ffi::force_field::ffi::MMFFMolProperties;
        type MMFFVerbosity = crate::ffi::force_field::ffi::MMFFVerbosity;
        type TermParams = crate::ffi::force_field::ffi::TermParams;

        fn uff_optimize_molecule(
            mol: Pin<&mut RWMol>,
            max_iters: u32,
            vdw_thresh: f64,
            conf_id: i32,
            ignore_interfrag_interactions: bool,
        ) -> Result<i32>;
        fn uff_optimize_molecule_confs(
            mol: Pin<&mut RWMol>,
            num_threads: i32,
            max_iters: u32,
            vdw_thresh: f64,
            ignore_interfrag_interactions: bool,
        ) -> Result<Vec<ConfOptimization>>;
        fn uff_get_molecule_force_field(
            mol: Pin<&mut RWMol>,
            vdw_thresh: f64,
            conf_id: i32,
            ignore_interfrag_interactions: bool,
        ) -> Result<UniquePtr<ForceField>>;
        fn uff_has_all_molecule_params(mol: &RWMol) -> Result<bool>;
        fn mmff_optimize_molecule(
            mol: Pin<&mut RWMol>,
            mmff_variant: &str,
            max_iters: u32,
            non_bonded_thresh: f64,
            conf_id: i32,
            ignore_interfrag_interactions: bool,
        ) -> Result<i32>;
        fn mmff_sanitize_molecule(mol: Pin<&mut RWMol>) -> Result<u32>;
        fn mmff_get_molecule_properties(
            mol: Pin<&mut RWMol>,
            mmff_variant: &str,
            verbosity: MMFFVerbosity,
        ) -> Result<UniquePtr<MMFFMolProperties>>;
        fn mmff_get_molecule_force_field(
            mol: Pin<&mut RWMol>,
            props: &MMFFMolProperties,
            non_bonded_thresh: f64,
            conf_id: i32,
            ignore_interfrag_interactions: bool,
        ) -> Result<UniquePtr<ForceField>>;
        fn create_empty_force_field_for_mol(
            mol: Pin<&mut RWMol>,
            conf_id: i32,
        ) -> Result<UniquePtr<ForceField>>;
        fn mmff_has_all_molecule_params(mol: &RWMol) -> Result<bool>;
        fn mmff_optimize_molecule_confs(
            mol: Pin<&mut RWMol>,
            num_threads: i32,
            max_iters: u32,
            mmff_variant: &str,
            non_bonded_thresh: f64,
            ignore_interfrag_interactions: bool,
        ) -> Result<Vec<ConfOptimization>>;
        fn get_uff_bond_stretch_params(mol: &RWMol, idx1: u32, idx2: u32) -> Result<TermParams>;
        fn get_uff_angle_bend_params(
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
        ) -> Result<TermParams>;
        fn get_uff_torsion_params(
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
            idx4: u32,
        ) -> Result<TermParams>;
        fn get_uff_inversion_params(
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
            idx4: u32,
        ) -> Result<TermParams>;
        fn get_uff_vdw_params(mol: &RWMol, idx1: u32, idx2: u32) -> Result<TermParams>;
        fn optimize_molecule(ff: Pin<&mut ForceField>, max_iters: u32) -> Result<i32>;
        fn optimize_molecule_confs(
            mol: Pin<&mut RWMol>,
            ff: Pin<&mut ForceField>,
            num_threads: i32,
            max_iters: u32,
        ) -> Result<Vec<ConfOptimization>>;
    }
}
