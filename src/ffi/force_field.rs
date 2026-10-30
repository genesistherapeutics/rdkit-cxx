// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MMFF`'s `MMFF_VERBOSITY` values, which RDKit 2025.09 declares
    /// in an unnamed enum and takes as `std::uint8_t`.
    #[repr(u8)]
    #[allow(clippy::enum_variant_names)]
    enum MMFFVerbosity {
        MMFF_VERBOSITY_NONE = 0,
        MMFF_VERBOSITY_LOW = 1,
        MMFF_VERBOSITY_HIGH = 2,
    }

    /// The MMFF terms a force field can include.
    #[repr(u8)]
    enum MMFFTerm {
        Bond,
        Angle,
        StretchBend,
        Oop,
        Torsion,
        VdW,
        Ele,
    }

    /// A parameter getter's result: `found`, then the term type and the
    /// values in Python's tuple order.
    struct TermParams {
        found: bool,
        term_type: u32,
        values: Vec<f64>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/force_field.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type ForceField;
        type MMFFMolProperties;

        fn force_field_calc_energy(ff: &ForceField) -> Result<f64>;
        fn force_field_calc_energy_with_pos(ff: Pin<&mut ForceField>, pos: &[f64]) -> Result<f64>;
        fn force_field_calc_grad(ff: &ForceField) -> Result<Vec<f64>>;
        fn force_field_calc_grad_with_pos(
            ff: Pin<&mut ForceField>,
            pos: &[f64],
        ) -> Result<Vec<f64>>;
        fn force_field_positions(ff: &ForceField) -> Result<Vec<f64>>;
        fn force_field_dimension(ff: &ForceField) -> Result<u32>;
        fn force_field_num_points(ff: &ForceField) -> Result<u32>;
        fn force_field_minimize(
            ff: Pin<&mut ForceField>,
            max_its: u32,
            force_tol: f64,
            energy_tol: f64,
        ) -> Result<i32>;
        fn force_field_minimize_trajectory(
            ff: Pin<&mut ForceField>,
            snapshot_freq: u32,
            max_its: u32,
            force_tol: f64,
            energy_tol: f64,
            energies: &mut Vec<f64>,
            positions: &mut Vec<f64>,
        ) -> Result<i32>;
        fn force_field_add_distance_constraint(
            ff: Pin<&mut ForceField>,
            idx1: u32,
            idx2: u32,
            min_len: f64,
            max_len: f64,
            force_constant: f64,
        ) -> Result<()>;
        fn force_field_add_fixed_point(ff: Pin<&mut ForceField>, idx: u32) -> Result<()>;
        fn force_field_add_relative_distance_constraint(
            ff: Pin<&mut ForceField>,
            idx1: u32,
            idx2: u32,
            relative: bool,
            min_len: f64,
            max_len: f64,
            force_constant: f64,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn force_field_add_angle_constraint(
            ff: Pin<&mut ForceField>,
            idx1: u32,
            idx2: u32,
            idx3: u32,
            relative: bool,
            min_angle_deg: f64,
            max_angle_deg: f64,
            force_constant: f64,
        ) -> Result<()>;
        #[allow(clippy::too_many_arguments)]
        fn force_field_add_torsion_constraint(
            ff: Pin<&mut ForceField>,
            idx1: u32,
            idx2: u32,
            idx3: u32,
            idx4: u32,
            relative: bool,
            min_dihedral_deg: f64,
            max_dihedral_deg: f64,
            force_constant: f64,
        ) -> Result<()>;
        fn force_field_add_position_constraint(
            ff: Pin<&mut ForceField>,
            idx: u32,
            max_displ: f64,
            force_constant: f64,
        ) -> Result<()>;
        fn force_field_initialize(ff: Pin<&mut ForceField>) -> Result<()>;
        fn force_field_add_extra_point(
            ff: Pin<&mut ForceField>,
            x: f64,
            y: f64,
            z: f64,
            fixed: bool,
        ) -> Result<u32>;
        fn force_field_get_extra_point_pos(ff: &ForceField, idx: u32) -> Result<Point3D>;

        fn mmff_mol_properties_new(
            mol: Pin<&mut RWMol>,
            variant: &str,
            verbosity: MMFFVerbosity,
        ) -> Result<UniquePtr<MMFFMolProperties>>;
        fn mmff_mol_properties_is_valid(props: &MMFFMolProperties) -> Result<bool>;
        fn mmff_mol_properties_get_atom_type(props: &MMFFMolProperties, idx: u32) -> Result<u32>;
        fn mmff_mol_properties_get_formal_charge(
            props: &MMFFMolProperties,
            idx: u32,
        ) -> Result<f64>;
        fn mmff_mol_properties_get_partial_charge(
            props: &MMFFMolProperties,
            idx: u32,
        ) -> Result<f64>;
        fn mmff_mol_properties_get_bond_stretch_params(
            props: &MMFFMolProperties,
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
        ) -> Result<TermParams>;
        fn mmff_mol_properties_get_angle_bend_params(
            props: &MMFFMolProperties,
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
        ) -> Result<TermParams>;
        fn mmff_mol_properties_get_stretch_bend_params(
            props: &MMFFMolProperties,
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
        ) -> Result<TermParams>;
        fn mmff_mol_properties_get_torsion_params(
            props: &MMFFMolProperties,
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
            idx4: u32,
        ) -> Result<TermParams>;
        fn mmff_mol_properties_get_oop_bend_params(
            props: &MMFFMolProperties,
            mol: &RWMol,
            idx1: u32,
            idx2: u32,
            idx3: u32,
            idx4: u32,
        ) -> Result<TermParams>;
        fn mmff_mol_properties_get_vdw_params(
            props: &MMFFMolProperties,
            idx1: u32,
            idx2: u32,
        ) -> Result<TermParams>;
        fn mmff_mol_properties_set_dielectric_model(
            props: Pin<&mut MMFFMolProperties>,
            model: u8,
        ) -> Result<()>;
        fn mmff_mol_properties_set_dielectric_constant(
            props: Pin<&mut MMFFMolProperties>,
            constant: f64,
        ) -> Result<()>;
        fn mmff_mol_properties_set_term(
            props: Pin<&mut MMFFMolProperties>,
            term: MMFFTerm,
            state: bool,
        ) -> Result<()>;
        fn mmff_mol_properties_set_variant(
            props: Pin<&mut MMFFMolProperties>,
            variant: &str,
        ) -> Result<()>;
        fn mmff_mol_properties_set_verbosity(
            props: Pin<&mut MMFFMolProperties>,
            verbosity: MMFFVerbosity,
        ) -> Result<()>;
    }

    // The force field helpers return these.
    impl UniquePtr<ForceField> {}
}

// SAFETY: a ForceField owns its contributions and extra points; its other
// points are in a molecule's conformer, which the Rust ForceField borrows
// mutably for its lifetime, and a molecule is Send (CONVENTIONS.md S1).
unsafe impl Send for ffi::ForceField {}

// SAFETY: MMFFMolProperties owns its data; its output stream is std::cout.
unsafe impl Send for ffi::MMFFMolProperties {}
