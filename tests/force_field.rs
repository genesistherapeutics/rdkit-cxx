// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{assert_close, assert_close_to, ethanol, ethanol_two_confs};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_force_field_helpers::*;
use rdkit_cxx::chem::rdmolops::SanitizeFlags;
use rdkit_cxx::chem::rdtrajectory::Snapshot;
use rdkit_cxx::chem::*;
use rdkit_cxx::force_field::*;
use rdkit_cxx::geometry::Point3D;

// Expected values below come from Python's rdForceField and
// rdForceFieldHelpers on the same molecules. Energies and gradients at given
// coordinates agree to the last bits between platforms; a minimization's
// steps amplify those bits, so its results are compared loosely.

fn assert_point(actual: Point3D, expected: (f64, f64, f64)) {
    assert_close(actual.x, expected.0);
    assert_close(actual.y, expected.1);
    assert_close(actual.z, expected.2);
}

/// Asserts an energy a minimization reached, relative to the larger.
fn assert_minimized(actual: f64, expected: f64) {
    assert_close_to(actual, expected, 1e-4);
}

/// Asserts coordinates a minimization reached, to 0.005 Å.
fn assert_minimized_point(actual: Point3D, expected: (f64, f64, f64)) {
    let expected = Point3D::new(expected.0, expected.1, expected.2);
    assert!(
        actual.distance(&expected) < 5e-3,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn params_defaults() {
    assert_eq!(
        MinimizeParams::default(),
        MinimizeParams {
            max_its: 200,
            force_tol: 1e-4,
            energy_tol: 1e-6,
        }
    );
    assert_eq!(
        MinimizeTrajectoryParams::default(),
        MinimizeTrajectoryParams {
            max_its: 200,
            force_tol: 1e-4,
            energy_tol: 1e-6,
        }
    );
    assert_eq!(
        MMFFMolPropertiesParams::default(),
        MMFFMolPropertiesParams {
            mmff_variant: "MMFF94".to_owned(),
            verbosity: MMFFVerbosity::MmffVerbosityNone,
        }
    );
    assert_eq!(
        UFFOptimizeMoleculeParams::default(),
        UFFOptimizeMoleculeParams {
            max_iters: 200,
            vdw_thresh: 10.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    );
    assert_eq!(
        UFFOptimizeMoleculeConfsParams::default(),
        UFFOptimizeMoleculeConfsParams {
            num_threads: 1,
            max_iters: 200,
            vdw_thresh: 10.0,
            ignore_interfrag_interactions: true,
        }
    );
    assert_eq!(
        UFFGetMoleculeForceFieldParams::default(),
        UFFGetMoleculeForceFieldParams {
            vdw_thresh: 10.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    );
    assert_eq!(
        MMFFOptimizeMoleculeParams::default(),
        MMFFOptimizeMoleculeParams {
            mmff_variant: "MMFF94".to_owned(),
            max_iters: 200,
            non_bonded_thresh: 100.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    );
    assert_eq!(
        MMFFGetMoleculePropertiesParams::default(),
        MMFFGetMoleculePropertiesParams {
            mmff_variant: "MMFF94".to_owned(),
            mmff_verbosity: MMFFVerbosity::MmffVerbosityNone,
        }
    );
    assert_eq!(
        MMFFGetMoleculeForceFieldParams::default(),
        MMFFGetMoleculeForceFieldParams {
            non_bonded_thresh: 100.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    );
    assert_eq!(
        MMFFOptimizeMoleculeConfsParams::default(),
        MMFFOptimizeMoleculeConfsParams {
            num_threads: 1,
            max_iters: 200,
            mmff_variant: "MMFF94".to_owned(),
            non_bonded_thresh: 100.0,
            ignore_interfrag_interactions: true,
        }
    );
    assert_eq!(
        OptimizeMoleculeConfsParams::default(),
        OptimizeMoleculeConfsParams {
            num_threads: 1,
            max_iters: 200,
        }
    );
}

#[test]
fn uff_force_field_minimizes_the_conformer() {
    let mut mol = ethanol();
    {
        let mut ff = uff_get_molecule_force_field(&mut mol).unwrap();
        assert_close(ff.calc_energy().unwrap(), 10.151428469229186);
        assert_eq!(ff.dimension(), 3);
        assert_eq!(ff.num_points(), 9);
        let positions = ff.positions();
        assert_eq!(positions.len(), 27);
        assert_close(positions[0], -0.964);
        let grad = ff.calc_grad().unwrap();
        assert_eq!(grad.len(), 27);
        assert_close(grad[0], -0.1886532881145184);
        assert_close(
            ff.calc_energy_with_pos(&positions).unwrap(),
            10.151428469229186,
        );
        assert_eq!(ff.calc_grad_with_pos(&positions).unwrap(), grad);
        assert_eq!(ff.minimize().unwrap(), 0);
        assert_minimized(ff.calc_energy().unwrap(), 1.5367465474722097);
    }
    // The force field moved the molecule's atoms.
    let pos = mol.get_conformer().unwrap().get_atom_position(0).unwrap();
    assert_minimized_point(
        pos,
        (
            -0.9258790873923073,
            0.07359893865426484,
            0.015828518902529374,
        ),
    );
}

#[test]
fn minimize_steps_and_trajectories() {
    let mut mol = ethanol();
    let mut ff = uff_get_molecule_force_field(&mut mol).unwrap();
    let one_step = MinimizeParams {
        max_its: 1,
        ..Default::default()
    };
    for _ in 0..30 {
        ff.minimize_with_params(&one_step).unwrap();
    }
    let positions = ff.positions();
    assert_minimized_point(
        Point3D::new(positions[0], positions[1], positions[2]),
        (
            -0.9705908182100972,
            0.04271027257846574,
            0.06469367054822592,
        ),
    );

    let mut mol = ethanol();
    let mut ff = uff_get_molecule_force_field(&mut mol).unwrap();
    let (res, snapshots) = ff.minimize_trajectory(10).unwrap();
    assert_eq!(res, 0);
    let energies: Vec<f64> = snapshots.iter().map(Snapshot::get_energy).collect();
    let expected = [
        5.576696893996235,
        2.940556528653744,
        1.5582474769637487,
        1.5367465474722097,
    ];
    assert_eq!(energies.len(), expected.len());
    for (actual, expected) in energies.iter().zip(expected) {
        assert_minimized(*actual, expected);
    }
    // The last snapshot holds the final coordinates.
    let last = snapshots.last().unwrap();
    let point = last.get_point_3d(0).unwrap();
    assert_close(point.x, ff.positions()[0]);
    assert_eq!(last.get_point_3d(9), Err(Error::Index { index: 9, len: 9 }));
    // RDKit requires a two-dimensional trajectory for points in 2D.
    assert!(matches!(last.get_point_2d(0), Err(Error::Runtime { .. })));
    let mut snapshot = last.clone();
    snapshot.set_energy(1.0);
    assert_eq!(snapshot.get_energy(), 1.0);
}

#[test]
fn constraints_and_extra_points() {
    let mut mol = ethanol();
    let mut ff = uff_get_molecule_force_field_with_params(
        &mut mol,
        &UFFGetMoleculeForceFieldParams {
            vdw_thresh: 2.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_close(ff.calc_energy().unwrap(), 10.151428469229186);
    ff.uff_add_torsion_constraint(3, 0, 1, 2, true, -10.0, 10.0, 5.0)
        .unwrap();
    ff.mmff_add_angle_constraint(0, 1, 2, false, 90.0, 90.0, 5.0)
        .unwrap();
    assert_close(ff.calc_energy().unwrap(), 1411.116130292972);
    ff.add_distance_constraint(0, 2, 1.0, 1.0, 10.0).unwrap();
    assert_close(ff.calc_energy().unwrap(), 1419.9853905040716);
    let five = MinimizeParams {
        max_its: 5,
        ..Default::default()
    };
    assert_eq!(ff.minimize_with_params(&five).unwrap(), 1);
    assert_minimized(ff.calc_energy().unwrap(), 47.95527901082138);

    let mut mol = ethanol();
    let mut ff = uff_get_molecule_force_field(&mut mol).unwrap();
    assert_eq!(ff.add_extra_point(1.0, 2.0, 3.0), 10);
    assert_eq!(ff.num_points(), 9);
    assert_eq!(ff.positions().len(), 30);
    assert_eq!(
        ff.get_extra_point_pos(0).unwrap(),
        Point3D::new(1.0, 2.0, 3.0)
    );
    assert_eq!(
        ff.get_extra_point_pos(1),
        Err(Error::Index { index: 1, len: 1 })
    );
    ff.initialize().unwrap();
    assert_eq!(ff.num_points(), 10);
    ff.uff_add_distance_constraint(0, 9, false, 1.0, 1.5, 100.0)
        .unwrap();
    assert_close(ff.calc_energy().unwrap(), 340.0146003197031);

    let mut mol = ethanol();
    let mut ff = create_empty_force_field_for_mol(&mut mol).unwrap();
    assert_eq!(ff.num_points(), 9);
    assert_eq!(ff.calc_energy().unwrap(), 0.0);
    assert_eq!(ff.minimize().unwrap(), 0);
    ff.uff_add_position_constraint(0, 0.0, 10.0).unwrap();
    ff.uff_add_distance_constraint(0, 1, false, 2.0, 2.0, 50.0)
        .unwrap();
    ff.add_fixed_point(1).unwrap();
    assert_close(ff.calc_energy().unwrap(), 5.685604480312075);
    assert_eq!(ff.minimize().unwrap(), 0);
    assert_minimized(ff.calc_energy().unwrap(), 0.9476007467186802);
}

#[test]
fn constraint_indices_are_checked() {
    let mut mol = ethanol();
    let mut ff = create_empty_force_field_for_mol(&mut mol).unwrap();
    let index = Err(Error::Index { index: 9, len: 9 });
    assert_eq!(ff.add_distance_constraint(0, 9, 1.0, 2.0, 1.0), index);
    assert_eq!(ff.add_fixed_point(9), index);
    assert_eq!(
        ff.uff_add_distance_constraint(9, 0, false, 1.0, 2.0, 1.0),
        index
    );
    assert_eq!(
        ff.mmff_add_angle_constraint(0, 1, 9, false, 10.0, 20.0, 1.0),
        index
    );
    assert_eq!(
        ff.uff_add_torsion_constraint(0, 1, 2, 9, false, 10.0, 20.0, 1.0),
        index
    );
    assert_eq!(ff.mmff_add_position_constraint(9, 1.0, 1.0), index);
    // RDKit's own checks of the other arguments are Error::Runtime.
    assert!(matches!(
        ff.uff_add_distance_constraint(0, 1, false, 2.0, 1.0, 1.0),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        ff.uff_add_angle_constraint(0, 1, 2, false, 10.0, 200.0, 1.0),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn points_added_since_initialize_are_refused() {
    // The contributions index every point, but the gradient and explicit
    // coordinates are sized by num_points until initialize: in Python,
    // CalcGrad here writes past the gradient's end and crashes.
    let mut mol = ethanol();
    let mut ff = create_empty_force_field_for_mol(&mut mol).unwrap();
    ff.add_extra_point_with_fixed(5.0, 1.0, 0.0, false);
    ff.uff_add_angle_constraint(0, 1, 9, false, 0.0, 0.0, 100.0)
        .unwrap();
    assert!(matches!(ff.calc_grad(), Err(Error::Value { .. })));
    let positions = vec![0.0; 27];
    assert!(matches!(
        ff.calc_energy_with_pos(&positions),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        ff.calc_grad_with_pos(&positions),
        Err(Error::Value { .. })
    ));
    assert!(matches!(ff.minimize(), Err(Error::Runtime { .. })));
    // calc_energy reads every point.
    assert!(ff.calc_energy().unwrap() > 0.0);
    ff.initialize().unwrap();
    assert_eq!(ff.calc_grad().unwrap().len(), 30);
    assert!(matches!(
        ff.calc_energy_with_pos(&positions),
        Err(Error::Value { .. })
    ));
    let positions = ff.positions();
    assert!(ff.calc_energy_with_pos(&positions).is_ok());
}

#[test]
fn minimizer_size_is_capped() {
    // BFGS sizes its inverse Hessian in unsigned int, which wraps above
    // 65535 coordinates.
    let mut mol = ethanol();
    let mut ff = create_empty_force_field_for_mol(&mut mol).unwrap();
    for i in 0..21837 {
        ff.add_extra_point(f64::from(i), 0.0, 0.0);
    }
    assert!(matches!(ff.minimize(), Err(Error::Value { .. })));
    assert!(matches!(
        ff.minimize_trajectory(1),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        optimize_molecule(&mut ff),
        Err(Error::Value { .. })
    ));
    let mut big = mol_from_smiles(&"C".repeat(21846)).unwrap();
    assert!(matches!(
        uff_optimize_molecule(&mut big),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        mmff_optimize_molecule(&mut big),
        Err(Error::Value { .. })
    ));
}

#[test]
fn empty_molecules_fail() {
    let mut mol = Mol::default();
    mol.add_conformer(&Conformer::new(0)).unwrap();
    assert!(matches!(
        uff_get_molecule_force_field(&mut mol),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        create_empty_force_field_for_mol(&mut mol),
        Err(Error::Runtime { .. })
    ));
    // RDKit's threads throw here, which calls std::terminate in Python.
    let two = UFFOptimizeMoleculeConfsParams {
        num_threads: 2,
        ..Default::default()
    };
    assert!(matches!(
        uff_optimize_molecule_confs_with_params(&mut mol, &two),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        uff_get_molecule_force_field(&mut Mol::default()),
        Err(Error::Value { .. })
    ));
}

#[test]
fn mmff_force_fields() {
    let mut mol = ethanol();
    let props = mmff_get_molecule_properties(&mut mol).unwrap().unwrap();
    let mut ff = mmff_get_molecule_force_field(&mut mol, &props).unwrap();
    assert_close(ff.calc_energy().unwrap(), 6.831185199278828);
    let params = MinimizeParams {
        max_its: 1000,
        ..Default::default()
    };
    assert_eq!(ff.minimize_with_params(&params).unwrap(), 0);
    assert_minimized(ff.calc_energy().unwrap(), -1.3368570621313476);

    let mut mol = ethanol();
    let mut props = MMFFMolProperties::new(&mut mol).unwrap();
    props.set_mmff_vdw_term_with_state(false);
    let ff = mmff_get_molecule_force_field(&mut mol, &props).unwrap();
    assert_close(ff.calc_energy().unwrap(), 4.818755395391008);

    let mut mol = ethanol();
    let props = MMFFMolProperties::new_with_params(
        &mut mol,
        &MMFFMolPropertiesParams {
            mmff_variant: "MMFF94s".to_owned(),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(props.get_mmff_variant(), "MMFF94s");
    let ff = mmff_get_molecule_force_field_with_params(
        &mut mol,
        &props,
        &MMFFGetMoleculeForceFieldParams {
            non_bonded_thresh: 2.0,
            ignore_interfrag_interactions: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_close(ff.calc_energy().unwrap(), 4.818755395391008);
}

#[test]
fn mmff_properties() {
    let mut mol = ethanol();
    let mut props = mmff_get_molecule_properties(&mut mol).unwrap().unwrap();
    let types: Vec<u32> = (0..9)
        .map(|i| props.get_mmff_atom_type(i).unwrap())
        .collect();
    assert_eq!(types, vec![1, 1, 6, 5, 5, 5, 5, 5, 21]);
    let charges: Vec<f64> = (0..9)
        .map(|i| props.get_mmff_partial_charge(i).unwrap())
        .collect();
    assert_eq!(
        charges,
        vec![0.0, 0.28, -0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4]
    );
    assert_eq!(props.get_mmff_formal_charge(2).unwrap(), 0.0);
    assert_eq!(
        props.get_mmff_atom_type(9),
        Err(Error::Index { index: 9, len: 9 })
    );
    assert_eq!(
        props.get_mmff_bond_stretch_params(&mol, 0, 1).unwrap(),
        Some((0, 4.258, 1.508))
    );
    assert_eq!(
        props.get_mmff_bond_stretch_params(&mol, 0, 2).unwrap(),
        None
    );
    assert_eq!(
        props.get_mmff_angle_bend_params(&mol, 0, 1, 2).unwrap(),
        Some((0, 0.992, 108.133))
    );
    assert_eq!(
        props.get_mmff_stretch_bend_params(&mol, 0, 1, 2).unwrap(),
        Some((0, 0.173, 0.417))
    );
    assert_eq!(
        props
            .get_mmff_torsion_params(&mol, 3, 0, 1, 2)
            .unwrap()
            .map(|t| (t.torsion_type, t.v1, t.v2, t.v3)),
        Some((0, -0.654, 1.072, 0.279))
    );
    assert_eq!(
        props.get_mmff_oop_bend_params(&mol, 0, 1, 2, 3).unwrap(),
        None
    );
    let vdw = props.get_mmff_vdw_params(0, 2).unwrap().unwrap();
    assert_close(vdw.r_ij_star_unscaled, 3.770659452485765);
    assert_close(vdw.epsilon_unscaled, 0.06793631574377197);
    assert_close(vdw.r_ij_star, 3.770659452485765);
    assert_close(vdw.epsilon, 0.06793631574377197);

    // An index past either the molecule or the properties.
    let small = mol_from_smiles("CC").unwrap();
    assert_eq!(
        props.get_mmff_bond_stretch_params(&small, 0, 2),
        Err(Error::Index { index: 2, len: 2 })
    );
    let big = mol_from_smiles("CCCCCCCCCCCC").unwrap();
    assert_eq!(
        props.get_mmff_angle_bend_params(&big, 9, 10, 11),
        Err(Error::Index { index: 9, len: 9 })
    );
    assert_eq!(
        props.get_mmff_vdw_params(0, 9),
        Err(Error::Index { index: 9, len: 9 })
    );

    assert_eq!(props.get_mmff_dielectric_model(), 1);
    assert_eq!(props.get_mmff_dielectric_constant(), 1.0);
    assert_eq!(props.get_mmff_variant(), "MMFF94");
    assert!(props.get_mmff_bond_term() && props.get_mmff_ele_term());
    props.set_mmff_dielectric_model_with_diel_model(2);
    props
        .set_mmff_dielectric_constant_with_diel_const(4.0)
        .unwrap();
    props.set_mmff_variant_with_mmff_variant("MMFF94s").unwrap();
    props.set_mmff_ele_term_with_state(false);
    props.set_mmff_verbosity_with_verbosity(MMFFVerbosity::MmffVerbosityLow);
    assert_eq!(props.get_mmff_dielectric_model(), 2);
    assert_eq!(props.get_mmff_dielectric_constant(), 4.0);
    assert_eq!(props.get_mmff_variant(), "MMFF94s");
    assert!(!props.get_mmff_ele_term());
    props.set_mmff_dielectric_model();
    props.set_mmff_dielectric_constant();
    props.set_mmff_variant();
    props.set_mmff_ele_term();
    props.set_mmff_verbosity();
    assert_eq!(props.get_mmff_dielectric_model(), 1);
    assert_eq!(props.get_mmff_dielectric_constant(), 1.0);
    assert_eq!(props.get_mmff_variant(), "MMFF94");
    assert!(props.get_mmff_ele_term());
    for (get, set_off, set_on) in [
        (
            MMFFMolProperties::get_mmff_angle_term as fn(&MMFFMolProperties) -> bool,
            MMFFMolProperties::set_mmff_angle_term_with_state as fn(&mut MMFFMolProperties, bool),
            MMFFMolProperties::set_mmff_angle_term as fn(&mut MMFFMolProperties),
        ),
        (
            MMFFMolProperties::get_mmff_bond_term,
            MMFFMolProperties::set_mmff_bond_term_with_state,
            MMFFMolProperties::set_mmff_bond_term,
        ),
        (
            MMFFMolProperties::get_mmff_stretch_bend_term,
            MMFFMolProperties::set_mmff_stretch_bend_term_with_state,
            MMFFMolProperties::set_mmff_stretch_bend_term,
        ),
        (
            MMFFMolProperties::get_mmff_oop_term,
            MMFFMolProperties::set_mmff_oop_term_with_state,
            MMFFMolProperties::set_mmff_oop_term,
        ),
        (
            MMFFMolProperties::get_mmff_torsion_term,
            MMFFMolProperties::set_mmff_torsion_term_with_state,
            MMFFMolProperties::set_mmff_torsion_term,
        ),
        (
            MMFFMolProperties::get_mmff_vdw_term,
            MMFFMolProperties::set_mmff_vdw_term_with_state,
            MMFFMolProperties::set_mmff_vdw_term,
        ),
    ] {
        set_off(&mut props, false);
        assert!(!get(&props));
        set_on(&mut props);
        assert!(get(&props));
    }
    assert!(matches!(
        props.set_mmff_variant_with_mmff_variant("foo"),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        props.set_mmff_dielectric_constant_with_diel_const(0.0),
        Err(Error::Runtime { .. })
    ));

    // The builder indexes the properties by the molecule's atoms.
    let mut other = mol_from_smiles("CC").unwrap();
    assert!(matches!(
        mmff_get_molecule_force_field(&mut other, &props),
        Err(Error::Value { .. })
    ));
}

#[test]
fn untypeable_molecules() {
    let mut plutonium = add_hs(&mol_from_smiles("[Pu]").unwrap()).unwrap();
    assert!(!uff_has_all_molecule_params(&plutonium).unwrap());
    assert!(!mmff_has_all_molecule_params(&plutonium).unwrap());
    assert!(
        mmff_get_molecule_properties(&mut plutonium)
            .unwrap()
            .is_none()
    );
    let props = MMFFMolProperties::new(&mut plutonium).unwrap();
    plutonium.add_conformer(&Conformer::new(1)).unwrap();
    assert!(matches!(
        mmff_get_molecule_force_field(&mut plutonium, &props),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(mmff_optimize_molecule(&mut plutonium).unwrap(), -1);
    assert_eq!(
        mmff_optimize_molecule_confs(&mut plutonium).unwrap(),
        vec![(-1, -1.0)]
    );
    let mol = ethanol();
    assert!(uff_has_all_molecule_params(&mol).unwrap());
    assert!(mmff_has_all_molecule_params(&mol).unwrap());
}

#[test]
fn optimize_molecules() {
    let mut mol = ethanol();
    assert_eq!(uff_optimize_molecule(&mut mol).unwrap(), 0);
    let pos = mol.get_conformer().unwrap().get_atom_position(0).unwrap();
    assert_minimized_point(
        pos,
        (
            -0.9258790873923073,
            0.07359893865426484,
            0.015828518902529374,
        ),
    );
    let mut mol = ethanol();
    assert_eq!(mmff_optimize_molecule(&mut mol).unwrap(), 0);
    let pos = mol.get_conformer().unwrap().get_atom_position(0).unwrap();
    assert_minimized_point(
        pos,
        (
            -0.8864157419367157,
            0.16352926331719358,
            -0.07264422730952377,
        ),
    );

    let uff = [(0, 1.5367465474722097), (0, 1.5367465492621029)];
    let mmff = [(0, -1.3368570621313476), (0, -1.336857064033163)];
    let assert_results = |actual: Vec<(i32, f64)>, expected: &[(i32, f64)]| {
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected) {
            assert_eq!(a.0, e.0);
            assert_minimized(a.1, e.1);
        }
    };
    for num_threads in [1, 2, 0, -1, i32::MIN] {
        let params = UFFOptimizeMoleculeConfsParams {
            num_threads,
            ..Default::default()
        };
        let mut mol = ethanol_two_confs();
        assert_results(
            uff_optimize_molecule_confs_with_params(&mut mol, &params).unwrap(),
            &uff,
        );
        let params = MMFFOptimizeMoleculeConfsParams {
            num_threads,
            ..Default::default()
        };
        let mut mol = ethanol_two_confs();
        assert_results(
            mmff_optimize_molecule_confs_with_params(&mut mol, &params).unwrap(),
            &mmff,
        );
    }
    let mut mol = ethanol_two_confs();
    assert_results(uff_optimize_molecule_confs(&mut mol).unwrap(), &uff);
    let mut mol = ethanol_two_confs();
    assert_results(mmff_optimize_molecule_confs(&mut mol).unwrap(), &mmff);
    let params = MMFFOptimizeMoleculeConfsParams {
        mmff_variant: "MMFF94s".to_owned(),
        max_iters: 10,
        ..Default::default()
    };
    let mut mol = ethanol_two_confs();
    assert_results(
        mmff_optimize_molecule_confs_with_params(&mut mol, &params).unwrap(),
        &[(1, 2.807953571390044), (1, -0.3135959826309159)],
    );

    let mut mol = ethanol();
    let mut ff = uff_get_molecule_force_field(&mut mol).unwrap();
    assert_eq!(optimize_molecule(&mut ff).unwrap(), 0);
    assert_minimized(ff.calc_energy().unwrap(), 1.5367465474722097);
    assert_eq!(optimize_molecule_with_max_iters(&mut ff, 10).unwrap(), 0);

    let mut benzene = mol_from_smiles("c1ccccc1").unwrap();
    assert_eq!(
        mmff_sanitize_molecule(&mut benzene).unwrap(),
        SanitizeFlags::SANITIZE_NONE
    );
}

#[test]
fn optimize_molecule_confs_restores_the_force_field() {
    let mut template = ethanol();
    let mut ff = uff_get_molecule_force_field(&mut template).unwrap();
    for num_threads in [1, 2] {
        let mut mol = ethanol_two_confs();
        let params = OptimizeMoleculeConfsParams {
            num_threads,
            ..Default::default()
        };
        let results = optimize_molecule_confs_with_params(&mut mol, &mut ff, &params).unwrap();
        assert_eq!(results.len(), 2);
        assert_minimized(results[0].1, 1.5367465474722097);
        assert_minimized(results[1].1, 1.5367465492621029);
        // In Python the force field now points into `mol`'s last conformer,
        // and reading it once `mol` is gone crashes.
        drop(mol);
        assert_close(ff.positions()[0], -0.964);
        assert_close(ff.calc_energy().unwrap(), 10.151428469229186);
    }
    let mut mol = ethanol_two_confs();
    assert_eq!(optimize_molecule_confs(&mut mol, &mut ff).unwrap().len(), 2);
    drop(ff);
    let pos = template
        .get_conformer()
        .unwrap()
        .get_atom_position(0)
        .unwrap();
    assert_point(pos, (-0.964, 0.0151, 0.041));
}

#[test]
fn optimize_molecule_confs_checks_the_force_field() {
    // RDKit writes one position per atom into the force field: past the end
    // of one for a smaller molecule, which corrupts Python's heap.
    let mut small = mol_from_smiles("CC").unwrap();
    small.add_conformer(&Conformer::new(2)).unwrap();
    let mut ff = create_empty_force_field_for_mol(&mut small).unwrap();
    let mut mol = ethanol_two_confs();
    assert!(matches!(
        optimize_molecule_confs(&mut mol, &mut ff),
        Err(Error::Value { .. })
    ));

    // RDKit's threads copy the force field without its extra points.
    let mut template = ethanol();
    let mut ff = create_empty_force_field_for_mol(&mut template).unwrap();
    ff.add_extra_point(0.0, 0.0, 0.0);
    ff.initialize().unwrap();
    let two = OptimizeMoleculeConfsParams {
        num_threads: 2,
        ..Default::default()
    };
    assert!(matches!(
        optimize_molecule_confs_with_params(&mut mol, &mut ff, &two),
        Err(Error::Value { .. })
    ));
    // One thread uses the force field itself, extra points and all.
    assert_eq!(optimize_molecule_confs(&mut mol, &mut ff).unwrap().len(), 2);

    // A satisfied position constraint has no gradient, so the minimizer
    // fails an invariant: in Python's threads, that calls std::terminate.
    let mut template = ethanol();
    let mut ff = create_empty_force_field_for_mol(&mut template).unwrap();
    ff.uff_add_position_constraint(0, 1.0, 10.0).unwrap();
    assert!(matches!(ff.minimize(), Err(Error::Runtime { .. })));
    assert!(matches!(
        optimize_molecule_confs_with_params(&mut mol, &mut ff, &two),
        Err(Error::Runtime { .. })
    ));
    // The force field is unchanged and still usable.
    assert_eq!(ff.num_points(), 9);
    assert_close(ff.positions()[0], -0.964);
}

#[test]
fn uff_parameters() {
    let mol = ethanol();
    let (kb, r0) = get_uff_bond_stretch_params(&mol, 0, 1).unwrap().unwrap();
    assert_close(kb, 699.591798712679);
    assert_close(r0, 1.514);
    let (ka, theta0) = get_uff_angle_bend_params(&mol, 0, 1, 2).unwrap().unwrap();
    assert_close(ka, 290.0327445933201);
    assert_close(theta0, 109.47);
    assert_close(
        get_uff_torsion_params(&mol, 3, 0, 1, 2).unwrap().unwrap(),
        2.119,
    );
    assert_eq!(get_uff_inversion_params(&mol, 0, 1, 2, 3).unwrap(), None);
    let (x_ij, d_ij) = get_uff_vdw_params(&mol, 0, 2).unwrap().unwrap();
    assert_close(x_ij, 3.671307668937595);
    assert_close(d_ij, 0.07937253933193772);
    assert_eq!(get_uff_bond_stretch_params(&mol, 0, 2).unwrap(), None);
    // Python raises RDKit's RuntimeError for these.
    let index = Error::Index { index: 20, len: 9 };
    assert_eq!(get_uff_bond_stretch_params(&mol, 0, 20), Err(index.clone()));
    assert_eq!(
        get_uff_angle_bend_params(&mol, 0, 1, 20),
        Err(index.clone())
    );
    assert_eq!(
        get_uff_torsion_params(&mol, 0, 1, 2, 20),
        Err(index.clone())
    );
    assert_eq!(
        get_uff_inversion_params(&mol, 0, 1, 2, 20),
        Err(index.clone())
    );
    assert_eq!(get_uff_vdw_params(&mol, 20, 0), Err(index));
}

#[test]
fn force_fields_are_send() {
    fn assert_send<T: Send>() {}
    assert_send::<ForceField<'_>>();
    assert_send::<MMFFMolProperties>();
}
