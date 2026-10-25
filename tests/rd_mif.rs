// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{assert_close, assert_close_to, ethanol, temp_path};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mif::*;
use rdkit_cxx::chem::rd_partial_charges::compute_gasteiger_charges;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::{Point3D, UniformRealValueGrid3D, UniformRealValueGrid3DParams};

// Expected values come from Python's rdMIF on the same molecule, ethanol
// with Gasteiger charges, at the point (3, 1.5, 1) with thresholds 100 and
// 9, and on the grid ConstructGrid(mol, margin=2.0, spacing=1.0): its first
// value and its least, without a threshold and with threshold 2.

const POINT: (f64, f64, f64) = (3.0, 1.5, 1.0);

fn charged_ethanol() -> Mol {
    let mut mol = ethanol();
    compute_gasteiger_charges(&mut mol).unwrap();
    mol
}

fn small_grid(mol: &Mol) -> UniformRealValueGrid3D {
    let params = ConstructGridParams {
        margin: 2.0,
        spacing: 1.0,
        ..Default::default()
    };
    construct_grid_with_params(mol, &params).unwrap()
}

fn values(grid: &UniformRealValueGrid3D) -> Vec<f64> {
    (0..grid.get_size())
        .map(|i| grid.get_val(i).unwrap())
        .collect()
}

fn min(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::INFINITY, f64::min)
}

struct Expected {
    at_100: f64,
    at_9: f64,
    first: f64,
    least: f64,
    least_within_2: f64,
}

fn expected(values: (f64, f64, f64, f64, f64)) -> Expected {
    Expected {
        at_100: values.0,
        at_9: values.1,
        first: values.2,
        least: values.3,
        least_within_2: values.4,
    }
}

// Energies near atoms are steep, so they are compared to 1e-6.
macro_rules! check_functor {
    ($mol:expr, $functor:expr, $calc:ident, $calc_threshold:ident, $expected:expr) => {{
        let functor = $functor;
        let expected = expected($expected);
        let (x, y, z) = POINT;
        assert_close_to(functor.call(x, y, z, 100.0), expected.at_100, 1e-6);
        assert_close_to(functor.call(x, y, z, 9.0), expected.at_9, 1e-6);
        assert_eq!(functor.call(x, y, z, -1.0), 0.0);
        let mut grid = small_grid($mol);
        $calc(&mut grid, &functor);
        let vals = values(&grid);
        assert_close_to(vals[0], expected.first, 1e-6);
        assert_close_to(min(&vals), expected.least, 1e-6);
        let mut grid = small_grid($mol);
        $calc_threshold(&mut grid, &functor, 2.0);
        assert_close_to(min(&values(&grid)), expected.least_within_2, 1e-6);
    }};
}

#[test]
fn param_defaults() {
    assert_eq!(
        ConstructGridParams::default(),
        ConstructGridParams {
            conf_id: -1,
            margin: 5.0,
            spacing: 0.5
        }
    );
    let coulomb = CoulombParams::default();
    assert_eq!(coulomb.conf_id, -1);
    assert_eq!(coulomb.probe_charge, 1.0);
    assert!(!coulomb.abs_val);
    assert_eq!(coulomb.charge_key, "_GasteigerCharge");
    assert_eq!(coulomb.softcore_param, 0.0);
    assert_eq!(coulomb.cutoff, 1.0);
    let from_charges = CoulombFromChargesParams::default();
    assert_eq!(
        (
            from_charges.probe_charge,
            from_charges.abs_val,
            from_charges.softcore_param,
            from_charges.cutoff
        ),
        (1.0, false, 0.0, 1.0)
    );
    let dielectric = CoulombDielectricParams::default();
    assert_eq!(dielectric.coulomb, coulomb);
    assert_eq!((dielectric.epsilon, dielectric.xi), (80.0, 4.0));
    let dielectric = CoulombDielectricFromChargesParams::default();
    assert_eq!(dielectric.coulomb, from_charges);
    assert_eq!((dielectric.epsilon, dielectric.xi), (80.0, 4.0));
    assert_eq!(
        MMFFVdWaalsParams::default(),
        MMFFVdWaalsParams {
            conf_id: -1,
            probe_atom_type: 6,
            scaling: false,
            cutoff: 1.0
        }
    );
    assert_eq!(
        UFFVdWaalsParams::default(),
        UFFVdWaalsParams {
            conf_id: -1,
            probe_atom_type: "O_3".to_string(),
            cutoff: 1.0
        }
    );
    assert_eq!(
        HBondParams::default(),
        HBondParams {
            conf_id: -1,
            probe_atom_type: "OH".to_string(),
            fixed: true,
            cutoff: 1.0
        }
    );
    assert_eq!(
        HydrophilicParams::default(),
        HydrophilicParams {
            conf_id: -1,
            fixed: true,
            cutoff: 1.0
        }
    );
    let cube = WriteToCubeFileParams::default();
    assert!(cube.mol.is_none());
    assert_eq!(cube.conf_id, -1);
}

#[test]
fn constructs_grids() {
    let mol = ethanol();
    let grid = construct_grid(&mol).unwrap();
    assert_eq!(
        (grid.get_num_x(), grid.get_num_y(), grid.get_num_z()),
        (27, 25, 25)
    );
    assert_close(grid.get_spacing(), 0.5);
    let offset = grid.get_offset();
    assert_close(offset.x, -6.5833);
    assert_close(offset.y, -6.4177);
    assert_close(offset.z, -6.23);
    let grid = small_grid(&mol);
    assert_eq!(
        (grid.get_num_x(), grid.get_num_y(), grid.get_num_z()),
        (7, 6, 6)
    );

    let missing = ConstructGridParams {
        conf_id: 4,
        ..Default::default()
    };
    assert!(matches!(
        construct_grid_with_params(&mol, &missing),
        Err(Error::Value { .. })
    ));
    let mut bare = mol.clone();
    bare.remove_all_conformers();
    assert!(matches!(construct_grid(&bare), Err(Error::Runtime { .. })));
    // RDKit casts the grid's size to unsigned int.
    let tiny = ConstructGridParams {
        spacing: 1e-5,
        ..Default::default()
    };
    assert!(matches!(
        construct_grid_with_params(&mol, &tiny),
        Err(Error::Value { .. })
    ));
    let mut far = mol.clone();
    far.get_conformer_mut()
        .unwrap()
        .set_atom_position(2, Point3D::new(f64::INFINITY, 0.0, 0.0))
        .unwrap();
    assert!(matches!(construct_grid(&far), Err(Error::Value { .. })));
}

#[test]
fn coulomb_fields() {
    let mol = charged_ethanol();
    check_functor!(
        &mol,
        Coulomb::new(&mol).unwrap(),
        calculate_descriptors,
        calculate_descriptors_with_threshold,
        (
            22.951017537779897,
            -45.14211502621587,
            7.166940911341788,
            -240.64474824058544,
            -400.6958479998257
        )
    );
    let params = CoulombParams {
        probe_charge: -1.0,
        abs_val: true,
        softcore_param: 0.5,
        ..Default::default()
    };
    check_functor!(
        &mol,
        Coulomb::new_with_params(&mol, &params).unwrap(),
        calculate_descriptors,
        calculate_descriptors_with_threshold,
        (
            -11.059202371605036,
            -55.823480570920275,
            -7.007557528833066,
            -308.0573114429,
            -386.4782227483537
        )
    );
    let charges: Vec<f64> = (0..9)
        .map(|i| {
            mol.get_atom_with_idx(i)
                .unwrap()
                .get_double_prop("_GasteigerCharge")
                .unwrap()
        })
        .collect();
    let positions = mol.get_conformer().unwrap().get_positions();
    check_functor!(
        &mol,
        Coulomb::from_charges(&charges, &positions).unwrap(),
        calculate_descriptors,
        calculate_descriptors_with_threshold,
        (
            22.951017537779897,
            -45.14211502621587,
            7.166940911341788,
            -240.64474824058544,
            -400.6958479998257
        )
    );
    // As in Python, the charges and positions must match.
    assert!(matches!(
        Coulomb::from_charges(&charges[..3], &positions),
        Err(Error::Value { .. })
    ));
    // Without charges, or the conformer.
    assert!(matches!(Coulomb::new(&ethanol()), Err(Error::Key { .. })));
    let missing = CoulombParams {
        conf_id: 2,
        ..Default::default()
    };
    assert!(matches!(
        Coulomb::new_with_params(&mol, &missing),
        Err(Error::Value { .. })
    ));
}

#[test]
fn coulomb_dielectric_fields() {
    let mol = charged_ethanol();
    check_functor!(
        &mol,
        CoulombDielectric::new(&mol).unwrap(),
        calculate_descriptors_coulomb_dielectric,
        calculate_descriptors_coulomb_dielectric_with_threshold,
        (
            0.5464527985185674,
            -1.074812262528951,
            0.17064145027004257,
            -50.220459070107985,
            -60.274526074343946
        )
    );
    let charges: Vec<f64> = (0..9)
        .map(|i| {
            mol.get_atom_with_idx(i)
                .unwrap()
                .get_double_prop("_GasteigerCharge")
                .unwrap()
        })
        .collect();
    let positions = mol.get_conformer().unwrap().get_positions();
    let params = CoulombDielectricFromChargesParams {
        epsilon: 60.0,
        xi: 2.0,
        ..Default::default()
    };
    check_functor!(
        &mol,
        CoulombDielectric::from_charges_with_params(&charges, &positions, &params).unwrap(),
        calculate_descriptors_coulomb_dielectric,
        calculate_descriptors_coulomb_dielectric_with_threshold,
        (
            0.7403554044444975,
            -1.4561972589102072,
            0.23119164230134606,
            -99.76582625634244,
            -117.83941303659198
        )
    );
    assert!(matches!(
        CoulombDielectric::from_charges(&charges, &positions[1..]),
        Err(Error::Value { .. })
    ));
    let mut bare = mol.clone();
    bare.remove_all_conformers();
    assert!(matches!(
        CoulombDielectric::new(&bare),
        Err(Error::Runtime { .. })
    ));
}

// RDKit's loop over pairs of atoms never ends for a molecule without atoms
// (Python hangs); the field is 0.
#[test]
fn coulomb_dielectric_without_atoms() {
    let mut mol = Mol::default();
    mol.add_conformer(&Conformer::new(0)).unwrap();
    let field = CoulombDielectric::new(&mol).unwrap();
    assert_eq!(field.call(1.0, 2.0, 3.0, 100.0), 0.0);
    let missing = CoulombDielectricParams {
        coulomb: CoulombParams {
            conf_id: 3,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(matches!(
        CoulombDielectric::new_with_params(&mol, &missing),
        Err(Error::Value { .. })
    ));
}

#[test]
fn van_der_waals_fields() {
    let mol = ethanol();
    check_functor!(
        &mol,
        MMFFVdWaals::new(&mol).unwrap(),
        calculate_descriptors_mmff_vd_waals,
        calculate_descriptors_mmff_vd_waals_with_threshold,
        (
            57.0093707772066,
            57.10409111926597,
            -0.03577634575222358,
            -0.2799143790454349,
            0.0
        )
    );
    let params = MMFFVdWaalsParams {
        probe_atom_type: 1,
        scaling: true,
        ..Default::default()
    };
    check_functor!(
        &mol,
        MMFFVdWaals::new_with_params(&mol, &params).unwrap(),
        calculate_descriptors_mmff_vd_waals,
        calculate_descriptors_mmff_vd_waals_with_threshold,
        (
            73.97256385737666,
            73.95894147521567,
            -0.046958897066036426,
            -0.24813946307312254,
            0.0
        )
    );
    check_functor!(
        &mol,
        UFFVdWaals::new(&mol).unwrap(),
        calculate_descriptors_uff_vd_waals,
        calculate_descriptors_uff_vd_waals_with_threshold,
        (
            1276.651228651169,
            1276.81672717869,
            -0.048315810195477316,
            -0.3520560988755732,
            0.0
        )
    );
    let params = UFFVdWaalsParams {
        probe_atom_type: "N_3".to_string(),
        ..Default::default()
    };
    check_functor!(
        &mol,
        UFFVdWaals::new_with_params(&mol, &params).unwrap(),
        calculate_descriptors_uff_vd_waals,
        calculate_descriptors_uff_vd_waals_with_threshold,
        (
            1793.53391235818,
            1793.6823416137017,
            -0.058914831785180075,
            -0.3901402279190317,
            0.0
        )
    );
}

// Python dereferences null for a probe type without parameters.
#[test]
fn probe_types_need_parameters() {
    let mol = ethanol();
    for probe_atom_type in [0, 200, u32::MAX] {
        let params = MMFFVdWaalsParams {
            probe_atom_type,
            ..Default::default()
        };
        assert!(matches!(
            MMFFVdWaals::new_with_params(&mol, &params),
            Err(Error::Value { .. })
        ));
    }
    for probe_atom_type in ["Xx", "", "O_3 "] {
        let params = UFFVdWaalsParams {
            probe_atom_type: probe_atom_type.to_string(),
            ..Default::default()
        };
        assert!(matches!(
            UFFVdWaals::new_with_params(&mol, &params),
            Err(Error::Value { .. })
        ));
    }
    let params = HBondParams {
        probe_atom_type: "S".to_string(),
        ..Default::default()
    };
    assert!(matches!(
        HBond::new_with_params(&mol, &params),
        Err(Error::Value { .. })
    ));
    // An atom MMFF cannot type.
    let mut odd = ethanol();
    odd.get_atom_with_idx_mut(0).unwrap().set_atomic_num(118);
    assert!(matches!(MMFFVdWaals::new(&odd), Err(Error::Value { .. })));
}

#[test]
fn hydrogen_bond_fields() {
    let mol = ethanol();
    let cases = [
        ("OH", true, (0.0, 0.0, 0.0, -7.773263854681051, 0.0)),
        (
            "OH",
            false,
            (
                664.7144912152708,
                664.7144912152708,
                0.0,
                -19.65232524570803,
                0.0,
            ),
        ),
        (
            "O",
            true,
            (
                452.461243286416,
                452.461243286416,
                0.0,
                -12.226170547686362,
                0.0,
            ),
        ),
        (
            "O",
            false,
            (
                400.9216257263478,
                400.9216257263478,
                0.0,
                -14.329300690521487,
                0.0,
            ),
        ),
        (
            "N",
            true,
            (
                1374.668896749454,
                1374.668896749454,
                0.0,
                -5.010398237450948,
                0.0,
            ),
        ),
        (
            "N",
            false,
            (
                1218.081099979116,
                1218.081099979116,
                0.0,
                -7.202150658728089,
                0.0,
            ),
        ),
        ("NH", true, (0.0, 0.0, 0.0, -9.776856134553752, 0.0)),
        (
            "NH",
            false,
            (
                2019.5372528600049,
                2019.5372528600049,
                0.0,
                -10.431501615997743,
                0.0,
            ),
        ),
    ];
    for (probe, fixed, values) in cases {
        let params = HBondParams {
            probe_atom_type: probe.to_string(),
            fixed,
            ..Default::default()
        };
        check_functor!(
            &mol,
            HBond::new_with_params(&mol, &params).unwrap(),
            calculate_descriptors_h_bond,
            calculate_descriptors_h_bond_with_threshold,
            values
        );
    }
    check_functor!(
        &mol,
        Hydrophilic::new(&mol).unwrap(),
        calculate_descriptors_hydrophilic,
        calculate_descriptors_hydrophilic_with_threshold,
        (0.0, 0.0, 0.0, -12.226170547686362, 0.0)
    );
    let params = HydrophilicParams {
        fixed: false,
        ..Default::default()
    };
    check_functor!(
        &mol,
        Hydrophilic::new_with_params(&mol, &params).unwrap(),
        calculate_descriptors_hydrophilic,
        calculate_descriptors_hydrophilic_with_threshold,
        (
            400.9216257263478,
            400.9216257263478,
            0.0,
            -19.65232524570803,
            0.0
        )
    );
}

// HBond writes each nitrogen's or oxygen's bonds into an array of 12, past
// its end for more neighbors: a segfault in Python.
#[test]
fn hydrogen_bonds_need_few_neighbors() {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    for smiles in [
        "N(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)C",
        "O(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)C",
    ] {
        let mut mol = mol_from_smiles_with_params(smiles, &params).unwrap();
        let mut conf = Conformer::new(mol.get_num_atoms());
        for i in 1..mol.get_num_atoms() {
            conf.set_atom_position(i, Point3D::new(f64::from(i), 1.0, 0.0))
                .unwrap();
        }
        mol.add_conformer(&conf).unwrap();
        assert!(matches!(HBond::new(&mol), Err(Error::Value { .. })));
        assert!(matches!(Hydrophilic::new(&mol), Err(Error::Value { .. })));
    }
}

#[test]
fn cube_files_round_trip() {
    let mol = ethanol();
    let params = UniformRealValueGrid3DParams {
        spacing: 1.0,
        off_set: Some(Point3D::new(0.0, 0.0, 0.0)),
    };
    let mut grid = UniformRealValueGrid3D::new_with_params(5.0, 5.0, 5.0, &params).unwrap();
    for i in 0..grid.get_size() {
        grid.set_val(i, f64::from(i) / 10.0).unwrap();
    }
    let path = temp_path("grid.cube");
    let with_mol = WriteToCubeFileParams {
        mol: Some(&mol),
        conf_id: -1,
    };
    write_to_cube_file_with_params(&grid, &path, &with_mol).unwrap();
    let (grid2, mol2) = read_from_cube_file(&path).unwrap();
    assert_eq!(grid2.get_size(), grid.get_size());
    for i in 0..grid.get_size() {
        assert_close_to(grid2.get_val(i).unwrap(), f64::from(i) / 10.0, 1e-5);
    }
    assert!(grid.compare_grids(&grid2));
    let mol2 = mol2.unwrap();
    assert_eq!(mol2.get_num_atoms(), 9);
    assert_eq!(mol2.get_num_bonds(), 0);
    let (pos, pos2) = (
        mol.get_conformer().unwrap().get_positions(),
        mol2.get_conformer().unwrap().get_positions(),
    );
    for (a, b) in pos.iter().zip(&pos2) {
        assert_close_to(a.x, b.x, 1e-4);
        assert_close_to(a.y, b.y, 1e-4);
        assert_close_to(a.z, b.z, 1e-4);
    }
    for i in 0..9 {
        assert_eq!(
            mol2.get_atom_with_idx(i).unwrap().get_atomic_num(),
            mol.get_atom_with_idx(i).unwrap().get_atomic_num()
        );
    }

    write_to_cube_file(&grid, &path).unwrap();
    let (grid3, mol3) = read_from_cube_file(&path).unwrap();
    assert!(grid.compare_grids(&grid3));
    assert!(mol3.is_none());

    let missing = WriteToCubeFileParams {
        mol: Some(&mol),
        conf_id: 7,
    };
    assert!(matches!(
        write_to_cube_file_with_params(&grid, &path, &missing),
        Err(Error::Value { .. })
    ));
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(read_from_cube_file(&path), Err(Error::Os { .. })));
    let dir = std::env::temp_dir();
    assert!(matches!(
        write_to_cube_file(&grid, &dir),
        Err(Error::Os { .. })
    ));
}

const CUBE_HEADER: &str = "title\ncomment\n";

fn read_cube(name: &str, text: &str) -> rdkit_cxx::Result<(UniformRealValueGrid3D, Option<Mol>)> {
    let path = temp_path(name);
    std::fs::write(&path, text).unwrap();
    let result = read_from_cube_file(&path);
    std::fs::remove_file(&path).unwrap();
    result
}

// RDKit reads cube files into uninitialized variables, so a value that does
// not parse, or is missing, has it read indeterminate values, and a negative
// atom count makes it allocate a conformer of about 2^32 atoms (bad_alloc in
// Python).
#[test]
fn cube_files_are_checked() {
    let grid = "2 1.0 0 0\n2 0 1.0 0\n2 0 0 1.0\n";
    let values = "1 2 3 4 5 6 7 8\n";
    let ok = format!("{CUBE_HEADER}0 0 0 0\n{grid}{values}");
    let (read, mol) = read_cube("ok.cube", &ok).unwrap();
    assert_eq!(read.get_size(), 8);
    assert!(mol.is_none());
    assert_close(read.get_val(7).unwrap(), 8.0);

    for (name, text) in [
        ("empty.cube", String::new()),
        ("header.cube", format!("{CUBE_HEADER}0 0 0 0\n2 1.0 0 0\n")),
        ("word.cube", format!("{CUBE_HEADER}0 0 0 x\n{grid}{values}")),
        (
            "negative.cube",
            format!("{CUBE_HEADER}-5 0 0 0\n{grid}{values}"),
        ),
        (
            "atoms.cube",
            format!("{CUBE_HEADER}2 0 0 0\n{grid}6 6 0 0 0\n{values}"),
        ),
        (
            "values.cube",
            format!("{CUBE_HEADER}0 0 0 0\n{grid}1 2 3 4 5 6 7\n"),
        ),
        (
            "nan.cube",
            format!("{CUBE_HEADER}0 0 0 0\n{grid}1 2 3 nan 5 6 7 8\n"),
        ),
        (
            "huge.cube",
            format!("{CUBE_HEADER}0 0 0 0\n2000 1.0 0 0\n2000 0 1.0 0\n2000 0 0 1.0\n{values}"),
        ),
    ] {
        assert!(
            matches!(read_cube(name, &text), Err(Error::Value { .. })),
            "{name}"
        );
    }
    // RDKit's own checks.
    let spacings = format!("{CUBE_HEADER}0 0 0 0\n2 1.0 0 0\n2 0 2.0 0\n2 0 0 1.0\n{values}");
    assert!(matches!(
        read_cube("spacings.cube", &spacings),
        Err(Error::Parse { .. })
    ));
    let flat = format!("{CUBE_HEADER}0 0 0 0\n0 1.0 0 0\n2 0 1.0 0\n2 0 0 1.0\n");
    assert!(matches!(
        read_cube("flat.cube", &flat),
        Err(Error::Runtime { .. })
    ));
}

// The functors that compute in buffers of their own while called are not
// Sync (rd_mif's module documentation).
#[test]
fn functors_move_between_threads() {
    fn send<T: Send>() {}
    fn sync<T: Sync>() {}
    send::<Coulomb>();
    send::<CoulombDielectric>();
    send::<MMFFVdWaals>();
    send::<UFFVdWaals>();
    send::<HBond>();
    send::<Hydrophilic>();
    sync::<Coulomb>();
    sync::<MMFFVdWaals>();
    sync::<UFFVdWaals>();
    let mol = ethanol();
    let field = HBond::new(&mol).unwrap();
    let moved = std::thread::spawn(move || field.call(3.0, 1.5, 1.0, 100.0))
        .join()
        .unwrap();
    assert_eq!(moved, 0.0);
}
