// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use std::path::PathBuf;

use common::{assert_close, ethanol};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_force_field_helpers::uff_get_molecule_force_field;
use rdkit_cxx::chem::rdtrajectory::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Ported from RDKit's testTrajectory.py; the trajectory files are RDKit's
// Code/GraphMol/test_data/water_coords*.

const AMBER: &str = "Coordinate of a water molecule
   0.1941767  -0.0136074  -0.2917594  -0.4915838  -0.4480844  -1.4087115
   0.1349280   1.0427680  -0.4088006
";

const AMBER_BAD: &str = "Coordinate of a water molecule
Coordinate of a water molecule
   0.1941767  -0.0136074  -0.2917594  -0.4915838  -0.4480844  -1.4087115
   0.1349280   1.0427680  -0.4088006
";

const AMBER_BAD2: &str = "Coordinate of a water molecule
   0.1349280   1.0427680  -0.4088006
";

const GROMOS_HEADER: &str = "TITLE
Coordinate of a water molecule
END
";

const GROMOS_FRAME: &str = "TIMESTEP
              0    0.000000000
END
POSITIONRED
    0.1941767  -0.0136074  -0.2917594
   -0.4915838  -0.4480844  -1.4087115
    0.1349280   1.0427680  -0.4088006
END
GENBOX
    1
    2.898752912    2.898752912    2.898752912
   90.000000000   90.000000000   90.000000000
    0.000000000    0.000000000    0.000000000
    0.000000000    0.000000000    0.000000000
END
";

// The positions block lacks its END.
const GROMOS_BAD: &str = "TITLE
Coordinate of a water molecule
END
POSITIONRED
    0.1941767  -0.0136074  -0.2917594
   -0.4915838  -0.4480844  -1.4087115
    0.1349280   1.0427680  -0.4088006
GENBOX
    1
END
";

// The positions block has one position of three.
const GROMOS_BAD2: &str = "TITLE
Coordinate of a water molecule
END
POSITIONRED
    0.1941767  -0.0136074  -0.2917594
END
";

fn write(name: &str, text: &str) -> PathBuf {
    let path = common::temp_path(name);
    std::fs::write(&path, text).unwrap();
    path
}

fn coords(n: usize) -> Vec<f64> {
    (0..n).map(|i| i as f64).collect()
}

fn runtime<T: std::fmt::Debug>(result: Result<T, Error>) {
    assert!(matches!(result, Err(Error::Runtime { .. })), "{result:?}");
}

fn value<T: std::fmt::Debug>(result: Result<T, Error>) {
    assert!(matches!(result, Err(Error::Value { .. })), "{result:?}");
}

#[test]
fn defaults_are_pythons() {
    assert_eq!(
        AddConformersToMolParams::default(),
        AddConformersToMolParams {
            from_cid: -1,
            to_cid: -1,
        }
    );
    assert_eq!(Snapshot::new(&[1.0]).get_energy(), 0.0);
    assert_eq!(
        Snapshot::new(&[1.0]),
        Snapshot::new_with_energy(&[1.0], 0.0)
    );
    let traj = Trajectory::new(3, 2);
    assert!(traj.is_empty());
    assert_eq!(traj, Trajectory::new_with_snapshot_list(3, 2, &[]).unwrap());
}

#[test]
fn snapshot() {
    // Snapshots in no trajectory cannot read their points.
    runtime(Snapshot::new(&[]).get_point_2d(12));
    runtime(Snapshot::new(&[0.0, 0.0, 0.0]).get_point_2d(0));
    runtime(Snapshot::new(&[0.0, 0.0, 0.0]).get_point_3d(0));
    let mut s = Snapshot::new_with_energy(&[0.0], 2.0);
    assert_eq!(s.get_energy(), 2.0);
    s.set_energy(3.0);
    assert_eq!(s.get_energy(), 3.0);
}

fn check_trajectory(dim: u32) {
    let np = 10;
    let ns = 5;
    let mut traj = Trajectory::new(dim, np);
    assert_eq!(traj.dimension(), dim);
    assert_eq!(traj.num_points(), np);
    let c = coords((np * dim) as usize);
    for i in 0..ns {
        assert_eq!(
            traj.add_snapshot(&Snapshot::new_with_energy(&c, i as f64))
                .unwrap(),
            i
        );
    }
    assert_eq!(traj.len(), ns as usize);
    assert_eq!(
        traj.get_snapshot(ns),
        Err(Error::Index {
            index: ns.into(),
            len: ns as usize,
        })
    );
    let first = traj.get_snapshot(0).unwrap();
    if dim == 2 {
        assert_eq!(
            first.get_point_2d(np),
            Err(Error::Index {
                index: np.into(),
                len: np as usize,
            })
        );
    } else {
        runtime(first.get_point_2d(np));
        runtime(first.get_point_2d(0));
    }
    for i in 0..np {
        let d = f64::from(dim);
        let point = first.get_point_3d(i).unwrap();
        assert_eq!(point.x, f64::from(i) * d);
        assert_eq!(point.y, f64::from(i) * d + 1.0);
        let z = if dim == 3 {
            f64::from(i) * d + 2.0
        } else {
            0.0
        };
        assert_eq!(point.z, z);
        if dim == 2 {
            let point = first.get_point_2d(i).unwrap();
            assert_eq!(
                (point.x, point.y),
                (f64::from(i) * d, f64::from(i) * d + 1.0)
            );
        }
    }
    for i in 0..ns {
        assert_eq!(traj.get_snapshot(i).unwrap().get_energy(), i as f64);
    }
    assert_eq!(traj.remove_snapshot(0).unwrap(), 0);
    assert_eq!(traj.len(), (ns - 1) as usize);
    for i in 0..ns - 1 {
        assert_eq!(traj.get_snapshot(i).unwrap().get_energy(), (i + 1) as f64);
    }
    assert_eq!(
        traj.insert_snapshot(0, &Snapshot::new_with_energy(&c, 999.0))
            .unwrap(),
        0
    );
    assert_eq!(traj.len(), ns as usize);
    let copy = traj.get_snapshot(0).unwrap().clone();
    // A copy keeps its trajectory's layout, as in Python.
    assert_eq!(
        copy.get_point_3d(1),
        traj.get_snapshot(0).unwrap().get_point_3d(1)
    );
    assert_eq!(traj.add_snapshot(&copy).unwrap(), ns);
    assert_eq!(traj.len(), (ns + 1) as usize);
    assert_eq!(traj.get_snapshot(0).unwrap().get_energy(), 999.0);
    assert_eq!(traj.get_snapshot(1).unwrap().get_energy(), 1.0);
    assert_eq!(traj.get_snapshot(ns).unwrap().get_energy(), 999.0);
    let traj2 = traj.clone();
    assert_eq!(traj.len(), traj2.len());
}

#[test]
fn trajectory_2d() {
    check_trajectory(2);
}

#[test]
fn trajectory_3d() {
    check_trajectory(3);
}

#[test]
fn trajectory_edits() {
    let mut traj = Trajectory::new(2, 1);
    let s = Snapshot::new(&[1.0, 2.0]);
    // Insertion takes any position up to the length.
    assert_eq!(
        traj.insert_snapshot(1, &s),
        Err(Error::Index { index: 1, len: 1 })
    );
    assert_eq!(traj.insert_snapshot(0, &s).unwrap(), 0);
    assert_eq!(traj.insert_snapshot(1, &s).unwrap(), 1);
    assert_eq!(
        traj.remove_snapshot(2),
        Err(Error::Index { index: 2, len: 2 })
    );
    // Removing the last returns the new length.
    assert_eq!(traj.remove_snapshot(1).unwrap(), 1);
    traj.clear();
    assert!(traj.is_empty());
    assert_eq!(traj.get_snapshot(0), Err(Error::Index { index: 0, len: 0 }));
    // Snapshots go into the trajectory as copies.
    let mut s = Snapshot::new_with_energy(&[1.0, 2.0], 1.0);
    traj.add_snapshot(&s).unwrap();
    s.set_energy(2.0);
    assert_eq!(traj.get_snapshot(0).unwrap().get_energy(), 1.0);
    runtime(s.get_point_2d(0));
    // A snapshot from one trajectory reads as another lays it out.
    let mut other = Trajectory::new(3, 1);
    other
        .add_snapshot(&Snapshot::new(&[1.0, 2.0, 3.0, 4.0]))
        .unwrap();
    traj.add_snapshot(other.get_snapshot(0).unwrap()).unwrap();
    let moved = traj.get_snapshot(1).unwrap();
    assert_eq!(moved.get_point_2d(0).unwrap().y, 2.0);
    // A dimension past 3 reads z as 0, as RDKit does.
    let mut wide = Trajectory::new(4, 1);
    wide.add_snapshot(&Snapshot::new(&[1.0, 2.0, 3.0, 4.0]))
        .unwrap();
    assert_eq!(
        wide.get_snapshot(0).unwrap().get_point_3d(0).unwrap(),
        Point3D::new(1.0, 2.0, 0.0)
    );
    // A dimension below 2 cannot read points; no points, neither.
    let mut line = Trajectory::new(1, 1);
    line.add_snapshot(&Snapshot::new(&[1.0])).unwrap();
    runtime(line.get_snapshot(0).unwrap().get_point_3d(0));
    let mut none = Trajectory::new(3, 0);
    none.add_snapshot(&Snapshot::new(&[1.0])).unwrap();
    runtime(none.get_snapshot(0).unwrap().get_point_3d(0));
}

#[test]
fn short_snapshots_are_refused() {
    // Python adds them and reads past their coordinates.
    let mut traj = Trajectory::new(3, 2);
    value(traj.add_snapshot(&Snapshot::new(&coords(5))));
    value(traj.insert_snapshot(0, &Snapshot::new(&coords(5))));
    value(Trajectory::new_with_snapshot_list(
        3,
        2,
        &[Snapshot::new(&coords(6)), Snapshot::new(&[])],
    ));
    assert!(traj.is_empty());
    // Extra coordinates are kept, as in Python.
    traj.add_snapshot(&Snapshot::new(&coords(7))).unwrap();
    assert_eq!(
        traj.get_snapshot(0).unwrap().get_point_3d(1).unwrap(),
        Point3D::new(3.0, 4.0, 5.0)
    );
    // With no points, a snapshot needs no coordinates.
    let mut empty = Trajectory::new(3, 0);
    empty.add_snapshot(&Snapshot::new(&[])).unwrap();
    runtime(empty.get_snapshot(0).unwrap().get_point_3d(0));
    let traj = Trajectory::new_with_snapshot_list(2, 1, &[Snapshot::new(&[1.0, 2.0])]).unwrap();
    assert_eq!(traj.len(), 1);
}

#[test]
fn read_amber() {
    let bad = write("water_coords_bad.trx", AMBER_BAD);
    // RDKit's readers need three dimensions.
    let mut traj = Trajectory::new(2, 0);
    runtime(read_amber_trajectory(&bad, &mut traj));
    let mut traj = Trajectory::new(3, 3);
    value(read_amber_trajectory(&bad, &mut traj));
    let bad2 = write("water_coords_bad2.trx", AMBER_BAD2);
    value(read_amber_trajectory(&bad2, &mut traj));
    assert!(traj.is_empty());
    let one = write("water_coords.trx", AMBER);
    assert_eq!(read_amber_trajectory(&one, &mut traj).unwrap(), 1);
    assert_eq!(traj.len(), 1);
    let snapshot = traj.get_snapshot(0).unwrap();
    assert_eq!(snapshot.get_energy(), 0.0);
    assert_eq!(
        snapshot.get_point_3d(2).unwrap(),
        Point3D::new(0.1349280, 1.0427680, -0.4088006)
    );
    let two = write("water_coords2.trx", &format!("{AMBER}{}", &AMBER[31..]));
    let mut traj = Trajectory::new(3, 3);
    assert_eq!(read_amber_trajectory(&two, &mut traj).unwrap(), 2);
    assert_eq!(traj.len(), 2);
    // Frames are appended.
    assert_eq!(read_amber_trajectory(&one, &mut traj).unwrap(), 1);
    assert_eq!(traj.len(), 3);
    for path in [bad, bad2, one, two] {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn read_amber_quirks() {
    let mut traj = Trajectory::new(3, 3);
    // RDKit keeps a frame only when text follows its last coordinate.
    let unterminated = write("unterminated.trx", AMBER.trim_end());
    assert_eq!(read_amber_trajectory(&unterminated, &mut traj).unwrap(), 0);
    // A last frame missing only its last coordinate is dropped.
    let short = AMBER.trim_end().strip_suffix("-0.4088006").unwrap();
    let short = write("short.trx", &format!("{short}\n"));
    assert_eq!(read_amber_trajectory(&short, &mut traj).unwrap(), 0);
    assert!(traj.is_empty());
    // With no points RDKit's reader never ends.
    let mut none = Trajectory::new(3, 0);
    value(read_amber_trajectory(&unterminated, &mut none));
    // Three times the points overflows RDKit's unsigned int.
    let mut huge = Trajectory::new(3, u32::MAX / 3 + 1);
    value(read_amber_trajectory(&unterminated, &mut huge));
    let mut traj = Trajectory::new(3, 3);
    assert!(matches!(
        read_amber_trajectory(common::temp_path("missing.trx"), &mut traj),
        Err(Error::Os { .. })
    ));
    for path in [unterminated, short] {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn read_gromos() {
    let bad = write("water_coords_bad.trc", GROMOS_BAD);
    let mut traj = Trajectory::new(2, 0);
    runtime(read_gromos_trajectory(&bad, &mut traj));
    let mut traj = Trajectory::new(3, 3);
    value(read_gromos_trajectory(&bad, &mut traj));
    let bad2 = write("water_coords_bad2.trc", GROMOS_BAD2);
    value(read_gromos_trajectory(&bad2, &mut traj));
    assert!(traj.is_empty());
    let one = write(
        "water_coords.trc",
        &format!("{GROMOS_HEADER}{GROMOS_FRAME}"),
    );
    assert_eq!(read_gromos_trajectory(&one, &mut traj).unwrap(), 1);
    assert_eq!(traj.len(), 1);
    // Nanometers become Angstrom.
    let point = traj.get_snapshot(0).unwrap().get_point_3d(0).unwrap();
    assert_close(point.x, 1.941767);
    let two = write(
        "water_coords2.trc",
        &format!("{GROMOS_HEADER}{GROMOS_FRAME}{GROMOS_FRAME}"),
    );
    let mut traj = Trajectory::new(3, 3);
    assert_eq!(read_gromos_trajectory(&two, &mut traj).unwrap(), 2);
    assert_eq!(traj.len(), 2);
    // Another block is refused.
    let other = write("other.trc", "TITLE\nx\nEND\nOTHER\nEND\n");
    value(read_gromos_trajectory(&other, &mut traj));
    // Comments in a positions block are skipped.
    let commented = write("commented.trc", "POSITION\n# atoms\n0.1 0.2 0.3\nEND\n");
    let mut traj = Trajectory::new(3, 1);
    assert_eq!(read_gromos_trajectory(&commented, &mut traj).unwrap(), 1);
    // With no points each positions block is an empty snapshot.
    let empty = write("empty.trc", "POSITION\nEND\nPOSITION\nEND\n");
    let mut none = Trajectory::new(3, 0);
    assert_eq!(read_gromos_trajectory(&empty, &mut none).unwrap(), 2);
    let mut huge = Trajectory::new(3, u32::MAX / 3 + 1);
    value(read_gromos_trajectory(&empty, &mut huge));
    assert!(matches!(
        read_gromos_trajectory(common::temp_path("missing.trc"), &mut traj),
        Err(Error::Os { .. })
    ));
    for path in [bad, bad2, one, two, other, commented, empty] {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn add_conformers_from_minimization() {
    let mut mol = ethanol();
    let snapshots = {
        let mut ff = uff_get_molecule_force_field(&mut mol).unwrap();
        let (res, snapshots) = ff.minimize_trajectory(10).unwrap();
        assert_eq!(res, 0);
        snapshots
    };
    let traj = Trajectory::new_with_snapshot_list(3, mol.get_num_atoms(), &snapshots).unwrap();
    let last = snapshots.last().unwrap().get_point_3d(0).unwrap();
    let position = mol.get_conformer().unwrap().get_atom_position(0).unwrap();
    mol.remove_conformer(0);
    assert_eq!(
        traj.add_conformers_to_mol(&mut mol).unwrap() as usize,
        snapshots.len()
    );
    assert_eq!(mol.get_num_conformers() as usize, snapshots.len());
    let conf = mol
        .get_conformer_with_id(snapshots.len() as i32 - 1)
        .unwrap();
    assert_eq!(conf.get_atom_position(0).unwrap(), last);
    // The last snapshot holds the final coordinates.
    assert_close(last.x, position.x);
    let mut traj = traj;
    traj.clear();
    let n1 = mol.get_num_conformers();
    assert_eq!(traj.add_conformers_to_mol(&mut mol).unwrap(), 0);
    assert_eq!(mol.get_num_conformers(), n1);
    assert_eq!(traj.get_snapshot(0), Err(Error::Index { index: 0, len: 0 }));
}

fn add_conformers_from_file(traj: &mut Trajectory, x: f64, z: f64) {
    let mut mol = mol_from_smiles("CCC").unwrap();
    for i in 0..2 {
        assert_eq!(traj.add_conformers_to_mol(&mut mol).unwrap(), 1);
        assert_eq!(mol.get_num_conformers(), i + 1);
        let conf = mol.get_conformer_with_id(i as i32).unwrap();
        assert_eq!(conf.get_num_atoms(), 3);
        assert_close(conf.get_atom_position(0).unwrap().x, x);
        assert_close(conf.get_atom_position(2).unwrap().z, z);
    }
    mol.remove_all_conformers();
    let from_1 = AddConformersToMolParams {
        from_cid: 1,
        ..Default::default()
    };
    runtime(traj.add_conformers_to_mol_with_params(&mut mol, &from_1));
    assert_eq!(mol.get_num_conformers(), 0);
}

#[test]
fn add_conformers_from_trajectory_files() {
    let amber = write("conformers.trx", AMBER);
    let mut traj = Trajectory::new(3, 3);
    read_amber_trajectory(&amber, &mut traj).unwrap();
    add_conformers_from_file(&mut traj, 0.1941767, -0.4088006);
    let gromos = write("conformers.trc", &format!("{GROMOS_HEADER}{GROMOS_FRAME}"));
    let mut traj = Trajectory::new(3, 3);
    read_gromos_trajectory(&gromos, &mut traj).unwrap();
    add_conformers_from_file(&mut traj, 1.941767, -4.088006);

    let two = write(
        "conformers2.trc",
        &format!("{GROMOS_HEADER}{GROMOS_FRAME}{GROMOS_FRAME}"),
    );
    let mut traj = Trajectory::new(3, 3);
    read_gromos_trajectory(&two, &mut traj).unwrap();
    let mut mol = mol_from_smiles("CCC").unwrap();
    assert_eq!(traj.add_conformers_to_mol(&mut mol).unwrap(), 2);
    assert_eq!(mol.get_num_conformers(), 2);
    mol.remove_all_conformers();
    let first = AddConformersToMolParams {
        from_cid: 0,
        to_cid: 0,
    };
    assert_eq!(
        traj.add_conformers_to_mol_with_params(&mut mol, &first)
            .unwrap(),
        1
    );
    assert_eq!(mol.get_num_conformers(), 1);
    let from_1 = AddConformersToMolParams {
        from_cid: 1,
        ..Default::default()
    };
    assert_eq!(
        traj.add_conformers_to_mol_with_params(&mut mol, &from_1)
            .unwrap(),
        1
    );
    assert_eq!(mol.get_num_conformers(), 2);
    for path in [amber, gromos, two] {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn add_conformers_refusals() {
    let mut traj = Trajectory::new(3, 3);
    for i in 0..3 {
        traj.add_snapshot(&Snapshot::new_with_energy(&coords(9), f64::from(i)))
            .unwrap();
    }
    // The point count must be the atom count.
    let mut ethane = mol_from_smiles("CC").unwrap();
    runtime(traj.add_conformers_to_mol(&mut ethane));
    let mut mol = mol_from_smiles("CCC").unwrap();
    for (from_cid, to_cid) in [(3, -1), (-1, 3), (2, 1)] {
        runtime(traj.add_conformers_to_mol_with_params(
            &mut mol,
            &AddConformersToMolParams { from_cid, to_cid },
        ));
    }
    // Any negative bound means the first or the last snapshot.
    let all = AddConformersToMolParams {
        from_cid: -5,
        to_cid: -7,
    };
    assert_eq!(
        traj.add_conformers_to_mol_with_params(&mut mol, &all)
            .unwrap(),
        3
    );
    // New conformers take the next ids; none past i32::MAX.
    let mut mol = mol_from_smiles("CCC").unwrap();
    let mut conf = Conformer::new(3);
    conf.set_id(i32::MAX as u32 - 2).unwrap();
    mol.add_conformer(&conf).unwrap();
    value(traj.add_conformers_to_mol(&mut mol));
    assert_eq!(mol.get_num_conformers(), 1);
    let last_two = AddConformersToMolParams {
        from_cid: 1,
        to_cid: 2,
    };
    assert_eq!(
        traj.add_conformers_to_mol_with_params(&mut mol, &last_two)
            .unwrap(),
        2
    );
    assert!(mol.get_conformer_with_id(i32::MAX).is_ok());
    // Points that cannot be read add no conformer.
    let mut line = Trajectory::new(1, 3);
    line.add_snapshot(&Snapshot::new(&coords(3))).unwrap();
    let mut mol = mol_from_smiles("CCC").unwrap();
    runtime(line.add_conformers_to_mol(&mut mol));
    assert_eq!(mol.get_num_conformers(), 0);
    // A molecule with no atoms takes conformers from any snapshot.
    let mut none = Trajectory::new(0, 0);
    none.add_snapshot(&Snapshot::new(&[])).unwrap();
    let mut empty = Mol::default();
    assert_eq!(none.add_conformers_to_mol(&mut empty).unwrap(), 1);
    // An empty trajectory adds nothing, whatever the negative bounds.
    let mut mol = mol_from_smiles("CCC").unwrap();
    assert_eq!(
        Trajectory::new(3, 3)
            .add_conformers_to_mol_with_params(&mut mol, &all)
            .unwrap(),
        0
    );
}
