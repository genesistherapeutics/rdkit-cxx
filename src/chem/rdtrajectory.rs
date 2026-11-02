// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdtrajectory`: trajectories of coordinate snapshots,
//! and RDKit's AMBER and GROMOS trajectory readers.
//!
//! [`Trajectory`] and [`Snapshot`] are Rust values. RDKit's snapshots point
//! back at the trajectory that holds them, which Python lets dangle; a Rust
//! snapshot copies the trajectory's dimension and point count instead.

use std::path::Path;

use crate::chem::{Conformer, Mol};
use crate::ffi::trajectory::ffi;
use crate::geometry::{Point2D, Point3D};
use crate::{Error, Result};

fn runtime(message: &str) -> Error {
    Error::Runtime {
        message: message.to_owned(),
    }
}

/// The coordinates and energy of one frame of a trajectory. Mirrors
/// `rdtrajectory.Snapshot`.
///
/// A snapshot reads its points as the trajectory it was taken from lays them
/// out: [`Trajectory::get_snapshot`] returns one that does, and so do its
/// clones, as in Python. A snapshot made with [`Snapshot::new`] belongs to no
/// trajectory, so its point getters are `Error::Runtime`, as in Python, until
/// a trajectory holds it. Snapshots from
/// [`ForceField::minimize_trajectory`](crate::force_field::rd_force_field::ForceField::minimize_trajectory)
/// read as the force field lays out its points; Python's belong to no
/// trajectory, so their point getters always raise.
#[derive(Clone, Debug, PartialEq)]
pub struct Snapshot {
    coords: Vec<f64>,
    energy: f64,
    // The dimension and point count of the trajectory the snapshot was taken
    // from (RDKit's `d_trajectory`), whose layout `coords` is long enough
    // for.
    layout: Option<(u32, u32)>,
}

impl Snapshot {
    /// A snapshot of `coord_list`, with energy 0, in no trajectory. Mirrors
    /// `Snapshot(coordList)`.
    pub fn new(coord_list: &[f64]) -> Snapshot {
        Snapshot::new_with_energy(coord_list, 0.0)
    }

    /// Mirrors `Snapshot(coordList, energy)`.
    pub fn new_with_energy(coord_list: &[f64], energy: f64) -> Snapshot {
        Snapshot {
            coords: coord_list.to_vec(),
            energy,
            layout: None,
        }
    }

    /// A snapshot laid out as a trajectory of `dimension` and `num_points`,
    /// whose `dimension * num_points` coordinates `coords` holds.
    pub(crate) fn with_layout(
        coords: Vec<f64>,
        energy: f64,
        dimension: u32,
        num_points: u32,
    ) -> Snapshot {
        debug_assert!(coords.len() as u64 >= u64::from(dimension) * u64::from(num_points));
        Snapshot {
            coords,
            energy,
            layout: Some((dimension, num_points)),
        }
    }

    // RDKit's checks in `Snapshot::getPoint2D` and `getPoint3D`, in its
    // order; returns the dimension and the first coordinate's index.
    fn point_start(
        &self,
        point_num: u32,
        dimension_ok: fn(u32) -> bool,
        dimension_message: &str,
    ) -> Result<(u32, usize)> {
        if self.coords.is_empty() {
            return Err(runtime("d_pos must not be NULL"));
        }
        let Some((dimension, num_points)) = self.layout else {
            return Err(runtime("d_trajectory must not be NULL"));
        };
        if !dimension_ok(dimension) {
            return Err(runtime(dimension_message));
        }
        if num_points == 0 {
            return Err(runtime("d_numPoints must be > 0"));
        }
        if point_num >= num_points {
            return Err(Error::index(point_num, num_points as usize));
        }
        Ok((dimension, point_num as usize * dimension as usize))
    }

    /// The coordinates of point `point_num`. A snapshot with no coordinates
    /// or in no trajectory, or from a trajectory whose dimension is not 2, is
    /// `Error::Runtime`; a point past the trajectory's is `Error::Index`.
    /// Mirrors `Snapshot.GetPoint2D`.
    pub fn get_point_2d(&self, point_num: u32) -> Result<Point2D> {
        let (_, i) = self.point_start(point_num, |d| d == 2, "d_dimension must be == 2")?;
        Ok(Point2D::new(self.coords[i], self.coords[i + 1]))
    }

    /// The coordinates of point `point_num`; z is 0 unless the trajectory's
    /// dimension is exactly 3 (RDKit reads only x and y of a point with more
    /// dimensions). A snapshot with no coordinates or in no trajectory, or
    /// from a trajectory whose dimension is below 2, is `Error::Runtime`; a
    /// point past the trajectory's is `Error::Index`. Mirrors
    /// `Snapshot.GetPoint3D`.
    pub fn get_point_3d(&self, point_num: u32) -> Result<Point3D> {
        let (dimension, i) = self.point_start(point_num, |d| d >= 2, "d_dimension must be >= 2")?;
        let z = if dimension == 3 {
            self.coords[i + 2]
        } else {
            0.0
        };
        Ok(Point3D::new(self.coords[i], self.coords[i + 1], z))
    }

    /// Mirrors `Snapshot.GetEnergy`.
    pub fn get_energy(&self) -> f64 {
        self.energy
    }

    /// Mirrors `Snapshot.SetEnergy`.
    pub fn set_energy(&mut self, energy: f64) {
        self.energy = energy;
    }
}

/// Snapshots of `num_points` coordinate tuples of `dimension` values each.
/// Mirrors `rdtrajectory.Trajectory`; `Clone` mirrors `Trajectory(other)`.
#[derive(Clone, Debug, PartialEq)]
pub struct Trajectory {
    dimension: u32,
    num_points: u32,
    snapshots: Vec<Snapshot>,
}

/// Parameters for `Trajectory::add_conformers_to_mol_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AddConformersToMolParams {
    /// The first snapshot to add; negative for the first.
    pub from_cid: i32,
    /// The last snapshot to add; negative for the last.
    pub to_cid: i32,
}

impl Default for AddConformersToMolParams {
    fn default() -> AddConformersToMolParams {
        AddConformersToMolParams {
            from_cid: -1,
            to_cid: -1,
        }
    }
}

impl Trajectory {
    /// An empty trajectory. Mirrors `Trajectory(dimension, numPoints)`.
    pub fn new(dimension: u32, num_points: u32) -> Trajectory {
        Trajectory {
            dimension,
            num_points,
            snapshots: Vec::new(),
        }
    }

    /// A trajectory holding copies of `snapshot_list`, which
    /// [`add_snapshot`](Trajectory::add_snapshot) checks. Mirrors
    /// `Trajectory(dimension, numPoints, snapshotList)`.
    pub fn new_with_snapshot_list(
        dimension: u32,
        num_points: u32,
        snapshot_list: &[Snapshot],
    ) -> Result<Trajectory> {
        let mut trajectory = Trajectory::new(dimension, num_points);
        for snapshot in snapshot_list {
            trajectory.add_snapshot(snapshot)?;
        }
        Ok(trajectory)
    }

    /// Mirrors `Trajectory.Dimension`.
    pub fn dimension(&self) -> u32 {
        self.dimension
    }

    /// Mirrors `Trajectory.NumPoints`.
    pub fn num_points(&self) -> u32 {
        self.num_points
    }

    /// The number of snapshots. Mirrors `len(Trajectory)`.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Whether the trajectory holds no snapshots. A Rust addition.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    // A copy of `s` laid out as this trajectory. RDKit reads a snapshot's
    // points from its array without a length, so Python reads past the end
    // of one with fewer than dimension * numPoints coordinates.
    fn adopt(&self, s: &Snapshot) -> Result<Snapshot> {
        let needed = u64::from(self.dimension) * u64::from(self.num_points);
        if (s.coords.len() as u64) < needed {
            return Err(Error::Value {
                message: format!(
                    "the snapshot has {} coordinates; the trajectory needs {needed}",
                    s.coords.len()
                ),
            });
        }
        // Python returns positions as unsigned int.
        if self.snapshots.len() >= u32::MAX as usize {
            return Err(Error::Value {
                message: "the trajectory holds as many snapshots as an index can name".to_owned(),
            });
        }
        Ok(Snapshot {
            layout: Some((self.dimension, self.num_points)),
            ..s.clone()
        })
    }

    /// Appends a copy of `s`, returning its index. A snapshot with fewer than
    /// `dimension * num_points` coordinates is `Error::Value`; Python adds it
    /// and reads past its end. Mirrors `Trajectory.AddSnapshot`.
    pub fn add_snapshot(&mut self, s: &Snapshot) -> Result<u32> {
        let snapshot = self.adopt(s)?;
        self.snapshots.push(snapshot);
        Ok((self.snapshots.len() - 1) as u32)
    }

    /// The snapshot at `snapshot_num`, or `Error::Index`. Mirrors
    /// `Trajectory.GetSnapshot`, which returns a copy.
    pub fn get_snapshot(&self, snapshot_num: u32) -> Result<&Snapshot> {
        self.snapshots
            .get(snapshot_num as usize)
            .ok_or_else(|| Error::index(snapshot_num, self.snapshots.len()))
    }

    /// Inserts a copy of `s` before the snapshot at `snapshot_num` (at the
    /// end for the length), returning its index. A position past the end is
    /// `Error::Index`, with the length plus one, RDKit's bound; a snapshot
    /// with too few coordinates is `Error::Value`, as in `add_snapshot`.
    /// Mirrors `Trajectory.InsertSnapshot`.
    pub fn insert_snapshot(&mut self, snapshot_num: u32, s: &Snapshot) -> Result<u32> {
        let index = snapshot_num as usize;
        if index > self.snapshots.len() {
            return Err(Error::index(snapshot_num, self.snapshots.len() + 1));
        }
        let snapshot = self.adopt(s)?;
        self.snapshots.insert(index, snapshot);
        Ok(snapshot_num)
    }

    /// Removes the snapshot at `snapshot_num`, or is `Error::Index`, and
    /// returns the index of the snapshot after it (the new length when it
    /// was the last). Mirrors `Trajectory.RemoveSnapshot`.
    pub fn remove_snapshot(&mut self, snapshot_num: u32) -> Result<u32> {
        let index = snapshot_num as usize;
        if index >= self.snapshots.len() {
            return Err(Error::index(snapshot_num, self.snapshots.len()));
        }
        self.snapshots.remove(index);
        Ok(snapshot_num)
    }

    /// Removes every snapshot. Mirrors `Trajectory.Clear`.
    pub fn clear(&mut self) {
        self.snapshots.clear();
    }

    /// Adds a conformer to `mol` for each snapshot, returning how many.
    /// Mirrors `Trajectory.AddConformersToMol(mol)`; see
    /// [`add_conformers_to_mol_with_params`](Trajectory::add_conformers_to_mol_with_params).
    pub fn add_conformers_to_mol(&self, mol: &mut Mol) -> Result<u32> {
        self.add_conformers_to_mol_with_params(mol, &AddConformersToMolParams::default())
    }

    /// Adds a conformer to `mol`, with a new id, for each snapshot from
    /// `from_cid` to `to_cid`, both included, returning how many. Points are
    /// read as [`Snapshot::get_point_3d`] reads them.
    ///
    /// `Error::Runtime` when the trajectory's point count is not the
    /// molecule's atom count, when either bound is not below the number of
    /// snapshots, when `from_cid` is past `to_cid` (once negative bounds
    /// become the first and last snapshots), or when the snapshots' points
    /// cannot be read (a dimension below 2). `Error::Value` when the
    /// molecule has too few conformer ids left (see
    /// [`Mol::add_conformer_with_assign_id`]). Everything is checked before a
    /// conformer is added, so an error leaves the molecule unchanged, as
    /// Python's does for every error it raises here. Mirrors
    /// `Trajectory.AddConformersToMol(mol, fromCid, toCid)`.
    pub fn add_conformers_to_mol_with_params(
        &self,
        mol: &mut Mol,
        params: &AddConformersToMolParams,
    ) -> Result<u32> {
        let num_atoms = mol.get_num_atoms();
        if self.num_points != num_atoms {
            return Err(runtime(
                "Number of atom mismatch between ROMol and Trajectory",
            ));
        }
        let size = self.snapshots.len() as i64;
        let (from, to) = (i64::from(params.from_cid), i64::from(params.to_cid));
        if from >= size {
            return Err(runtime("from must be < size()"));
        }
        if to >= size {
            return Err(runtime("to must be < size()"));
        }
        let from = from.max(0);
        let to = if to < 0 { size - 1 } else { to };
        if size == 0 {
            return Ok(0);
        }
        if from > to {
            return Err(runtime("from must be <= to"));
        }
        let snapshots = &self.snapshots[from as usize..=to as usize];
        // Every snapshot is laid out as the trajectory, so the first point
        // reads as every other does.
        if num_atoms > 0 {
            snapshots[0].get_point_3d(0)?;
        }
        // RDKit gives each new conformer one more than the largest id.
        let max_id = mol
            .get_conformers()
            .map(|conf| i64::from(conf.get_id()))
            .max()
            .unwrap_or(-1);
        if max_id + snapshots.len() as i64 > i64::from(i32::MAX) {
            return Err(Error::Value {
                message: "no conformer id is left to assign".to_owned(),
            });
        }
        for snapshot in snapshots {
            let mut conf = Conformer::new(num_atoms);
            for i in 0..num_atoms {
                conf.set_atom_position(i, snapshot.get_point_3d(i)?)?;
            }
            mol.add_conformer_with_assign_id(&conf, true)?;
        }
        Ok(snapshots.len() as u32)
    }
}

// The checks Python's readers make before reading; RDKit's readers fill a
// trajectory of dimension 3.
fn read_trajectory(
    traj: &mut Trajectory,
    read: impl FnOnce(u32, &mut Vec<f64>) -> std::result::Result<u32, cxx::Exception>,
) -> Result<u32> {
    if traj.dimension != 3 {
        return Err(runtime("The trajectory must have dimension == 3"));
    }
    let mut coords = Vec::new();
    let count = read(traj.num_points, &mut coords)?;
    // Python returns positions as unsigned int.
    if traj.snapshots.len() as u64 + u64::from(count) > u64::from(u32::MAX) {
        return Err(Error::Value {
            message: "the trajectory holds as many snapshots as an index can name".to_owned(),
        });
    }
    let size = 3 * traj.num_points as usize;
    traj.snapshots.extend((0..count as usize).map(|i| {
        Snapshot::with_layout(
            coords[i * size..(i + 1) * size].to_vec(),
            0.0,
            3,
            traj.num_points,
        )
    }));
    Ok(count)
}

/// Appends the frames of the AMBER trajectory file `f_name` (a title line,
/// then each frame's `num_points * 3` coordinates, separated by whitespace)
/// to `traj`, returning how many. Mirrors `rdtrajectory.ReadAmberTrajectory`,
/// quirks included: a frame is kept only when more text (a newline will do)
/// follows its last coordinate, and a last frame missing only its last
/// coordinate is dropped without an error.
///
/// `Error::Runtime` for a trajectory whose dimension is not 3. `Error::Value`
/// for one with no points (RDKit's reader never ends) or more than
/// `u32::MAX / 3` (its arrays' sizes overflow), for text that is not a
/// number, and for a file that ends inside a frame. `Error::Os` for a file
/// that cannot be opened. On an error `traj` is unchanged; Python keeps the
/// frames read before it.
pub fn read_amber_trajectory(f_name: impl AsRef<Path>, traj: &mut Trajectory) -> Result<u32> {
    let path = f_name.as_ref().as_os_str().as_encoded_bytes();
    read_trajectory(traj, |num_points, coords| {
        ffi::read_amber_trajectory(path, num_points, coords)
    })
}

/// Appends the `POSITION` and `POSITIONRED` blocks of the GROMOS trajectory
/// file `f_name` to `traj`, in Angstrom (GROMOS writes nanometers), and
/// returns how many. Lines holding `#` inside a block are comments, and
/// `TITLE`, `TIMESTEP`, `VELOCITY`, `VELOCITYRED`, `GENBOX`, and `BOX`
/// blocks are skipped. Mirrors `rdtrajectory.ReadGromosTrajectory`.
///
/// `Error::Runtime` for a trajectory whose dimension is not 3. `Error::Value`
/// for one with more than `u32::MAX / 3` points (RDKit's arrays' sizes
/// overflow), for another block, for a position line without three numbers,
/// and for a block without `num_points` positions followed by `END`.
/// `Error::Os` for a file that cannot be opened or read. On an error `traj`
/// is unchanged; Python keeps the frames read before it.
pub fn read_gromos_trajectory(f_name: impl AsRef<Path>, traj: &mut Trajectory) -> Result<u32> {
    let path = f_name.as_ref().as_os_str().as_encoded_bytes();
    read_trajectory(traj, |num_points, coords| {
        ffi::read_gromos_trajectory(path, num_points, coords)
    })
}
