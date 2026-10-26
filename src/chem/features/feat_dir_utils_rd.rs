// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Features.FeatDirUtilsRD`, a pure-Python module: the
//! direction vectors of donors, acceptors, and aromatic rings in a
//! conformer.
//!
//! Each `get_*_feat_vects` function returns the feature's directions, each
//! as a start (the feature's atom or center) and an end `scale` along it,
//! with how to draw them: `"linear"`, or `"cone"` where the direction is
//! uncertain. A feature whose atoms do not have the neighbors the function
//! assumes is `Error::Value` (Python's `AssertionError`).

use crate::chem::rdchem::{AtomRef, BondType, ConformerRef};
use crate::geometry::rd_geometry::Point3D;
use crate::{Error, Result};

/// A feature's directions, each a start and an end, and how to draw them.
pub type FeatVects = (Vec<(Point3D, Point3D)>, &'static str);

fn assertion(message: &str) -> Error {
    Error::Value {
        message: message.to_owned(),
    }
}

/// The cross product of two vectors, as Python's NumPy expression computes
/// it. Mirrors `FeatDirUtilsRD.cross`.
pub fn cross(v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    [
        v1[1] * v2[2] - v1[2] * v2[1],
        -v1[0] * v2[2] + v1[2] * v2[0],
        v1[0] * v2[1] - v1[1] * v2[0],
    ]
}

/// The atoms bonded to `atom_id` in an adjacency matrix. Mirrors
/// `FeatDirUtilsRD.findNeighbors`.
pub fn find_neighbors(atom_id: usize, adj_mat: &[Vec<i32>]) -> Result<Vec<usize>> {
    let row = adj_mat.get(atom_id).ok_or(Error::Index {
        index: atom_id as u32,
        len: adj_mat.len() as u32,
    })?;
    Ok(row
        .iter()
        .enumerate()
        .filter(|&(_, &nid)| nid >= 1)
        .map(|(i, _)| i)
        .collect())
}

/// The normalized mean of the unit vectors from `center` to the atoms.
fn find_avg_vec(conf: &ConformerRef<'_>, center: Point3D, nbrs: &[AtomRef<'_>]) -> Result<Point3D> {
    let mut avg_vec: Option<Point3D> = None;
    for nbr in nbrs {
        let mut pt = conf.get_atom_position(nbr.get_idx())?;
        pt -= center;
        pt.normalize()?;
        avg_vec = Some(match avg_vec {
            None => pt,
            Some(v) => v + pt,
        });
    }
    let mut avg_vec = avg_vec.ok_or_else(|| assertion("no neighbors to average"))?;
    avg_vec.normalize()?;
    Ok(avg_vec)
}

/// The aromatic ring's normals at `feat_loc`, each way. Mirrors
/// `FeatDirUtilsRD.GetAromaticFeatVects`.
pub fn get_aromatic_feat_vects(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    feat_loc: Point3D,
) -> Result<FeatVects> {
    get_aromatic_feat_vects_with_scale(conf, feat_atoms, feat_loc, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetAromaticFeatVects(conf, featAtoms, featLoc, scale)`.
pub fn get_aromatic_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    feat_loc: Point3D,
    scale: f64,
) -> Result<FeatVects> {
    let head = feat_loc;
    if feat_atoms.len() < 2 {
        return Err(Error::Index {
            index: feat_atoms.len() as u32,
            len: feat_atoms.len() as u32,
        });
    }
    let v1 = conf.get_atom_position(feat_atoms[0])? - head;
    let v2 = conf.get_atom_position(feat_atoms[1])? - head;
    let mut norm1 = v1.cross_product(&v2);
    norm1.normalize()?;
    norm1 *= scale;
    let norm2 = head - norm1;
    norm1 += head;
    Ok((vec![(head, norm1), (head, norm2)], "linear"))
}

/// The rotation matrix for `theta` degrees about the unit axis `ax`.
fn rotation_matrix(theta: f64, ax: Point3D) -> [[f64; 3]; 3] {
    let theta = std::f64::consts::PI * theta / 180.0;
    let (c, s) = (theta.cos(), theta.sin());
    let t = 1.0 - c;
    let (x, y, z) = (ax.x, ax.y, ax.z);
    [
        [t * x * x + c, t * x * y + s * z, t * x * z - s * y],
        [t * x * y - s * z, t * y * y + c, t * y * z + s * x],
        [t * x * z + s * y, t * y * z - s * x, t * z * z + c],
    ]
}

fn rotate(mat: &[[f64; 3]; 3], p: Point3D) -> Point3D {
    let row = |r: &[f64; 3]| r[0] * p.x + r[1] * p.y + r[2] * p.z;
    Point3D::new(row(&mat[0]), row(&mat[1]), row(&mat[2]))
}

/// `pt` rotated `theta` degrees about the unit axis `ax`. NumPy's matrix
/// product may round differently in the last place. Mirrors
/// `FeatDirUtilsRD.ArbAxisRotation`.
pub fn arb_axis_rotation(theta: f64, ax: Point3D, pt: Point3D) -> Point3D {
    rotate(&rotation_matrix(theta, ax), pt)
}

/// Each of the points, rotated. Mirrors `FeatDirUtilsRD.ArbAxisRotation`
/// with a list of points.
pub fn arb_axis_rotation_list(theta: f64, ax: Point3D, pts: &[Point3D]) -> Vec<Point3D> {
    let mat = rotation_matrix(theta, ax);
    pts.iter().map(|&p| rotate(&mat, p)).collect()
}

/// The feature's one atom.
fn single_atom(feat_atoms: &[u32]) -> Result<u32> {
    match feat_atoms {
        [aid] => Ok(*aid),
        _ => Err(assertion("the feature must have one atom")),
    }
}

/// An acceptor with two heavy neighbors: the direction opposite their
/// mean, or for an oxygen, two directions tilted 54.5 degrees from it,
/// as for sp3. Mirrors `FeatDirUtilsRD.GetAcceptor2FeatVects`.
pub fn get_acceptor2_feat_vects(conf: &ConformerRef<'_>, feat_atoms: &[u32]) -> Result<FeatVects> {
    get_acceptor2_feat_vects_with_scale(conf, feat_atoms, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetAcceptor2FeatVects(conf, featAtoms, scale)`.
pub fn get_acceptor2_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    scale: f64,
) -> Result<FeatVects> {
    let aid = single_atom(feat_atoms)?;
    let cpt = conf.get_atom_position(aid)?;
    let mol = conf.get_owning_mol();
    let atom = mol.get_atom_with_idx(aid)?;
    // Python pops the neighbors from the end.
    let nbrs: Vec<AtomRef<'_>> = atom
        .get_neighbors()
        .into_iter()
        .rev()
        .filter(|n| n.get_atomic_num() != 1)
        .collect();
    if nbrs.len() != 2 {
        return Err(assertion("the acceptor must have two heavy neighbors"));
    }
    let mut bvec = find_avg_vec(conf, cpt, &nbrs)?;
    bvec *= -scale;
    if atom.get_atomic_num() == 8 {
        let v1 = conf.get_atom_position(nbrs[0].get_idx())? - cpt;
        let v2 = conf.get_atom_position(nbrs[1].get_idx())? - cpt;
        let mut rot_axis = v1 - v2;
        rot_axis.normalize()?;
        let bv1 = arb_axis_rotation(54.5, rot_axis, bvec) + cpt;
        let bv2 = arb_axis_rotation(-54.5, rot_axis, bvec) + cpt;
        return Ok((vec![(cpt, bv1), (cpt, bv2)], "linear"));
    }
    bvec += cpt;
    Ok((vec![(cpt, bvec)], "linear"))
}

/// Whether the vectors to three neighbors are coplanar within `tol`.
fn check_planarity(
    conf: &ConformerRef<'_>,
    cpt: Point3D,
    nbrs: &[AtomRef<'_>],
    tol: f64,
) -> Result<bool> {
    if nbrs.len() != 3 {
        return Err(assertion("the atom must have three neighbors"));
    }
    let v = |i: usize| -> Result<Point3D> { Ok(conf.get_atom_position(nbrs[i].get_idx())? - cpt) };
    let normal = v(0)?.cross_product(&v(1)?);
    Ok(v(2)?.dot_product(&normal).abs() <= tol)
}

/// The fourth arm of a tetrahedral atom with three neighbors, unless they
/// are coplanar.
fn tetrahedral_feat_vect(
    conf: &ConformerRef<'_>,
    aid: u32,
    scale: f64,
) -> Result<Vec<(Point3D, Point3D)>> {
    let mol = conf.get_owning_mol();
    let cpt = conf.get_atom_position(aid)?;
    let nbrs = mol.get_atom_with_idx(aid)?.get_neighbors();
    if check_planarity(conf, cpt, &nbrs, 0.1)? {
        return Ok(Vec::new());
    }
    let mut bvec = find_avg_vec(conf, cpt, &nbrs)?;
    bvec *= -scale;
    bvec += cpt;
    Ok(vec![(cpt, bvec)])
}

/// A donor with three heavy neighbors: the fourth arm of their tetrahedron.
/// Mirrors `FeatDirUtilsRD.GetDonor3FeatVects`.
pub fn get_donor3_feat_vects(conf: &ConformerRef<'_>, feat_atoms: &[u32]) -> Result<FeatVects> {
    get_donor3_feat_vects_with_scale(conf, feat_atoms, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetDonor3FeatVects(conf, featAtoms, scale)`.
pub fn get_donor3_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    scale: f64,
) -> Result<FeatVects> {
    Ok((
        tetrahedral_feat_vect(conf, single_atom(feat_atoms)?, scale)?,
        "linear",
    ))
}

/// An acceptor with three heavy neighbors: the fourth arm of their
/// tetrahedron. Mirrors `FeatDirUtilsRD.GetAcceptor3FeatVects`.
pub fn get_acceptor3_feat_vects(conf: &ConformerRef<'_>, feat_atoms: &[u32]) -> Result<FeatVects> {
    get_acceptor3_feat_vects_with_scale(conf, feat_atoms, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetAcceptor3FeatVects(conf, featAtoms, scale)`.
pub fn get_acceptor3_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    scale: f64,
) -> Result<FeatVects> {
    Ok((
        tetrahedral_feat_vect(conf, single_atom(feat_atoms)?, scale)?,
        "linear",
    ))
}

fn toward(conf: &ConformerRef<'_>, cpt: Point3D, aid: u32, scale: f64) -> Result<Point3D> {
    let mut v = conf.get_atom_position(aid)? - cpt;
    v.normalize()?;
    v *= scale;
    Ok(v + cpt)
}

/// A donor with two heavy neighbors: with no hydrogens, the direction
/// opposite their mean; with one, toward it, and, unless the atom is
/// planar, opposite the heavy atoms too; with two or more, toward each.
/// Mirrors `FeatDirUtilsRD.GetDonor2FeatVects`.
pub fn get_donor2_feat_vects(conf: &ConformerRef<'_>, feat_atoms: &[u32]) -> Result<FeatVects> {
    get_donor2_feat_vects_with_scale(conf, feat_atoms, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetDonor2FeatVects(conf, featAtoms, scale)`.
pub fn get_donor2_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    scale: f64,
) -> Result<FeatVects> {
    let aid = single_atom(feat_atoms)?;
    let mol = conf.get_owning_mol();
    let cpt = conf.get_atom_position(aid)?;
    let nbrs = mol.get_atom_with_idx(aid)?.get_neighbors();
    if nbrs.len() < 2 {
        return Err(assertion("the donor must have two neighbors"));
    }
    let (hydrogens, heavy): (Vec<AtomRef<'_>>, Vec<AtomRef<'_>>) =
        nbrs.iter().partition(|n| n.get_atomic_num() == 1);
    match nbrs.len() {
        2 => {
            if !hydrogens.is_empty() {
                return Err(assertion("the donor must have no hydrogens"));
            }
            let mut bvec = find_avg_vec(conf, cpt, &heavy)?;
            bvec *= -scale;
            bvec += cpt;
            Ok((vec![(cpt, bvec)], "linear"))
        }
        3 => {
            if hydrogens.len() != 1 {
                return Err(assertion("the donor must have one hydrogen"));
            }
            let bvec = toward(conf, cpt, hydrogens[0].get_idx(), scale)?;
            if check_planarity(conf, cpt, &nbrs, 1.0e-2)? {
                return Ok((vec![(cpt, bvec)], "linear"));
            }
            let mut ovec = find_avg_vec(conf, cpt, &heavy)?;
            ovec *= -scale;
            ovec += cpt;
            Ok((vec![(cpt, bvec), (cpt, ovec)], "linear"))
        }
        _ => {
            let res = hydrogens
                .iter()
                .map(|h| Ok((cpt, toward(conf, cpt, h.get_idx(), scale)?)))
                .collect::<Result<Vec<_>>>()?;
            Ok((res, "linear"))
        }
    }
}

/// The first heavy neighbor of the atom, as Python looks for it.
fn first_heavy_neighbor<'m>(conf: &ConformerRef<'m>, aid: u32) -> Result<AtomRef<'m>> {
    conf.get_owning_mol()
        .get_atom_with_idx(aid)?
        .get_neighbors()
        .into_iter()
        .find(|n| n.get_atomic_num() != 1)
        .ok_or_else(|| assertion("the atom has no heavy neighbor"))
}

/// A donor with one heavy neighbor: the direction away from it, as a cone.
/// Mirrors `FeatDirUtilsRD.GetDonor1FeatVects`.
pub fn get_donor1_feat_vects(conf: &ConformerRef<'_>, feat_atoms: &[u32]) -> Result<FeatVects> {
    get_donor1_feat_vects_with_scale(conf, feat_atoms, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetDonor1FeatVects(conf, featAtoms, scale)`.
pub fn get_donor1_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    scale: f64,
) -> Result<FeatVects> {
    let aid = single_atom(feat_atoms)?;
    let hnbr = first_heavy_neighbor(conf, aid)?.get_idx();
    let cpt = conf.get_atom_position(aid)?;
    let v1 = toward(conf, cpt, hnbr, -scale)?;
    Ok((vec![(cpt, v1)], "cone"))
}

/// An acceptor with one heavy neighbor: bonded by more than a single bond,
/// or to a sulfur, the direction away from it, as a cone; otherwise two
/// directions 120 degrees from the bond, in the plane of the neighbor's
/// next neighbor. Mirrors `FeatDirUtilsRD.GetAcceptor1FeatVects`.
pub fn get_acceptor1_feat_vects(conf: &ConformerRef<'_>, feat_atoms: &[u32]) -> Result<FeatVects> {
    get_acceptor1_feat_vects_with_scale(conf, feat_atoms, 1.5)
}

/// With directions `scale` long. Mirrors
/// `FeatDirUtilsRD.GetAcceptor1FeatVects(conf, featAtoms, scale)`.
pub fn get_acceptor1_feat_vects_with_scale(
    conf: &ConformerRef<'_>,
    feat_atoms: &[u32],
    scale: f64,
) -> Result<FeatVects> {
    let aid = single_atom(feat_atoms)?;
    let mol = conf.get_owning_mol();
    let cpt = conf.get_atom_position(aid)?;
    let heavy_at = first_heavy_neighbor(conf, aid)?;
    let bond = mol
        .get_bond_between_atoms(aid, heavy_at.get_idx())?
        .expect("neighbors are bonded");
    // Python's `GetBondType() > BondType.SINGLE`: any type but those two.
    let multiple = !matches!(
        bond.get_bond_type(),
        BondType::Unspecified | BondType::Single
    );
    let sulfur = heavy_at.get_atomic_num() == 16;
    if multiple || sulfur {
        let v1 = toward(conf, cpt, heavy_at.get_idx(), -scale)?;
        return Ok((vec![(cpt, v1)], "cone"));
    }
    let hv_nbr = heavy_at
        .get_neighbors()
        .into_iter()
        .find(|n| n.get_idx() != aid)
        .ok_or_else(|| assertion("the neighbor has no other neighbor"))?;
    let v_heavy = conf.get_atom_position(heavy_at.get_idx())?;
    let pt1 = conf.get_atom_position(hv_nbr.get_idx())? - v_heavy;
    let v1 = v_heavy - cpt;
    let mut rot_axis = v1.cross_product(&pt1);
    rot_axis.normalize()?;
    let end = |theta: f64| -> Result<Point3D> {
        let mut bv = arb_axis_rotation(theta, rot_axis, v1);
        bv.normalize()?;
        bv *= scale;
        Ok(bv + cpt)
    };
    Ok((vec![(cpt, end(120.0)?), (cpt, end(-120.0)?)], "linear"))
}
