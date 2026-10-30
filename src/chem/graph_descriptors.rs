// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.GraphDescriptors`, a Python module: topological
//! descriptors (connectivity and shape indices, information content,
//! Balaban's J, and Bertz's complexity index).
//!
//! The connectivity indices of valence and Hall and Kier's kappas call
//! `rdkit_cxx::chem::rd_mol_descriptors`, as Python's do; `Chi0`, `Chi1`,
//! `Ipc`, `BalabanJ`, and `BertzCT` are computed here, as Python computes
//! them. Python keeps the distance and adjacency matrices it computes on
//! the molecule and reads them back unless forced; these never do, so each
//! call computes them afresh, as Python does for a new molecule. Python's
//! Python reference implementations (`_pyChi0v` and the like) are not
//! ported.

use std::collections::{BTreeMap, HashMap};
use std::hash::Hash;

use crate::chem::rd_mol_descriptors as descriptors;
use crate::chem::{
    BondType, GetDistanceMatrixParams, Mol, get_adjacency_matrix, get_distance_matrix_with_params,
};
use crate::ml::info_theory::rd_info_theory::info_entropy;
use crate::py::{np_sum, py_sum};
use crate::{Error, Result};

/// Hall and Kier's alpha values by element symbol, for sp, sp2, and sp3
/// atoms, with `None` where an element has none for that hybridization.
/// Mirrors `GraphDescriptors.hallKierAlphas`.
pub fn hall_kier_alphas() -> BTreeMap<String, [Option<f64>; 3]> {
    [
        ("Br", [None, None, Some(0.48)]),
        ("C", [Some(-0.22), Some(-0.13), Some(0.0)]),
        ("Cl", [None, None, Some(0.29)]),
        ("F", [None, None, Some(-0.07)]),
        ("H", [Some(0.0), Some(0.0), Some(0.0)]),
        ("I", [None, None, Some(0.73)]),
        ("N", [Some(-0.29), Some(-0.2), Some(-0.04)]),
        ("O", [None, Some(-0.2), Some(-0.04)]),
        ("P", [None, Some(0.3), Some(0.43)]),
        ("S", [None, Some(0.22), Some(0.35)]),
    ]
    .into_iter()
    .map(|(symbol, alphas)| (symbol.to_owned(), alphas))
    .collect()
}

/// A square matrix with one row per atom; `Error::Value` for another shape,
/// where Python's NumPy arithmetic would fail.
fn check_square(mat: &[Vec<f64>], n: usize) -> Result<()> {
    if mat.len() != n || mat.iter().any(|row| row.len() != n) {
        return Err(Error::Value {
            message: format!("the matrix must have {n} rows of {n} values, one per atom"),
        });
    }
    Ok(())
}

fn distance_matrix(mol: &Mol, use_bo: bool) -> Result<Vec<Vec<f64>>> {
    get_distance_matrix_with_params(
        mol,
        &GetDistanceMatrixParams {
            use_bo,
            use_atom_wts: false,
        },
    )
}

/// The characteristic polynomial's coefficients for the adjacency matrix
/// `a`, by the Le Verrier-Faddeev-Frame method, as Python's
/// `Graphs.CharacteristicPolynomial` computes them.
fn characteristic_polynomial(a: &[Vec<f64>]) -> Vec<f64> {
    let n = a.len();
    let mut res = Vec::with_capacity(n + 1);
    res.push(1.0);
    let mut an: Vec<Vec<f64>> = a.to_vec();
    for k in 1..=n {
        let diagonal: Vec<f64> = an.iter().enumerate().map(|(i, row)| row[i]).collect();
        let trace = np_sum(&diagonal);
        let coefficient = 1.0 / k as f64 * trace;
        res.push(coefficient);
        let mut bn = an;
        for (i, row) in bn.iter_mut().enumerate() {
            row[i] -= coefficient;
        }
        an = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| (0..n).map(|l| a[i][l] * bn[l][j]).sum())
                    .collect()
            })
            .collect();
    }
    for value in &mut res[1..] {
        *value *= -1.0;
    }
    res
}

/// Parameters for `ipc_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IpcParams {
    /// Returns the information content alone, not times the total.
    pub avg: bool,
    /// A distance matrix to use, whose entries of 1 are the adjacencies;
    /// one is computed when `None`.
    pub d_mat: Option<Vec<Vec<f64>>>,
    /// Computes the distance matrix even when `d_mat` is given.
    pub force_d_mat: bool,
}

/// The information content of the coefficients of the characteristic
/// polynomial of the hydrogen-suppressed graph's adjacency matrix, times
/// their total (D. Bonchev and N. Trinajstic, J. Chem. Phys. 67, 4517-4533
/// (1977), eq. 6). Mirrors `GraphDescriptors.Ipc`.
pub fn ipc(mol: &Mol) -> Result<f64> {
    ipc_with_params(mol, &IpcParams::default())
}

/// Mirrors `GraphDescriptors.Ipc(mol, avg, dMat, forceDMat)`. `Error::Value`
/// for a `d_mat` that is not one row and column per atom.
pub fn ipc_with_params(mol: &Mol, params: &IpcParams) -> Result<f64> {
    let n = mol.get_num_atoms() as usize;
    let computed;
    let d_mat = match &params.d_mat {
        Some(d_mat) if !params.force_d_mat => d_mat,
        _ => {
            computed = distance_matrix(mol, false)?;
            &computed
        }
    };
    check_square(d_mat, n)?;
    let adj_mat: Vec<Vec<f64>> = d_mat
        .iter()
        .map(|row| row.iter().map(|&d| f64::from(u8::from(d == 1.0))).collect())
        .collect();
    let c_poly: Vec<f64> = characteristic_polynomial(&adj_mat)
        .into_iter()
        .map(f64::abs)
        .collect();
    let entropy = info_entropy(&c_poly)?;
    if params.avg {
        Ok(entropy)
    } else {
        Ok(c_poly.iter().sum::<f64>() * entropy)
    }
}

/// Parameters for `avg_ipc_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AvgIpcParams {
    /// A distance matrix to use, as in [`IpcParams::d_mat`].
    pub d_mat: Option<Vec<Vec<f64>>>,
    /// Computes the distance matrix even when `d_mat` is given.
    pub force_d_mat: bool,
}

/// The average information content of the coefficients of the
/// characteristic polynomial (eq. 7 of Bonchev and Trinajstic): [`ipc`]
/// without the total. Mirrors `GraphDescriptors.AvgIpc`.
pub fn avg_ipc(mol: &Mol) -> Result<f64> {
    avg_ipc_with_params(mol, &AvgIpcParams::default())
}

/// Mirrors `GraphDescriptors.AvgIpc(mol, dMat, forceDMat)`.
pub fn avg_ipc_with_params(mol: &Mol, params: &AvgIpcParams) -> Result<f64> {
    ipc_with_params(
        mol,
        &IpcParams {
            avg: true,
            d_mat: params.d_mat.clone(),
            force_d_mat: params.force_d_mat,
        },
    )
}

/// Hall and Kier's alpha. Mirrors `GraphDescriptors.HallKierAlpha`, which
/// calls `rdMolDescriptors.CalcHallKierAlpha`.
pub fn hall_kier_alpha(mol: &Mol) -> Result<f64> {
    descriptors::calc_hall_kier_alpha(mol)
}

/// Hall and Kier's first kappa shape index. Mirrors
/// `GraphDescriptors.Kappa1`, which calls `rdMolDescriptors.CalcKappa1`.
pub fn kappa1(mol: &Mol) -> Result<f64> {
    descriptors::calc_kappa1(mol)
}

/// The second. Mirrors `GraphDescriptors.Kappa2`.
pub fn kappa2(mol: &Mol) -> Result<f64> {
    descriptors::calc_kappa2(mol)
}

/// The third. Mirrors `GraphDescriptors.Kappa3`.
pub fn kappa3(mol: &Mol) -> Result<f64> {
    descriptors::calc_kappa3(mol)
}

/// The connectivity index of order 0, from the atoms' degrees (Rev. Comp.
/// Chem. 2, 367-422 (1991), eqs. 1, 9, and 10); atoms without neighbors
/// add nothing. Mirrors `GraphDescriptors.Chi0`.
pub fn chi0(mol: &Mol) -> f64 {
    mol.get_atoms()
        .map(|atom| atom.get_degree())
        .filter(|&degree| degree != 0)
        .map(|degree| (1.0 / f64::from(degree)).sqrt())
        .fold(0.0, |sum, x| sum + x)
}

/// The connectivity index of order 1, from the products of bonded atoms'
/// degrees (eqs. 1, 11, and 12). Mirrors `GraphDescriptors.Chi1`.
pub fn chi1(mol: &Mol) -> f64 {
    mol.get_bonds()
        .map(|bond| bond.get_begin_atom().get_degree() * bond.get_end_atom().get_degree())
        .filter(|&product| product != 0)
        .map(|product| (1.0 / f64::from(product)).sqrt())
        .fold(0.0, |sum, x| sum + x)
}

/// Defines functions that call an `rdMolDescriptors` function of a
/// molecule.
macro_rules! calc {
    ($($(#[$meta:meta])* $name:ident = $calc:ident;)*) => {
        $(
            $(#[$meta])*
            pub fn $name(mol: &Mol) -> Result<f64> {
                descriptors::$calc(mol)
            }
        )*
    };
}

calc! {
    /// Kier and Hall's valence connectivity index of order 0. Mirrors
    /// `GraphDescriptors.Chi0v`, which calls `rdMolDescriptors.CalcChi0v`.
    chi0v = calc_chi0v;
    /// Of order 1. Mirrors `GraphDescriptors.Chi1v`.
    chi1v = calc_chi1v;
    /// Of order 2. Mirrors `GraphDescriptors.Chi2v`.
    chi2v = calc_chi2v;
    /// Of order 3. Mirrors `GraphDescriptors.Chi3v`.
    chi3v = calc_chi3v;
    /// Of order 4. Mirrors `GraphDescriptors.Chi4v`.
    chi4v = calc_chi4v;
    /// As `chi0v`, from valence electrons rather than valence. Mirrors
    /// `GraphDescriptors.Chi0n`, which calls `rdMolDescriptors.CalcChi0n`.
    chi0n = calc_chi0n;
    /// Of order 1. Mirrors `GraphDescriptors.Chi1n`.
    chi1n = calc_chi1n;
    /// Of order 2. Mirrors `GraphDescriptors.Chi2n`.
    chi2n = calc_chi2n;
    /// Of order 3. Mirrors `GraphDescriptors.Chi3n`.
    chi3n = calc_chi3n;
    /// Of order 4. Mirrors `GraphDescriptors.Chi4n`.
    chi4n = calc_chi4n;
}

/// The valence connectivity index of order `n`. Mirrors
/// `GraphDescriptors.ChiNv_`, which calls `rdMolDescriptors.CalcChiNv`.
pub fn chi_nv(mol: &Mol, n: u32) -> Result<f64> {
    descriptors::calc_chi_nv(mol, n)
}

/// As `chi_nv`, from valence electrons. Mirrors `GraphDescriptors.ChiNn_`,
/// which calls `rdMolDescriptors.CalcChiNn`.
pub fn chi_nn(mol: &Mol, n: u32) -> Result<f64> {
    descriptors::calc_chi_nn(mol, n)
}

/// Parameters for `balaban_j_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BalabanJParams {
    /// A distance matrix to use; one weighted by bond order is computed
    /// when `None`.
    pub d_mat: Option<Vec<Vec<f64>>>,
    /// Computes the distance matrix even when `d_mat` is given.
    pub force_d_mat: bool,
}

/// Balaban's J (Chem. Phys. Lett. 89, 399-404 (1982)), from the distance
/// matrix weighted by bond order. `Error::Value` for a molecule without
/// atoms, for which Python raises. Mirrors `GraphDescriptors.BalabanJ`.
pub fn balaban_j(mol: &Mol) -> Result<f64> {
    balaban_j_with_params(mol, &BalabanJParams::default())
}

/// Mirrors `GraphDescriptors.BalabanJ(mol, dMat, forceDMat)`. `Error::Value`
/// for a `d_mat` that is not one row and column per atom.
pub fn balaban_j_with_params(mol: &Mol, params: &BalabanJParams) -> Result<f64> {
    let n = mol.get_num_atoms() as usize;
    let computed;
    let d_mat = match &params.d_mat {
        Some(d_mat) if !params.force_d_mat => d_mat,
        _ => {
            computed = distance_matrix(mol, true)?;
            &computed
        }
    };
    if n == 0 {
        return Err(Error::Value {
            message: "BalabanJ needs at least one atom".to_owned(),
        });
    }
    check_square(d_mat, n)?;
    let adj_mat = get_adjacency_matrix(mol)?;
    // Python's sum over the matrix adds its rows in order.
    let mut s = vec![0.0; n];
    for row in d_mat {
        for (total, value) in s.iter_mut().zip(row) {
            *total += value;
        }
    }
    let q = mol.get_num_bonds();
    let mu = i64::from(q) - n as i64 + 1;
    let mut sum = 0.0;
    for i in 0..n {
        for j in i..n {
            if adj_mat[i][j] == 1 {
                sum += 1.0 / (s[i] * s[j]).sqrt();
            }
        }
    }
    if mu + 1 != 0 {
        Ok(f64::from(q) / (mu + 1) as f64 * sum)
    } else {
        Ok(0.0)
    }
}

/// Parameters for `bertz_ct_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct BertzCTParams {
    /// How many of each atom's nearest distances decide which atoms are
    /// equivalent; a negative cutoff counts from the end, as a Python slice
    /// does.
    pub cutoff: i32,
    /// A distance matrix to use when `force_d_mat` is false; one without
    /// bond orders is computed when `None`.
    pub d_mat: Option<Vec<Vec<f64>>>,
    /// Computes the distance matrix, weighted by bond order, whatever
    /// `d_mat` holds.
    pub force_d_mat: bool,
}

impl Default for BertzCTParams {
    fn default() -> BertzCTParams {
        BertzCTParams {
            cutoff: 100,
            d_mat: None,
            force_d_mat: true,
        }
    }
}

/// Values keyed in the order their keys first appear, as a Python `dict`
/// keeps them, so sums over them add in Python's order.
struct OrderedCounts<K> {
    index: HashMap<K, usize>,
    values: Vec<f64>,
}

impl<K: Eq + Hash> OrderedCounts<K> {
    fn new() -> OrderedCounts<K> {
        OrderedCounts {
            index: HashMap::new(),
            values: Vec::new(),
        }
    }

    fn add(&mut self, key: K, value: f64) {
        let next = self.values.len();
        let i = *self.index.entry(key).or_insert(next);
        if i == next {
            self.values.push(0.0);
        }
        self.values[i] += value;
    }
}

/// The classes of atoms whose sorted rows of the distance matrix agree, to
/// four decimals, in their first `cutoff` entries, numbered from 1 in the
/// order they first appear. Python's `_AssignSymmetryClasses`.
fn assign_symmetry_classes(bd_mat: &[Vec<f64>], num_atoms: usize, cutoff: i32) -> Vec<usize> {
    let mut keys_seen: HashMap<String, usize> = HashMap::new();
    bd_mat[..num_atoms]
        .iter()
        .map(|row| {
            let mut sorted = row.clone();
            // NaN last, as NumPy sorts. RDKit 2025.03's Python sorts each row
            // as a list, whose order around a NaN (which compares false with
            // everything) depends on the row; only a given `d_mat` holds NaN.
            sorted.sort_by(|a, b| match (a.is_nan(), b.is_nan()) {
                (false, false) => a.partial_cmp(b).expect("neither is NaN"),
                (a_nan, b_nan) => a_nan.cmp(&b_nan),
            });
            let len = sorted.len() as i64;
            let end = if cutoff < 0 {
                (len + i64::from(cutoff)).max(0)
            } else {
                i64::from(cutoff).min(len)
            };
            let key: String = sorted[..end as usize]
                .iter()
                .map(|value| format!("{value:.4},"))
                .collect();
            let next = keys_seen.len();
            *keys_seen.entry(key).or_insert(next) + 1
        })
        .collect()
}

/// Bertz's complexity index (S. H. Bertz, J. Am. Chem. Soc. 103,
/// 3599-3601 (1981)), from the connections between bonds and the
/// distribution of elements, with aromatic bonds of order 1.5 and the others
/// of their `BondType`'s number, as RDKit 2025.03's Python takes them (so a
/// dative bond is of order 17). 0 for fewer than two atoms. Mirrors
/// `GraphDescriptors.BertzCT`.
pub fn bertz_ct(mol: &Mol) -> Result<f64> {
    bertz_ct_with_params(mol, &BertzCTParams::default())
}

/// Mirrors `GraphDescriptors.BertzCT(mol, cutoff, dMat, forceDMat)`. Python
/// without `forceDMat` or `dMat` reads a matrix an earlier `Ipc` or
/// `BalabanJ` call left on the molecule, if any; this computes the distance
/// matrix without bond orders, as Python does for a new molecule.
/// `Error::Value` for a `d_mat` with fewer rows than atoms, or rows of
/// another length, and where Python's `math.log` raises: when the bonds'
/// connections total 0 or less (zero-order bonds).
pub fn bertz_ct_with_params(mol: &Mol, params: &BertzCTParams) -> Result<f64> {
    let num_atoms = mol.get_num_atoms() as usize;
    if num_atoms < 2 {
        return Ok(0.0);
    }
    let computed;
    let d_mat = if params.force_d_mat {
        computed = distance_matrix(mol, true)?;
        &computed
    } else if let Some(d_mat) = &params.d_mat {
        if d_mat.len() < num_atoms || d_mat.iter().any(|row| row.len() != d_mat[0].len()) {
            return Err(Error::Value {
                message: format!("the matrix needs at least {num_atoms} rows of one length"),
            });
        }
        d_mat
    } else {
        computed = distance_matrix(mol, false)?;
        &computed
    };

    // Python's _CreateBondDictEtc.
    let mut bond_dict: HashMap<(usize, usize), f64> = HashMap::new();
    let mut neighbor_list: Vec<Vec<usize>> = vec![Vec::new(); num_atoms];
    for bond in mol.get_bonds() {
        let (mut atom1, mut atom2) = (
            bond.get_begin_atom_idx() as usize,
            bond.get_end_atom_idx() as usize,
        );
        if atom1 > atom2 {
            (atom1, atom2) = (atom2, atom1);
        }
        // Python's float(bond.GetBondType()), the enum's number, which is the
        // bond order only from single to hextuple bonds.
        let order = if bond.get_is_aromatic() || bond.get_bond_type() == BondType::Aromatic {
            1.5
        } else {
            f64::from(bond.get_bond_type().to_ffi().repr)
        };
        bond_dict.insert((atom1, atom2), order);
        if !neighbor_list[atom1].contains(&atom2) {
            neighbor_list[atom1].push(atom2);
        }
        if !neighbor_list[atom2].contains(&atom1) {
            neighbor_list[atom2].push(atom1);
        }
    }
    for neighbors in &mut neighbor_list {
        neighbors.sort_unstable();
    }
    let bond_order = |a: usize, b: usize| bond_dict[&(a.min(b), a.max(b))];

    let symmetry_classes = assign_symmetry_classes(d_mat, num_atoms, params.cutoff);
    let mut atom_types: OrderedCounts<i32> = OrderedCounts::new();
    let mut connections: OrderedCounts<(usize, usize, usize)> = OrderedCounts::new();
    for (atom_idx, atom) in mol.get_atoms().enumerate() {
        atom_types.add(atom.get_atomic_num(), 1.0);
        let hinge_class = symmetry_classes[atom_idx];
        let neighbors = &neighbor_list[atom_idx];
        for (i, &neighbor_i) in neighbors.iter().enumerate() {
            let ni_class = symmetry_classes[neighbor_i];
            let bond_i_order = bond_order(atom_idx, neighbor_i);
            if bond_i_order > 1.0 && neighbor_i > atom_idx {
                let num_connections = bond_i_order * (bond_i_order - 1.0) / 2.0;
                // Python keys these by a pair; 0 never names a class.
                let key = (hinge_class.min(ni_class), hinge_class.max(ni_class), 0);
                connections.add(key, num_connections);
            }
            for &neighbor_j in &neighbors[i + 1..] {
                let nj_class = symmetry_classes[neighbor_j];
                let bond_j_order = bond_order(atom_idx, neighbor_j);
                let key = (ni_class.min(nj_class), hinge_class, ni_class.max(nj_class));
                connections.add(key, bond_i_order * bond_j_order);
            }
        }
    }
    let connection_list = if connections.values.is_empty() {
        vec![1.0]
    } else {
        connections.values
    };

    // Python's _CalculateEntropies.
    let tot_connections = py_sum(connection_list.iter().copied());
    if tot_connections <= 0.0 {
        return Err(Error::Value {
            message: "math domain error: the connections total 0 or less".to_owned(),
        });
    }
    let connection_ie =
        tot_connections * (info_entropy(&connection_list)? + tot_connections.ln() / 2f64.ln());
    let atom_type_ie = num_atoms as f64 * info_entropy(&atom_types.values)?;
    Ok(atom_type_ie + connection_ie)
}
