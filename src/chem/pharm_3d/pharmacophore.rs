// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm3D.Pharmacophore`, a Python module: 3D
//! pharmacophores, as features with distance bounds between them or with
//! explicit positions and radii.
//!
//! Python also takes molecules' features, copying their family, type, and
//! position into free features; here the features are free ones
//! (`FreeChemicalFeature::new(family, type, pos)` copies one).

use std::fmt;
use std::io::Read;

use crate::chem::rd_chemical_features::FreeChemicalFeature;
use crate::geometry::rd_geometry::Point3D;
use crate::{Error, Result};

/// Python's `'%W.Pf'`, and `'% W.Pf'` with `space_sign`.
pub(crate) fn py_fixed(x: f64, width: usize, precision: usize, space_sign: bool) -> String {
    let body = if x.is_nan() {
        "nan".to_owned()
    } else if x.is_infinite() {
        "inf".to_owned()
    } else {
        format!("{:.precision$}", x.abs())
    };
    let sign = if x.is_sign_negative() && !x.is_nan() {
        "-"
    } else if space_sign {
        " "
    } else {
        ""
    };
    format!("{:>width$}", format!("{sign}{body}"))
}

fn index_error(i: usize, len: usize) -> Error {
    Error::index(i, len)
}

/// A pharmacophore: features with lower and upper bounds on the distances
/// between them, in 3D and topologically. Mirrors
/// `Pharmacophore.Pharmacophore`.
#[derive(Clone, Debug)]
pub struct Pharmacophore {
    feats: Vec<FreeChemicalFeature>,
    bounds_mat: Vec<Vec<f64>>,
    bounds_mat_2d: Vec<Vec<i64>>,
}

impl Pharmacophore {
    /// A pharmacophore whose 3D bounds are the features' distances and
    /// whose topological upper bounds are 1000. Mirrors
    /// `Pharmacophore(feats)`.
    pub fn new(feats: &[FreeChemicalFeature]) -> Pharmacophore {
        Pharmacophore::new_with_init_mats(feats, true)
    }

    /// With all bounds 0 unless `init_mats`. Mirrors `Pharmacophore(feats,
    /// initMats)`.
    pub fn new_with_init_mats(feats: &[FreeChemicalFeature], init_mats: bool) -> Pharmacophore {
        let nf = feats.len();
        let mut res = Pharmacophore {
            feats: feats.to_vec(),
            bounds_mat: vec![vec![0.0; nf]; nf],
            bounds_mat_2d: vec![vec![0; nf]; nf],
        };
        if init_mats {
            for i in 1..nf {
                let loci = res.feats[i].get_pos();
                for j in 0..i {
                    let dist = loci.distance(&res.feats[j].get_pos());
                    res.bounds_mat[i][j] = dist;
                    res.bounds_mat[j][i] = dist;
                }
            }
            for i in 0..nf {
                for j in i + 1..nf {
                    res.bounds_mat_2d[i][j] = 1000;
                }
            }
        }
        res
    }

    /// Mirrors `Pharmacophore.getFeatures`.
    pub fn get_features(&self) -> &[FreeChemicalFeature] {
        &self.feats
    }

    /// The feature at `i`; `Error::Index` past the end. Mirrors
    /// `Pharmacophore.getFeature`.
    pub fn get_feature(&self, i: usize) -> Result<&FreeChemicalFeature> {
        self.feats
            .get(i)
            .ok_or_else(|| index_error(i, self.feats.len()))
    }

    fn check(&self, i: usize, j: usize) -> Result<()> {
        let nf = self.feats.len();
        for k in [i, j] {
            if k >= nf {
                return Err(index_error(k, nf));
            }
        }
        Ok(())
    }

    /// Python's `checkBounds`, which raises ValueError.
    fn check_bounds(&self, i: usize, j: usize) -> Result<()> {
        self.check(i, j).map_err(|_| Error::Value {
            message: "Index out of bound".to_owned(),
        })
    }

    /// The upper bound on the distance between features `i` and `j`.
    /// Indices past the features are `Error::Index`. Mirrors
    /// `Pharmacophore.getUpperBound`.
    pub fn get_upper_bound(&self, i: usize, j: usize) -> Result<f64> {
        self.check(i, j)?;
        Ok(self.bounds_mat[i.min(j)][i.max(j)])
    }

    /// The lower bound. Mirrors `Pharmacophore.getLowerBound`.
    pub fn get_lower_bound(&self, i: usize, j: usize) -> Result<f64> {
        self.check(i, j)?;
        Ok(self.bounds_mat[i.max(j)][i.min(j)])
    }

    /// Sets the upper bound; indices past the features are `Error::Index`.
    /// Mirrors `Pharmacophore.setUpperBound`.
    pub fn set_upper_bound(&mut self, i: usize, j: usize, val: f64) -> Result<()> {
        self.set_upper_bound_with_check_bounds(i, j, val, false)
    }

    /// With `check_bounds`, such indices are `Error::Value` instead, as
    /// Python's check raises. Mirrors `Pharmacophore.setUpperBound(i, j,
    /// val, checkBounds)`.
    pub fn set_upper_bound_with_check_bounds(
        &mut self,
        i: usize,
        j: usize,
        val: f64,
        check_bounds: bool,
    ) -> Result<()> {
        if check_bounds {
            self.check_bounds(i, j)?;
        }
        self.check(i, j)?;
        self.bounds_mat[i.min(j)][i.max(j)] = val;
        Ok(())
    }

    /// Sets the lower bound. Mirrors `Pharmacophore.setLowerBound`.
    pub fn set_lower_bound(&mut self, i: usize, j: usize, val: f64) -> Result<()> {
        self.set_lower_bound_with_check_bounds(i, j, val, false)
    }

    /// Mirrors `Pharmacophore.setLowerBound(i, j, val, checkBounds)`.
    pub fn set_lower_bound_with_check_bounds(
        &mut self,
        i: usize,
        j: usize,
        val: f64,
        check_bounds: bool,
    ) -> Result<()> {
        if check_bounds {
            self.check_bounds(i, j)?;
        }
        self.check(i, j)?;
        self.bounds_mat[i.max(j)][i.min(j)] = val;
        Ok(())
    }

    /// The upper bound on the topological distance. Mirrors
    /// `Pharmacophore.getUpperBound2D`.
    pub fn get_upper_bound_2d(&self, i: usize, j: usize) -> Result<i64> {
        self.check(i, j)?;
        Ok(self.bounds_mat_2d[i.min(j)][i.max(j)])
    }

    /// The lower bound on the topological distance. Mirrors
    /// `Pharmacophore.getLowerBound2D`.
    pub fn get_lower_bound_2d(&self, i: usize, j: usize) -> Result<i64> {
        self.check(i, j)?;
        Ok(self.bounds_mat_2d[i.max(j)][i.min(j)])
    }

    /// Mirrors `Pharmacophore.setUpperBound2D`.
    pub fn set_upper_bound_2d(&mut self, i: usize, j: usize, val: i64) -> Result<()> {
        self.set_upper_bound_2d_with_check_bounds(i, j, val, false)
    }

    /// Mirrors `Pharmacophore.setUpperBound2D(i, j, val, checkBounds)`.
    pub fn set_upper_bound_2d_with_check_bounds(
        &mut self,
        i: usize,
        j: usize,
        val: i64,
        check_bounds: bool,
    ) -> Result<()> {
        if check_bounds {
            self.check_bounds(i, j)?;
        }
        self.check(i, j)?;
        self.bounds_mat_2d[i.min(j)][i.max(j)] = val;
        Ok(())
    }

    /// Mirrors `Pharmacophore.setLowerBound2D`.
    pub fn set_lower_bound_2d(&mut self, i: usize, j: usize, val: i64) -> Result<()> {
        self.set_lower_bound_2d_with_check_bounds(i, j, val, false)
    }

    /// Mirrors `Pharmacophore.setLowerBound2D(i, j, val, checkBounds)`.
    pub fn set_lower_bound_2d_with_check_bounds(
        &mut self,
        i: usize,
        j: usize,
        val: i64,
        check_bounds: bool,
    ) -> Result<()> {
        if check_bounds {
            self.check_bounds(i, j)?;
        }
        self.check(i, j)?;
        self.bounds_mat_2d[i.max(j)][i.min(j)] = val;
        Ok(())
    }
}

/// A table of the features' families with the lower bounds below the
/// diagonal and the upper bounds above it. Mirrors
/// `Pharmacophore.__str__`.
impl fmt::Display for Pharmacophore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut res = " ".repeat(14);
        for feat in &self.feats {
            res += &format!("{:>13} ", feat.get_family());
        }
        res += "\n";
        let nf = self.feats.len();
        for (i, feat) in self.feats.iter().enumerate() {
            res += &format!("{:>13} ", feat.get_family());
            for j in 0..nf {
                let (val, space) = match j.cmp(&i) {
                    std::cmp::Ordering::Less => (self.bounds_mat[i][j], false),
                    std::cmp::Ordering::Greater => (self.bounds_mat[i][j], false),
                    std::cmp::Ordering::Equal => (0.0, true),
                };
                res += &py_fixed(val, 13, 3, space);
                res += " ";
            }
            res += "\n";
        }
        f.write_str(&res)
    }
}

/// A pharmacophore of features at explicit positions with radii. Python's
/// is empty until it is given features. Mirrors
/// `Pharmacophore.ExplicitPharmacophore`.
#[derive(Clone, Debug, Default)]
pub struct ExplicitPharmacophore {
    feats: Vec<FreeChemicalFeature>,
    radii: Vec<f64>,
}

impl ExplicitPharmacophore {
    /// A pharmacophore of these features and radii; lists of different
    /// lengths are `Error::Value`. Python makes an empty one if either list
    /// is empty. Mirrors `ExplicitPharmacophore(feats, radii)`; its
    /// defaults, `None`, make the empty pharmacophore `default` makes.
    pub fn new(feats: &[FreeChemicalFeature], radii: &[f64]) -> Result<ExplicitPharmacophore> {
        let mut res = ExplicitPharmacophore::default();
        if !feats.is_empty() && !radii.is_empty() {
            res.initialize_feats(feats.to_vec(), radii.to_vec())?;
        }
        Ok(res)
    }

    fn initialize_feats(&mut self, feats: Vec<FreeChemicalFeature>, radii: Vec<f64>) -> Result<()> {
        if feats.len() != radii.len() {
            return Err(Error::Value {
                message: "len(feats)!=len(radii)".to_owned(),
            });
        }
        self.feats = feats;
        self.radii = radii;
        Ok(())
    }

    /// Mirrors `ExplicitPharmacophore.getFeatures`.
    pub fn get_features(&self) -> &[FreeChemicalFeature] {
        &self.feats
    }

    /// Mirrors `ExplicitPharmacophore.getRadii`.
    pub fn get_radii(&self) -> &[f64] {
        &self.radii
    }

    /// The feature at `i`; `Error::Index` past the end. Mirrors
    /// `ExplicitPharmacophore.getFeature`.
    pub fn get_feature(&self, i: usize) -> Result<&FreeChemicalFeature> {
        self.feats
            .get(i)
            .ok_or_else(|| index_error(i, self.feats.len()))
    }

    /// The radius at `i`; `Error::Index` past the end. Mirrors
    /// `ExplicitPharmacophore.getRadius`.
    pub fn get_radius(&self, i: usize) -> Result<f64> {
        self.radii
            .get(i)
            .copied()
            .ok_or_else(|| index_error(i, self.radii.len()))
    }

    /// Mirrors `ExplicitPharmacophore.setRadius`.
    pub fn set_radius(&mut self, i: usize, rad: f64) -> Result<()> {
        let len = self.radii.len();
        *self.radii.get_mut(i).ok_or_else(|| index_error(i, len))? = rad;
        Ok(())
    }

    /// Reads the pharmacophore from text. Python splits the text at the
    /// two characters `\n`, not at line ends, and this does too. Mirrors
    /// `ExplicitPharmacophore.initFromString`.
    pub fn init_from_string(&mut self, text: &str) -> Result<()> {
        let lines: Vec<&str> = text.split("\\n").collect();
        self.init_from_lines(&lines)
    }

    /// Reads the pharmacophore from a file's lines; `Error::Os` if it
    /// cannot be read and `Error::Value` if it is not UTF-8. Mirrors
    /// `ExplicitPharmacophore.initFromFile`.
    pub fn init_from_file<R: Read + ?Sized>(&mut self, in_f: &mut R) -> Result<()> {
        let mut text = String::new();
        in_f.read_to_string(&mut text)
            .map_err(|err| match err.kind() {
                std::io::ErrorKind::InvalidData => Error::Value {
                    message: format!("the pharmacophore is not UTF-8: {err}"),
                },
                _ => Error::Os {
                    message: err.to_string(),
                },
            })?;
        let lines: Vec<&str> = text.split_inclusive('\n').collect();
        self.init_from_lines(&lines)
    }

    /// Reads the pharmacophore from lines of a family, x, y, z, and radius,
    /// each with any comment after `#`. A line with fewer fields or a
    /// number that does not parse is `Error::Parse`, leaving the
    /// pharmacophore as it was (Python logs the error and returns). Mirrors
    /// `ExplicitPharmacophore.initFromLines`.
    pub fn init_from_lines(&mut self, lines: &[&str]) -> Result<()> {
        let mut feats = Vec::new();
        let mut rads = Vec::new();
        for (line_num, line) in lines.iter().enumerate() {
            let txt = line.split('#').next().unwrap_or_default().trim();
            if txt.is_empty() {
                continue;
            }
            let split_l: Vec<&str> = txt.split([' ', '\t']).filter(|s| !s.is_empty()).collect();
            if split_l.len() < 5 {
                return Err(Error::Parse {
                    message: format!(
                        "Input line {line_num} only contains {} fields, 5 are required. Read failed.",
                        split_l.len()
                    ),
                });
            }
            let number = |s: &str| -> Result<f64> {
                s.parse::<f64>().map_err(|_| Error::Parse {
                    message: format!("Error parsing a number of line {line_num}. Read failed."),
                })
            };
            let (x, y, z, rad) = (
                number(split_l[1])?,
                number(split_l[2])?,
                number(split_l[3])?,
                number(split_l[4])?,
            );
            let name = split_l[0];
            feats.push(FreeChemicalFeature::new(name, name, Point3D::new(x, y, z)));
            rads.push(rad);
        }
        self.initialize_feats(feats, rads)
    }
}

/// One line per feature: its family, position, and radius. Mirrors
/// `ExplicitPharmacophore.__str__`.
impl fmt::Display for ExplicitPharmacophore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (feat, rad) in self.feats.iter().zip(&self.radii) {
            let p = feat.get_pos();
            writeln!(
                f,
                "{:>12}    {} {} {}    {}",
                feat.get_family(),
                py_fixed(p.x, 8, 4, true),
                py_fixed(p.y, 8, 4, true),
                py_fixed(p.z, 8, 4, true),
                py_fixed(*rad, 5, 2, true),
            )?;
        }
        Ok(())
    }
}
