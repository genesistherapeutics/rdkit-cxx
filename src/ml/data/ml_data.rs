// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Data.MLData`, a Python module: data sets for
//! machine learning, rows of variables followed by results.
//!
//! Python's general data sets hold any objects; RDKit reads integers,
//! floats, and strings into them, which [`Value`] holds.

use std::fmt;

use crate::{Error, Result};

/// A value in an [`MLDataSet`]: what RDKit's readers make of a field.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Value {
    /// An integer.
    Int(i64),
    /// A float.
    Float(f64),
    /// Text.
    Str(String),
}

impl Value {
    /// The value as a float, for a number.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Int(i) => Some(*i as f64),
            Value::Float(x) => Some(*x),
            Value::Str(_) => None,
        }
    }
}

/// Python's `repr` of a float: the shortest text that reads back as it,
/// in positional notation for exponents from -4 to 15.
pub(crate) fn float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_owned();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_owned();
    }
    let sci = format!("{x:e}");
    let (mantissa, exp) = sci.split_once('e').expect("{:e} has an exponent");
    let exp: i32 = exp.parse().expect("{:e} has an integer exponent");
    let (sign, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => ("-", m),
        None => ("", mantissa),
    };
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if (-4..16).contains(&exp) {
        let n = digits.len() as i32;
        let text = if exp < 0 {
            format!("0.{}{digits}", "0".repeat((-exp - 1) as usize))
        } else if exp + 1 >= n {
            format!("{digits}{}.0", "0".repeat((exp + 1 - n) as usize))
        } else {
            let (a, b) = digits.split_at((exp + 1) as usize);
            format!("{a}.{b}")
        };
        format!("{sign}{text}")
    } else {
        let (first, rest) = digits.split_at(1);
        let mant = if rest.is_empty() {
            first.to_owned()
        } else {
            format!("{first}.{rest}")
        };
        let esign = if exp < 0 { '-' } else { '+' };
        format!("{sign}{mant}e{esign}{:02}", exp.abs())
    }
}

/// Python's `str` of the value.
impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Value::Int(i) => write!(f, "{i}"),
            Value::Float(x) => f.write_str(&float_repr(*x)),
            Value::Str(s) => f.write_str(s),
        }
    }
}

/// Python's `math.floor(d) == d` test and the floor, for a number; `Err`
/// for a float with no integer floor (infinities and NaN), which Python
/// raises for.
fn integral(d: &Value) -> Result<Option<Option<i64>>> {
    Ok(match d {
        Value::Int(i) => Some(Some(*i)),
        Value::Float(x) => {
            if !x.is_finite() {
                return Err(Error::Value {
                    message: format!("cannot convert float {} to integer", float_repr(*x)),
                });
            }
            let fl = x.floor();
            Some(if fl == *x { Some(fl as i64) } else { None })
        }
        Value::Str(_) => None,
    })
}

/// Parameters for [`MLDataSet::new_with_params`] and
/// [`MLQuantDataSet::new_with_params`], with Python's defaults (A2): `None`
/// computes a value from the data.
#[derive(Clone, Debug, PartialEq)]
pub struct MLDataSetParams {
    /// The number of variables; the row length less the results by default.
    pub n_vars: Option<usize>,
    /// The number of points; the number of rows by default.
    pub n_pts: Option<usize>,
    /// The number of values each column takes, -1 plus 1 (that is, 0) for
    /// columns that are not integers.
    pub n_possible_vals: Option<Vec<i64>>,
    /// The quantization bounds of each column, empty for none.
    pub q_bounds: Option<Vec<Vec<f64>>>,
    /// The variables' names; empty names by default.
    pub var_names: Option<Vec<String>>,
    /// The points' names; empty names by default.
    pub pt_names: Option<Vec<String>>,
    /// The number of result columns, at the end of each row.
    pub n_results: usize,
}

impl Default for MLDataSetParams {
    fn default() -> MLDataSetParams {
        MLDataSetParams {
            n_vars: None,
            n_pts: None,
            n_possible_vals: None,
            q_bounds: None,
            var_names: None,
            pt_names: None,
            n_results: 1,
        }
    }
}

fn first_row_len<T>(data: &[Vec<T>]) -> Result<usize> {
    data.first().map(Vec::len).ok_or(Error::index(0, 0))
}

/// The getters both data sets share.
macro_rules! data_set_getters {
    () => {
        /// The number of result columns. Mirrors `MLDataSet.GetNResults`.
        pub fn get_n_results(&self) -> usize {
            self.n_results
        }

        /// The number of variables. Mirrors `MLDataSet.GetNVars`.
        pub fn get_n_vars(&self) -> usize {
            self.n_vars
        }

        /// The number of points. Mirrors `MLDataSet.GetNPts`.
        pub fn get_n_pts(&self) -> usize {
            self.n_pts
        }

        /// The number of values each column takes, 0 for columns that are
        /// not integers. Mirrors `MLDataSet.GetNPossibleVals`.
        pub fn get_n_possible_vals(&self) -> &[i64] {
            &self.n_possible_vals
        }

        /// The quantization bounds of each column. Mirrors
        /// `MLDataSet.GetQuantBounds`.
        pub fn get_quant_bounds(&self) -> &[Vec<f64>] {
            &self.q_bounds
        }

        /// The variables' names. Mirrors `MLDataSet.GetVarNames`.
        pub fn get_var_names(&self) -> &[String] {
            &self.var_names
        }

        /// The points' names. Mirrors `MLDataSet.GetPtNames`.
        pub fn get_pt_names(&self) -> &[String] {
            &self.pt_names
        }
    };
}

/// A data set of general values. Mirrors `MLData.MLDataSet`.
#[derive(Clone, Debug, PartialEq)]
pub struct MLDataSet {
    data: Vec<Vec<Value>>,
    n_results: usize,
    n_vars: usize,
    n_pts: usize,
    q_bounds: Vec<Vec<f64>>,
    n_possible_vals: Vec<i64>,
    var_names: Vec<String>,
    pt_names: Vec<String>,
}

impl MLDataSet {
    /// A data set of the rows `data`, with one result column. Mirrors
    /// `MLDataSet(data)`.
    ///
    /// Errors as [`new_with_params`](Self::new_with_params).
    pub fn new(data: Vec<Vec<Value>>) -> Result<MLDataSet> {
        MLDataSet::new_with_params(data, &MLDataSetParams::default())
    }

    /// Mirrors `MLDataSet(data, nVars, nPts, nPossibleVals, qBounds,
    /// varNames, ptNames, nResults)`.
    ///
    /// `Error::Index` for no rows unless both the number of variables and
    /// the quantization bounds are given, or for a row shorter than the
    /// first, as Python raises; `Error::Value` for an infinite or NaN float
    /// in a column whose values are counted.
    pub fn new_with_params(data: Vec<Vec<Value>>, params: &MLDataSetParams) -> Result<MLDataSet> {
        MLDataSet::from_params(data, params.clone())
    }

    /// As [`new_with_params`](MLDataSet::new_with_params), taking the
    /// parameters.
    pub(crate) fn from_params(data: Vec<Vec<Value>>, params: MLDataSetParams) -> Result<MLDataSet> {
        let n_results = params.n_results;
        let n_vars = match params.n_vars {
            Some(n) => n,
            None => first_row_len(&data)?
                .checked_sub(n_results)
                .ok_or_else(|| Error::Value {
                    message: "the rows are shorter than the results".to_owned(),
                })?,
        };
        let n_pts = params.n_pts.unwrap_or(data.len());
        let q_bounds = match params.q_bounds {
            Some(q) => q,
            None => vec![Vec::new(); first_row_len(&data)?],
        };
        let mut set = MLDataSet {
            data,
            n_results,
            n_vars,
            n_pts,
            q_bounds,
            n_possible_vals: Vec::new(),
            var_names: params
                .var_names
                .unwrap_or_else(|| vec![String::new(); n_vars]),
            pt_names: params
                .pt_names
                .unwrap_or_else(|| vec![String::new(); n_pts]),
        };
        set.n_possible_vals = match params.n_possible_vals {
            Some(n) => n,
            None => set.calc_n_possible()?,
        };
        Ok(set)
    }

    /// Python's `_CalcNPossible`.
    fn calc_n_possible(&self) -> Result<Vec<i64>> {
        let n_vars = self.n_vars + self.n_results;
        let mut n_possible = vec![-1i64; n_vars];
        let mut cols: Vec<usize> = (0..n_vars).collect();
        for (i, bounds) in self.q_bounds.iter().enumerate() {
            if !bounds.is_empty() {
                let slot = n_possible.get_mut(i).ok_or(Error::index(i, n_vars))?;
                *slot = bounds.len() as i64;
                cols.retain(|&c| c != i);
            }
        }
        for pt in &self.data {
            for col in cols.clone() {
                let d = pt.get(col).ok_or(Error::index(col, pt.len()))?;
                match integral(d)? {
                    Some(Some(fl)) => n_possible[col] = n_possible[col].max(fl),
                    _ => {
                        n_possible[col] = -1;
                        cols.retain(|&c| c != col);
                    }
                }
            }
        }
        Ok(n_possible.into_iter().map(|x| x + 1).collect())
    }

    data_set_getters!();

    /// The point's name followed by its values. Mirrors `dataSet[idx]`.
    ///
    /// `Error::Index` past the points.
    pub fn get(&self, idx: usize) -> Result<Vec<Value>> {
        let (name, row) = self
            .pt_names
            .get(idx)
            .zip(self.data.get(idx))
            .ok_or(Error::index(idx, self.data.len()))?;
        let mut res = vec![Value::Str(name.clone())];
        res.extend(row.iter().cloned());
        Ok(res)
    }

    /// Replaces a point with `val`, its name (a string) followed by its
    /// values. Mirrors `dataSet[idx] = val`.
    ///
    /// `Error::Value` for a value of the wrong length, as Python raises, or
    /// a name that is not a string; `Error::Index` past the points.
    pub fn set(&mut self, idx: usize, val: Vec<Value>) -> Result<()> {
        if val.len() != self.n_vars + self.n_results + 1 {
            return Err(Error::Value {
                message: "bad value in assignment".to_owned(),
            });
        }
        if idx >= self.data.len() || idx >= self.pt_names.len() {
            return Err(Error::index(idx, self.data.len()));
        }
        let mut val = val.into_iter();
        let Some(Value::Str(name)) = val.next() else {
            return Err(Error::Value {
                message: "a point's name is a string".to_owned(),
            });
        };
        self.pt_names[idx] = name;
        self.data[idx] = val.collect();
        Ok(())
    }

    /// Each point's name followed by its values. Mirrors
    /// `MLDataSet.GetNamedData`.
    ///
    /// `Error::Index` when there are fewer rows or names than points.
    pub fn get_named_data(&self) -> Result<Vec<Vec<Value>>> {
        (0..self.n_pts).map(|i| self.get(i)).collect()
    }

    /// A copy of the rows. Mirrors `MLDataSet.GetAllData`.
    pub fn get_all_data(&self) -> Vec<Vec<Value>> {
        self.data.clone()
    }

    /// The rows without their result columns. Mirrors
    /// `MLDataSet.GetInputData`.
    pub fn get_input_data(&self) -> Vec<Vec<Value>> {
        // Python slices `x[:-v]`, which is empty for no results.
        let v = self.n_results;
        self.data
            .iter()
            .map(|x| {
                let end = if v == 0 { 0 } else { x.len().saturating_sub(v) };
                x[..end].to_vec()
            })
            .collect()
    }

    /// Each row's result columns. Mirrors `MLDataSet.GetResults`, which
    /// gives the results themselves, not one-item lists, when there is one
    /// result column.
    pub fn get_results(&self) -> Vec<Vec<Value>> {
        let v = self.n_results.max(1);
        self.data
            .iter()
            .map(|x| x[x.len().saturating_sub(v)..].to_vec())
            .collect()
    }

    /// Adds a point, its name (a string) followed by its values. Mirrors
    /// `MLDataSet.AddPoint`.
    ///
    /// `Error::Value` for a point without a name that is a string.
    pub fn add_point(&mut self, pt: Vec<Value>) -> Result<()> {
        let mut pt = pt.into_iter();
        let Some(Value::Str(name)) = pt.next() else {
            return Err(Error::Value {
                message: "a point's name is a string".to_owned(),
            });
        };
        self.data.push(pt.collect());
        self.pt_names.push(name);
        self.n_pts += 1;
        Ok(())
    }

    /// Adds points and their names. Mirrors `MLDataSet.AddPoints`.
    ///
    /// `Error::Value` for different numbers of points and names, as Python
    /// raises.
    pub fn add_points(&mut self, pts: Vec<Vec<Value>>, names: Vec<String>) -> Result<()> {
        if pts.len() != names.len() {
            return Err(Error::Value {
                message: "input length mismatch".to_owned(),
            });
        }
        self.data.extend(pts);
        self.pt_names.extend(names);
        self.n_pts = self.data.len();
        Ok(())
    }
}

/// A data set of quantized (integer) values; the results need no
/// quantization bounds. Mirrors `MLData.MLQuantDataSet`.
///
/// Python's version inherits `MLDataSet`'s indexing and point-adding
/// methods, which fail on its NumPy array; they are not ported.
#[derive(Clone, Debug, PartialEq)]
pub struct MLQuantDataSet {
    data: Vec<Vec<i64>>,
    n_results: usize,
    n_vars: usize,
    n_pts: usize,
    q_bounds: Vec<Vec<f64>>,
    n_possible_vals: Vec<i64>,
    var_names: Vec<String>,
    pt_names: Vec<String>,
}

impl MLQuantDataSet {
    /// A data set of the rows `data`, with one result column. Mirrors
    /// `MLQuantDataSet(data)`.
    ///
    /// Errors as [`new_with_params`](Self::new_with_params).
    pub fn new(data: Vec<Vec<i64>>) -> Result<MLQuantDataSet> {
        MLQuantDataSet::new_with_params(data, &MLDataSetParams::default())
    }

    /// Mirrors `MLQuantDataSet(data, nVars, nPts, nPossibleVals, qBounds,
    /// varNames, ptNames, nResults)`. Each column's number of possible
    /// values is its largest value plus 1.
    ///
    /// `Error::Index` for no rows unless the number of variables is given;
    /// `Error::Value` for rows of different lengths (which NumPy cannot
    /// make an integer array of) or no rows and no possible values.
    pub fn new_with_params(
        data: Vec<Vec<i64>>,
        params: &MLDataSetParams,
    ) -> Result<MLQuantDataSet> {
        MLQuantDataSet::from_params(data, params.clone())
    }

    /// As [`new_with_params`](MLQuantDataSet::new_with_params), taking the
    /// parameters.
    pub(crate) fn from_params(
        data: Vec<Vec<i64>>,
        params: MLDataSetParams,
    ) -> Result<MLQuantDataSet> {
        let ncols = data.first().map_or(0, Vec::len);
        if data.iter().any(|row| row.len() != ncols) {
            return Err(Error::Value {
                message: "the rows have different lengths".to_owned(),
            });
        }
        let n_results = params.n_results;
        let n_vars = match params.n_vars {
            Some(n) => n,
            None => first_row_len(&data)?
                .checked_sub(n_results)
                .ok_or_else(|| Error::Value {
                    message: "the rows are shorter than the results".to_owned(),
                })?,
        };
        let n_pts = params.n_pts.unwrap_or(data.len());
        let n_possible_vals = match params.n_possible_vals {
            Some(n) => n,
            None => (0..ncols)
                .map(|c| data.iter().map(|row| row[c]).max().map(|m| m + 1))
                .collect::<Option<Vec<i64>>>()
                .ok_or_else(|| Error::Value {
                    message: "no rows to count values in".to_owned(),
                })?,
        };
        Ok(MLQuantDataSet {
            data,
            n_results,
            n_vars,
            n_pts,
            q_bounds: params.q_bounds.unwrap_or_else(|| vec![Vec::new(); n_vars]),
            n_possible_vals,
            var_names: params
                .var_names
                .unwrap_or_else(|| vec![String::new(); n_vars]),
            pt_names: params
                .pt_names
                .unwrap_or_else(|| vec![String::new(); n_pts]),
        })
    }

    data_set_getters!();

    /// Each point's name followed by its values. Mirrors
    /// `MLQuantDataSet.GetNamedData`.
    ///
    /// `Error::Index` when there are fewer rows or names than points.
    pub fn get_named_data(&self) -> Result<Vec<(String, Vec<i64>)>> {
        (0..self.n_pts)
            .map(|i| {
                self.pt_names
                    .get(i)
                    .zip(self.data.get(i))
                    .map(|(n, r)| (n.clone(), r.clone()))
                    .ok_or(Error::index(i, self.data.len()))
            })
            .collect()
    }

    /// A copy of the rows. Mirrors `MLQuantDataSet.GetAllData`.
    pub fn get_all_data(&self) -> Vec<Vec<i64>> {
        self.data.clone()
    }

    /// The rows without their result columns. Mirrors
    /// `MLQuantDataSet.GetInputData`.
    pub fn get_input_data(&self) -> Vec<Vec<i64>> {
        // Python slices `x[:-v]`, which is empty for no results.
        let v = self.n_results;
        self.data
            .iter()
            .map(|x| {
                let end = if v == 0 { 0 } else { x.len().saturating_sub(v) };
                x[..end].to_vec()
            })
            .collect()
    }

    /// Each row's result columns. Mirrors `MLQuantDataSet.GetResults`,
    /// which gives the results themselves, not one-item lists, when there
    /// is one result column.
    pub fn get_results(&self) -> Vec<Vec<i64>> {
        let v = self.n_results.max(1);
        self.data
            .iter()
            .map(|x| x[x.len().saturating_sub(v)..].to_vec())
            .collect()
    }
}
