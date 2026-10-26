// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Data.DataUtils`, a Python module: reading, writing,
//! and sampling data sets.
//!
//! A `.qdat` file holds quantized data: any lines, a line with `Variable
//! Table`, a line `# name [bounds]` per variable (the last is the result),
//! a line starting `# ----`, then one line per point, its name and integer
//! values separated by spaces or tabs. A `.dat` file is the same with
//! values of any kind (integers, floats, or text), each line ending in a
//! separator, as RDKit's writer leaves it.
//!
//! The random choices draw from the generator `rdkit_cxx::rd_random`
//! shares, as Python's draw from `RDRandom`, so with the same seed
//! ([`init_random_numbers`]) they are Python's.

use std::cmp::Ordering;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::ml::data::ml_data::{MLDataSet, MLDataSetParams, MLQuantDataSet, Value, float_repr};
use crate::rd_random;
use crate::utils::fileutils::{move_to_matching_line, read_line};
use crate::utils::{py_float, py_int, py_strip, read_line_universal};
use crate::{Error, Result};

fn index_error(index: usize, len: usize) -> Error {
    Error::index(index, len)
}

fn os_error(e: std::io::Error) -> Error {
    Error::Os {
        message: e.to_string(),
    }
}

/// A random permutation of `0..n_to_do`. Mirrors `DataUtils.permutation`.
pub fn permutation(n_to_do: usize) -> Vec<usize> {
    let mut res: Vec<usize> = (0..n_to_do).collect();
    rd_random::shuffle_with_random(&mut res, rd_random::random)
        .expect("the generator's floats are in [0, 1)");
    res
}

/// Python's `str` of a list of floats.
fn list_repr(values: &[f64]) -> String {
    let items: Vec<String> = values.iter().map(|&x| float_repr(x)).collect();
    format!("[{}]", items.join(", "))
}

/// Writes a `.qdat` file: the variable table, then each example's values
/// separated by spaces. Mirrors `DataUtils.WriteData`.
///
/// `Error::Index` for fewer bounds than names; `Error::Os` when writing
/// fails.
pub fn write_data<W: Write>(
    out_file: &mut W,
    var_names: &[String],
    q_bounds: &[Vec<f64>],
    examples: &[Vec<Value>],
) -> Result<()> {
    let mut text =
        String::from("# Quantized data from DataUtils\n# ----------\n# Variable Table\n");
    for (i, name) in var_names.iter().enumerate() {
        let bounds = q_bounds.get(i).ok_or(index_error(i, q_bounds.len()))?;
        text.push_str(&format!("# {name} {}\n", list_repr(bounds)));
    }
    text.push_str("# ----------\n");
    for example in examples {
        let items: Vec<String> = example.iter().map(Value::to_string).collect();
        text.push_str(&items.join(" "));
        text.push('\n');
    }
    out_file.write_all(text.as_bytes()).map_err(os_error)
}

/// The variables' names and quantization bounds from a `.qdat` or `.dat`
/// file's variable table, leaving the file after it. Mirrors
/// `DataUtils.ReadVars`.
///
/// `Error::Runtime` with no variable table; `Error::Index` for a variable
/// line without `[`, or a table without its closing line; `Error::Value`
/// for a bound that is not a number; `Error::Os` when reading fails.
pub fn read_vars<R: BufRead>(in_file: &mut R) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    move_to_matching_line(in_file, "Variable Table")?;
    let mut var_names = Vec::new();
    let mut q_bounds = Vec::new();
    loop {
        let in_line = read_line(in_file)?;
        if in_line.contains("# ----") {
            break;
        }
        let rest: String = in_line.chars().skip(2).collect();
        let mut split_line = rest.split('[');
        let name = split_line.next().unwrap_or("");
        let bounds = split_line.next().ok_or(index_error(1, 1))?;
        var_names.push(py_strip(name).to_owned());
        // Python drops the last two characters, `]` and the newline.
        let n = bounds.chars().count();
        let text: String = bounds.chars().take(n.saturating_sub(2)).collect();
        let parsed = if text.is_empty() {
            Vec::new()
        } else {
            text.split(',')
                .map(|item| {
                    py_float(item).ok_or_else(|| Error::Value {
                        message: format!("could not convert string to float: '{item}'"),
                    })
                })
                .collect::<Result<Vec<f64>>>()?
        };
        q_bounds.push(parsed);
    }
    Ok((var_names, q_bounds))
}

/// Python's `re.split('[ ]+|[\t]+', line)`: splits at each run of spaces
/// and each run of tabs.
fn split_fields(line: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0;
    let bytes = line.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        if c == b' ' || c == b'\t' {
            fields.push(&line[start..i]);
            while i < bytes.len() && bytes[i] == c {
                i += 1;
            }
            start = i;
        } else {
            i += 1;
        }
    }
    fields.push(&line[start..]);
    fields
}

/// The points' names and integer values from the rest of a `.qdat` file,
/// skipping lines that start with `#`. Mirrors
/// `DataUtils.ReadQuantExamples`.
///
/// `Error::Value` for a value that is not an integer; `Error::Os` when
/// reading fails.
pub fn read_quant_examples<R: BufRead>(in_file: &mut R) -> Result<(Vec<String>, Vec<Vec<i64>>)> {
    let mut examples = Vec::new();
    let mut names = Vec::new();
    loop {
        let in_line = read_line(in_file)?;
        if in_line.is_empty() {
            break;
        }
        if in_line.starts_with('#') {
            continue;
        }
        let res_arr = split_fields(&in_line);
        if res_arr.len() > 1 {
            let values = res_arr[1..]
                .iter()
                .map(|x| {
                    py_int(x).ok_or_else(|| Error::Value {
                        message: format!("invalid literal for int() with base 10: {x:?}"),
                    })
                })
                .collect::<Result<Vec<i64>>>()?;
            examples.push(values);
            names.push(res_arr[0].to_owned());
        }
    }
    Ok((names, examples))
}

/// A field as RDKit's readers take it: an integer, else a float, else the
/// text.
fn general_value(d: &str) -> Value {
    if let Some(i) = py_int(d) {
        Value::Int(i)
    } else if let Some(x) = py_float(d) {
        Value::Float(x)
    } else {
        Value::Str(d.to_owned())
    }
}

/// The points' names and values (integers, floats, or text) from the rest
/// of a `.dat` file, skipping lines that start with `#`. As in Python, the
/// last field of each line, after its last separator, is dropped. Mirrors
/// `DataUtils.ReadGeneralExamples`.
///
/// `Error::Os` when reading fails.
pub fn read_general_examples<R: BufRead>(
    in_file: &mut R,
) -> Result<(Vec<String>, Vec<Vec<Value>>)> {
    let mut examples = Vec::new();
    let mut names = Vec::new();
    loop {
        let in_line = read_line(in_file)?;
        if in_line.is_empty() {
            break;
        }
        if in_line.starts_with('#') {
            continue;
        }
        let mut res_arr = split_fields(&in_line);
        res_arr.pop();
        if res_arr.len() > 1 {
            examples.push(res_arr[1..].iter().map(|d| general_value(d)).collect());
            names.push(res_arr[0].to_owned());
        }
    }
    Ok((names, examples))
}

fn open(file_name: &Path) -> Result<BufReader<fs::File>> {
    Ok(BufReader::new(fs::File::open(file_name).map_err(os_error)?))
}

/// A quantized data set from a `.qdat` file. Mirrors
/// `DataUtils.BuildQuantDataSet`.
///
/// Errors as [`read_vars`], [`read_quant_examples`], and
/// [`MLQuantDataSet::new_with_params`].
pub fn build_quant_data_set(file_name: impl AsRef<Path>) -> Result<MLQuantDataSet> {
    let mut in_file = open(file_name.as_ref())?;
    let (var_names, q_bounds) = read_vars(&mut in_file)?;
    let (pt_names, examples) = read_quant_examples(&mut in_file)?;
    MLQuantDataSet::new_with_params(
        examples,
        MLDataSetParams {
            q_bounds: Some(q_bounds),
            var_names: Some(var_names),
            pt_names: Some(pt_names),
            ..Default::default()
        },
    )
}

/// A data set from a `.dat` file. Mirrors `DataUtils.BuildDataSet`.
///
/// Errors as [`read_vars`], [`read_general_examples`], and
/// [`MLDataSet::new_with_params`].
pub fn build_data_set(file_name: impl AsRef<Path>) -> Result<MLDataSet> {
    let mut in_file = open(file_name.as_ref())?;
    let (var_names, q_bounds) = read_vars(&mut in_file)?;
    let (pt_names, examples) = read_general_examples(&mut in_file)?;
    MLDataSet::new_with_params(
        examples,
        MLDataSetParams {
            q_bounds: Some(q_bounds),
            var_names: Some(var_names),
            pt_names: Some(pt_names),
            ..Default::default()
        },
    )
}

/// The number of values each variable takes in the examples: the number
/// of its quantization bounds when it has some, else its largest value
/// plus 1 when every value is an integer, else 0. Variable `i` is column
/// `order[i]` of the examples. Mirrors `DataUtils.CalcNPossibleUsingMap(data,
/// order, qBounds)`.
///
/// Errors as [`calc_n_possible_using_map_with_n_q_bounds`].
pub fn calc_n_possible_using_map(
    data: &[Vec<Value>],
    order: &[usize],
    q_bounds: &[Vec<f64>],
) -> Result<Vec<i64>> {
    calc_n_possible_using_map_with_n_q_bounds(data, order, q_bounds, None)
}

/// As [`calc_n_possible_using_map`]; a variable whose `n_q_bounds` entry is
/// nonzero takes 0. Mirrors `DataUtils.CalcNPossibleUsingMap(data, order,
/// qBounds, nQBounds)`, whose `silent` argument, when false, prints the
/// arguments.
///
/// `Error::Value` unless `order` is as long as the bounds or `n_q_bounds`
/// (Python's assertion), or for an infinite or NaN float; `Error::Index`
/// for a column past an example.
pub fn calc_n_possible_using_map_with_n_q_bounds(
    data: &[Vec<Value>],
    order: &[usize],
    q_bounds: &[Vec<f64>],
    n_q_bounds: Option<&[i64]>,
) -> Result<Vec<i64>> {
    let n_q_bounds = n_q_bounds.filter(|n| !n.is_empty());
    let ok = (!q_bounds.is_empty() && order.len() == q_bounds.len())
        || n_q_bounds.is_some_and(|n| order.len() == n.len());
    if !ok {
        return Err(Error::Value {
            message: "order/qBounds mismatch".to_owned(),
        });
    }
    let n_vars = order.len();
    let mut n_possible = vec![-1i64; n_vars];
    let mut cols: Vec<usize> = (0..n_vars).collect();
    for i in 0..n_vars {
        if n_q_bounds.is_some_and(|n| n[i] != 0) {
            n_possible[i] = -1;
            cols.retain(|&c| c != i);
        } else {
            let bounds = q_bounds.get(i).ok_or(index_error(i, q_bounds.len()))?;
            if !bounds.is_empty() {
                n_possible[i] = bounds.len() as i64;
                cols.retain(|&c| c != i);
            }
        }
    }
    for pt in data {
        for col in cols.clone() {
            let d = pt
                .get(order[col])
                .ok_or(index_error(order[col], pt.len()))?;
            let as_int = match d {
                Value::Int(i) => Some(*i),
                Value::Float(x) => {
                    if !x.is_finite() {
                        return Err(Error::Value {
                            message: format!("cannot convert float {} to integer", float_repr(*x)),
                        });
                    }
                    (x.trunc() == *x).then_some(x.trunc() as i64)
                }
                Value::Str(_) => None,
            };
            match as_int {
                Some(v) => n_possible[col] = n_possible[col].max(v),
                None => {
                    n_possible[col] = -1;
                    cols.retain(|&c| c != col);
                }
            }
        }
    }
    Ok(n_possible.into_iter().map(|x| x + 1).collect())
}

/// The items of `vect` at `ensemble_ids`. Mirrors
/// `DataUtils.TakeEnsemble(vect, ensembleIds)`.
///
/// `Error::Index` for an id past the end.
pub fn take_ensemble<T: Clone>(vect: &[T], ensemble_ids: &[usize]) -> Result<Vec<T>> {
    take_ensemble_with_is_data_vect(vect, ensemble_ids, false)
}

/// As [`take_ensemble`]; with `is_data_vect`, `vect` is a data row, its
/// first item (the name) and last (the result) kept around the items at
/// `ensemble_ids` counted after the first. Mirrors
/// `DataUtils.TakeEnsemble(vect, ensembleIds, isDataVect)`.
pub fn take_ensemble_with_is_data_vect<T: Clone>(
    vect: &[T],
    ensemble_ids: &[usize],
    is_data_vect: bool,
) -> Result<Vec<T>> {
    let at = |x: usize| vect.get(x).cloned().ok_or(index_error(x, vect.len()));
    if is_data_vect {
        let mut res = vec![at(0)?];
        for &x in ensemble_ids {
            res.push(at(x + 1)?);
        }
        res.push(at(vect.len().wrapping_sub(1))?);
        Ok(res)
    } else {
        ensemble_ids.iter().map(|&x| at(x)).collect()
    }
}

/// Parameters for [`text_to_data_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TextToDataParams {
    /// Columns, by name, to leave out.
    pub ignore_cols: Vec<String>,
    /// The only columns, by name, to keep, in this order; when given and
    /// not empty, `ignore_cols` is not used. A name the header lacks takes
    /// the last column, as Python's -1 index does.
    pub only_cols: Option<Vec<String>>,
}

/// A data set from rows of text, such as a CSV reader gives: a header of
/// column names, then a row per point, its name first. Values become
/// integers, else floats, else text; empty rows are skipped. Mirrors
/// `DataUtils.TextToData(reader)`.
///
/// Errors as [`text_to_data_with_params`].
pub fn text_to_data<I: IntoIterator<Item = Vec<String>>>(reader: I) -> Result<MLDataSet> {
    text_to_data_with_params(reader, &TextToDataParams::default())
}

/// Mirrors `DataUtils.TextToData(reader, ignoreCols, onlyCols)`.
///
/// `Error::Value` for no header, or a row not as long as the header, as
/// Python raises; `Error::Index` with no columns kept or no rows; errors
/// from [`MLDataSet::new_with_params`] pass through.
pub fn text_to_data_with_params<I: IntoIterator<Item = Vec<String>>>(
    reader: I,
    params: &TextToDataParams,
) -> Result<MLDataSet> {
    let mut reader = reader.into_iter();
    let var_names = reader.next().ok_or_else(|| Error::Value {
        message: "there is no header".to_owned(),
    })?;
    let n_cols = var_names.len();
    let keep_cols: Vec<usize> = match params.only_cols.as_ref().filter(|o| !o.is_empty()) {
        None => (0..n_cols)
            .filter(|&i| !params.ignore_cols.contains(&var_names[i]))
            .collect(),
        Some(only_cols) => {
            let mut keep: Vec<Option<usize>> = vec![None; only_cols.len()];
            for (i, name) in var_names.iter().enumerate() {
                if let Some(at) = only_cols.iter().position(|c| c == name) {
                    keep[at] = Some(i);
                }
            }
            keep.into_iter()
                .map(|k| k.or(n_cols.checked_sub(1)).ok_or(index_error(0, 0)))
                .collect::<Result<_>>()?
        }
    };
    let kept_names: Vec<String> = keep_cols.iter().map(|&x| var_names[x].clone()).collect();
    let n_vars = kept_names.len();
    let mut vals = Vec::new();
    let mut pt_names = Vec::new();
    for split_line in reader {
        if split_line.is_empty() {
            continue;
        }
        if split_line.len() != n_cols {
            return Err(Error::Value {
                message: "unequal line lengths".to_owned(),
            });
        }
        let tmp: Vec<&String> = keep_cols.iter().map(|&x| &split_line[x]).collect();
        pt_names.push(tmp.first().ok_or(index_error(0, 0))?.to_string());
        vals.push(
            tmp[1..n_vars.max(1)]
                .iter()
                .map(|s| general_value(s))
                .collect(),
        );
    }
    MLDataSet::new_with_params(
        vals,
        MLDataSetParams {
            var_names: Some(kept_names),
            pt_names: Some(pt_names),
            ..Default::default()
        },
    )
}

/// Python's `csv.reader` with the default dialect and `delimiter`, over
/// text whose line ends are all `\n`.
fn read_csv(text: &str, delimiter: char) -> Vec<Vec<String>> {
    #[derive(PartialEq)]
    enum State {
        StartRecord,
        StartField,
        InField,
        InQuotedField,
        QuoteInQuotedField,
    }
    let mut rows = Vec::new();
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut state = State::StartRecord;
    for c in text.chars() {
        if state == State::StartRecord {
            if c == '\n' {
                rows.push(Vec::new());
                continue;
            }
            state = State::StartField;
        }
        match state {
            State::StartField => {
                if c == '\n' {
                    fields.push(std::mem::take(&mut field));
                    rows.push(std::mem::take(&mut fields));
                    state = State::StartRecord;
                } else if c == '"' {
                    state = State::InQuotedField;
                } else if c == delimiter {
                    fields.push(std::mem::take(&mut field));
                } else {
                    field.push(c);
                    state = State::InField;
                }
            }
            State::InField => {
                if c == '\n' {
                    fields.push(std::mem::take(&mut field));
                    rows.push(std::mem::take(&mut fields));
                    state = State::StartRecord;
                } else if c == delimiter {
                    fields.push(std::mem::take(&mut field));
                    state = State::StartField;
                } else {
                    field.push(c);
                }
            }
            State::InQuotedField => {
                if c == '"' {
                    state = State::QuoteInQuotedField;
                } else {
                    field.push(c);
                }
            }
            State::QuoteInQuotedField => {
                if c == '"' {
                    field.push('"');
                    state = State::InQuotedField;
                } else if c == delimiter {
                    fields.push(std::mem::take(&mut field));
                    state = State::StartField;
                } else if c == '\n' {
                    fields.push(std::mem::take(&mut field));
                    rows.push(std::mem::take(&mut fields));
                    state = State::StartRecord;
                } else {
                    field.push(c);
                    state = State::InField;
                }
            }
            State::StartRecord => unreachable!("handled above"),
        }
    }
    if state != State::StartRecord {
        fields.push(field);
        rows.push(fields);
    }
    rows
}

/// A data set from a CSV file (a name ending in `.csv`, in any case) or a
/// tab-separated one, as [`text_to_data`] reads its rows. Mirrors
/// `DataUtils.TextFileToData(fName)`.
///
/// `Error::Os` when the file cannot be read, `Error::Value` when it is not
/// UTF-8; otherwise as [`text_to_data_with_params`].
pub fn text_file_to_data(f_name: impl AsRef<Path>) -> Result<MLDataSet> {
    text_file_to_data_with_only_cols(f_name, None)
}

/// As [`text_file_to_data`], keeping only the columns `only_cols` as
/// [`TextToDataParams`] describes. Mirrors `DataUtils.TextFileToData(fName,
/// onlyCols)`.
pub fn text_file_to_data_with_only_cols(
    f_name: impl AsRef<Path>,
    only_cols: Option<Vec<String>>,
) -> Result<MLDataSet> {
    let path = f_name.as_ref();
    let name = path.to_string_lossy();
    let ext = name.rsplit('.').next().unwrap_or("");
    let delimiter = if ext.to_uppercase() == "CSV" {
        ','
    } else {
        '\t'
    };
    // Python's text files turn `\r\n` and `\r` into `\n` before the CSV
    // reader sees them.
    let mut in_file = open(path)?;
    let mut text = String::new();
    loop {
        let line = read_line_universal(&mut in_file).map_err(|e| {
            if e.kind() == std::io::ErrorKind::InvalidData {
                Error::Value {
                    message: e.to_string(),
                }
            } else {
                os_error(e)
            }
        })?;
        if line.is_empty() {
            break;
        }
        text.push_str(&line);
    }
    text_to_data_with_params(
        read_csv(&text, delimiter),
        &TextToDataParams {
            only_cols,
            ..Default::default()
        },
    )
}

/// Seeds the random numbers, from the first of the two seeds. Mirrors
/// `DataUtils.InitRandomNumbers`, which seeds `RDRandom` and Python's
/// `random`, the same generator.
pub fn init_random_numbers(seed: (i64, i64)) {
    rd_random::seed(seed.0);
}

/// Parameters for [`filter_data_with_params`] and
/// [`filter_data_indices_only_with_params`], with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FilterDataParams {
    /// The column of the target value, counting from the end when
    /// negative.
    pub col: i64,
    /// The rows to filter, by index; `None` or empty for all.
    pub indices_to_use: Option<Vec<usize>>,
}

impl Default for FilterDataParams {
    fn default() -> FilterDataParams {
        FilterDataParams {
            col: -1,
            indices_to_use: None,
        }
    }
}

/// Python's index into a row.
fn column(row: &[f64], col: i64) -> Result<f64> {
    let n = row.len();
    let idx = if col < 0 { n as i64 + col } else { col };
    usize::try_from(idx)
        .ok()
        .and_then(|i| row.get(i).copied())
        .ok_or_else(|| Error::Value {
            message: "target column index out of range".to_owned(),
        })
}

fn zero_division() -> Error {
    Error::Value {
        message: "float division by zero".to_owned(),
    }
}

/// Python's `seq[:n]` for a possibly negative `n`.
fn head<T: Copy>(seq: &[T], n: i64) -> Vec<T> {
    let len = seq.len() as i64;
    let end = if n < 0 { (len + n).max(0) } else { n.min(len) };
    seq[..end as usize].to_vec()
}

/// The shared work of the filters: the rows (sorted by the column), the
/// sort order, and, in a random order, whether to keep each.
#[allow(clippy::type_complexity)]
fn filter(
    in_data: &[Vec<f64>],
    val: f64,
    frac: f64,
    params: &FilterDataParams,
) -> Result<(Vec<Vec<f64>>, Vec<usize>, Vec<(usize, bool)>)> {
    if !(0.0..=1.0).contains(&frac) {
        return Err(Error::Value {
            message: "filter fraction out of bounds".to_owned(),
        });
    }
    column(
        in_data.first().ok_or_else(|| Error::Value {
            message: "target column index out of range".to_owned(),
        })?,
        params.col,
    )?;
    let tmp: Vec<Vec<f64>> = match params.indices_to_use.as_ref().filter(|i| !i.is_empty()) {
        Some(indices) => indices
            .iter()
            .map(|&x| in_data.get(x).cloned().ok_or(index_error(x, in_data.len())))
            .collect::<Result<_>>()?,
        None => in_data.to_vec(),
    };
    let n_orig = tmp.len();
    let keys: Vec<f64> = tmp
        .iter()
        .map(|row| column(row, params.col))
        .collect::<Result<_>>()?;
    let mut sort_order: Vec<usize> = (0..n_orig).collect();
    sort_order.sort_by(|&a, &b| keys[a].partial_cmp(&keys[b]).unwrap_or(Ordering::Equal));
    let tmp: Vec<Vec<f64>> = sort_order.iter().map(|&x| tmp[x].clone()).collect();
    let keys: Vec<f64> = sort_order.iter().map(|&x| keys[x]).collect();
    let mut start = 0;
    while start < n_orig && keys[start] != val {
        start += 1;
    }
    if start >= n_orig {
        return Err(Error::Value {
            message: format!("target value ({}) not found in data", float_repr(val)),
        });
    }
    let mut finish = start + 1;
    while finish < n_orig && keys[finish] == val {
        finish += 1;
    }
    let n_with_val = (finish - start) as i64;
    let n_others = n_orig as i64 - n_with_val;
    let curr_frac = n_with_val as f64 / n_orig as f64;
    let (n_tgt_final, n_others_final);
    if curr_frac < frac {
        let mut n_tgt = n_with_val;
        let mut n_final = (n_with_val as f64 / frac).round_ties_even() as i64;
        n_others_final = n_final - n_tgt;
        loop {
            if n_final == 0 {
                return Err(zero_division());
            }
            if (n_tgt as f64 / n_final as f64).partial_cmp(&frac) != Some(Ordering::Greater) {
                break;
            }
            n_tgt -= 1;
            n_final -= 1;
        }
        n_tgt_final = n_tgt;
    } else {
        let mut n_oth = n_others;
        if frac == 1.0 {
            return Err(zero_division());
        }
        let mut n_final = (n_others as f64 / (1.0 - frac)).round_ties_even() as i64;
        n_tgt_final = n_final - n_oth;
        loop {
            if n_final == 0 {
                return Err(zero_division());
            }
            if (n_tgt_final as f64 / n_final as f64).partial_cmp(&frac) != Some(Ordering::Less) {
                break;
            }
            n_oth -= 1;
            n_final -= 1;
        }
        n_others_final = n_oth;
    }
    let others: Vec<usize> = (0..start).chain(finish..n_orig).collect();
    let others_take = permutation(n_others as usize);
    let others: Vec<usize> = head(&others_take, n_others_final)
        .into_iter()
        .map(|x| others[x])
        .collect();
    let targets: Vec<usize> = (start..finish).collect();
    let targets_take = permutation(n_with_val as usize);
    let targets: Vec<usize> = head(&targets_take, n_tgt_final)
        .into_iter()
        .map(|x| targets[x])
        .collect();
    let mut to_keep = vec![false; n_orig];
    for i in targets.into_iter().chain(others) {
        to_keep[i] = true;
    }
    let order = permutation(n_orig)
        .into_iter()
        .map(|i| (i, to_keep[i]))
        .collect();
    Ok((tmp, sort_order, order))
}

/// Samples the rows so that rows whose last column is `val` make up
/// `frac` of them, as nearly as can be: returns the rows kept and those
/// rejected, each in a random order. Mirrors `DataUtils.FilterData(inData,
/// val, frac)`.
///
/// Errors as [`filter_data_with_params`].
#[allow(clippy::type_complexity)]
pub fn filter_data(
    in_data: &[Vec<f64>],
    val: f64,
    frac: f64,
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    filter_data_with_params(in_data, val, frac, &FilterDataParams::default())
}

/// Mirrors `DataUtils.FilterData(inData, val, frac, col, indicesToUse)`.
///
/// `Error::Value` for a fraction outside 0 to 1, a column past the first
/// row, or a value no row has, as Python raises, and where Python divides
/// by zero (a fraction no sample can meet); `Error::Index` for an index to
/// use past the rows.
#[allow(clippy::type_complexity)]
pub fn filter_data_with_params(
    in_data: &[Vec<f64>],
    val: f64,
    frac: f64,
    params: &FilterDataParams,
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    let (tmp, _, order) = filter(in_data, val, frac, params)?;
    let (mut res, mut rej) = (Vec::new(), Vec::new());
    for (i, keep) in order {
        if keep {
            res.push(tmp[i].clone());
        } else {
            rej.push(tmp[i].clone());
        }
    }
    Ok((res, rej))
}

/// As [`filter_data`], returning the rows' indices. Mirrors
/// `DataUtils.FilterData(inData, val, frac, indicesOnly=1)`.
///
/// Errors as [`filter_data_with_params`].
pub fn filter_data_indices_only(
    in_data: &[Vec<f64>],
    val: f64,
    frac: f64,
) -> Result<(Vec<usize>, Vec<usize>)> {
    filter_data_indices_only_with_params(in_data, val, frac, &FilterDataParams::default())
}

/// As [`filter_data_with_params`], returning the rows' indices (into
/// `in_data`). Mirrors `DataUtils.FilterData(inData, val, frac, col,
/// indicesToUse, indicesOnly=1)`.
pub fn filter_data_indices_only_with_params(
    in_data: &[Vec<f64>],
    val: f64,
    frac: f64,
    params: &FilterDataParams,
) -> Result<(Vec<usize>, Vec<usize>)> {
    let (_, sort_order, order) = filter(in_data, val, frac, params)?;
    let indices = params.indices_to_use.as_ref().filter(|i| !i.is_empty());
    let (mut res, mut rej) = (Vec::new(), Vec::new());
    for (i, keep) in order {
        let idx = match indices {
            Some(indices) => indices[sort_order[i]],
            None => sort_order[i],
        };
        if keep {
            res.push(idx);
        } else {
            rej.push(idx);
        }
    }
    Ok((res, rej))
}

/// Parameters for [`count_results_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CountResultsParams {
    /// The column counted, counting from the end when negative.
    pub col: i64,
    /// Bounds that bin the values: a value goes in the first bin whose
    /// bound it is below, or past the last. `None` or empty counts the
    /// values themselves.
    pub bounds: Option<Vec<f64>>,
}

impl Default for CountResultsParams {
    fn default() -> CountResultsParams {
        CountResultsParams {
            col: -1,
            bounds: None,
        }
    }
}

/// How often each value of the last column occurs, in the order the values
/// first occur. Mirrors `DataUtils.CountResults(inData)`.
///
/// Errors as [`count_results_with_params`].
pub fn count_results(in_data: &[Vec<f64>]) -> Result<Vec<(f64, usize)>> {
    count_results_with_params(in_data, &CountResultsParams::default())
}

/// Mirrors `DataUtils.CountResults(inData, col, bounds)`; with bounds, the
/// keys are the bins' indices.
///
/// `Error::Value` for a column past a row.
pub fn count_results_with_params(
    in_data: &[Vec<f64>],
    params: &CountResultsParams,
) -> Result<Vec<(f64, usize)>> {
    let bounds = params.bounds.as_ref().filter(|b| !b.is_empty());
    let mut counts: Vec<(f64, usize)> = Vec::new();
    for p in in_data {
        let act = column(p, params.col)?;
        let r = match bounds {
            None => act,
            Some(bounds) => bounds.iter().position(|&b| act < b).unwrap_or(bounds.len()) as f64,
        };
        match counts.iter_mut().find(|(k, _)| *k == r) {
            Some(entry) => entry.1 += 1,
            None => counts.push((r, 1)),
        }
    }
    Ok(counts)
}

/// Replaces the results of a data set with random values between their
/// smallest and largest. Mirrors `DataUtils.RandomizeActivities(dataSet)`.
///
/// Errors as [`randomize_activities_with_shuffle`].
pub fn randomize_activities(data_set: &mut MLDataSet) -> Result<()> {
    randomize_activities_with_shuffle(data_set, false)
}

/// As [`randomize_activities`]; with `shuffle`, shuffles the results
/// instead, which keeps the number of each. Mirrors
/// `DataUtils.RandomizeActivities(dataSet, shuffle)`, whose `runDetails`
/// argument, a composite model's run record, is only marked.
///
/// `Error::Value` for a data set with more than one result column, or,
/// without `shuffle`, no points or a result that is not a number;
/// `Error::Index` when the data set has fewer rows than points.
pub fn randomize_activities_with_shuffle(data_set: &mut MLDataSet, shuffle: bool) -> Result<()> {
    if data_set.get_n_results() > 1 {
        return Err(Error::Value {
            message: "the data set has more than one result column".to_owned(),
        });
    }
    let n_pts = data_set.get_n_pts();
    let mut acts: Vec<Value> = data_set
        .get_results()
        .into_iter()
        .filter_map(|r| r.into_iter().next_back())
        .collect();
    if shuffle {
        rd_random::shuffle_with_random(&mut acts, rd_random::random)?;
    } else {
        let nums = acts
            .iter()
            .map(|a| {
                a.as_f64().ok_or_else(|| Error::Value {
                    message: "the results must be numbers".to_owned(),
                })
            })
            .collect::<Result<Vec<f64>>>()?;
        let empty = || Error::Value {
            message: "the data set has no results".to_owned(),
        };
        let max = nums.iter().copied().reduce(f64::max).ok_or_else(empty)?;
        let min = nums.iter().copied().reduce(f64::min).ok_or_else(empty)?;
        let rng = max - min;
        acts = (0..n_pts)
            .map(|_| Value::Float(min + rd_random::random() * rng))
            .collect();
    }
    for (i, act) in acts.into_iter().enumerate().take(n_pts) {
        let mut tmp = data_set.get(i)?;
        let last = tmp.len() - 1;
        tmp[last] = act;
        data_set.set(i, tmp)?;
    }
    Ok(())
}
