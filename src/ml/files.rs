// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.files`, a Python module: reading data files with
//! comments.
//!
//! Python's comment and line-trailer arguments are regular expressions;
//! here they are plain text.

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::Path;

use crate::utils::fileutils::read_line;
use crate::utils::{py_float, py_int, py_split, py_strip};
use crate::{Error, Result};

/// Parameters for [`ReFile::new_with_params`], with Python's defaults (A2).
/// Python's `mode` argument, which must be `'r'` for the methods to work,
/// is not taken.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReFileParams {
    /// The text that starts a comment, which runs to the end of the line.
    pub comment: String,
    /// More text that ends a line's content; `None` or empty for none.
    pub trailer: Option<String>,
}

impl Default for ReFileParams {
    fn default() -> ReFileParams {
        ReFileParams {
            comment: "#".to_owned(),
            trailer: Some("\n".to_owned()),
        }
    }
}

/// A file read line by line with comments and surrounding whitespace
/// stripped and blank lines skipped. Mirrors `files.ReFile`.
#[derive(Debug)]
pub struct ReFile {
    in_file: BufReader<File>,
    markers: Vec<String>,
}

fn os_error(e: std::io::Error) -> Error {
    Error::Os {
        message: e.to_string(),
    }
}

impl ReFile {
    /// Opens a file, with `#` comments. Mirrors `ReFile(fileName)`.
    ///
    /// `Error::Os` when the file cannot be opened.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<ReFile> {
        ReFile::new_with_params(file_name, &ReFileParams::default())
    }

    /// Mirrors `ReFile(fileName, mode, comment, trailer)`.
    pub fn new_with_params<P: AsRef<Path>>(file_name: P, params: &ReFileParams) -> Result<ReFile> {
        let mut markers = vec![params.comment.clone()];
        if let Some(trailer) = params.trailer.as_ref().filter(|t| !t.is_empty()) {
            markers.push(trailer.clone());
        }
        Ok(ReFile {
            in_file: BufReader::new(File::open(file_name).map_err(os_error)?),
            markers,
        })
    }

    /// A line's content: the text before the first comment or trailer,
    /// stripped of whitespace.
    fn content(&self, line: &str) -> String {
        let end = self
            .markers
            .iter()
            .filter(|m| !m.is_empty())
            .filter_map(|m| line.find(m.as_str()))
            .min();
        // An empty comment, a regular expression matching everywhere,
        // splits before the first character.
        let end = if self.markers.iter().any(String::is_empty) {
            Some(0)
        } else {
            end
        };
        py_strip(&line[..end.unwrap_or(line.len())]).to_owned()
    }

    /// The next line with content, or an empty string at the end of the
    /// file. Mirrors `ReFile.readline`.
    ///
    /// `Error::Os` when reading fails, and `Error::Value` for a line that is
    /// not UTF-8.
    pub fn readline(&mut self) -> Result<String> {
        loop {
            let line = read_line(&mut self.in_file)?;
            if line.is_empty() {
                return Ok(String::new());
            }
            let result = self.content(&line);
            if !result.is_empty() {
                return Ok(result);
            }
        }
    }

    /// The remaining lines with content. Mirrors `ReFile.readlines`.
    ///
    /// Errors as [`readline`](Self::readline).
    pub fn readlines(&mut self) -> Result<Vec<String>> {
        let mut res = Vec::new();
        loop {
            let line = read_line(&mut self.in_file)?;
            if line.is_empty() {
                return Ok(res);
            }
            let result = self.content(&line);
            if !result.is_empty() {
                res.push(result);
            }
        }
    }

    /// Goes back to the start of the file. Mirrors `ReFile.rewind`.
    ///
    /// `Error::Os` when seeking fails.
    pub fn rewind(&mut self) -> Result<()> {
        self.in_file.seek(SeekFrom::Start(0)).map_err(os_error)?;
        Ok(())
    }
}

/// How [`read_data_file_with_params`] reads values. Mirrors the types
/// `files.ReadDataFile`'s `dataType` takes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum DataType {
    /// Floats, as Python's `float`, `numpy.float32`, and `numpy.float64`
    /// read them.
    #[default]
    Float,
    /// Integers, as every other type reads them (Python's `int`).
    Int,
}

/// Parameters for [`read_data_file_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadDataFileParams {
    /// The text that starts a comment.
    pub comment: String,
    /// The column of the dependent variable.
    pub dep_var_col: usize,
    /// How values are read.
    pub data_type: DataType,
}

impl Default for ReadDataFileParams {
    fn default() -> ReadDataFileParams {
        ReadDataFileParams {
            comment: "#".to_owned(),
            dep_var_col: 0,
            data_type: DataType::Float,
        }
    }
}

/// The independent variables (rows) and the dependent variable (the first
/// column) of a file of whitespace-separated numbers with `#` comments.
/// Mirrors `files.ReadDataFile(fileName)`.
///
/// Python's raises `TypeError` in RDKit 2026.09.1 (it assigns a `map`
/// object to a NumPy row); this reads what it means to.
///
/// `Error::Os` when the file cannot be read; `Error::Value` for a value
/// that is not a number, a line not as long as the first, or no lines;
/// `Error::Index` for a line without the dependent column.
pub fn read_data_file(file_name: impl AsRef<Path>) -> Result<(Vec<Vec<f64>>, Vec<f64>)> {
    read_data_file_with_params(file_name, &ReadDataFileParams::default())
}

/// Mirrors `files.ReadDataFile(fileName, comment, depVarCol, dataType)`,
/// except that Python's ignores `comment` (its file always takes `#`).
pub fn read_data_file_with_params(
    file_name: impl AsRef<Path>,
    params: &ReadDataFileParams,
) -> Result<(Vec<Vec<f64>>, Vec<f64>)> {
    let mut in_file = ReFile::new_with_params(
        file_name,
        &ReFileParams {
            comment: params.comment.clone(),
            ..Default::default()
        },
    )?;
    let data_lines = in_file.readlines()?;
    let first = data_lines.first().ok_or_else(|| Error::Value {
        message: "the file has no data".to_owned(),
    })?;
    let n_ind_vars = py_split(first).len().saturating_sub(1);
    let conv = |s: &str| -> Result<f64> {
        let v = match params.data_type {
            DataType::Float => py_float(s),
            DataType::Int => py_int(s).map(|i| i as f64),
        };
        v.ok_or_else(|| Error::Value {
            message: format!("could not convert {s:?} to a number"),
        })
    };
    let mut ind_var_mat = Vec::with_capacity(data_lines.len());
    let mut dep_var_vect = Vec::with_capacity(data_lines.len());
    for line in &data_lines {
        let mut split_line = py_split(line);
        if params.dep_var_col >= split_line.len() {
            return Err(Error::index(params.dep_var_col, split_line.len()));
        }
        dep_var_vect.push(conv(split_line.remove(params.dep_var_col))?);
        if split_line.len() != n_ind_vars {
            return Err(Error::Value {
                message: "the lines have different numbers of values".to_owned(),
            });
        }
        ind_var_mat.push(split_line.into_iter().map(conv).collect::<Result<_>>()?);
    }
    Ok((ind_var_mat, dep_var_vect))
}
