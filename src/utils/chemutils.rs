// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.utils.chemutils`, a Python module: utilities with
//! chemical know-how.

use crate::utils::{py_float, py_int};
use crate::{Error, Result};

/// The elements of a simple composition and their counts, in order:
/// `Fe3Al` gives `[("Fe", 3.0), ("Al", 1.0)]`. Each element is a capital
/// letter and an optional lowercase one, and its count the digits and
/// periods after it (1 when there are none). Mirrors
/// `chemutils.SplitComposition`, whose counts of 1 are integers.
///
/// `Error::Value` for a count that is not a number, such as `1.2.3`, as
/// Python raises.
pub fn split_composition(comp_str: &str) -> Result<Vec<(String, f64)>> {
    let chars: Vec<char> = comp_str.chars().collect();
    let mut res = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        if !chars[i].is_ascii_uppercase() {
            i += 1;
            continue;
        }
        let start = i;
        i += 1;
        if i < chars.len() && chars[i].is_ascii_lowercase() {
            i += 1;
        }
        let element: String = chars[start..i].iter().collect();
        let count_start = i;
        while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
            i += 1;
        }
        let count: String = chars[count_start..i].iter().collect();
        let count = if count.is_empty() {
            1.0
        } else {
            py_float(&count).ok_or_else(|| Error::Value {
                message: format!("could not convert string to float: '{count}'"),
            })?
        };
        res.push((element, count));
    }
    Ok(res)
}

/// Parameters for [`config_to_num_electrons_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ConfigToNumElectronsParams {
    /// Leaves out full d shells, unless the configuration has only one
    /// shell after its core.
    pub ignore_full_d: bool,
    /// Leaves out full f shells, likewise.
    pub ignore_full_f: bool,
}

/// The number of electrons in an electron configuration such as
/// `[Ne] 3s^2 3p^4`, the shells after the first space-separated item.
/// Mirrors `chemutils.ConfigToNumElectrons(config)`.
///
/// `Error::Index` for a shell without `^`, and `Error::Value` for one whose
/// count is not an integer, as Python raises.
pub fn config_to_num_electrons(config: &str) -> Result<i64> {
    config_to_num_electrons_with_params(config, &ConfigToNumElectronsParams::default())
}

/// Mirrors `chemutils.ConfigToNumElectrons(config, ignoreFullD,
/// ignoreFullF)`.
pub fn config_to_num_electrons_with_params(
    config: &str,
    params: &ConfigToNumElectronsParams,
) -> Result<i64> {
    let arr: Vec<&str> = config.split(' ').collect();
    let mut n_el = 0;
    for shell in &arr[1..] {
        let l: Vec<&str> = shell.split('^').collect();
        let count = l.get(1).ok_or(Error::Index { index: 1, len: 1 })?;
        let mut incr = py_int(count).ok_or_else(|| Error::Value {
            message: format!("invalid literal for int() with base 10: '{count}'"),
        })?;
        if params.ignore_full_f && incr == 14 && l[0].contains('f') && arr.len() > 2 {
            incr = 0;
        }
        if params.ignore_full_d && incr == 10 && l[0].contains('d') && arr.len() > 2 {
            incr = 0;
        }
        n_el += incr;
    }
    Ok(n_el)
}
