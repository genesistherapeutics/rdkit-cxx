// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Data.Transforms`, a Python module: the data
//! transformations offered to users, which RDKit never implemented.

use crate::ml::data::ml_data::MLDataSet;

/// A data transformation, which changes a data set.
/// A Rust addition, naming the functions of `GetAvailTransforms`.
pub type Transform = fn(&mut MLDataSet);

fn center_t_form(_data_set: &mut MLDataSet) {}

fn normalize_t_form(_data_set: &mut MLDataSet) {}

fn standard_t_form(_data_set: &mut MLDataSet) {}

/// The available transformations, as `(name, function, description)`.
/// Their functions, as in Python, do nothing. Mirrors
/// `Transforms.GetAvailTransforms`.
pub fn get_avail_transforms() -> Vec<(&'static str, Transform, &'static str)> {
    vec![
        ("Center", center_t_form, "translates so that mean(x)=0"),
        ("Normalize", normalize_t_form, "scales so that dot(x,x)=1"),
        ("Standardize", standard_t_form, "scales so that dev(x)=0"),
    ]
}
