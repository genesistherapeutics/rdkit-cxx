// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Descriptors.MoleculeDescriptors`, a Python module:
//! a calculator of a chosen list of the descriptors in
//! `rdkit.Chem.Descriptors`, for model building.

use crate::chem::Mol;
use crate::chem::descriptors::{DescriptorFn, lookup};

/// What `calc_descriptors` gives a descriptor that fails, as Python's does.
const FAILED: f64 = -666.0;

/// What a name `rdkit.Chem.Descriptors` lacks gives, as in Python.
const UNKNOWN: f64 = 777.0;

/// A calculator of the descriptors named in its list, by their names in
/// `rdkit.Chem.Descriptors`: those of its `descList` and the
/// `AUTOCORR2D_*` descriptors. Mirrors
/// `MoleculeDescriptors.MolecularDescriptorCalculator`.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularDescriptorCalculator {
    simple_list: Vec<String>,
}

impl MolecularDescriptorCalculator {
    /// A calculator of the descriptors `simple_list` names. Mirrors
    /// `MolecularDescriptorCalculator(simpleList)`.
    pub fn new(simple_list: &[&str]) -> MolecularDescriptorCalculator {
        MolecularDescriptorCalculator {
            simple_list: simple_list.iter().map(|&name| name.to_owned()).collect(),
        }
    }

    /// Each descriptor's value for `mol`, in the list's order: −666 for one
    /// that fails, after its error is printed to standard error, as Python
    /// prints the traceback, and 777 for a name `rdkit.Chem.Descriptors`
    /// lacks. Mirrors `MolecularDescriptorCalculator.CalcDescriptors`.
    pub fn calc_descriptors(&self, mol: &Mol) -> Vec<f64> {
        self.get_descriptor_funcs()
            .into_iter()
            .map(|function| {
                function(mol).unwrap_or_else(|error| {
                    eprintln!("{error}");
                    FAILED
                })
            })
            .collect()
    }

    /// The descriptors' names. Mirrors
    /// `MolecularDescriptorCalculator.GetDescriptorNames`.
    pub fn get_descriptor_names(&self) -> Vec<String> {
        self.simple_list.clone()
    }

    /// Each descriptor's summary: the first paragraph of Python's docstring
    /// for it, on one line, or `N/A`. Mirrors
    /// `MolecularDescriptorCalculator.GetDescriptorSummaries`.
    pub fn get_descriptor_summaries(&self) -> Vec<String> {
        self.simple_list
            .iter()
            .map(|name| lookup(name).map_or("N/A", |entry| entry.2).to_owned())
            .collect()
    }

    /// Each descriptor's function; for a name `rdkit.Chem.Descriptors`
    /// lacks, one that returns 777. Mirrors
    /// `MolecularDescriptorCalculator.GetDescriptorFuncs`.
    pub fn get_descriptor_funcs(&self) -> Vec<DescriptorFn> {
        self.simple_list
            .iter()
            .map(|name| lookup(name).map_or(unknown as DescriptorFn, |entry| entry.3))
            .collect()
    }

    /// Each descriptor's version, or `N/A`. Mirrors
    /// `MolecularDescriptorCalculator.GetDescriptorVersions`.
    pub fn get_descriptor_versions(&self) -> Vec<String> {
        self.simple_list
            .iter()
            .map(|name| lookup(name).map_or("N/A", |entry| entry.1).to_owned())
            .collect()
    }

    /// Prints the descriptors' names to standard output. Mirrors
    /// `Descriptors.DescriptorCalculator.ShowDescriptors`, the base class's
    /// method.
    pub fn show_descriptors(&self) {
        print!("{}", self.descriptors_text());
    }

    /// What `show_descriptors` prints.
    fn descriptors_text(&self) -> String {
        let mut text = String::from("#---------\nSimple:\n");
        for desc in &self.simple_list {
            text.push_str(desc);
            text.push('\n');
        }
        text
    }
}

fn unknown(_: &Mol) -> crate::Result<f64> {
    Ok(UNKNOWN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn show_descriptors_prints_pythons_text() {
        let calc = MolecularDescriptorCalculator::new(&["MolLogP", "TPSA"]);
        assert_eq!(
            calc.descriptors_text(),
            "#---------\nSimple:\nMolLogP\nTPSA\n"
        );
    }
}
