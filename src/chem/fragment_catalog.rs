// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FragmentCatalog`, a Python module over
//! `rdfragcatalog`, which it re-exports: information gains read back from
//! files, the links between a catalog's fragments, and the molecules that
//! have a fragment.

use std::collections::BTreeSet;
use std::path::Path;

use crate::chem::Mol;
use crate::data_structs::ExplicitBitVect;
use crate::py::{parse_float, parse_int};
use crate::{Error, Result};

pub use crate::chem::rdfragcatalog::*;

/// A fragment's bit and its information gain, as a gains file lists them.
/// Mirrors `FragmentCatalog.BitGainsInfo`; `default` has Python's class
/// attributes.
#[derive(Clone, Debug, PartialEq)]
pub struct BitGainsInfo {
    /// The fragment's bit; -1 by default.
    pub id: i64,
    /// The fragment's description.
    pub description: String,
    /// The bit's information gain.
    pub gain: f64,
    /// The number of molecules of each class with the bit; `None` by
    /// default.
    pub n_per_class: Option<Vec<i64>>,
}

impl Default for BitGainsInfo {
    fn default() -> BitGainsInfo {
        BitGainsInfo {
            id: -1,
            description: String::new(),
            gain: 0.0,
            n_per_class: None,
        }
    }
}

/// Parameters for `process_gains_file_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ProcessGainsFileParams {
    /// The number of bits after which to stop; `None` for all. As in
    /// Python, `Some(0)` reads them all, since it is checked after each bit.
    pub n_to_do: Option<usize>,
    /// What separates a line's fields.
    pub delim: String,
    /// The second field is the description.
    pub have_descriptions: bool,
}

impl Default for ProcessGainsFileParams {
    fn default() -> ProcessGainsFileParams {
        ProcessGainsFileParams {
            n_to_do: None,
            delim: ",".to_owned(),
            have_descriptions: true,
        }
    }
}

/// Lines as Python reads a text file: `\r\n` and `\r` read as `\n`, each
/// line keeping its end.
fn universal_lines(text: &str) -> Vec<String> {
    let text = text.replace("\r\n", "\n").replace('\r', "\n");
    text.split_inclusive('\n').map(str::to_owned).collect()
}

/// The bits of a gains file, as `BuildFragmentCatalog.OutputGainsData`
/// writes them: a header line, then per bit its id, description, gain, and
/// counts per class. Mirrors `FragmentCatalog.ProcessGainsFile(fileName)`.
pub fn process_gains_file(file_name: impl AsRef<Path>) -> Result<Vec<BitGainsInfo>> {
    process_gains_file_with_params(file_name, &ProcessGainsFileParams::default())
}

/// `process_gains_file` with parameters. Mirrors
/// `FragmentCatalog.ProcessGainsFile(fileName, nToDo, delim,
/// haveDescriptions)`, which fails on Python 3, calling the Python 2
/// `xreadlines`; this reads the file as that meant to. `Error::Os` when the
/// file can't be read, and `Error::Value` for a field that isn't a number
/// or a line with too few fields.
pub fn process_gains_file_with_params(
    file_name: impl AsRef<Path>,
    params: &ProcessGainsFileParams,
) -> Result<Vec<BitGainsInfo>> {
    let text = std::fs::read_to_string(file_name).map_err(|e| Error::Os {
        message: e.to_string(),
    })?;
    let missing = || Error::Value {
        message: "list index out of range".to_owned(),
    };
    let mut res = Vec::new();
    for (i, line) in universal_lines(&text).iter().enumerate() {
        if i == 0 {
            continue;
        }
        let split: Vec<&str> = line.split(params.delim.as_str()).map(str::trim).collect();
        let mut bit = BitGainsInfo {
            id: parse_int(split[0])?,
            ..BitGainsInfo::default()
        };
        let mut col = 1;
        if params.have_descriptions {
            bit.description = split.get(col).ok_or_else(missing)?.to_string();
            col += 1;
        }
        bit.gain = parse_float(split.get(col).ok_or_else(missing)?)?;
        col += 1;
        bit.n_per_class = Some(
            split[col.min(split.len())..]
                .iter()
                .map(|s| parse_int(s))
                .collect::<Result<_>>()?,
        );
        res.push(bit);
        if Some(res.len()) == params.n_to_do {
            break;
        }
    }
    Ok(res)
}

/// The links between bits of a catalog: from `build_adjacency_list`.
/// A Rust addition, naming Python's adjacency list.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AdjacencyList {
    /// Per bit, in the order first given, the bits of the fragments one
    /// bond larger that contain its fragment. Python's dict `adjs`.
    pub adjs: Vec<(u32, Vec<u32>)>,
    /// Per fragment size in bonds, in the order first seen, its bits.
    /// Python's dict `levels`.
    pub levels: Vec<(u32, BTreeSet<u32>)>,
}

impl AdjacencyList {
    fn add_level(&mut self, order: u32, bit: u32) {
        match self.levels.iter_mut().find(|(o, _)| *o == order) {
            Some((_, bits)) => {
                bits.insert(bit);
            }
            None => self.levels.push((order, BTreeSet::from([bit]))),
        }
    }
}

/// The links between bits of a catalog, limited to the bits given. Mirrors
/// `FragmentCatalog.BuildAdjacencyList(catalog, bits)`.
pub fn build_adjacency_list(catalog: &FragCatalog, bits: &[BitGainsInfo]) -> Result<AdjacencyList> {
    build_adjacency_list_with_limit_inclusion(catalog, bits, true)
}

/// `build_adjacency_list`, with `limit_inclusion` false linking all bits.
/// Mirrors `FragmentCatalog.BuildAdjacencyList(catalog, bits,
/// limitInclusion)`. Fragments without bits are left out, where Python
/// links 2^32 - 1. `Error::Value` for a negative id, and the catalog's
/// errors for an id past its bits.
pub fn build_adjacency_list_with_limit_inclusion(
    catalog: &FragCatalog,
    bits: &[BitGainsInfo],
    limit_inclusion: bool,
) -> Result<AdjacencyList> {
    let bit_ids: Vec<u32> = bits
        .iter()
        .map(|b| {
            u32::try_from(b.id).map_err(|_| Error::Value {
                message: format!("invalid bit id {}", b.id),
            })
        })
        .collect::<Result<_>>()?;
    let no_entry = |bit: u32| Error::Value {
        message: format!("no fragment has bit {bit}"),
    };
    let mut res = AdjacencyList::default();
    for &bit_id in &bit_ids {
        let entry = catalog
            .get_bit_entry_id(bit_id)?
            .ok_or_else(|| no_entry(bit_id))?;
        let mut tmp = Vec::new();
        res.add_level(catalog.get_entry_order(entry)?, bit_id);
        for down in catalog.get_entry_down_ids(entry)? {
            let Some(id) = catalog.get_entry_bit_id(down)? else {
                continue;
            };
            if !limit_inclusion || bit_ids.contains(&id) {
                tmp.push(id);
                res.add_level(catalog.get_entry_order(down)?, id);
            }
        }
        match res.adjs.iter_mut().find(|(b, _)| *b == bit_id) {
            Some((_, adj)) => *adj = tmp,
            None => res.adjs.push((bit_id, tmp)),
        }
    }
    Ok(res)
}

/// The molecules whose fingerprints, paired by position, set a bit.
/// Mirrors `FragmentCatalog.GetMolsMatchingBit(mols, bit, fps)`; for a
/// `BitGainsInfo`, pass its id. `Error::Index` for a bit past a
/// fingerprint or fewer fingerprints than molecules.
pub fn get_mols_matching_bit<'m>(
    mols: &'m [Mol],
    bit: u32,
    fps: &[ExplicitBitVect],
) -> Result<Vec<&'m Mol>> {
    let mut res = Vec::new();
    for (i, mol) in mols.iter().enumerate() {
        let fp = fps.get(i).ok_or(Error::index(i, fps.len()))?;
        if fp.get_bit(bit)? {
            res.push(mol);
        }
    }
    Ok(res)
}
