// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.MolKey.InchiInfo`, a Python module: the layers of
//! an InChI, and their sp3 stereo and mobile hydrogens.
//!
//! Python parses InChIs with regular expressions, which its module
//! attributes hold compiled; the crate has none, so those attributes are
//! the patterns' text and the parsing matches what each pattern matches.
//! Python's `logging` handler (`console`) is not ported.

use std::collections::BTreeMap;

use crate::chem::{
    FindMolChiralCentersParams, find_mol_chiral_centers_with_params, mol_from_inchi,
};
use crate::{Error, Result};

/// Splits off the version. Mirrors `InchiInfo.version_re`.
pub fn version_re() -> &'static str {
    "(.*?)/(.*)"
}

/// Splits off a reconnected layer. Mirrors `InchiInfo.reconnected_re`.
pub fn reconnected_re() -> &'static str {
    "(.*?)/r(.*)"
}

/// Splits off a fixed-hydrogen layer. Mirrors `InchiInfo.fixed_h_re`.
pub fn fixed_h_re() -> &'static str {
    "(.*?)/f(.*)"
}

/// Splits off an isotope layer. Mirrors `InchiInfo.isotope_re`.
pub fn isotope_re() -> &'static str {
    "(.*?)/i(.*)"
}

/// Finds a stereo layer followed by another layer. Mirrors
/// `InchiInfo.stereo_re`.
pub fn stereo_re() -> &'static str {
    r".*/t(.*?)/.*"
}

/// Finds a stereo layer. Mirrors `InchiInfo.stereo_all_re`.
pub fn stereo_all_re() -> &'static str {
    r".*/t([^/]+)"
}

/// Finds an undefined stereocenter. Mirrors `InchiInfo.undef_stereo_re`.
pub fn undef_stereo_re() -> &'static str {
    r"(\d+)\?"
}

/// Finds any stereocenter. Mirrors `InchiInfo.all_stereo_re`.
pub fn all_stereo_re() -> &'static str {
    r"(\d+)[?+-]"
}

/// Finds a defined stereocenter. Mirrors `InchiInfo.defined_stereo_re`.
pub fn defined_stereo_re() -> &'static str {
    r"(\d+)[+-]"
}

/// Finds the hydrogen layer. Mirrors `InchiInfo.h_layer_re`.
pub fn h_layer_re() -> &'static str {
    r".*/h(.*)/?"
}

/// Finds a mobile hydrogen group. Mirrors `InchiInfo.mobile_h_group_re`.
pub fn mobile_h_group_re() -> &'static str {
    r"(\(H.+?\))"
}

/// Finds a mobile hydrogen group's atoms; the module does not use it.
/// Mirrors `InchiInfo.mobile_h_atoms_re`.
pub fn mobile_h_atoms_re() -> &'static str {
    r",(\d+)"
}

/// The text up to the first newline, which `.` in Python's patterns does
/// not cross.
fn line(s: &str) -> &str {
    s.split('\n').next().unwrap_or_default()
}

/// `re.match('(.*?)' + sep + '(.*)', s)`: the text around the first `sep`.
fn split_first<'a>(s: &'a str, sep: &str) -> Option<(&'a str, &'a str)> {
    let i = line(s).find(sep)?;
    Some((&s[..i], line(&s[i + sep.len()..])))
}

/// `re.match(r'.*/t(.*?)/.*', s)`: the text between the last `/t` that
/// another `/` follows and that `/`.
fn stereo_match(s: &str) -> Option<&str> {
    let s = line(s);
    let mut found = None;
    for (i, _) in s.match_indices("/t") {
        if let Some(j) = s[i + 2..].find('/') {
            found = Some(&s[i + 2..i + 2 + j]);
        }
    }
    found
}

/// `re.match(r'.*/t([^/]+)', s)`: the text after the last `/t` that is
/// followed by something other than `/`, up to the next `/`.
fn stereo_all_match(s: &str) -> Option<&str> {
    let s = line(s);
    let mut found = None;
    for (i, _) in s.match_indices("/t") {
        let rest = &s[i + 2..];
        let end = rest.find('/').unwrap_or(rest.len());
        if end > 0 {
            found = Some(&rest[..end]);
        }
    }
    found
}

/// The number of `(\d+)[chars]` matches: characters of `chars` after a
/// digit.
fn count_after_digit(s: &str, chars: &[char]) -> usize {
    let mut prev_digit = false;
    let mut n = 0;
    for c in s.chars() {
        if prev_digit && chars.contains(&c) {
            n += 1;
        }
        prev_digit = c.is_ascii_digit();
    }
    n
}

/// `re.findall(r'(\(H.+?\))', s)`.
fn mobile_h_groups(s: &str) -> Vec<&str> {
    let mut groups = Vec::new();
    let mut pos = 0;
    while let Some(i) = s[pos..].find("(H") {
        let start = pos + i;
        // At least one character, not a newline, then the first ")".
        let body = &s[start + 2..];
        let mut end = None;
        for (j, c) in body.char_indices() {
            if c == '\n' {
                break;
            }
            if c == ')' && j > 0 {
                end = Some(start + 2 + j + 1);
                break;
            }
        }
        match end {
            Some(end) => {
                groups.push(&s[start..end]);
                pos = end;
            }
            None => pos = start + 1,
        }
    }
    groups
}

/// Whether the InChI's molecule has no chiral centers, assigned or not; a
/// layer that is not an InChI is tried as one after `InChI=1/`. False when
/// neither reads. Mirrors `InchiInfo._is_achiral_by_symmetry`.
fn is_achiral_by_symmetry(inchi: &str) -> Result<bool> {
    let mol = match mol_from_inchi(inchi)? {
        Some(mol) => Some(mol),
        None => mol_from_inchi(&format!("InChI=1/{inchi}"))?,
    };
    let Some(mut mol) = mol else {
        return Ok(false);
    };
    let params = FindMolChiralCentersParams {
        force: true,
        include_unassigned: true,
        ..Default::default()
    };
    Ok(find_mol_chiral_centers_with_params(&mut mol, &params).is_ok_and(|c| c.is_empty()))
}

/// A layer's sp3 stereo: Python's 4-tuple of the numbers of stereocenters
/// and of undefined ones, whether the molecule is a meso form, and the
/// stereo layer's text.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct Sp3Stereo {
    /// The number of stereocenters.
    pub num_stereo: u32,
    /// The number of undefined stereocenters, at most `num_stereo`.
    pub num_undef_stereo: u32,
    /// Whether the molecule is a meso form (with chiral centers and a plane
    /// of symmetry).
    pub is_meso: bool,
    /// The comma-separated atom numbers with sp3 stereochemistry, as the
    /// stereo layer gives them.
    pub stereo: String,
}

/// A layer's mobile hydrogens: the number of groups, and the groups,
/// separated by commas.
pub type MobileH = (u32, String);

/// An InChI's layers: by connection layer (`id`, or `id_disconnected` and
/// `id_reconnected`), by fixed-hydrogen layer (`main`, and `fixed_h`), the
/// non-isotopic and isotopic layers' text. Mirrors
/// `InchiInfo.InchiInfo`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InchiInfo {
    parsed_inchi: BTreeMap<String, BTreeMap<String, BTreeMap<String, String>>>,
}

impl InchiInfo {
    /// The layers of an InChI. A string without `/` is `Error::Value`, where
    /// Python fails on its regular expression's `None`. Mirrors
    /// `InchiInfo(inchi_str)`.
    pub fn new(inchi_str: &str) -> Result<InchiInfo> {
        let (_, rest) = split_first(inchi_str, "/").ok_or_else(|| Error::Value {
            message: format!("not an InChI: {inchi_str}"),
        })?;
        let mut connection_layers = BTreeMap::new();
        match split_first(rest, "/r") {
            Some((disconnected, reconnected)) => {
                connection_layers.insert("id_disconnected", disconnected);
                connection_layers.insert("id_reconnected", reconnected);
            }
            None => {
                connection_layers.insert("id", rest);
            }
        }
        let mut parsed_inchi = BTreeMap::new();
        for (conn_layer, text) in connection_layers {
            let mut fixed_h_layers = BTreeMap::new();
            match split_first(text, "/f") {
                Some((main, fixed_h)) => {
                    fixed_h_layers.insert("main", main);
                    fixed_h_layers.insert("fixed_h", fixed_h);
                }
                None => {
                    fixed_h_layers.insert("main", text);
                }
            }
            let mut inchi = BTreeMap::new();
            for (fixed_layer, text) in fixed_h_layers {
                let mut layers = BTreeMap::new();
                match split_first(text, "/i") {
                    Some((non_isotopic, isotopic)) => {
                        layers.insert("non-isotopic".to_owned(), non_isotopic.to_owned());
                        layers.insert("isotopic".to_owned(), isotopic.to_owned());
                    }
                    None => {
                        layers.insert("non-isotopic".to_owned(), text.to_owned());
                    }
                }
                inchi.insert(fixed_layer.to_owned(), layers);
            }
            parsed_inchi.insert(conn_layer.to_owned(), inchi);
        }
        Ok(InchiInfo { parsed_inchi })
    }

    /// The layers' text, by connection, fixed-hydrogen, and isotope layer.
    /// Mirrors the `parsed_inchi` attribute.
    pub fn parsed_inchi(&self) -> &BTreeMap<String, BTreeMap<String, BTreeMap<String, String>>> {
        &self.parsed_inchi
    }

    /// Each layer's sp3 stereo, by fixed-hydrogen and isotope layer; with a
    /// reconnected layer, the reconnected one's, as in Python, which writes
    /// over the disconnected one's. A molecule is meso when its stereo
    /// layer, with no following layer, defines more than one center, or when
    /// it has more than one undefined center and no chiral centers by
    /// symmetry. Mirrors `InchiInfo.get_sp3_stereo`.
    pub fn get_sp3_stereo(&self) -> Result<BTreeMap<String, BTreeMap<String, Sp3Stereo>>> {
        let mut sp3_stereo = BTreeMap::new();
        for fixed in self.parsed_inchi.values() {
            for (fixed_layer, isos) in fixed {
                let mut by_iso = BTreeMap::new();
                for (iso_layer, inchi_layer) in isos {
                    let mut is_meso = false;
                    let stereo = if let Some(s) = stereo_match(inchi_layer) {
                        s
                    } else if let Some(s) = stereo_all_match(inchi_layer) {
                        is_meso = count_after_digit(s, &['+', '-']) > 1;
                        s
                    } else {
                        ""
                    };
                    let num_stereo = count_after_digit(stereo, &['?', '+', '-']) as u32;
                    let num_undef_stereo = count_after_digit(stereo, &['?']) as u32;
                    is_meso =
                        is_meso || (num_undef_stereo > 1 && is_achiral_by_symmetry(inchi_layer)?);
                    by_iso.insert(
                        iso_layer.clone(),
                        Sp3Stereo {
                            num_stereo,
                            num_undef_stereo,
                            is_meso,
                            stereo: stereo.to_owned(),
                        },
                    );
                }
                sp3_stereo.insert(fixed_layer.clone(), by_iso);
            }
        }
        Ok(sp3_stereo)
    }

    /// Each layer's mobile hydrogen groups, by fixed-hydrogen and isotope
    /// layer, as for `get_sp3_stereo`. Mirrors `InchiInfo.get_mobile_h`.
    pub fn get_mobile_h(&self) -> BTreeMap<String, BTreeMap<String, MobileH>> {
        let mut mobile_h = BTreeMap::new();
        for fixed in self.parsed_inchi.values() {
            for (fixed_layer, isos) in fixed {
                let mut by_iso = BTreeMap::new();
                for (iso_layer, inchi_layer) in isos {
                    // `.*/h(.*)/?`: everything after the last "/h".
                    let line = line(inchi_layer);
                    let (num_groups, groups) = match line.rfind("/h") {
                        Some(i) => {
                            let groups = mobile_h_groups(&line[i + 2..]);
                            (groups.len() as u32, groups.join(","))
                        }
                        None => (0, String::new()),
                    };
                    by_iso.insert(iso_layer.clone(), (num_groups, groups));
                }
                mobile_h.insert(fixed_layer.clone(), by_iso);
            }
        }
        mobile_h
    }
}
