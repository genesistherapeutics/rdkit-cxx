// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.EState.AtomTypes`, a pure-Python module: SMARTS
//! definitions of Hall and Kier's EState atom types (Hall and Kier, JCICS
//! 35, 1039-1045 (1995), table 1).
//!
//! Python builds the patterns once, into the module's `esPatterns`; here
//! each call builds them anew, since a `Mol` cannot be shared between
//! threads (CONVENTIONS.md S1).

use crate::Result;
use crate::chem::{Mol, SubstructMatchParameters, mol_from_smarts};

/// The atom types' names and SMARTS, in Python's order (`_rawD`).
pub(crate) const RAW_D: [(&str, &str); 79] = [
    ("sLi", "[LiD1]-*"),
    ("ssBe", "[BeD2](-*)-*"),
    ("ssssBe", "[BeD4](-*)(-*)(-*)-*"),
    ("ssBH", "[BD2H](-*)-*"),
    ("sssB", "[BD3](-*)(-*)-*"),
    ("ssssB", "[BD4](-*)(-*)(-*)-*"),
    ("sCH3", "[CD1H3]-*"),
    ("dCH2", "[CD1H2]=*"),
    ("ssCH2", "[CD2H2](-*)-*"),
    ("tCH", "[CD1H]#*"),
    ("dsCH", "[CD2H](=*)-*"),
    ("aaCH", "[C,c;D2H](:*):*"),
    ("sssCH", "[CD3H](-*)(-*)-*"),
    ("ddC", "[CD2H0](=*)=*"),
    ("tsC", "[CD2H0](#*)-*"),
    ("dssC", "[CD3H0](=*)(-*)-*"),
    ("aasC", "[C,c;D3H0](:*)(:*)-*"),
    ("aaaC", "[C,c;D3H0](:*)(:*):*"),
    ("ssssC", "[CD4H0](-*)(-*)(-*)-*"),
    ("sNH3", "[ND1H3]-*"),
    ("sNH2", "[ND1H2]-*"),
    ("ssNH2", "[ND2H2](-*)-*"),
    ("dNH", "[ND1H]=*"),
    ("ssNH", "[ND2H](-*)-*"),
    ("aaNH", "[N,nD2H](:*):*"),
    ("tN", "[ND1H0]#*"),
    ("sssNH", "[ND3H](-*)(-*)-*"),
    ("dsN", "[ND2H0](=*)-*"),
    ("aaN", "[N,nD2H0](:*):*"),
    ("sssN", "[ND3H0](-*)(-*)-*"),
    ("ddsN", "[ND3H0](~[OD1H0])(~[OD1H0])-,:*"),
    ("aasN", "[N,nD3H0](:*)(:*)-,:*"),
    ("ssssN", "[ND4H0](-*)(-*)(-*)-*"),
    ("sOH", "[OD1H]-*"),
    ("dO", "[OD1H0]=*"),
    ("ssO", "[OD2H0](-*)-*"),
    ("aaO", "[O,oD2H0](:*):*"),
    ("sF", "[FD1]-*"),
    ("sSiH3", "[SiD1H3]-*"),
    ("ssSiH2", "[SiD2H2](-*)-*"),
    ("sssSiH", "[SiD3H1](-*)(-*)-*"),
    ("ssssSi", "[SiD4H0](-*)(-*)(-*)-*"),
    ("sPH2", "[PD1H2]-*"),
    ("ssPH", "[PD2H1](-*)-*"),
    ("sssP", "[PD3H0](-*)(-*)-*"),
    ("dsssP", "[PD4H0](=*)(-*)(-*)-*"),
    ("sssssP", "[PD5H0](-*)(-*)(-*)(-*)-*"),
    ("sSH", "[SD1H1]-*"),
    ("dS", "[SD1H0]=*"),
    ("ssS", "[SD2H0](-*)-*"),
    ("aaS", "[S,sD2H0](:*):*"),
    ("dssS", "[SD3H0](=*)(-*)-*"),
    ("ddssS", "[SD4H0](~[OD1H0])(~[OD1H0])(-*)-*"),
    ("sCl", "[ClD1]-*"),
    ("sGeH3", "[GeD1H3](-*)"),
    ("ssGeH2", "[GeD2H2](-*)-*"),
    ("sssGeH", "[GeD3H1](-*)(-*)-*"),
    ("ssssGe", "[GeD4H0](-*)(-*)(-*)-*"),
    ("sAsH2", "[AsD1H2]-*"),
    ("ssAsH", "[AsD2H1](-*)-*"),
    ("sssAs", "[AsD3H0](-*)(-*)-*"),
    ("sssdAs", "[AsD4H0](=*)(-*)(-*)-*"),
    ("sssssAs", "[AsD5H0](-*)(-*)(-*)(-*)-*"),
    ("sSeH", "[SeD1H1]-*"),
    ("dSe", "[SeD1H0]=*"),
    ("ssSe", "[SeD2H0](-*)-*"),
    ("aaSe", "[SeD2H0](:*):*"),
    ("dssSe", "[SeD3H0](=*)(-*)-*"),
    ("ddssSe", "[SeD4H0](=*)(=*)(-*)-*"),
    ("sBr", "[BrD1]-*"),
    ("sSnH3", "[SnD1H3]-*"),
    ("ssSnH2", "[SnD2H2](-*)-*"),
    ("sssSnH", "[SnD3H1](-*)(-*)-*"),
    ("ssssSn", "[SnD4H0](-*)(-*)(-*)-*"),
    ("sI", "[ID1]-*"),
    ("sPbH3", "[PbD1H3]-*"),
    ("ssPbH2", "[PbD2H2](-*)-*"),
    ("sssPbH", "[PbD3H1](-*)(-*)-*"),
    ("ssssPb", "[PbD4H0](-*)(-*)(-*)-*"),
];

/// The atom types' names and patterns, in Python's order. Mirrors
/// `AtomTypes.esPatterns`, after `BuildPatts()`.
pub fn es_patterns() -> Vec<(&'static str, Mol)> {
    RAW_D
        .iter()
        .map(|&(name, sma)| (name, mol_from_smarts(sma).expect("the EState SMARTS parse")))
        .collect()
}

/// The atom types' names and patterns. Python stores them in
/// `esPatterns`; this returns them, as [`es_patterns`] does. Mirrors
/// `AtomTypes.BuildPatts`.
pub fn build_patts() -> Vec<(&'static str, Mol)> {
    es_patterns()
}

/// From other `(name, SMARTS)` pairs. A SMARTS that does not parse leaves
/// `None` in its place, after a warning to standard error, as Python's
/// does. Mirrors `AtomTypes.BuildPatts(rawV)`.
pub fn build_patts_with_raw_v(raw_v: &[(&str, &str)]) -> Vec<Option<(String, Mol)>> {
    raw_v
        .iter()
        .map(|&(name, sma)| match mol_from_smarts(sma) {
            Ok(patt) => Some((name.to_owned(), patt)),
            Err(_) => {
                eprintln!("WARNING: problems with pattern {sma} (name: {name}), skipped.");
                None
            }
        })
        .collect()
}

/// Each atom's EState types: the names of the patterns that match with it
/// first, in the patterns' order. Mirrors `AtomTypes.TypeAtoms`.
pub fn type_atoms(mol: &Mol) -> Result<Vec<Vec<String>>> {
    let params = SubstructMatchParameters {
        uniquify: false,
        ..SubstructMatchParameters::default()
    };
    let mut res: Vec<Vec<String>> = vec![Vec::new(); mol.get_num_atoms() as usize];
    for (name, patt) in es_patterns() {
        for m in mol.get_substruct_matches_with_params(&patt, &params)? {
            let types = &mut res[m[0] as usize];
            if !types.iter().any(|t| t == name) {
                types.push(name.to_owned());
            }
        }
    }
    Ok(res)
}
