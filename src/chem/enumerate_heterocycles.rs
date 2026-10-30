// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.EnumerateHeterocycles`, a Python module: the
//! aromatic rings a molecule's could become by swapping their atoms between
//! carbon, nitrogen, oxygen, and sulfur, one atom at a time.
//!
//! Python parses its reactions once, into the module's `REACTION_CACHE`;
//! here each enumeration parses its own, since a `ChemicalReaction` cannot be
//! shared between threads (CONVENTIONS.md S4).

use std::collections::{HashSet, VecDeque};

use crate::Result;
use crate::chem::rd_chem_reactions::{ChemicalReaction, reaction_from_smarts};
use crate::chem::{Mol, SmilesWriteParams, mol_to_smiles_with_params, sanitize_mol};

/// A pattern for mutating an atom of an aromatic ring, with examples for
/// testing it. Mirrors the `HeteroAtomReaction` named tuples
/// `GetHeterocycleReactionSmarts` returns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HeteroAtomReaction {
    /// The left side of the reaction SMARTS, matching the atom to mutate;
    /// empty for a kind of atom with no mutations.
    pub smarts: &'static str,
    /// The element mutated: `c`, `n`, `o`, or `s` (`o,s` for either).
    pub convert_from: &'static str,
    /// Right sides of the reaction SMARTS, separated by commas, each a
    /// reaction.
    pub convert_to: &'static str,
    /// Ring systems, separated by commas, that `smarts` matches.
    pub example: &'static str,
    /// Ring systems, separated by commas, that `smarts` must not match.
    pub negative_example: &'static str,
    /// What `smarts` matches, for people.
    pub description: &'static str,
}

/// The patterns for mutating atoms of aromatic rings. Mirrors
/// `EnumerateHeterocycles.GetHeterocycleReactionSmarts`.
pub fn get_heterocycle_reaction_smarts() -> Vec<HeteroAtomReaction> {
    vec![
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([n;x2;r5:2])[n;x2;r5:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([n:2])[n:3]"#,
            example: r#"c1nc[nH]c1,c1ncn([*])c1,c1nncn([*])1"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 2 nitrogens"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([n;h1;!+;!-;x2;r5;!$(*=*);!$(*[*;r5]=*):2])[n;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[*;r5;x3]):3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[o:1]([n;H0:2])[n;H0:3],[s:1]([n;H0:2])[n;H0:3]"#,
            example: r#"c1nc[nH]c1,[*]c1nc[nH]c1"#,
            negative_example: r#"[*]c1ncno1,[*]n1[nH]cnc1=O,[*]c1nc2[nH]cnn2c(=O)c1[*],c1nc2cscc2[nH]1,c1c[nH]c2c1[nH]cn2,c1ccc2c(c1)[nH]cn2"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 nitrogen not bonded to another atom that is in 2 rings, bonded to 1 nitrogen with H and no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([n;h1;!+;!-;x2;r5;!$(*=*);!$(*[*;r5]=*):2])[n;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[n;r5;x3]);$(*[c;x3;$(*[*;r6])]):3] "#,
            convert_from: r#"c"#,
            convert_to: r#"[o:1]([n;H0:2])[n;H0:3],[s:1]([n;H0:2])[n;H0:3]"#,
            example: r#"c1ccc2c(c1)[nH]cn2"#,
            negative_example: r#"[*]c1ncno1,[*]n1[nH]cnc1=O,[*]c1nc2[nH]cnn2c(=O)c1[*],c1nc2cscc2[nH]1,c1c[nH]c2c1[nH]cn2"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 nitrogen bonded to a bridgehead carbon in 6-member ring, bonded to 1 nitrogen with H and no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;x2;r5:2])[o,s;x2;r5:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([c:2])[*:3]"#,
            example: r#"[*]c1ccoc1,c1ccoc1[*],n1ccoc1[*]"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon and 1 oxygen or sulfur"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;r5:2])[n;r5;h0;!-;!+:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([c:2])[n:3]"#,
            example: r#"c1nc[nH]c1[*],n1cc[nH]c1[*],n1cc[nH][n+]1[*],[*]n1cc2c(c1)CCCC2,[*]N1Cc2cn([*])cc2C1,[*]c1c2c(cn1[*])CCCC2,[*]c1cn2c(c1[*])CCC2"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with no H no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;x2;r5:2])[n;x2;r5;h1;+:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([c:2])[n:3]"#,
            example: r#"c1c([*])sc[nH+]1,c1csc([*])[nH+]1"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with H and charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c:2])[n;x2;r5;h1;!+;!-:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([c:2])[n:3]"#,
            example: r#"[*]n1[nH]ccc1=O,[*]n1cc[nH]c1=O,[*]c1nc2[nH]ccn2c(=O)c1[*],[*]c1ccc[nH]1,[*]c1ncc[nH]1,[*][n+]1ccc[nH]1,[*][n+]1cc[nH]c1,c1ccc2c(c1)cc[nH]2"#,
            negative_example: r#"[*]c1cn2ccnc2nc1[*]"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with H and no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[*;r5;x3]):2])[n;x2;r5;!$(*=*);!$(*[*;r5]=*);!$(*[*;r5;x3]);h1;!+;!-:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[o:1]([c:2])[n;H0:3],[s:1]([c:2])[n;H0:3]"#,
            example: r#"[*]c1ccc[nH]1,[*]c1ncc[nH]1,[*][n+]1ccc[nH]1,[*][n+]1cc[nH]c1"#,
            negative_example: r#"c1ccc2c(c1)cc[nH]2,[*]n1[nH]ccc1=O,[*]n1cc[nH]c1=O,[*]c1nc2[nH]ccn2c(=O)c1[*],[*]c1cn2ccnc2nc1[*]"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with H and no charge, no exocyclic double bonds or fused rings"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;x2;r5:2])[n;h0;+;x2;r5:3]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([c:2])[n:3]"#,
            example: r#"CC[n+]1ccco1,[*][n+]1ccco1,[*][n+]1ccc[nH]1"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with no H and with charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;x2;r5:2])[n;h0;+;x2;r5:3][n;x2;r5;h1:4]"#,
            convert_from: r#"c"#,
            convert_to: r#"[o:1]([c:2])[n:3][n;H0:4],[s:1]([c:2])[n:3][n;H0:4]"#,
            example: r#"[*][n+]1ccc[nH]1"#,
            negative_example: r#"CC[n+]1ccco1,[*][n+]1ccco1"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 1 carbon, bonded to 1 nitrogen with no H and with charge, that nitrogen is then bonded to a nitrogen with a hydrogen"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;x2;D2;r5:1]([c;x2;r5:2])[c;x2;r5:3][o,s;x2;r5:4]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]([c:2])[c:3][*:4]"#,
            example: r#"s1cccn1,o1c([*])ccn1"#,
            negative_example: r#"s1cc([*])cn1"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 2 carbons, oxygen or sulfur in ring"#,
        },
        HeteroAtomReaction {
            smarts: r#"[!o;!s;r5:1]1~[c;r5:2]~[c;x2;D2;r5:3]~[c;r5:4]~[!o;!s;r5:5]~1"#,
            convert_from: r#"c"#,
            convert_to: r#"[*:1]1~[c:2]~[n:3]~[c:4]~[*:5]~1"#,
            example: r#"[*]c1cc[nH]c1,[*]n1cccc1,[*]n1ccc2nccc-2c1"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 2 carbons, no oxygen or sulfur in ring"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;h1;!+;!-;r5:1]1~[c;!$(*=*);r5:2]~[c;x2;D2;r5:3]~[c;!$(*=*);r5:4]~[!o;!s;r5:5]~1"#,
            convert_from: r#"c"#,
            convert_to: r#"[n;H0:1]1~[c:2]~[o:3]~[c:4]~[*:5]1,[n;H0:1]1~[c:2]~[s:3]~[c:4]~[*:5]~1"#,
            example: r#"c1cc[nH]c1,c1ccc2c(c1)cc[nH]2"#,
            negative_example: r#"[*]n1cccc1,[*]n1cc([*])cc1,[*]n1ccc(=O)[nH]1"#,
            description: r#"aromatic carbon in 5 membered ring, bonded to 2 carbons without carbonyls, no oxygen or sulfur in ring, and a nitrogen with a hydrogen in the ring"#,
        },
        HeteroAtomReaction {
            smarts: r#"[c;h1;D2;r6;$(a1ccc[c,n,o][c,n,o]1),$(a1cc[c,n,o][c,n,o]c1),$(a1c[c,n,o][c,n,o]cc1),$(a1[c,n,o][c,n,o]ccc1),$(a1cc[n,o]c[n,o]1),$(a1c[n,o]c[n,o]c1),$(a1[n,o]c[n,o]cc1),$(a1[n,o]ccc[n,o]1),$(a1c[n,o]cc[n,o]1),$(a1[n,o]cc[n,o]c1):1]"#,
            convert_from: r#"c"#,
            convert_to: r#"[n:1]"#,
            example: r#"c1ccccc1,c1cccnc1,c1nccnc1,[*]c1cc(=O)cco1,c1c([*])[nH]c(=O)[nH]c1=O,[*]c1c[nH]c(=O)[nH]c1=O"#,
            negative_example: r#""#,
            description: r#"aromatic carbon in 6 membered ring, has 1H, ring has less than 3 nitrogens"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;x2;r5;h1;!+;!-;!$(*[o,s]):1]"#,
            convert_from: r#"n"#,
            convert_to: r#"[o:1],[s:1]"#,
            example: r#"c1n[nH]nn1,c1cc[nH]c1"#,
            negative_example: r#"c1ccn([*])c1,[*]c1o[nH]c(=O)c1[*]"#,
            description: r#"aromatic nitrogen in 5 membered ring, 1H, no charge, not next to sulfur or oxygen"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;x2;X2;r5;h0;!+;!-;$(*[a;r5][a;r5]):1]"#,
            convert_from: r#"n"#,
            convert_to: r#"[c:1]"#,
            example: r#"c1cnc[nH]1"#,
            negative_example: r#"c1ccn([*])c1"#,
            description: r#"aromatic nitrogen in 5 membered ring, bonded to 2 heavy atoms, no H, no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#""#,
            convert_from: r#"n"#,
            convert_to: r#"No mutations"#,
            example: r#""#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 5 membered ring, bonded to 3 heavy atoms, no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;x2;X3;r5;+;$(*[a;r5][a;r5]):1]"#,
            convert_from: r#"n"#,
            convert_to: r#"[c:1]"#,
            example: r#"c1coc[n+]1C"#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 5 membered ring, bonded to 3 heavy atoms, charged"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;x2;X2;r6;h0;!+;!-;$(*[a;r6][a;r6]):1]"#,
            convert_from: r#"n"#,
            convert_to: r#"[c:1]"#,
            example: r#"c1cccnc1,c1nccnc1"#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 6 membered ring, bonded to 2 heavy atoms, no hydrogens, no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#""#,
            convert_from: r#"n"#,
            convert_to: r#"No mutations"#,
            example: r#""#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 6 membered ring, bonded to 2 heavy atoms, 1H, no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;x2;r6;h1;+;$(*[a;r6][a;r6]):1]"#,
            convert_from: r#"n"#,
            convert_to: r#"[c:1]"#,
            example: r#"c1ccc2c(c1)ccc[nH+]2"#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 6 membered ring, bonded to 2 heavy atoms, 1H, charged"#,
        },
        HeteroAtomReaction {
            smarts: r#"[n;x2;r6;X3;+;$(*[a;r6][a;r6]):1]"#,
            convert_from: r#"n"#,
            convert_to: r#"[c:1]"#,
            example: r#"[*][n+]1cccc2c1cccc2"#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 6 membered ring, bonded to 3 heavy atoms, charged"#,
        },
        HeteroAtomReaction {
            smarts: r#""#,
            convert_from: r#"n"#,
            convert_to: r#"No mutations"#,
            example: r#""#,
            negative_example: r#""#,
            description: r#"aromatic nitrogen in 6 membered ring, bonded to 3 heavy atoms, no charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[s;x2;r5;!+;!-;$(*([c;r5])[c;r5]):1]"#,
            convert_from: r#"s"#,
            convert_to: r#"[n;H1:1],[o:1]"#,
            example: r#"c1ccsc1,c1ncsc1"#,
            negative_example: r#""#,
            description: r#"aromatic sulfur in 5 membered ring, bonded to 2 carbons"#,
        },
        HeteroAtomReaction {
            smarts: r#"[o;x2;r5;!+;!-;$(*([c;r5])[c;r5]):1]"#,
            convert_from: r#"o"#,
            convert_to: r#"[n;H1:1],[s:1]"#,
            example: r#"c1ccoc1,c1ncoc1"#,
            negative_example: r#""#,
            description: r#"aromatic oxygen in 5 membered ring, bonded to 2 carbons"#,
        },
        HeteroAtomReaction {
            smarts: r#"[o,s;x2;r5;!+;!-:1]([c;x2;r5;!$(*=*):2])[n;x2;r5;h0;D2;!+;!-:3]"#,
            convert_from: r#"o,s"#,
            convert_to: r#"[c:1]([c:2])[n;H1:3]"#,
            example: r#"c1cnoc1,c1csnc1,c1c([*])snc1"#,
            negative_example: r#"[*]n1occc1=O,[*][n+]1ccco1,[*]n1c(=O)[nH]sc1=O,[*]c1nsc(=O)o1"#,
            description: r#"aromatic oxygen or sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen with no charge with no hydrogen and only bonded to 2 atoms"#,
        },
        HeteroAtomReaction {
            smarts: r#"[o,s;x2;r5;!+;!-:1]([c;x2;r5:2])[n;x2;r5;+:3]"#,
            convert_from: r#"o,s"#,
            convert_to: r#"[c:1]([c:2])[n;+0:3]"#,
            example: r#"[*][n+]1ccco1"#,
            negative_example: r#""#,
            description: r#"aromatic oxygen or sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen in ring with a charge"#,
        },
        HeteroAtomReaction {
            smarts: r#"[s;x2;r5;!+;!-;$(*([c;r5])[n;h1;r5]):1]"#,
            convert_from: r#"s"#,
            convert_to: r#"[n;H1:1],[o:1]"#,
            example: r#"[*]c1s[nH]c(=O)c1[*]"#,
            negative_example: r#""#,
            description: r#"aromatic sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen with a hydrogen"#,
        },
        HeteroAtomReaction {
            smarts: r#"[o;x2;r5;!+;!-;$(*([c;r5])[n;h1;r5]):1]"#,
            convert_from: r#"o"#,
            convert_to: r#"[n;H1:1],[s:1]"#,
            example: r#"[*]c1o[nH]c(=O)c1[*]"#,
            negative_example: r#""#,
            description: r#"aromatic oxygen in 5 membered ring, bonded to 1 carbon, 1 nitrogen with a hydrogen"#,
        },
        HeteroAtomReaction {
            smarts: r#"[s;x2;r5;!+;!-:1]([c;x2;r5:2])[n;x2;r5;h0;D3;!+;!-:3]"#,
            convert_from: r#"s"#,
            convert_to: r#"[n;H1:1]([c:2])[n:3],[o:1]([c:2])[n:3]"#,
            example: r#"[*]n1sc(=O)n([*])c1=O"#,
            negative_example: r#""#,
            description: r#"aromatic sulfur in 5 membered ring, bonded to 1 carbon, 1 nitrogen with no charge bonded to 3 atoms "#,
        },
        HeteroAtomReaction {
            smarts: r#"[o;x2;r5;!+;!-:1]([c;x2;r5:2])[n;x2;r5;h0;D3;!+;!-:3]"#,
            convert_from: r#"o"#,
            convert_to: r#"[n;H1:1]([c:2])[n:3],[s:1]([c:2])[n:3]"#,
            example: r#"[*]n1oc(=O)n([*])c1=O"#,
            negative_example: r#""#,
            description: r#"aromatic oxygen in 5 membered ring, bonded to 1 carbon, 1 nitrogen with no charge bonded to 3 atoms "#,
        },
        HeteroAtomReaction {
            smarts: r#"[o;r6:1]"#,
            convert_from: r#"o"#,
            convert_to: r#"[n;H1:1]"#,
            example: r#"[*]c1cccc(=O)o1"#,
            negative_example: r#""#,
            description: r#"aromatic oxygen in 6 membered ring"#,
        },
    ]
}

/// The reactions of `get_heterocycle_reaction_smarts`: one for each right
/// side of each pattern that has a left side. Mirrors
/// `EnumerateHeterocycles.GetHeterocycleReactions`.
pub fn get_heterocycle_reactions() -> Vec<ChemicalReaction> {
    let mut res = Vec::new();
    for row in get_heterocycle_reaction_smarts() {
        if row.smarts.is_empty() {
            continue;
        }
        for product in row.convert_to.split(',') {
            res.push(
                reaction_from_smarts(&format!("{}>>{product}", row.smarts))
                    .expect("the heterocycle reaction SMARTS parse"),
            );
        }
    }
    res
}

/// The heterocycles of a molecule, made as they are asked for. Built by
/// [`enumerate_heterocycles`].
/// A Rust addition: the iterator of Python's generator.
pub struct Heterocycles {
    reactions: Vec<ChemicalReaction>,
    depth: Option<u32>,
    starting_points: VecDeque<(Mol, u32)>,
    /// The molecule being mutated, its depth, and the next reaction.
    current: Option<(Mol, u32, usize)>,
    pending: VecDeque<Mol>,
    seen: HashSet<String>,
    done: bool,
}

impl std::fmt::Debug for Heterocycles {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Heterocycles").finish_non_exhaustive()
    }
}

impl Heterocycles {
    fn next_mol(&mut self) -> Result<Option<Mol>> {
        loop {
            if let Some(mut newmol) = self.pending.pop_front() {
                let (_, curdepth, _) = self
                    .current
                    .as_ref()
                    .expect("products come from a molecule");
                let curdepth = *curdepth;
                sanitize_mol(&mut newmol)?;
                let smiles = mol_to_smiles_with_params(
                    &newmol,
                    &SmilesWriteParams {
                        do_isomeric_smiles: true,
                        ..Default::default()
                    },
                )?;
                if !self.seen.insert(smiles) {
                    continue;
                }
                self.starting_points
                    .push_back((newmol.clone(), curdepth + 1));
                return Ok(Some(newmol));
            }
            match &mut self.current {
                Some((curmol, _, rxn_idx)) if *rxn_idx < self.reactions.len() => {
                    let rxn = &mut self.reactions[*rxn_idx];
                    *rxn_idx += 1;
                    for products in rxn.run_reactants(std::slice::from_ref(curmol))? {
                        if let Some(first) = products.into_iter().next() {
                            self.pending.push_back(first);
                        }
                    }
                }
                _ => {
                    let Some((curmol, curdepth)) = self.starting_points.pop_front() else {
                        return Ok(None);
                    };
                    if self.depth.is_some_and(|depth| curdepth >= depth) {
                        return Ok(None);
                    }
                    self.current = Some((curmol, curdepth, 0));
                }
            }
        }
    }
}

/// The next heterocycle.
impl Iterator for Heterocycles {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match self.next_mol() {
            Ok(Some(mol)) => Some(Ok(mol)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Every aromatic ring system the molecule's could become, mutating one
/// atom at a time and again on each result, each once, by its isomeric
/// SMILES; the molecule itself comes back when it is one of them. Mirrors
/// `EnumerateHeterocycles.EnumerateHeterocycles(inputmol)`.
pub fn enumerate_heterocycles(inputmol: &Mol) -> Heterocycles {
    enumerate_heterocycles_with_depth(inputmol, None)
}

/// As [`enumerate_heterocycles`], mutating at most `depth` atoms; `None`
/// for no limit. Mirrors
/// `EnumerateHeterocycles.EnumerateHeterocycles(inputmol, depth)`.
pub fn enumerate_heterocycles_with_depth(inputmol: &Mol, depth: Option<u32>) -> Heterocycles {
    Heterocycles {
        reactions: get_heterocycle_reactions(),
        depth,
        starting_points: VecDeque::from([(inputmol.clone(), 0)]),
        current: None,
        pending: VecDeque::new(),
        seen: HashSet::new(),
        done: false,
    }
}
