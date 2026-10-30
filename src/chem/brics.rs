// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.BRICS`, a Python module: the BRICS algorithm of
//! Degen et al., _ChemMedChem_ 3, 1503-1507 (2008), which breaks molecules
//! at retrosynthetically interesting bonds and builds new ones from the
//! pieces.
//!
//! Python builds its patterns and reactions once, on import; each function
//! here that returns them makes new ones, since neither a `Mol` nor a
//! `ChemicalReaction` can be shared between threads (CONVENTIONS.md S1,
//! S4), and the functions that use them build them once per call. Where
//! Python draws from the `random` module's global generator, the ports take
//! a [`Random`], which draws what Python's would from the same seed.
//! `BRICSBuild` is an iterator, made as it is asked for, as Python's
//! generator is. Python's `silent` arguments, which only print progress,
//! are not ported.

use std::collections::{BTreeSet, HashSet, VecDeque};

use crate::chem::rd_chem_reactions::{ChemicalReaction, reaction_from_smarts};
use crate::chem::{
    Atom, EditableMol, Mol, fragment_on_brics_bonds, mol_from_smarts, mol_from_smiles,
    mol_to_smiles, sanitize_mol,
};
use crate::rd_random::{Random, shuffle_with_random};
use crate::{Error, Result};

/// The BRICS environments, as SMARTS by label, in Python's order. `#L8` is
/// a definition Python keeps but no reaction uses. Mirrors `BRICS.environs`.
pub fn environs() -> Vec<(&'static str, &'static str)> {
    vec![
        ("L1", "[C;D3]([#0,#6,#7,#8])(=O)"),
        ("L3", "[O;D2]-;!@[#0,#6,#1]"),
        ("L4", "[C;!D1;!$(C=*)]-;!@[#6]"),
        (
            "L5",
            "[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]",
        ),
        ("L6", "[C;D3;!R](=O)-;!@[#0,#6,#7,#8]"),
        ("L7a", "[C;D2,D3]-[#6]"),
        ("L7b", "[C;D2,D3]-[#6]"),
        ("#L8", "[C;!R;!D1]-;!@[#6]"),
        ("L8", "[C;!R;!D1;!$(C!-*)]"),
        ("L9", "[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]"),
        ("L10", "[N;R;$(N(@C(=O))@[C,N,O,S])]"),
        ("L11", "[S;D2](-;!@[#0,#6])"),
        ("L12", "[S;D4]([#6,#0])(=O)(=O)"),
        ("L13", "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]"),
        ("L14", "[c;$(c(:[c,n,o,s]):[n,o,s])]"),
        ("L14b", "[c;$(c(:[c,n,o,s]):[n,o,s])]"),
        ("L15", "[C;$(C(-;@C)-;@C)]"),
        ("L16", "[c;$(c(:c):c)]"),
        ("L16b", "[c;$(c(:c):c)]"),
    ]
}

fn environ(label: &str) -> &'static str {
    environs()
        .into_iter()
        .find(|(l, _)| *l == label)
        .map(|(_, sma)| sma)
        .expect("every BRICS reaction names an environment")
}

/// The pairs of environments BRICS cleaves between, with the bond between
/// them, grouped by the first environment. Mirrors `BRICS.reactionDefs`.
pub fn reaction_defs() -> Vec<Vec<(&'static str, &'static str, &'static str)>> {
    vec![
        // L1
        vec![("1", "3", "-"), ("1", "5", "-"), ("1", "10", "-")],
        // L3
        vec![
            ("3", "4", "-"),
            ("3", "13", "-"),
            ("3", "14", "-"),
            ("3", "15", "-"),
            ("3", "16", "-"),
        ],
        // L4
        vec![("4", "5", "-"), ("4", "11", "-")],
        // L5
        vec![
            ("5", "12", "-"),
            ("5", "14", "-"),
            ("5", "16", "-"),
            ("5", "13", "-"),
            ("5", "15", "-"),
        ],
        // L6
        vec![
            ("6", "13", "-"),
            ("6", "14", "-"),
            ("6", "15", "-"),
            ("6", "16", "-"),
        ],
        // L7
        vec![("7a", "7b", "=")],
        // L8
        vec![
            ("8", "9", "-"),
            ("8", "10", "-"),
            ("8", "13", "-"),
            ("8", "14", "-"),
            ("8", "15", "-"),
            ("8", "16", "-"),
        ],
        // L9
        vec![
            ("9", "13", "-"),
            ("9", "14", "-"),
            ("9", "15", "-"),
            ("9", "16", "-"),
        ],
        // L10
        vec![
            ("10", "13", "-"),
            ("10", "14", "-"),
            ("10", "15", "-"),
            ("10", "16", "-"),
        ],
        // L11
        vec![
            ("11", "13", "-"),
            ("11", "14", "-"),
            ("11", "15", "-"),
            ("11", "16", "-"),
        ],
        // L13
        vec![("13", "14", "-"), ("13", "15", "-"), ("13", "16", "-")],
        // L14
        vec![("14", "14", "-"), ("14", "15", "-"), ("14", "16", "-")],
        // L15
        vec![("15", "16", "-")],
        // L16
        vec![("16", "16", "-")],
    ]
}

/// A label without its letters, as `re.sub('[a-z,A-Z]', '', label)`.
fn strip_letters(label: &str) -> String {
    label
        .chars()
        .filter(|c| !c.is_ascii_alphabetic() && *c != ',')
        .collect()
}

/// The reaction SMARTS that break each pair of `reaction_defs`, capping each
/// side with a dummy labeled by its environment. Mirrors `BRICS.smartsGps`.
pub fn smarts_gps() -> Vec<Vec<String>> {
    reaction_defs()
        .into_iter()
        .map(|gp| {
            gp.into_iter()
                .map(|(g1, g2, bnd)| {
                    let r1 = environ(&format!("L{g1}"));
                    let r2 = environ(&format!("L{g2}"));
                    let (g1, g2) = (strip_letters(g1), strip_letters(g2));
                    format!("[$({r1}):1]{bnd};!@[$({r2}):2]>>[{g1}*]-[*:1].[{g2}*]-[*:2]")
                })
                .collect()
        })
        .collect()
}

/// Each environment's pattern, by label, in `environs`' order. Mirrors
/// `BRICS.environMatchers`.
pub fn environ_matchers() -> Vec<(&'static str, Mol)> {
    environs()
        .into_iter()
        .map(|(env, sma)| (env, mol_from_smarts(sma).expect("the BRICS SMARTS parse")))
        .collect()
}

/// One of `bond_matchers`: Python's tuple `(i1, i2, bType, patt)`.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct BRICSBondMatcher {
    /// The first environment's label, without its `L`.
    pub i1: &'static str,
    /// The second environment's label, without its `L`.
    pub i2: &'static str,
    /// The bond between them, as SMARTS.
    pub b_type: &'static str,
    /// A pattern for the bond between the two environments.
    pub patt: Mol,
}

/// For each group of `reaction_defs`, each pair's labels, its bond, and a
/// pattern for the bond between the two environments. Mirrors
/// `BRICS.bondMatchers`.
pub fn bond_matchers() -> Vec<Vec<BRICSBondMatcher>> {
    reaction_defs()
        .into_iter()
        .map(|compats| {
            compats
                .into_iter()
                .map(|(i1, i2, b_type)| {
                    let e1 = environ(&format!("L{i1}"));
                    let e2 = environ(&format!("L{i2}"));
                    let patt = mol_from_smarts(&format!("[$({e1})]{b_type};!@[$({e2})]"))
                        .expect("the BRICS SMARTS parse");
                    BRICSBondMatcher {
                        i1,
                        i2,
                        b_type,
                        patt,
                    }
                })
                .collect()
        })
        .collect()
}

/// The reactions of `smarts_gps`. Mirrors `BRICS.reactions`.
pub fn reactions() -> Vec<Vec<ChemicalReaction>> {
    smarts_gps()
        .iter()
        .map(|gp| gp.iter().map(|sma| reaction(sma)).collect())
        .collect()
}

fn reaction(smarts: &str) -> ChemicalReaction {
    reaction_from_smarts(smarts).expect("the BRICS reaction SMARTS parse")
}

/// A reverse reaction, and patterns for the dummies of its two reactants.
struct ReverseReaction {
    rxn: ChemicalReaction,
    matchers: [Mol; 2],
}

fn reverse_reactions_with_matchers() -> Vec<ReverseReaction> {
    let mut res = Vec::new();
    for gp in reaction_defs() {
        for (g1, g2, bnd) in gp {
            let r1 = environ(&format!("L{g1}"));
            let r2 = environ(&format!("L{g2}"));
            let (g1, g2) = (strip_letters(g1), strip_letters(g2));
            let rs = format!("[$({r1}):1]{bnd};!@[$({r2}):2]");
            let ps = format!("[{g1}*]-[*:1].[{g2}*]-[*:2]");
            let dummy = |label: &str| {
                mol_from_smiles(&format!("[{label}*]")).expect("a labeled dummy parses")
            };
            res.push(ReverseReaction {
                rxn: reaction(&format!("{ps}>>{rs}")),
                matchers: [dummy(&g1), dummy(&g2)],
            });
        }
    }
    res
}

/// The reactions that join the fragments of `reactions` back together, in
/// one list. Python gives each a private `_matchers` attribute, patterns for
/// its reactants' dummies, which `brics_build` uses. Mirrors
/// `BRICS.reverseReactions`.
pub fn reverse_reactions() -> Vec<ChemicalReaction> {
    reverse_reactions_with_matchers()
        .into_iter()
        .map(|r| r.rxn)
        .collect()
}

/// A dummy atom, which marks a fragment's attachment points. Mirrors
/// `BRICS.dummyPattern`.
pub fn dummy_pattern() -> Mol {
    mol_from_smiles("[*]").expect("a dummy parses")
}

/// A bond BRICS would cleave: its atoms, and the labels of their
/// environments, without letters.
/// A Rust addition, naming the tuples `FindBRICSBonds` yields.
pub type BRICSBond = ((u32, u32), (String, String));

/// The bonds BRICS would cleave, each once, in the order of
/// `bond_matchers`. Mirrors `BRICS.FindBRICSBonds(mol)`, whose generator
/// this collects.
pub fn find_brics_bonds(mol: &Mol) -> Result<Vec<BRICSBond>> {
    find_bonds(mol, None)
}

/// As [`find_brics_bonds`], with the groups of `bond_matchers` and the
/// pairs in each in a random order. Mirrors
/// `BRICS.FindBRICSBonds(mol, randomizeOrder=True)`, which draws from the
/// `random` module's generator.
pub fn find_brics_bonds_with_random_order(
    mol: &Mol,
    random: &mut Random,
) -> Result<Vec<BRICSBond>> {
    find_bonds(mol, Some(random))
}

fn find_bonds(mol: &Mol, mut random: Option<&mut Random>) -> Result<Vec<BRICSBond>> {
    let bond_matchers = bond_matchers();
    let mut indices: Vec<usize> = (0..bond_matchers.len()).collect();
    if let Some(random) = random.as_deref_mut() {
        shuffle_with_random(&mut indices, || random.random()).expect("random() is in [0, 1)");
    }
    let mut env_matches = std::collections::HashMap::new();
    for (env, patt) in environ_matchers() {
        env_matches.insert(env.to_owned(), mol.has_substruct_match(&patt)?);
    }
    let mut bonds_done = HashSet::new();
    let mut res = Vec::new();
    for gp_idx in indices {
        let mut compats: Vec<&BRICSBondMatcher> = bond_matchers[gp_idx].iter().collect();
        if let Some(random) = random.as_deref_mut() {
            shuffle_with_random(&mut compats, || random.random()).expect("random() is in [0, 1)");
        }
        for BRICSBondMatcher { i1, i2, patt, .. } in compats {
            if !env_matches[&format!("L{i1}")] || !env_matches[&format!("L{i2}")] {
                continue;
            }
            let matches = mol.get_substruct_matches(patt)?;
            let (i1, i2) = (strip_letters(i1), strip_letters(i2));
            for m in matches {
                let pair = (m[0], m[1]);
                if !bonds_done.contains(&pair) && !bonds_done.contains(&(pair.1, pair.0)) {
                    bonds_done.insert(pair);
                    res.push((pair, (i1.clone(), i2.clone())));
                }
            }
        }
    }
    Ok(res)
}

/// Parameters for `break_brics_bonds_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct BreakBRICSBondsParams {
    /// The bonds to break, as `find_brics_bonds` returns them; `None` (or
    /// empty, as in Python) for all of them, which RDKit's
    /// `FragmentOnBRICSBonds` finds and breaks.
    pub bonds: Option<Vec<BRICSBond>>,
    /// Sanitizes the result.
    pub sanitize: bool,
}

impl Default for BreakBRICSBondsParams {
    fn default() -> BreakBRICSBondsParams {
        BreakBRICSBondsParams {
            bonds: None,
            sanitize: true,
        }
    }
}

/// The molecule broken at its BRICS bonds, each end capped with a dummy
/// whose isotope is its environment's label, as one molecule of the
/// fragments. Mirrors `BRICS.BreakBRICSBonds(mol)`.
pub fn break_brics_bonds(mol: &Mol) -> Result<Mol> {
    break_brics_bonds_with_params(mol, &BreakBRICSBondsParams::default())
}

/// As [`break_brics_bonds`], for the bonds `params.bonds`, which keeps the
/// molecule's conformers, with each dummy where the atom it replaces was. A
/// bond between atoms that are not bonded is `Error::Value`, and a label
/// that is not a number is `Error::Value`, as in Python. Mirrors
/// `BRICS.BreakBRICSBonds(mol, bonds, sanitize)`.
pub fn break_brics_bonds_with_params(mol: &Mol, params: &BreakBRICSBondsParams) -> Result<Mol> {
    let bonds = match &params.bonds {
        Some(bonds) if !bonds.is_empty() => bonds,
        _ => {
            let mut res = fragment_on_brics_bonds(mol)?;
            if params.sanitize {
                sanitize_mol(&mut res)?;
            }
            return Ok(res);
        }
    };
    let mut e_mol = EditableMol::new(mol);
    let mut n_ats = mol.get_num_atoms();
    let mut dummy_positions = Vec::new();
    let dummy = |label: &str| -> Result<Atom> {
        let isotope: u32 = label.trim().parse().map_err(|_| Error::Value {
            message: format!("invalid literal for int() with base 10: '{label}'"),
        })?;
        let mut atom = Atom::from_num(0);
        atom.set_isotope(isotope);
        atom.set_no_implicit(true);
        Ok(atom)
    };
    for ((ia, ib), (da, db)) in bonds {
        let (ia, ib) = (*ia, *ib);
        let bond_type = mol
            .get_bond_between_atoms(ia, ib)?
            .ok_or_else(|| Error::Value {
                message: format!("atoms {ia} and {ib} are not bonded"),
            })?
            .get_bond_type();
        e_mol.remove_bond(ia, ib)?;
        let idxa = n_ats;
        n_ats += 1;
        e_mol.add_atom(&dummy(da)?);
        e_mol.add_bond_with_order(ia, idxa, bond_type)?;
        let idxb = n_ats;
        n_ats += 1;
        e_mol.add_atom(&dummy(db)?);
        e_mol.add_bond_with_order(ib, idxb, bond_type)?;
        if mol.get_num_conformers() > 0 {
            dummy_positions.push((idxa, ib));
            dummy_positions.push((idxb, ia));
        }
    }
    let mut res = e_mol.get_mol();
    if params.sanitize {
        sanitize_mol(&mut res)?;
    }
    for conf in mol.get_conformers() {
        let id = conf.get_id() as i32;
        let mut res_conf = res.get_conformer_with_id_mut(id)?;
        for &(ia, pa) in &dummy_positions {
            res_conf.set_atom_position(ia, conf.get_atom_position(pa)?)?;
        }
    }
    Ok(res)
}

/// Parameters for `brics_decompose_with_params` and its forms, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct BRICSDecomposeParams {
    /// The SMILES of the nodes already found, from earlier decompositions;
    /// a molecule among them decomposes to nothing. Python's set grows with
    /// the decomposition's nodes, which `brics_decompose_with_outputs`
    /// returns (A9).
    pub all_nodes: Option<BTreeSet<String>>,
    /// Rejects a cleavage that leaves a fragment with fewer atoms, not
    /// counting dummies.
    pub min_fragment_size: u32,
    /// Only these reactions, as (group, reaction) indices into `reactions`;
    /// `None` (or empty, as in Python) for all.
    pub only_use_reactions: Option<Vec<(usize, usize)>>,
    /// Also returns the fragments that were broken further.
    pub keep_non_leaf_nodes: bool,
    /// Breaks the molecule once with each reaction, without breaking the
    /// pieces further.
    pub single_pass: bool,
}

impl Default for BRICSDecomposeParams {
    fn default() -> BRICSDecomposeParams {
        BRICSDecomposeParams {
            all_nodes: None,
            min_fragment_size: 1,
            only_use_reactions: None,
            keep_non_leaf_nodes: false,
            single_pass: false,
        }
    }
}

/// What `brics_decompose_with_outputs` returns besides the fragments.
/// A Rust addition, for Python's output arguments (A9).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct BRICSDecomposeOutputs {
    /// `params.all_nodes` with this decomposition's nodes added.
    pub all_nodes: BTreeSet<String>,
}

/// The SMILES of the molecule's BRICS fragments: what remains after every
/// cleavage, each once. Python returns a set; the SMILES here are sorted.
/// Mirrors `BRICS.BRICSDecompose(mol)`.
pub fn brics_decompose(mol: &Mol) -> Result<BTreeSet<String>> {
    brics_decompose_with_params(mol, &BRICSDecomposeParams::default())
}

/// As [`brics_decompose`]; with `keep_non_leaf_nodes` or `single_pass`,
/// every node found, `params.all_nodes` included, as Python returns that
/// set. Mirrors `BRICS.BRICSDecompose(mol, allNodes, minFragmentSize,
/// onlyUseReactions, silent, keepNonLeafNodes, singlePass)`.
pub fn brics_decompose_with_params(
    mol: &Mol,
    params: &BRICSDecomposeParams,
) -> Result<BTreeSet<String>> {
    Ok(brics_decompose_with_outputs(mol, params)?.0)
}

/// As [`brics_decompose_with_params`], with the nodes found added to
/// `params.all_nodes`, as Python adds them to its `allNodes` (A9).
pub fn brics_decompose_with_outputs(
    mol: &Mol,
    params: &BRICSDecomposeParams,
) -> Result<(BTreeSet<String>, BRICSDecomposeOutputs)> {
    let mut all_nodes = params.all_nodes.clone().unwrap_or_default();
    let res = match decompose(mol, params, &mut all_nodes)? {
        None => BTreeSet::new(),
        Some((leaves, _)) if !(params.single_pass || params.keep_non_leaf_nodes) => {
            leaves.into_iter().map(|(smi, _)| smi).collect()
        }
        Some(_) => all_nodes.clone(),
    };
    Ok((res, BRICSDecomposeOutputs { all_nodes }))
}

/// The molecule's BRICS fragments as molecules, in the order Python's
/// dictionaries hold them; with `keep_non_leaf_nodes` or `single_pass`,
/// every node found in this call, the molecule first. A molecule in
/// `params.all_nodes` gives none. Mirrors `BRICS.BRICSDecompose(mol, ...,
/// returnMols=True)` (A4).
pub fn brics_decompose_return_mols(mol: &Mol, params: &BRICSDecomposeParams) -> Result<Vec<Mol>> {
    let mut all_nodes = params.all_nodes.clone().unwrap_or_default();
    Ok(match decompose(mol, params, &mut all_nodes)? {
        None => Vec::new(),
        Some((leaves, _)) if !(params.single_pass || params.keep_non_leaf_nodes) => {
            leaves.into_iter().map(|(_, m)| m).collect()
        }
        Some((_, found)) => found.into_iter().map(|(_, m)| m).collect(),
    })
}

/// A list of (SMILES, molecule) pairs with distinct SMILES, in the order a
/// Python dictionary keeps them.
type Pool = Vec<(String, Mol)>;

/// The decomposition: the final pool and every node found, or `None` for a
/// molecule already among `all_nodes`.
fn decompose(
    mol: &Mol,
    params: &BRICSDecomposeParams,
    all_nodes: &mut BTreeSet<String>,
) -> Result<Option<(Pool, Pool)>> {
    let m_smi = mol_to_smiles(mol)?;
    if all_nodes.contains(&m_smi) {
        return Ok(None);
    }
    let only_use = params
        .only_use_reactions
        .as_ref()
        .filter(|only| !only.is_empty());
    let mut active_pool: VecDeque<(String, Mol)> = VecDeque::from([(m_smi.clone(), mol.clone())]);
    all_nodes.insert(m_smi.clone());
    let mut found_mols: Pool = vec![(m_smi, mol.clone())];
    for (gp_idx, mut reaction_gp) in reactions().into_iter().enumerate() {
        let mut new_pool = VecDeque::new();
        while let Some((n_smi, mol)) = active_pool.pop_front() {
            let mut matched = false;
            for (rxn_idx, reaction) in reaction_gp.iter_mut().enumerate() {
                if only_use.is_some_and(|only| !only.contains(&(gp_idx, rxn_idx))) {
                    continue;
                }
                let ps = reaction.run_reactants(std::slice::from_ref(&mol))?;
                for prod_seq in ps {
                    let mut seq_ok = true;
                    // Python sorts the products by size, so that small
                    // fragments disqualify a cleavage first.
                    let mut t_seq = Vec::with_capacity(prod_seq.len());
                    for (idx, prod) in prod_seq.iter().enumerate() {
                        t_seq.push((prod.get_num_atoms_with_only_explicit(true)?, idx));
                    }
                    t_seq.sort();
                    let mut prods: Vec<Option<Mol>> = prod_seq.into_iter().map(Some).collect();
                    // Each product with its SMILES, or the error that
                    // sanitizing it raised, which Python raises when it
                    // reads the SMILES it never set.
                    let mut sized: Vec<(Mol, Result<String>)> = Vec::new();
                    for &(nats, idx) in &t_seq {
                        let mut prod = prods[idx].take().expect("each product once");
                        if let Err(e) = sanitize_mol(&mut prod) {
                            sized.push((prod, Err(e)));
                            continue;
                        }
                        let p_smi = mol_to_smiles(&prod)?;
                        if params.min_fragment_size > 0 {
                            let n_dummies = p_smi.matches('*').count() as u32;
                            if nats < params.min_fragment_size + n_dummies {
                                seq_ok = false;
                                break;
                            }
                        }
                        sized.push((prod, Ok(p_smi)));
                    }
                    if seq_ok {
                        matched = true;
                        for (prod, p_smi) in sized {
                            let p_smi = p_smi?;
                            if !all_nodes.contains(&p_smi) {
                                if !params.single_pass {
                                    active_pool.push_back((p_smi.clone(), prod.clone()));
                                }
                                all_nodes.insert(p_smi.clone());
                                found_mols.push((p_smi, prod));
                            }
                        }
                    }
                }
            }
            if params.single_pass || params.keep_non_leaf_nodes || !matched {
                new_pool.push_back((n_smi, mol));
            }
        }
        active_pool = new_pool;
    }
    Ok(Some((active_pool.into_iter().collect(), found_mols)))
}

/// Parameters for `brics_build_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug)]
pub struct BRICSBuildParams {
    /// Yields only molecules without attachment points.
    pub only_complete_mols: bool,
    /// The molecules to build from; `None` (or empty, as in Python) for the
    /// fragments themselves.
    pub seeds: Option<Vec<Mol>>,
    /// Yields each molecule, by its canonical SMILES, once.
    pub uniquify: bool,
    /// Shuffles the seeds and reactions at each step, drawing from the
    /// `Random` the build is given.
    pub scramble_reagents: bool,
    /// How many times a product with attachment points is built on again.
    pub max_depth: u32,
}

impl Default for BRICSBuildParams {
    fn default() -> BRICSBuildParams {
        BRICSBuildParams {
            only_complete_mols: true,
            seeds: None,
            uniquify: true,
            scramble_reagents: true,
            max_depth: 3,
        }
    }
}

/// One level of Python's recursive generator.
struct Frame {
    seeds: Vec<Mol>,
    seed_idx: usize,
    rxn_order: Vec<usize>,
    rxn_pos: usize,
    rxn_checked: bool,
    frag_idx: usize,
    seed_is_r1: bool,
    seed_is_r2: bool,
    next_steps: Vec<Mol>,
    pending: VecDeque<Mol>,
    seen: HashSet<String>,
    max_depth: u32,
}

enum Step {
    Yield(Mol),
    Child(Vec<Mol>),
    Done,
}

/// New molecules built from BRICS fragments, made as they are asked for.
/// Built by [`brics_build`].
/// A Rust addition: the iterator of Python's generator.
pub struct BRICSBuild<'a> {
    fragments: &'a [Mol],
    random: &'a mut Random,
    params: BRICSBuildParams,
    reverse: Vec<ReverseReaction>,
    dummy: Mol,
    stack: Vec<Frame>,
}

impl std::fmt::Debug for BRICSBuild<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BRICSBuild").finish_non_exhaustive()
    }
}

impl BRICSBuild<'_> {
    fn frame(&mut self, seeds: Vec<Mol>, max_depth: u32) -> Frame {
        let mut seeds = seeds;
        let mut rxn_order: Vec<usize> = (0..self.reverse.len()).collect();
        if self.params.scramble_reagents {
            shuffle_with_random(&mut seeds, || self.random.random())
                .expect("random() is in [0, 1)");
            shuffle_with_random(&mut rxn_order, || self.random.random())
                .expect("random() is in [0, 1)");
        }
        Frame {
            seeds,
            seed_idx: 0,
            rxn_order,
            rxn_pos: 0,
            rxn_checked: false,
            frag_idx: 0,
            seed_is_r1: false,
            seed_is_r2: false,
            next_steps: Vec::new(),
            pending: VecDeque::new(),
            seen: HashSet::new(),
            max_depth,
        }
    }

    /// The next thing the top frame does.
    fn step(&mut self) -> Result<Step> {
        let uniquify = self.params.uniquify;
        let only_complete = self.params.only_complete_mols;
        let Some(f) = self.stack.last_mut() else {
            return Ok(Step::Done);
        };
        loop {
            if let Some(p) = f.pending.pop_front() {
                if uniquify && !f.seen.insert(mol_to_smiles(&p)?) {
                    continue;
                }
                if p.has_substruct_match(&self.dummy)? {
                    f.next_steps.push(p.clone());
                    if !only_complete {
                        return Ok(Step::Yield(p));
                    }
                } else {
                    return Ok(Step::Yield(p));
                }
                continue;
            }
            if f.seed_idx >= f.seeds.len() {
                return Ok(Step::Done);
            }
            if f.rxn_pos >= f.rxn_order.len() {
                // This seed is done; build on its products.
                let steps = std::mem::take(&mut f.next_steps);
                f.seed_idx += 1;
                f.rxn_pos = 0;
                f.rxn_checked = false;
                f.frag_idx = 0;
                f.seed_is_r1 = false;
                f.seed_is_r2 = false;
                if !steps.is_empty() && f.max_depth > 0 {
                    return Ok(Step::Child(steps));
                }
                continue;
            }
            let rxn = &mut self.reverse[f.rxn_order[f.rxn_pos]];
            let seed = &f.seeds[f.seed_idx];
            if !f.rxn_checked {
                // Python never clears these within a seed.
                if seed.has_substruct_match(&rxn.matchers[0])? {
                    f.seed_is_r1 = true;
                }
                if seed.has_substruct_match(&rxn.matchers[1])? {
                    f.seed_is_r2 = true;
                }
                f.rxn_checked = true;
            }
            let Some(fragment) = self.fragments.get(f.frag_idx) else {
                f.rxn_pos += 1;
                f.rxn_checked = false;
                f.frag_idx = 0;
                continue;
            };
            f.frag_idx += 1;
            let mut ps = None;
            if fragment.has_substruct_match(&rxn.matchers[0])? && f.seed_is_r2 {
                ps = Some(rxn.rxn.run_reactants(&[fragment.clone(), seed.clone()])?);
            }
            if fragment.has_substruct_match(&rxn.matchers[1])? && f.seed_is_r1 {
                ps = Some(rxn.rxn.run_reactants(&[seed.clone(), fragment.clone()])?);
            }
            for p in ps.into_iter().flatten() {
                if let Some(first) = p.into_iter().next() {
                    f.pending.push_back(first);
                }
            }
        }
    }

    fn next_mol(&mut self) -> Result<Option<Mol>> {
        loop {
            match self.step()? {
                Step::Done => {
                    if self.stack.pop().is_none() || self.stack.is_empty() {
                        return Ok(None);
                    }
                }
                Step::Child(seeds) => {
                    let max_depth = self.stack.last().map_or(0, |f| f.max_depth) - 1;
                    let child = self.frame(seeds, max_depth);
                    self.stack.push(child);
                }
                Step::Yield(mol) => {
                    // Each level above the one that built the molecule
                    // passes it on only if it has not yielded it.
                    let mut ok = true;
                    if self.params.uniquify && self.stack.len() > 1 {
                        let smi = mol_to_smiles(&mol)?;
                        let n = self.stack.len();
                        for frame in self.stack[..n - 1].iter_mut().rev() {
                            if !frame.seen.insert(smi.clone()) {
                                ok = false;
                                break;
                            }
                        }
                    }
                    if ok {
                        return Ok(Some(mol));
                    }
                }
            }
        }
    }
}

/// The next molecule built.
impl Iterator for BRICSBuild<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        match self.next_mol() {
            Ok(mol) => mol.map(Ok),
            Err(e) => {
                self.stack.clear();
                Some(Err(e))
            }
        }
    }
}

/// New molecules built from BRICS fragments by joining their attachment
/// points with BRICS's reverse reactions, three levels deep, each once and
/// each complete, with the seeds and reactions in an order drawn from
/// `random`. Mirrors `BRICS.BRICSBuild(fragments)`, which draws from the
/// `random` module's generator.
pub fn brics_build<'a>(fragments: &'a [Mol], random: &'a mut Random) -> BRICSBuild<'a> {
    brics_build_with_params(fragments, random, &BRICSBuildParams::default())
}

/// As [`brics_build`]; `random` is drawn from only with
/// `params.scramble_reagents`. Mirrors `BRICS.BRICSBuild(fragments,
/// onlyCompleteMols, seeds, uniquify, scrambleReagents, maxDepth)`.
pub fn brics_build_with_params<'a>(
    fragments: &'a [Mol],
    random: &'a mut Random,
    params: &BRICSBuildParams,
) -> BRICSBuild<'a> {
    let mut build = BRICSBuild {
        fragments,
        random,
        params: params.clone(),
        reverse: reverse_reactions_with_matchers(),
        dummy: dummy_pattern(),
        stack: Vec::new(),
    };
    let seeds = match &params.seeds {
        Some(seeds) if !seeds.is_empty() => seeds.clone(),
        _ => fragments.to_vec(),
    };
    let root = build.frame(seeds, params.max_depth);
    build.stack.push(root);
    build
}
