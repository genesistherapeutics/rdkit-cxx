// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Recap`, a Python module: the RECAP algorithm of
//! Lewell et al., _J. Chem. Inf. Comput. Sci._ 38, 511-522 (1998), which
//! breaks molecules at bonds chemical reactions make, into a hierarchy of
//! fragments.
//!
//! Python's nodes hold their children, and weak references to their
//! parents, in dictionaries by SMILES, and its `allNodes` dictionary holds
//! every node of a decomposition. Here a [`RecapHierarchy`] owns the nodes,
//! as that dictionary does, and a [`RecapHierarchyNode`] borrows one; maps
//! by SMILES are lists of (SMILES, node) pairs in the order Python's
//! dictionaries keep them.

use std::collections::HashMap;

use crate::Result;
use crate::chem::rd_chem_reactions::{ChemicalReaction, reaction_from_smarts};
use crate::chem::{Mol, mol_to_smiles, sanitize_mol};

/// The RECAP reactions, as reaction SMARTS: urea, amide, ester, amines,
/// cyclic amines, ether, olefin, aromatic nitrogen to aliphatic carbon,
/// lactam nitrogen to aliphatic carbon, aromatic carbon to aromatic carbon,
/// aromatic nitrogen to aromatic carbon (not in the standard set), and
/// sulphonamide. Mirrors `Recap.reactionDefs`.
pub fn reaction_defs() -> Vec<&'static str> {
    vec![
        "[#7;+0;D2,D3:1]!@C(!@=O)!@[#7;+0;D2,D3:2]>>*[#7:1].[#7:2]*",
        "[C;!$(C([#7])[#7]):1](=!@[O:2])!@[#7;+0;!D1:3]>>*[C:1]=[O:2].*[#7:3]",
        "[C:1](=!@[O:2])!@[O;+0:3]>>*[C:1]=[O:2].[O:3]*",
        "[N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*:1])-!@[*:2]>>*[*:1].[*:2]*",
        "[#7;R;D3;+0:1]-!@[*:2]>>*[#7:1].[*:2]*",
        "[#6:1]-!@[O;+0]-!@[#6:2]>>[#6:1]*.*[#6:2]",
        "[C:1]=!@[C:2]>>[C:1]*.*[C:2]",
        "[n;+0:1]-!@[C:2]>>[n:1]*.[C:2]*",
        "[O:3]=[C:4]-@[N;+0:1]-!@[C:2]>>[O:3]=[C:4]-[N:1]*.[C:2]*",
        "[c:1]-!@[c:2]>>[c:1]*.*[c:2]",
        "[n;+0:1]-!@[c:2]>>[n:1]*.*[c:2]",
        "[#7;+0;D2,D3:1]-!@[S:2](=[O:3])=[O:4]>>[#7:1]*.*[S:2](=[O:3])=[O:4]",
    ]
}

/// The reactions of `reaction_defs`. Mirrors `Recap.reactions`.
pub fn reactions() -> Vec<ChemicalReaction> {
    reaction_defs()
        .into_iter()
        .map(|sma| reaction_from_smarts(sma).expect("the RECAP reaction SMARTS parse"))
        .collect()
}

/// A list of (SMILES, node index) pairs that keeps the first position of
/// each SMILES, as assigning to a Python dictionary does.
#[derive(Clone, Debug, Default)]
struct NodeMap(Vec<(String, usize)>);

impl NodeMap {
    fn insert(&mut self, smiles: &str, idx: usize) {
        match self.0.iter_mut().find(|(s, _)| s == smiles) {
            Some(entry) => entry.1 = idx,
            None => self.0.push((smiles.to_owned(), idx)),
        }
    }
}

#[derive(Clone, Debug)]
struct NodeData {
    mol: Mol,
    smiles: String,
    children: NodeMap,
    parents: NodeMap,
}

/// The nodes of RECAP decompositions, by SMILES: what Python's `allNodes`
/// dictionary holds. Decomposing more molecules into one hierarchy reuses
/// the nodes it has. `Default` is an empty hierarchy.
#[derive(Clone, Debug, Default)]
pub struct RecapHierarchy {
    nodes: Vec<NodeData>,
    by_smiles: HashMap<String, usize>,
    root: Option<usize>,
}

/// Parameters for `recap_decompose_with_params` and
/// `RecapHierarchy::recap_decompose`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RecapDecomposeParams {
    /// Rejects a cleavage that leaves a fragment with fewer atoms, not
    /// counting dummies; 0 instead rejects fragments of at most three
    /// aliphatic carbons.
    pub min_fragment_size: u32,
    /// Only these reactions, as indices into `reactions`; `None` (or empty,
    /// as in Python) for all.
    pub only_use_reactions: Option<Vec<usize>>,
}

/// The molecule's RECAP decomposition, whose root is the molecule. Mirrors
/// `Recap.RecapDecompose(mol)`, which returns the root node.
pub fn recap_decompose(mol: &Mol) -> Result<RecapHierarchy> {
    recap_decompose_with_params(mol, &RecapDecomposeParams::default())
}

/// As [`recap_decompose`]. Mirrors `Recap.RecapDecompose(mol,
/// minFragmentSize=..., onlyUseReactions=...)`.
pub fn recap_decompose_with_params(
    mol: &Mol,
    params: &RecapDecomposeParams,
) -> Result<RecapHierarchy> {
    let mut hierarchy = RecapHierarchy::default();
    hierarchy.recap_decompose(mol, params)?;
    Ok(hierarchy)
}

impl RecapHierarchy {
    /// An empty hierarchy, as Python's `allNodes={}`.
    pub fn new() -> RecapHierarchy {
        RecapHierarchy::default()
    }

    /// Decomposes the molecule into this hierarchy, adding the nodes it
    /// lacks, and returns the molecule's node, which becomes the root; a
    /// molecule already in the hierarchy is not decomposed again. Mirrors
    /// `Recap.RecapDecompose(mol, allNodes, minFragmentSize,
    /// onlyUseReactions)` with `allNodes` this hierarchy.
    ///
    /// A product of a cleavage RDKit cannot sanitize is its sanitization
    /// error, as in Python, which fails when it reads the product's SMILES.
    pub fn recap_decompose(
        &mut self,
        mol: &Mol,
        params: &RecapDecomposeParams,
    ) -> Result<RecapHierarchyNode<'_>> {
        let m_smi = mol_to_smiles(mol)?;
        if let Some(&idx) = self.by_smiles.get(&m_smi) {
            self.root = Some(idx);
            return Ok(self.node(idx));
        }
        let only_use = params
            .only_use_reactions
            .as_ref()
            .filter(|only| !only.is_empty());
        let res = self.add(mol.clone(), m_smi);
        let mut active_pool = std::collections::VecDeque::from([res]);
        let mut reactions = reactions();
        while let Some(node) = active_pool.pop_front() {
            let n_smi = self.nodes[node].smiles.clone();
            for (rxn_idx, reaction) in reactions.iter_mut().enumerate() {
                if only_use.is_some_and(|only| !only.contains(&rxn_idx)) {
                    continue;
                }
                let ps = reaction.run_reactants(std::slice::from_ref(&self.nodes[node].mol))?;
                for prod_seq in ps {
                    let mut t_seq = Vec::with_capacity(prod_seq.len());
                    for (idx, prod) in prod_seq.iter().enumerate() {
                        t_seq.push((prod.get_num_atoms_with_only_explicit(true)?, idx));
                    }
                    t_seq.sort();
                    let mut prods: Vec<Option<Mol>> = prod_seq.into_iter().map(Some).collect();
                    let mut seq_ok = true;
                    let mut sized: Vec<(Mol, Result<String>)> = Vec::new();
                    for &(nats, idx) in &t_seq {
                        let mut prod = prods[idx].take().expect("each product once");
                        if let Err(e) = sanitize_mol(&mut prod) {
                            sized.push((prod, Err(e)));
                            continue;
                        }
                        let p_smi = mol_to_smiles(&prod)?;
                        if params.min_fragment_size > 0 {
                            let n_dummies = p_smi.matches('*').count() as i64;
                            if i64::from(nats) - n_dummies < i64::from(params.min_fragment_size) {
                                seq_ok = false;
                                break;
                            }
                        } else if ["", "C", "CC", "CCC"]
                            .contains(&p_smi.replace('*', "").replace("()", "").as_str())
                        {
                            // Too small, once the dummies and the empty
                            // branches they leave are gone.
                            seq_ok = false;
                            break;
                        }
                        sized.push((prod, Ok(p_smi)));
                    }
                    if !seq_ok {
                        continue;
                    }
                    for (prod, p_smi) in sized {
                        let p_smi = p_smi?;
                        let p_node = match self.by_smiles.get(&p_smi) {
                            Some(&p_node) => p_node,
                            None => {
                                let p_node = self.add(prod, p_smi.clone());
                                active_pool.push_back(p_node);
                                p_node
                            }
                        };
                        self.nodes[p_node].parents.insert(&n_smi, node);
                        self.nodes[node].children.insert(&p_smi, p_node);
                    }
                }
            }
        }
        self.root = Some(res);
        Ok(self.node(res))
    }

    fn add(&mut self, mol: Mol, smiles: String) -> usize {
        let idx = self.nodes.len();
        self.by_smiles.insert(smiles.clone(), idx);
        self.nodes.push(NodeData {
            mol,
            smiles,
            children: NodeMap::default(),
            parents: NodeMap::default(),
        });
        idx
    }

    fn node(&self, idx: usize) -> RecapHierarchyNode<'_> {
        RecapHierarchyNode {
            hierarchy: self,
            idx,
        }
    }

    /// The node of the molecule most recently decomposed into the
    /// hierarchy: what Python's `RecapDecompose` returns. `None` for an
    /// empty hierarchy.
    pub fn root(&self) -> Option<RecapHierarchyNode<'_>> {
        self.root.map(|idx| self.node(idx))
    }

    /// The node with this SMILES, as `allNodes[smiles]`.
    pub fn get(&self, smiles: &str) -> Option<RecapHierarchyNode<'_>> {
        self.by_smiles.get(smiles).map(|&idx| self.node(idx))
    }

    /// The SMILES of every node, in the order the nodes were found, as
    /// `allNodes.keys()`.
    pub fn keys(&self) -> Vec<&str> {
        self.nodes.iter().map(|n| n.smiles.as_str()).collect()
    }

    /// The number of nodes, as `len(allNodes)`.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether the hierarchy has no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

/// A node of a RECAP hierarchy: a molecule or fragment, and the fragments
/// it breaks into. It borrows the hierarchy that owns it. Mirrors
/// `Recap.RecapHierarchyNode`.
#[derive(Clone, Copy)]
pub struct RecapHierarchyNode<'h> {
    hierarchy: &'h RecapHierarchy,
    idx: usize,
}

impl std::fmt::Debug for RecapHierarchyNode<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("RecapHierarchyNode")
            .field(&self.smiles())
            .finish()
    }
}

impl PartialEq for RecapHierarchyNode<'_> {
    /// The same node of the same hierarchy, as Python compares nodes.
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.hierarchy, other.hierarchy) && self.idx == other.idx
    }
}

impl Eq for RecapHierarchyNode<'_> {}

/// (SMILES, node) pairs, in the order a Python dictionary keeps them.
pub type RecapNodes<'h> = Vec<(&'h str, RecapHierarchyNode<'h>)>;

impl<'h> RecapHierarchyNode<'h> {
    fn data(&self) -> &'h NodeData {
        &self.hierarchy.nodes[self.idx]
    }

    fn nodes(&self, map: &'h NodeMap) -> RecapNodes<'h> {
        map.0
            .iter()
            .map(|(s, idx)| (s.as_str(), self.hierarchy.node(*idx)))
            .collect()
    }

    /// The node's molecule, its dummies marking where it was cut. Mirrors
    /// `RecapHierarchyNode.mol`.
    pub fn mol(&self) -> &'h Mol {
        &self.data().mol
    }

    /// The molecule's canonical SMILES. Mirrors `RecapHierarchyNode.smiles`.
    pub fn smiles(&self) -> &'h str {
        &self.data().smiles
    }

    /// The fragments the molecule breaks into, by SMILES. Mirrors
    /// `RecapHierarchyNode.children`.
    pub fn children(&self) -> RecapNodes<'h> {
        self.nodes(&self.data().children)
    }

    /// The nodes that break into this one, by SMILES. Mirrors
    /// `RecapHierarchyNode.parents`, whose values are weak references.
    pub fn parents(&self) -> RecapNodes<'h> {
        self.nodes(&self.data().parents)
    }

    /// Every node below this one, by SMILES. Mirrors
    /// `RecapHierarchyNode.GetAllChildren`.
    pub fn get_all_children(&self) -> RecapNodes<'h> {
        let mut res = NodeMap::default();
        for (smi, child) in &self.data().children.0 {
            res.insert(smi, *child);
            self.hierarchy.node(*child).gac_recurse(&mut res, false);
        }
        self.nodes_of(res)
    }

    /// The nodes below this one that break no further, by SMILES: the
    /// standard RECAP fragments. Mirrors `RecapHierarchyNode.GetLeaves`.
    pub fn get_leaves(&self) -> RecapNodes<'h> {
        let mut res = NodeMap::default();
        for (smi, child) in &self.data().children.0 {
            let child_node = self.hierarchy.node(*child);
            if child_node.data().children.0.is_empty() {
                res.insert(smi, *child);
            } else {
                child_node.gac_recurse(&mut res, true);
            }
        }
        self.nodes_of(res)
    }

    /// The nodes of a map, by their SMILES, which are its keys.
    fn nodes_of(&self, map: NodeMap) -> RecapNodes<'h> {
        map.0
            .into_iter()
            .map(|(_, idx)| {
                let node = self.hierarchy.node(idx);
                (node.smiles(), node)
            })
            .collect()
    }

    fn gac_recurse(&self, res: &mut NodeMap, terminal_only: bool) {
        for (smi, child) in &self.data().children.0 {
            let child_node = self.hierarchy.node(*child);
            if !terminal_only || child_node.data().children.0.is_empty() {
                res.insert(smi, *child);
            }
            child_node.gac_recurse(res, terminal_only);
        }
    }

    /// The nodes above this one that have no parents, this one when it has
    /// none. Mirrors `RecapHierarchyNode.getUltimateParents`.
    pub fn get_ultimate_parents(&self) -> Vec<RecapHierarchyNode<'h>> {
        let parents = &self.data().parents.0;
        if parents.is_empty() {
            return vec![*self];
        }
        let mut res: Vec<RecapHierarchyNode<'h>> = Vec::new();
        for (_, p) in parents {
            for u_p in self.hierarchy.node(*p).get_ultimate_parents() {
                if !res.contains(&u_p) {
                    res.push(u_p);
                }
            }
        }
        res
    }
}
