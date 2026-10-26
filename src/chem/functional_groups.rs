// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FunctionalGroups`, a Python module: a hierarchy of
//! functional groups, read from a file of labels, SMARTS, and names, and
//! fingerprints that count a molecule's groups.
//!
//! Python reads its default hierarchy from `Functional_Group_Hierarchy.txt`
//! in RDKit's data directory; the crate embeds that file
//! (`wrapper/third_party/rdkit/Data/Functional_Group_Hierarchy.txt`). Python
//! also keeps the last hierarchy it built in module globals, and returns it
//! again unless asked to `force` a new one; here each call builds a new one.

use std::path::{Path, PathBuf};

use crate::chem::py_compat::{splitlines, str_repr};
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::{Mol, mol_from_smarts};
use crate::{Error, Result};

/// RDKit's functional group hierarchy, as `BuildFuncGroupHierarchy` reads it
/// by default.
pub(crate) const DEFAULT_HIERARCHY: &str =
    include_str!("../../wrapper/third_party/rdkit/Data/Functional_Group_Hierarchy.txt");

/// A functional group, and the more specific groups under it. Mirrors
/// `FunctionalGroups.FGHierarchyNode`.
#[derive(Clone, Debug)]
pub struct FGHierarchyNode {
    /// The group's name, such as `Aromatic`.
    pub name: String,
    /// The group's pattern, whose first atom is the group's distinguishing
    /// feature.
    pub pattern: Mol,
    /// The pattern's SMARTS.
    pub smarts: String,
    /// The group's label, its parents' labels before it, separated by dots,
    /// such as `AcidChloride.Aromatic`.
    pub label: String,
    /// A reaction SMARTS that removes the group; empty for none.
    pub rxn_smarts: String,
    /// The label of the group above this one; `None` at the top. Python
    /// holds a weak reference to the parent node.
    pub parent: Option<String>,
    /// The groups under this one.
    pub children: Vec<FGHierarchyNode>,
    /// A reaction that removes the group; Python's module never sets it.
    pub removal_reaction: Option<ChemicalReaction>,
}

/// Parameters for `FGHierarchyNode::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FGHierarchyNodeParams {
    /// The pattern's SMARTS.
    pub smarts: String,
    /// The group's label.
    pub label: String,
    /// A reaction SMARTS that removes the group.
    pub rxn_smarts: String,
    /// The label of the group above this one; Python takes the parent node.
    pub parent: Option<String>,
}

impl FGHierarchyNode {
    /// A group named `name` with the pattern `patt`, and no children.
    /// Mirrors `FGHierarchyNode(name, patt)`.
    pub fn new(name: &str, patt: Mol) -> FGHierarchyNode {
        FGHierarchyNode::new_with_params(name, patt, &FGHierarchyNodeParams::default())
    }

    /// Mirrors `FGHierarchyNode(name, patt, smarts, label, rxnSmarts,
    /// parent)`.
    pub fn new_with_params(
        name: &str,
        patt: Mol,
        params: &FGHierarchyNodeParams,
    ) -> FGHierarchyNode {
        FGHierarchyNode {
            name: name.to_owned(),
            pattern: patt,
            smarts: params.smarts.clone(),
            label: params.label.clone(),
            rxn_smarts: params.rxn_smarts.clone(),
            parent: params.parent.clone(),
            children: Vec::new(),
            removal_reaction: None,
        }
    }

    /// The number of groups in the node's subtree, the node's own included.
    /// Mirrors `len(FGHierarchyNode)`.
    pub fn len(&self) -> usize {
        1 + self
            .children
            .iter()
            .map(FGHierarchyNode::len)
            .sum::<usize>()
    }

    /// Always false: a node counts itself.
    pub fn is_empty(&self) -> bool {
        false
    }
}

/// Parameters for `build_func_group_hierarchy_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BuildFuncGroupHierarchyParams {
    /// The file to read the hierarchy from; it wins over `data`.
    pub file_nm: Option<PathBuf>,
    /// The hierarchy's text; without it or a file, RDKit's own hierarchy.
    pub data: Option<String>,
}

/// RDKit's functional group hierarchy: its top-level groups, each with the
/// groups under it. Mirrors `FunctionalGroups.BuildFuncGroupHierarchy()`.
pub fn build_func_group_hierarchy() -> Result<Vec<FGHierarchyNode>> {
    build_func_group_hierarchy_with_params(&BuildFuncGroupHierarchyParams::default())
}

/// The hierarchy a file or text defines: one group per line, as its label,
/// SMARTS, name, and optionally a removal reaction SMARTS, separated by
/// tabs, with `//` comments. A malformed line is `Error::Value` with
/// Python's message (Python's `FuncGroupFileParseError` is a `ValueError`),
/// and a file that cannot be read is `Error::Os`. Mirrors
/// `FunctionalGroups.BuildFuncGroupHierarchy(fileNm, data)`.
pub fn build_func_group_hierarchy_with_params(
    params: &BuildFuncGroupHierarchyParams,
) -> Result<Vec<FGHierarchyNode>> {
    let text = match (&params.file_nm, &params.data) {
        (Some(file), _) => read_text(file)?,
        (None, Some(data)) if !data.is_empty() => data.clone(),
        _ => DEFAULT_HIERARCHY.to_owned(),
    };
    parse_hierarchy(&text)
}

fn read_text(file: &Path) -> Result<String> {
    let text = std::fs::read_to_string(file).map_err(|e| Error::Os {
        message: format!("{}: {e}", file.display()),
    })?;
    // Python reads text files with universal newlines.
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// The fields of a line, split at runs of tabs, as `re.split('\t+', line)`.
fn split_tabs(line: &str) -> Vec<&str> {
    let mut parts = Vec::new();
    let mut start = 0;
    let mut in_tab = false;
    for (i, c) in line.char_indices() {
        if c == '\t' {
            if !in_tab {
                parts.push(&line[start..i]);
                in_tab = true;
            }
            start = i + 1;
        } else {
            in_tab = false;
        }
    }
    parts.push(&line[start..]);
    parts
}

fn parse_error(message: String) -> Error {
    Error::Value { message }
}

fn parse_hierarchy(text: &str) -> Result<Vec<FGHierarchyNode>> {
    let mut labels = std::collections::HashSet::new();
    let mut res: Vec<FGHierarchyNode> = Vec::new();
    for (i, line) in splitlines(text).into_iter().enumerate() {
        let line_no = i + 1;
        let line = line.trim().split("//").next().unwrap_or_default();
        if line.is_empty() {
            continue;
        }
        let split_l = split_tabs(line);
        if split_l.len() < 3 {
            return Err(parse_error(format!(
                "Input line {line_no} ({}) is not long enough.",
                str_repr(line)
            )));
        }
        let label = split_l[0].trim();
        if labels.contains(label) {
            return Err(parse_error(format!("Duplicate label on line {line_no}.")));
        }
        let label_hierarchy: Vec<&str> = label.split('.').collect();
        let parent = if label_hierarchy.len() > 1 {
            for i in 0..label_hierarchy.len() - 1 {
                let tmp = label_hierarchy[..=i].join(".");
                if !labels.contains(tmp.as_str()) {
                    return Err(parse_error(format!(
                        "Hierarchy member {tmp} (line {line_no}) not found."
                    )));
                }
            }
            Some(label_hierarchy[..label_hierarchy.len() - 1].join("."))
        } else {
            None
        };
        let smarts = split_l[1];
        let patt = mol_from_smarts(smarts).map_err(|_| {
            parse_error(format!(
                "Smarts \"{smarts}\" (line {line_no}) could not be parsed."
            ))
        })?;
        let name = split_l[2].trim();
        let rxn_smarts = split_l.get(3).copied().unwrap_or_default();
        let node = FGHierarchyNode::new_with_params(
            name,
            patt,
            &FGHierarchyNodeParams {
                smarts: smarts.to_owned(),
                label: label.to_owned(),
                rxn_smarts: rxn_smarts.to_owned(),
                parent: parent.clone(),
            },
        );
        labels.insert(label.to_owned());
        match parent {
            Some(parent) => find_mut(&mut res, &parent)
                .expect("a parent's label leads to it")
                .children
                .push(node),
            None => res.push(node),
        }
    }
    Ok(res)
}

/// The node with this label, which its parents' labels lead to.
fn find_mut<'a>(nodes: &'a mut [FGHierarchyNode], label: &str) -> Option<&'a mut FGHierarchyNode> {
    let parts: Vec<&str> = label.split('.').collect();
    let mut nodes = nodes;
    for i in 0..parts.len() {
        let prefix = parts[..=i].join(".");
        let idx = nodes.iter().position(|n| n.label == prefix)?;
        if i + 1 == parts.len() {
            return Some(&mut nodes[idx]);
        }
        nodes = &mut nodes[idx].children;
    }
    None
}

/// Counts the group's matches (by their first atoms) into `res` from
/// `idx`, then its children's when it matches, and returns the next index.
fn set_node_bits(
    mol: &Mol,
    node: &FGHierarchyNode,
    res: &mut [u32],
    mut idx: usize,
) -> Result<usize> {
    let ms = mol.get_substruct_matches(&node.pattern)?;
    let mut seen = std::collections::HashSet::new();
    let count = ms.iter().filter(|m| seen.insert(m[0])).count() as u32;
    if count > 0 {
        res[idx] = count;
        idx += 1;
        for child in &node.children {
            idx = set_node_bits(mol, child, res, idx)?;
        }
    } else {
        idx += node.len();
    }
    Ok(idx)
}

/// How many of each group the molecule has, one count per node of the
/// hierarchy, in depth-first order; a group whose parent the molecule
/// lacks counts 0. Mirrors `FunctionalGroups.CreateMolFingerprint`.
pub fn create_mol_fingerprint(mol: &Mol, hierarchy: &[FGHierarchyNode]) -> Result<Vec<u32>> {
    let tot_l = hierarchy.iter().map(FGHierarchyNode::len).sum();
    let mut res = vec![0; tot_l];
    let mut idx = 0;
    for entry in hierarchy {
        idx = set_node_bits(mol, entry, &mut res, idx)?;
    }
    Ok(res)
}
