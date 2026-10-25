// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.Scaffolds.rdScaffoldNetwork`: networks of the
//! scaffolds of molecules, linked by the fragmentations, generalizations,
//! and removals of attachment points that lead from one to the next.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smiles;
//! use rdkit_cxx::chem::scaffolds::rd_scaffold_network::{
//!     ScaffoldNetworkParams, create_scaffold_network,
//! };
//!
//! let mols = vec![mol_from_smiles("c1ccccc1CC1NC(=O)CCC1")?];
//! let network = create_scaffold_network(&mols, &ScaffoldNetworkParams::default())?;
//! for edge in &network.edges {
//!     println!("{} -> {}", network.nodes[edge.begin_idx], network.nodes[edge.end_idx]);
//! }
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

use std::fmt;

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::error::expect_ok;
use crate::ffi::scaffold_network::ffi;
use crate::{Error, Result};

rdkit_enum! {
    /// How an edge's end follows from its beginning. Mirrors
    /// `rdScaffoldNetwork.EdgeType`.
    EdgeType: ffi::EdgeType {
        Fragment = Fragment,
        Generic = Generic,
        GenericBond = GenericBond,
        RemoveAttachment = RemoveAttachment,
        Initialize = Initialize,
    }
}

impl fmt::Display for EdgeType {
    /// The name RDKit's `operator<<` writes.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            EdgeType::Fragment => "Fragment",
            EdgeType::Generic => "Generic",
            EdgeType::GenericBond => "GenericBond",
            EdgeType::RemoveAttachment => "RemoveAttachment",
            EdgeType::Initialize => "Initialize",
        })
    }
}

/// An edge between two nodes of a network, by their indices in its
/// `nodes`. Mirrors `rdScaffoldNetwork.NetworkEdge`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct NetworkEdge {
    /// The node the edge leaves. Mirrors `NetworkEdge.beginIdx`.
    pub begin_idx: usize,
    /// The node the edge reaches. Mirrors `NetworkEdge.endIdx`.
    pub end_idx: usize,
    /// Mirrors `NetworkEdge.type`.
    pub type_: EdgeType,
}

/// `NetworkEdge( 0->1, type:Fragment )`, as RDKit writes edges. Mirrors
/// `NetworkEdge.__str__`.
impl fmt::Display for NetworkEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "NetworkEdge( {}->{}, type:{} )",
            self.begin_idx, self.end_idx, self.type_
        )
    }
}

/// The parameters of a network: which scaffolds to include, how to flatten
/// molecules before fragmenting them, and the bond-breaker reactions that
/// fragment them. Mirrors `rdScaffoldNetwork.ScaffoldNetworkParams`.
///
/// The bond breakers are set when the parameters are made (`new`, or
/// `brics_scaffold_params`), as in Python, which does not expose them.
/// Clones share them. `ScaffoldNetworkParams` is `Send` but not `Sync`, as
/// RDKit's reactions are (CONVENTIONS.md S4).
pub struct ScaffoldNetworkParams {
    /// Includes scaffolds with every atom made a dummy. Mirrors
    /// `ScaffoldNetworkParams.includeGenericScaffolds`.
    pub include_generic_scaffolds: bool,
    /// Includes scaffolds with every bond made single. Mirrors
    /// `ScaffoldNetworkParams.includeGenericBondScaffolds`.
    pub include_generic_bond_scaffolds: bool,
    /// Includes scaffolds with their attachment points removed. Mirrors
    /// `ScaffoldNetworkParams.includeScaffoldsWithoutAttachments`.
    pub include_scaffolds_without_attachments: bool,
    /// Includes scaffolds with their attachment points. Mirrors
    /// `ScaffoldNetworkParams.includeScaffoldsWithAttachments`.
    pub include_scaffolds_with_attachments: bool,
    /// Appends the molecules' names to their nodes. Mirrors
    /// `ScaffoldNetworkParams.includeNames`.
    pub include_names: bool,
    /// Keeps only the first fragment of each bond break. Mirrors
    /// `ScaffoldNetworkParams.keepOnlyFirstFragment`.
    pub keep_only_first_fragment: bool,
    /// Flattens and prunes molecules before fragmenting them. Mirrors
    /// `ScaffoldNetworkParams.pruneBeforeFragmenting`.
    pub prune_before_fragmenting: bool,
    /// Removes isotopes when flattening. Mirrors
    /// `ScaffoldNetworkParams.flattenIsotopes`.
    pub flatten_isotopes: bool,
    /// Removes chirality and bond stereo when flattening. Mirrors
    /// `ScaffoldNetworkParams.flattenChirality`.
    pub flatten_chirality: bool,
    /// Keeps only the largest fragment when flattening. Mirrors
    /// `ScaffoldNetworkParams.flattenKeepLargest`.
    pub flatten_keep_largest: bool,
    /// Counts the molecules each scaffold was found in. Mirrors
    /// `ScaffoldNetworkParams.collectMolCounts`.
    pub collect_mol_counts: bool,
    /// Never null; only its reactions are read.
    bond_breakers: UniquePtr<ffi::ScaffoldNetworkParams>,
}

impl ScaffoldNetworkParams {
    fn from_ffi(ptr: UniquePtr<ffi::ScaffoldNetworkParams>) -> ScaffoldNetworkParams {
        let raw = ptr.as_ref().expect("RDKit returned null parameters");
        let flags = expect_ok(ffi::scaffold_network_params_flags(raw));
        ScaffoldNetworkParams {
            include_generic_scaffolds: flags.include_generic_scaffolds,
            include_generic_bond_scaffolds: flags.include_generic_bond_scaffolds,
            include_scaffolds_without_attachments: flags.include_scaffolds_without_attachments,
            include_scaffolds_with_attachments: flags.include_scaffolds_with_attachments,
            include_names: flags.include_names,
            keep_only_first_fragment: flags.keep_only_first_fragment,
            prune_before_fragmenting: flags.prune_before_fragmenting,
            flatten_isotopes: flags.flatten_isotopes,
            flatten_chirality: flags.flatten_chirality,
            flatten_keep_largest: flags.flatten_keep_largest,
            collect_mol_counts: flags.collect_mol_counts,
            bond_breakers: ptr,
        }
    }

    fn flags(&self) -> ffi::ScaffoldNetworkFlags {
        ffi::ScaffoldNetworkFlags {
            include_generic_scaffolds: self.include_generic_scaffolds,
            include_generic_bond_scaffolds: self.include_generic_bond_scaffolds,
            include_scaffolds_without_attachments: self.include_scaffolds_without_attachments,
            include_scaffolds_with_attachments: self.include_scaffolds_with_attachments,
            include_names: self.include_names,
            keep_only_first_fragment: self.keep_only_first_fragment,
            prune_before_fragmenting: self.prune_before_fragmenting,
            flatten_isotopes: self.flatten_isotopes,
            flatten_chirality: self.flatten_chirality,
            flatten_keep_largest: self.flatten_keep_largest,
            collect_mol_counts: self.collect_mol_counts,
        }
    }

    fn raw(&self) -> &ffi::ScaffoldNetworkParams {
        self.bond_breakers.as_ref().unwrap()
    }

    /// Parameters whose bond breakers are the given reaction SMARTS, each
    /// with one reactant template and two product templates, and the
    /// default flags. SMARTS that do not parse are `Error::Parse`, and
    /// other reactions `Error::Value`. Mirrors
    /// `ScaffoldNetworkParams(bondBreakerSmartsList)`.
    pub fn new(bond_breaker_smarts_list: &[&str]) -> Result<ScaffoldNetworkParams> {
        Ok(ScaffoldNetworkParams::from_ffi(
            ffi::scaffold_network_params_from_smarts(bond_breaker_smarts_list)?,
        ))
    }

    /// The number of bond-breaker reactions.
    pub fn num_bond_breakers(&self) -> usize {
        expect_ok(ffi::scaffold_network_params_num_bond_breakers(self.raw()))
    }
}

/// RDKit's defaults, whose one bond breaker cuts acyclic bonds from ring
/// atoms. Mirrors `ScaffoldNetworkParams()`.
impl Default for ScaffoldNetworkParams {
    fn default() -> ScaffoldNetworkParams {
        ScaffoldNetworkParams::from_ffi(expect_ok(ffi::scaffold_network_params_new()))
    }
}

/// A copy that shares the bond breakers.
impl Clone for ScaffoldNetworkParams {
    fn clone(&self) -> ScaffoldNetworkParams {
        ScaffoldNetworkParams {
            bond_breakers: expect_ok(ffi::scaffold_network_params_clone(self.raw())),
            ..*self
        }
    }
}

impl fmt::Debug for ScaffoldNetworkParams {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ScaffoldNetworkParams")
            .field("include_generic_scaffolds", &self.include_generic_scaffolds)
            .field(
                "include_generic_bond_scaffolds",
                &self.include_generic_bond_scaffolds,
            )
            .field(
                "include_scaffolds_without_attachments",
                &self.include_scaffolds_without_attachments,
            )
            .field(
                "include_scaffolds_with_attachments",
                &self.include_scaffolds_with_attachments,
            )
            .field("include_names", &self.include_names)
            .field("keep_only_first_fragment", &self.keep_only_first_fragment)
            .field("prune_before_fragmenting", &self.prune_before_fragmenting)
            .field("flatten_isotopes", &self.flatten_isotopes)
            .field("flatten_chirality", &self.flatten_chirality)
            .field("flatten_keep_largest", &self.flatten_keep_largest)
            .field("collect_mol_counts", &self.collect_mol_counts)
            .field("num_bond_breakers", &self.num_bond_breakers())
            .finish()
    }
}

/// A scaffold network: the scaffolds' SMILES as nodes, and the edges
/// between them. Mirrors `rdScaffoldNetwork.ScaffoldNetwork`, whose fields
/// Python reads only.
///
/// `update_scaffold_network` needs one count per node (`Error::Value`
/// otherwise); RDKit's networks always have one.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScaffoldNetwork {
    /// The scaffolds' SMILES, with the molecules' names when
    /// `include_names` was set. Mirrors `ScaffoldNetwork.nodes`.
    pub nodes: Vec<String>,
    /// How often each node was reached while the network was built.
    /// Mirrors `ScaffoldNetwork.counts`.
    pub counts: Vec<u32>,
    /// How many molecules reached each node, for the nodes up to the last
    /// one counted (`collect_mol_counts`). Mirrors `ScaffoldNetwork.molCounts`.
    pub mol_counts: Vec<u32>,
    /// Mirrors `ScaffoldNetwork.edges`.
    pub edges: Vec<NetworkEdge>,
}

impl ScaffoldNetwork {
    fn to_ffi(&self) -> UniquePtr<ffi::ScaffoldNetwork> {
        let edges: Vec<ffi::NetworkEdgeData> = self
            .edges
            .iter()
            .map(|edge| ffi::NetworkEdgeData {
                begin_idx: edge.begin_idx,
                end_idx: edge.end_idx,
                edge_type: edge.type_.to_ffi(),
            })
            .collect();
        expect_ok(ffi::scaffold_network_new(
            &self.nodes,
            &self.counts,
            &self.mol_counts,
            &edges,
        ))
    }

    /// Invalid UTF-8 in a node (from a molecule's name) is `Error::Value`.
    fn from_ffi(network: &ffi::ScaffoldNetwork) -> Result<ScaffoldNetwork> {
        let nodes = expect_ok(ffi::scaffold_network_nodes(network))
            .iter()
            .map(|node| {
                node.to_str()
                    .map(str::to_owned)
                    .map_err(|err| Error::Value {
                        message: format!("invalid UTF-8 in a string from RDKit: {err}"),
                    })
            })
            .collect::<Result<Vec<String>>>()?;
        Ok(ScaffoldNetwork {
            nodes,
            counts: expect_ok(ffi::scaffold_network_counts(network)),
            mol_counts: expect_ok(ffi::scaffold_network_mol_counts(network)),
            edges: expect_ok(ffi::scaffold_network_edges(network))
                .into_iter()
                .map(|edge| NetworkEdge {
                    begin_idx: edge.begin_idx,
                    end_idx: edge.end_idx,
                    type_: EdgeType::from_ffi(edge.edge_type),
                })
                .collect(),
        })
    }

    /// The network a pickle holds, from `to_binary` here or in Python. A
    /// network without one count per node, or with an edge type RDKit
    /// lacks, is `Error::Value`. Mirrors `ScaffoldNetwork(pkl)`.
    ///
    /// # Safety
    ///
    /// `pkl` must be bytes `to_binary` (this crate's or RDKit's) returned.
    /// RDKit reads it as a Boost.Serialization archive, which it does not
    /// validate (CONVENTIONS.md U2).
    pub unsafe fn from_pkl(pkl: &[u8]) -> Result<ScaffoldNetwork> {
        let network = ffi::scaffold_network_from_binary(pkl)?;
        ScaffoldNetwork::from_ffi(network.as_ref().expect("RDKit returned a null network"))
    }

    /// The network as RDKit pickles it, a Boost.Serialization text archive.
    /// `Error::Value` if RDKit was built without Boost.Serialization.
    pub fn to_binary(&self) -> Result<Vec<u8>> {
        Ok(ffi::scaffold_network_to_binary(&self.to_ffi())?)
    }
}

/// A new network for the molecules. Mirrors
/// `rdScaffoldNetwork.CreateScaffoldNetwork`.
pub fn create_scaffold_network(
    mols: &[Mol],
    params: &ScaffoldNetworkParams,
) -> Result<ScaffoldNetwork> {
    let mut network = ScaffoldNetwork::default();
    update_scaffold_network(mols, &mut network, params)?;
    Ok(network)
}

/// Adds the molecules to the network. Parameters that include neither
/// scaffolds with attachment points nor scaffolds without them are
/// `Error::Value`, as is a network without one count per node (which
/// RDKit's update would write past). A molecule that fails (a fragment
/// that does not sanitize) is an error that leaves the molecules before it
/// added, as in Python. Mirrors `rdScaffoldNetwork.UpdateScaffoldNetwork`.
pub fn update_scaffold_network(
    mols: &[Mol],
    network: &mut ScaffoldNetwork,
    params: &ScaffoldNetworkParams,
) -> Result<()> {
    let mut updated = network.to_ffi();
    let result = with_borrowed_mols_shared(mols.iter().map(|mol| ("", mol)), |list| {
        ffi::update_scaffold_network(list, updated.pin_mut(), &params.flags(), params.raw())
    });
    *network = ScaffoldNetwork::from_ffi(&updated)?;
    Ok(result?)
}

/// Parameters that fragment molecules with the BRICS rules, keeping both
/// fragments of each break. Mirrors `rdScaffoldNetwork.BRICSScaffoldParams`.
pub fn brics_scaffold_params() -> ScaffoldNetworkParams {
    ScaffoldNetworkParams::from_ffi(expect_ok(ffi::brics_scaffold_params()))
}
