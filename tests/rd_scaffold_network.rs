// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::chem::scaffolds::rd_scaffold_network::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdScaffoldNetwork.

fn mols(smiles: &[&str]) -> Vec<Mol> {
    smiles.iter().map(|s| mol_from_smiles(s).unwrap()).collect()
}

fn edge_strings(network: &ScaffoldNetwork) -> Vec<String> {
    network.edges.iter().map(ToString::to_string).collect()
}

const TWO: [&str; 2] = ["c1ccccc1CC1NC(=O)CCC1", "c1cccnc1CC1NC(=O)CCC1"];

const TWO_NODES: [&str; 12] = [
    "O=C1CCCC(Cc2ccccc2)N1",
    "*c1ccccc1",
    "**1:*:*:*:*:*:1",
    "c1ccccc1",
    "*1:*:*:*:*:*:1",
    "*C1CCCC(=O)N1",
    "**1****(=*)*1",
    "O=C1CCCCN1",
    "*=*1*****1",
    "O=C1CCCC(Cc2ccccn2)N1",
    "*c1ccccn1",
    "c1ccncc1",
];

const TWO_COUNTS: [u32; 12] = [1, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1];

#[test]
fn params_defaults() {
    let params = ScaffoldNetworkParams::default();
    assert!(params.include_generic_scaffolds);
    assert!(!params.include_generic_bond_scaffolds);
    assert!(params.include_scaffolds_without_attachments);
    assert!(params.include_scaffolds_with_attachments);
    assert!(!params.include_names);
    assert!(params.keep_only_first_fragment);
    assert!(params.prune_before_fragmenting);
    assert!(params.flatten_isotopes);
    assert!(params.flatten_chirality);
    assert!(params.flatten_keep_largest);
    assert!(params.collect_mol_counts);
    assert_eq!(params.num_bond_breakers(), 1);

    let brics = brics_scaffold_params();
    assert!(!brics.keep_only_first_fragment);
    assert!(brics.include_generic_scaffolds);
    assert_eq!(brics.num_bond_breakers(), 46);

    let mut changed = ScaffoldNetworkParams::default();
    changed.include_names = true;
    let copy = changed.clone();
    assert!(copy.include_names);
    assert_eq!(copy.num_bond_breakers(), 1);
    assert!(format!("{copy:?}").contains("include_names: true"));
}

#[test]
fn bond_breakers() {
    let params =
        ScaffoldNetworkParams::new(&["[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]"]).unwrap();
    assert_eq!(params.num_bond_breakers(), 1);
    assert_eq!(
        ScaffoldNetworkParams::new(&[]).unwrap().num_bond_breakers(),
        0
    );
    // RDKit's constructor refuses SMARTS that do not parse and reactions
    // that are not one reactant to two products.
    assert!(matches!(
        ScaffoldNetworkParams::new(&["foo"]),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        ScaffoldNetworkParams::new(&["[C:1]>>[X"]),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        ScaffoldNetworkParams::new(&["[C:1]>>[C:1]"]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        ScaffoldNetworkParams::new(&["[C:1].[N:2]>>[C:1].[N:2]"]),
        Err(Error::Value { .. })
    ));

    // No bond breakers: only the flattened molecules.
    let network =
        create_scaffold_network(&mols(&TWO), &ScaffoldNetworkParams::new(&[]).unwrap()).unwrap();
    assert_eq!(
        network.nodes,
        ["O=C1CCCC(Cc2ccccc2)N1", "O=C1CCCC(Cc2ccccn2)N1"]
    );
}

#[test]
fn create_and_update() {
    let params = ScaffoldNetworkParams::default();
    let network = create_scaffold_network(&mols(&TWO), &params).unwrap();
    assert_eq!(network.nodes, TWO_NODES);
    assert_eq!(network.counts, TWO_COUNTS);
    assert_eq!(network.mol_counts, TWO_COUNTS);
    assert_eq!(
        edge_strings(&network),
        [
            "NetworkEdge( 0->1, type:Fragment )",
            "NetworkEdge( 1->2, type:Generic )",
            "NetworkEdge( 1->3, type:RemoveAttachment )",
            "NetworkEdge( 3->4, type:Generic )",
            "NetworkEdge( 0->5, type:Fragment )",
            "NetworkEdge( 5->6, type:Generic )",
            "NetworkEdge( 5->7, type:RemoveAttachment )",
            "NetworkEdge( 7->8, type:Generic )",
            "NetworkEdge( 9->10, type:Fragment )",
            "NetworkEdge( 10->2, type:Generic )",
            "NetworkEdge( 10->11, type:RemoveAttachment )",
            "NetworkEdge( 11->4, type:Generic )",
            "NetworkEdge( 9->5, type:Fragment )",
        ]
    );
    assert_eq!(
        network.edges[2],
        NetworkEdge {
            begin_idx: 1,
            end_idx: 3,
            type_: EdgeType::RemoveAttachment
        }
    );

    let mut updated = create_scaffold_network(&mols(&TWO[..1]), &params).unwrap();
    update_scaffold_network(&mols(&TWO[1..]), &mut updated, &params).unwrap();
    assert_eq!(updated, network);

    // Neither kind of scaffold is RDKit's ValueErrorException.
    let mut neither = ScaffoldNetworkParams::default();
    neither.include_scaffolds_with_attachments = false;
    neither.include_scaffolds_without_attachments = false;
    assert!(matches!(
        create_scaffold_network(&mols(&TWO), &neither),
        Err(Error::Value { .. })
    ));
}

#[test]
fn flags() {
    let mut mol = mol_from_smiles("c1ccccc1C1CC1").unwrap();
    mol.set_prop("_Name", "cpd1");
    let mut params = ScaffoldNetworkParams::default();
    params.include_names = true;
    let network = create_scaffold_network(std::slice::from_ref(&mol), &params).unwrap();
    assert_eq!(
        network.nodes[..3],
        [
            "c1ccc(C2CC2)cc1 cpd1",
            "c1ccc(C2CC2)cc1",
            "*1**1*1:*:*:*:*:*:1"
        ]
    );
    assert_eq!(
        network.edges[0].to_string(),
        "NetworkEdge( 0->1, type:Initialize )"
    );
    assert_eq!(network.counts, [1; 11]);

    let mut params = ScaffoldNetworkParams::default();
    params.include_generic_scaffolds = false;
    params.collect_mol_counts = false;
    let network = create_scaffold_network(&mols(&["c1ccccc1C1CC1"]), &params).unwrap();
    assert_eq!(
        network.nodes,
        [
            "c1ccc(C2CC2)cc1",
            "*c1ccccc1",
            "c1ccccc1",
            "*C1CC1",
            "C1CC1"
        ]
    );
    assert_eq!(network.counts, [1; 5]);
    assert!(network.mol_counts.is_empty());

    let network =
        create_scaffold_network(&mols(&["c1ccccc1CC1NC(=O)CCC1"]), &brics_scaffold_params())
            .unwrap();
    assert_eq!(network.nodes.len(), 19);
    assert_eq!(network.nodes[1], "[8*]Cc1ccccc1");
    assert_eq!(network.edges.len(), 23);
}

#[test]
fn invalid_utf8_names() {
    // A name read from a Latin-1 file reaches a node.
    let mol = mol_from_mol_file(common::latin1_mol_file()).unwrap();
    let mut params = ScaffoldNetworkParams::default();
    params.include_names = true;
    assert!(matches!(
        create_scaffold_network(&[mol], &params),
        Err(Error::Value { .. })
    ));
}

#[test]
fn pickles() {
    let mut mol = mol_from_smiles("c1ccccc1C1CC1").unwrap();
    mol.set_prop("_Name", "cpd1");
    let mut params = ScaffoldNetworkParams::default();
    params.include_names = true;
    let network = create_scaffold_network(&[mol], &params).unwrap();
    // Python's pickle of the same network.
    let python = "22 serialization::archive 20 0 1 0 0 11 0 20 c1ccc(C2CC2)cc1 cpd1 15 \
                  c1ccc(C2CC2)cc1 19 *1**1*1:*:*:*:*:*:1 9 *c1ccccc1 15 **1:*:*:*:*:*:1 8 \
                  c1ccccc1 14 *1:*:*:*:*:*:1 6 *C1CC1 6 **1**1 5 C1CC1 5 *1**1 11 0 1 1 1 1 \
                  1 1 1 1 1 1 1 11 0 1 1 1 1 1 1 1 1 1 1 1 0 0 10 0 0 0 0 1 5 1 2 2 1 3 1 3 4 \
                  2 3 5 4 5 6 2 1 7 1 7 8 2 7 9 4 9 10 2";
    // Boost writes its archive version after the signature, so only the
    // signature is compared.
    let pkl = network.to_binary().unwrap();
    assert!(pkl.starts_with(b"22 serialization::archive "));
    // SAFETY: RDKit wrote both.
    unsafe {
        assert_eq!(
            ScaffoldNetwork::from_pkl(python.as_bytes()).unwrap(),
            network
        );
        assert_eq!(ScaffoldNetwork::from_pkl(&pkl).unwrap(), network);
        let empty = ScaffoldNetwork::default().to_binary().unwrap();
        assert_eq!(
            ScaffoldNetwork::from_pkl(&empty).unwrap(),
            ScaffoldNetwork::default()
        );
    }
}

// RDKit's update indexes the counts by node unchecked: Python segfaults
// updating a network whose pickle has a node and no counts.
#[test]
fn networks_need_a_count_per_node() {
    let params = ScaffoldNetworkParams::default();
    let mut network = ScaffoldNetwork {
        nodes: vec!["O=C1CCCC(Cc2ccccc2)N1".to_owned()],
        ..ScaffoldNetwork::default()
    };
    assert!(matches!(
        update_scaffold_network(&mols(&TWO[..1]), &mut network, &params),
        Err(Error::Value { .. })
    ));
    assert_eq!(network.nodes.len(), 1);

    let unchecked =
        "22 serialization::archive 20 0 1 0 0 1 0 21 O=C1CCCC(Cc2ccccc2)N1 0 0 0 0 0 0 0 0";
    // SAFETY: the archive is well formed; only its counts are short.
    assert!(matches!(
        unsafe { ScaffoldNetwork::from_pkl(unchecked.as_bytes()) },
        Err(Error::Value { .. })
    ));
    // An edge type RDKit lacks.
    let edge = "22 serialization::archive 20 0 1 0 0 1 0 1 C 1 0 1 1 0 1 0 0 1 0 0 0 0 0 9";
    // SAFETY: as above.
    assert!(matches!(
        unsafe { ScaffoldNetwork::from_pkl(edge.as_bytes()) },
        Err(Error::Value { .. })
    ));

    // Other counts are RDKit's to resize, and edges are not read.
    let mut network = ScaffoldNetwork {
        nodes: vec!["C".to_owned()],
        counts: vec![7],
        mol_counts: vec![1, 2, 3],
        edges: vec![NetworkEdge {
            begin_idx: 40,
            end_idx: 50,
            type_: EdgeType::Generic,
        }],
    };
    update_scaffold_network(&mols(&["CC"]), &mut network, &params).unwrap();
    assert_eq!(network.nodes, ["C", "CC", ""]);
    assert_eq!(network.counts, [7, 1, 2]);
    assert_eq!(network.mol_counts, [1, 3, 4]);
    assert_eq!(
        edge_strings(&network),
        [
            "NetworkEdge( 40->50, type:Generic )",
            "NetworkEdge( 1->2, type:Initialize )",
            "NetworkEdge( 2->2, type:Generic )",
        ]
    );
}

// The fragmentation runs reactions, which match substructures under
// control_c_mutex; networks built on several threads at once agree.
#[test]
fn threads() {
    let params = ScaffoldNetworkParams::default();
    let handles: Vec<_> = (0..4)
        .map(|_| {
            let params = params.clone();
            std::thread::spawn(move || {
                (0..5)
                    .map(|_| create_scaffold_network(&mols(&TWO), &params).unwrap())
                    .collect::<Vec<_>>()
            })
        })
        .collect();
    for handle in handles {
        for network in handle.join().unwrap() {
            assert_eq!(network.nodes, TWO_NODES);
            assert_eq!(network.counts, TWO_COUNTS);
        }
    }
}
