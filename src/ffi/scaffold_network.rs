// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::ScaffoldNetwork::EdgeType`; cxx checks the values.
    #[repr(i32)]
    enum EdgeType {
        Fragment = 1,
        Generic = 2,
        GenericBond = 3,
        RemoveAttachment = 4,
        Initialize = 5,
    }

    /// `ScaffoldNetworkParams`' flags.
    struct ScaffoldNetworkFlags {
        include_generic_scaffolds: bool,
        include_generic_bond_scaffolds: bool,
        include_scaffolds_without_attachments: bool,
        include_scaffolds_with_attachments: bool,
        keep_only_first_fragment: bool,
        prune_before_fragmenting: bool,
        flatten_isotopes: bool,
        flatten_chirality: bool,
        flatten_keep_largest: bool,
        collect_mol_counts: bool,
    }

    /// `NetworkEdge`.
    struct NetworkEdgeData {
        begin_idx: usize,
        end_idx: usize,
        edge_type: EdgeType,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/scaffold_network.h");

        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type EdgeType;
        type ScaffoldNetworkParams;
        type ScaffoldNetwork;

        fn scaffold_network_params_new() -> Result<UniquePtr<ScaffoldNetworkParams>>;
        fn scaffold_network_params_from_smarts(
            smarts: &[&str],
        ) -> Result<UniquePtr<ScaffoldNetworkParams>>;
        fn brics_scaffold_params() -> Result<UniquePtr<ScaffoldNetworkParams>>;
        fn scaffold_network_params_clone(
            params: &ScaffoldNetworkParams,
        ) -> Result<UniquePtr<ScaffoldNetworkParams>>;
        fn scaffold_network_params_num_bond_breakers(
            params: &ScaffoldNetworkParams,
        ) -> Result<usize>;
        fn scaffold_network_params_flags(
            params: &ScaffoldNetworkParams,
        ) -> Result<ScaffoldNetworkFlags>;

        fn scaffold_network_new(
            nodes: &[String],
            counts: &[u32],
            mol_counts: &[u32],
            edges: &[NetworkEdgeData],
        ) -> Result<UniquePtr<ScaffoldNetwork>>;
        fn scaffold_network_nodes(network: &ScaffoldNetwork) -> Result<&CxxVector<CxxString>>;
        fn scaffold_network_counts(network: &ScaffoldNetwork) -> Result<Vec<u32>>;
        fn scaffold_network_mol_counts(network: &ScaffoldNetwork) -> Result<Vec<u32>>;
        fn scaffold_network_edges(network: &ScaffoldNetwork) -> Result<Vec<NetworkEdgeData>>;
        fn update_scaffold_network(
            mols: Pin<&mut BorrowedMols>,
            network: Pin<&mut ScaffoldNetwork>,
            flags: &ScaffoldNetworkFlags,
            params: &ScaffoldNetworkParams,
        ) -> Result<()>;
        fn scaffold_network_to_binary(network: &ScaffoldNetwork) -> Result<Vec<u8>>;
        fn scaffold_network_from_binary(data: &[u8]) -> Result<UniquePtr<ScaffoldNetwork>>;
    }
}

// SAFETY: parameters own their flags and share their bond breakers, which
// RDKit parses and initializes when it makes them, through shared_ptrs
// (atomic counts); copies on other threads only run the reactions, which
// RDKit's const runReactant reads, under control_c_mutex.
unsafe impl Send for ffi::ScaffoldNetworkParams {}
