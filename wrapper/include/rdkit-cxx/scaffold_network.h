/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdScaffoldNetwork: scaffold networks built by fragmenting molecules with
// bond-breaker reactions.
//
// RDKit's update indexes a network's counts by node without checking that
// there is a count per node, so a network the wrappers did not build (a
// Rust value, a pickle) is checked first.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>

#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using ScaffoldNetwork = RDKit::ScaffoldNetwork::ScaffoldNetwork;
using ScaffoldNetworkParams = RDKit::ScaffoldNetwork::ScaffoldNetworkParams;
using EdgeType = RDKit::ScaffoldNetwork::EdgeType;

// Shared structs, defined by the bridge (src/ffi/scaffold_network.rs).
struct ScaffoldNetworkFlags;
struct NetworkEdgeData;

// Parameters with RDKit's default bond breaker, and with the given reaction
// SMARTS, each with one reactant and two products (ValueErrorException
// otherwise). Only their reactions are read: the flags come with each call.
std::unique_ptr<ScaffoldNetworkParams> scaffold_network_params_new();
std::unique_ptr<ScaffoldNetworkParams> scaffold_network_params_from_smarts(
    rust::Slice<const rust::Str> smarts);
std::unique_ptr<ScaffoldNetworkParams> brics_scaffold_params();
// A copy, which shares the initialized reactions.
std::unique_ptr<ScaffoldNetworkParams> scaffold_network_params_clone(
    const ScaffoldNetworkParams &params);
std::size_t scaffold_network_params_num_bond_breakers(
    const ScaffoldNetworkParams &params);
ScaffoldNetworkFlags scaffold_network_params_flags(
    const ScaffoldNetworkParams &params);

std::unique_ptr<ScaffoldNetwork> scaffold_network_new(
    rust::Slice<const rust::String> nodes,
    rust::Slice<const std::uint32_t> counts,
    rust::Slice<const std::uint32_t> mol_counts,
    rust::Slice<const NetworkEdgeData> edges);
const std::vector<std::string> &scaffold_network_nodes(
    const ScaffoldNetwork &network);
rust::Vec<std::uint32_t> scaffold_network_counts(
    const ScaffoldNetwork &network);
rust::Vec<std::uint32_t> scaffold_network_mol_counts(
    const ScaffoldNetwork &network);
rust::Vec<NetworkEdgeData> scaffold_network_edges(
    const ScaffoldNetwork &network);

// Python's UpdateScaffoldNetwork, which adds the molecules to the network.
// A molecule that fails partway leaves the ones before it added.
void update_scaffold_network(BorrowedMols &mols, ScaffoldNetwork &network,
                             const ScaffoldNetworkFlags &flags,
                             const ScaffoldNetworkParams &params);

// The network's Boost.Serialization text archive, and a network read from
// one, checked as update_scaffold_network checks networks.
rust::Vec<std::uint8_t> scaffold_network_to_binary(
    const ScaffoldNetwork &network);
std::unique_ptr<ScaffoldNetwork> scaffold_network_from_binary(
    rust::Slice<const std::uint8_t> data);

}  // namespace rdkit_cxx
