// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/scaffold_network.h"

#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/scaffold_network.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

// RDKit's addEntryIfMissing increments a node's count at the node's index
// without checking the counts' size, so a network with fewer counts than
// nodes is written past the end of its counts. Its molCounts are resized as
// needed, and its edges are not read.
void check_network(const ScaffoldNetwork &network) {
  if (network.counts.size() != network.nodes.size()) {
    throw ValueErrorException("a scaffold network needs one count per node");
  }
}

bool valid_edge_type(EdgeType type) {
  switch (type) {
    case EdgeType::Fragment:
    case EdgeType::Generic:
    case EdgeType::GenericBond:
    case EdgeType::RemoveAttachment:
    case EdgeType::Initialize:
      return true;
  }
  return false;
}

ScaffoldNetworkParams with_flags(const ScaffoldNetworkParams &params,
                                 const ScaffoldNetworkFlags &flags) {
  ScaffoldNetworkParams res(params);
  res.includeGenericScaffolds = flags.include_generic_scaffolds;
  res.includeGenericBondScaffolds = flags.include_generic_bond_scaffolds;
  res.includeScaffoldsWithoutAttachments =
      flags.include_scaffolds_without_attachments;
  res.includeScaffoldsWithAttachments =
      flags.include_scaffolds_with_attachments;
  res.includeNames = flags.include_names;
  res.keepOnlyFirstFragment = flags.keep_only_first_fragment;
  res.pruneBeforeFragmenting = flags.prune_before_fragmenting;
  res.flattenIsotopes = flags.flatten_isotopes;
  res.flattenChirality = flags.flatten_chirality;
  res.flattenKeepLargest = flags.flatten_keep_largest;
  res.collectMolCounts = flags.collect_mol_counts;
  return res;
}

}  // namespace

std::unique_ptr<ScaffoldNetworkParams> scaffold_network_params_new() {
  return std::make_unique<ScaffoldNetworkParams>();
}

// RDKit's constructor throws for SMARTS that do not parse (a null
// reaction) and for reactions without one reactant and two products.
std::unique_ptr<ScaffoldNetworkParams> scaffold_network_params_from_smarts(
    rust::Slice<const rust::Str> smarts) {
  std::vector<std::string> list;
  for (auto s : smarts) {
    list.emplace_back(s);
    // ScaffoldNetworkParams parses each as reaction SMARTS (file_checks.h).
    check_reaction_cx(list.back());
  }
  return std::make_unique<ScaffoldNetworkParams>(list);
}

std::unique_ptr<ScaffoldNetworkParams> brics_scaffold_params() {
  return std::make_unique<ScaffoldNetworkParams>(
      RDKit::ScaffoldNetwork::getBRICSNetworkParams());
}

std::unique_ptr<ScaffoldNetworkParams> scaffold_network_params_clone(
    const ScaffoldNetworkParams &params) {
  return std::make_unique<ScaffoldNetworkParams>(params);
}

std::size_t scaffold_network_params_num_bond_breakers(
    const ScaffoldNetworkParams &params) {
  return params.bondBreakersRxns.size();
}

ScaffoldNetworkFlags scaffold_network_params_flags(
    const ScaffoldNetworkParams &p) {
  return ScaffoldNetworkFlags{p.includeGenericScaffolds,
                              p.includeGenericBondScaffolds,
                              p.includeScaffoldsWithoutAttachments,
                              p.includeScaffoldsWithAttachments,
                              p.includeNames,
                              p.keepOnlyFirstFragment,
                              p.pruneBeforeFragmenting,
                              p.flattenIsotopes,
                              p.flattenChirality,
                              p.flattenKeepLargest,
                              p.collectMolCounts};
}

std::unique_ptr<ScaffoldNetwork> scaffold_network_new(
    rust::Slice<const rust::String> nodes,
    rust::Slice<const std::uint32_t> counts,
    rust::Slice<const std::uint32_t> mol_counts,
    rust::Slice<const NetworkEdgeData> edges) {
  auto network = std::make_unique<ScaffoldNetwork>();
  for (const auto &node : nodes) {
    network->nodes.emplace_back(node);
  }
  network->counts.assign(counts.begin(), counts.end());
  network->molCounts.assign(mol_counts.begin(), mol_counts.end());
  for (const auto &edge : edges) {
    network->edges.emplace_back(edge.begin_idx, edge.end_idx, edge.edge_type);
  }
  return network;
}

const std::vector<std::string> &scaffold_network_nodes(
    const ScaffoldNetwork &network) {
  return network.nodes;
}

rust::Vec<std::uint32_t> scaffold_network_counts(
    const ScaffoldNetwork &network) {
  rust::Vec<std::uint32_t> res;
  for (auto count : network.counts) {
    res.push_back(count);
  }
  return res;
}

rust::Vec<std::uint32_t> scaffold_network_mol_counts(
    const ScaffoldNetwork &network) {
  rust::Vec<std::uint32_t> res;
  for (auto count : network.molCounts) {
    res.push_back(count);
  }
  return res;
}

rust::Vec<NetworkEdgeData> scaffold_network_edges(
    const ScaffoldNetwork &network) {
  rust::Vec<NetworkEdgeData> res;
  for (const auto &edge : network.edges) {
    res.push_back(NetworkEdgeData{edge.beginIdx, edge.endIdx, edge.type});
  }
  return res;
}

// Fragmenting runs the bond breakers, and flattening and pruning sanitize
// and write SMILES, all on the calling thread; the reactions match
// substructures (control_c_lock.h). RDKit computes distance matrices only on
// its own copies of the molecules, with force (MurckoDecompose).
void update_scaffold_network(BorrowedMols &mols, ScaffoldNetwork &network,
                             const ScaffoldNetworkFlags &flags,
                             const ScaffoldNetworkParams &params) {
  check_network(network);
  // RDKit removes atoms from copies of each molecule (MurckoDecompose, its
  // fragmenters) and writes them as SMILES (check_removals).
  for (const auto &mol : mols.mols()) {
    check_removals(*mol);
  }
  auto with = with_flags(params, flags);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  RDKit::ScaffoldNetwork::updateScaffoldNetwork(mols.mols(), network, with);
}

rust::Vec<std::uint8_t> scaffold_network_to_binary(
    const ScaffoldNetwork &network) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  std::stringstream oss;
  {
    boost::archive::text_oarchive oa(oss);
    oa << network;
  }
  auto text = oss.str();
  rust::Vec<std::uint8_t> res;
  res.reserve(text.size());
  for (char c : text) {
    res.push_back(static_cast<std::uint8_t>(c));
  }
  return res;
#else
  (void)network;
  throw ValueErrorException(
      "Pickling of ScaffoldNetwork instances is not enabled in this build");
#endif
}

std::unique_ptr<ScaffoldNetwork> scaffold_network_from_binary(
    rust::Slice<const std::uint8_t> data) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  auto network = std::make_unique<ScaffoldNetwork>(
      std::string(reinterpret_cast<const char *>(data.data()), data.size()));
  // RDKit's networks keep these; Rust's edge types are only RDKit's.
  check_network(*network);
  for (const auto &edge : network->edges) {
    if (!valid_edge_type(edge.type)) {
      throw ValueErrorException("a scaffold network edge has an unknown type");
    }
  }
  return network;
#else
  (void)data;
  throw ValueErrorException(
      "Pickling of ScaffoldNetwork instances is not enabled in this build");
#endif
}

}  // namespace rdkit_cxx
