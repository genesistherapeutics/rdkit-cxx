// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/rascal_mces.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <list>
#include <mutex>
#include <set>
#include <string>
#include <utility>

#include <GraphMol/MolOps.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/rascal_mces.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace Rascal = RDKit::RascalMCES;

RascalOptionsData to_data(const Rascal::RascalOptions &o) {
  return RascalOptionsData{o.similarityThreshold,
                           o.completeAromaticRings,
                           o.ringMatchesRingOnly,
                           o.completeSmallestRings,
                           o.exactConnectionsMatch,
                           o.singleLargestFrag,
                           o.minFragSize,
                           o.maxFragSeparation,
                           o.allBestMCESs,
                           o.maxBestMCESs,
                           o.timeout,
                           o.returnEmptyMCES,
                           o.maxBondMatchPairs,
                           rust::String::lossy(o.equivalentAtoms),
                           o.ignoreBondOrders,
                           o.ignoreAtomAromaticity,
                           o.minCliqueSize};
}

RascalClusterOptionsData to_data(const Rascal::RascalClusterOptions &o) {
  return RascalClusterOptionsData{o.similarityCutoff,
                                  o.a,
                                  o.b,
                                  o.minFragSize,
                                  o.minIntraClusterSim,
                                  o.clusterMergeSim,
                                  o.maxNumFrags,
                                  o.numThreads};
}

Rascal::RascalClusterOptions cluster_options_from(
    const RascalClusterOptionsData &d) {
  Rascal::RascalClusterOptions o;
  o.similarityCutoff = d.similarity_cutoff;
  o.a = d.a;
  o.b = d.b;
  o.minFragSize = d.min_frag_size;
  o.minIntraClusterSim = d.min_intra_cluster_sim;
  o.clusterMergeSim = d.cluster_merge_sim;
  o.maxNumFrags = d.max_num_frags;
  o.numThreads = d.num_threads;
  // The clustering's MCESs use this as their similarityThreshold.
  if (std::isnan(o.similarityCutoff)) {
    throw ValueErrorException("similarityCutoff is NaN");
  }
  return o;
}

rust::Vec<IndexPair> to_pairs(const std::vector<std::pair<int, int>> &pairs) {
  rust::Vec<IndexPair> res;
  res.reserve(pairs.size());
  for (const auto &[first, second] : pairs) {
    res.push_back(IndexPair{static_cast<std::uint32_t>(first),
                            static_cast<std::uint32_t>(second)});
  }
  return res;
}

rust::Vec<IndexList> to_lists(
    const std::vector<std::vector<unsigned int>> &clusters) {
  rust::Vec<IndexList> res;
  for (const auto &cluster : clusters) {
    res.push_back(IndexList{to_vec(cluster)});
  }
  return res;
}

// ---------------------------------------------------------------------------
// The clustering. RDKit's rascalCluster and rascalButinaCluster
// (RascalCluster.cpp, RascalButinaCluster.cpp) compute their proximity graph on
// std::threads that catch nothing, so an exception in one (a molecule without
// ring info, an atom without a computed valence) calls std::terminate, and each
// MCES they compute matches substructures (deltaYExchangePossible in
// RascalMCES.cpp), whose ControlCHandlers the threads would construct at once
// (control_c_lock.h). The wrapper computes the graph on the calling thread,
// holding control_c_mutex (as it does for the MCESs that splitting clusters
// computes), and runs the rest of RDKit's clustering, which those files keep to
// themselves, as ported from them below; numThreads only changes how RDKit
// would split the work, so it is not read.

// As details::ClusNode in RascalDetails.h.
struct ClusNode {
  std::shared_ptr<RascalResult> d_res;
  double d_sim;
  unsigned int d_mol1Num, d_mol2Num;
};

using ProxGraph = std::vector<std::vector<ClusNode>>;

// As details::calcMolMolSimilarity.
ClusNode calc_mol_mol_similarity(const RDKit::ROMol &mol1,
                                 const RDKit::ROMol &mol2, unsigned int i,
                                 unsigned int j,
                                 const Rascal::RascalOptions &opts,
                                 const Rascal::RascalClusterOptions &clusOpts) {
  auto res = Rascal::rascalMCES(mol1, mol2, opts);
  ClusNode cn{};
  cn.d_mol1Num = i;
  cn.d_mol2Num = j;
  if (res.empty()) {
    // tier1Sim and tier2Sim were above the threshold, but no MCES
    // was found.
    cn.d_sim = 0.0;
  } else {
    if (res.front().getBondMatches().empty()) {
      cn.d_sim = 0.0;
    } else {
      res.front().trimSmallFrags();
      res.front().largestFragsOnly(clusOpts.maxNumFrags);
      cn.d_sim = res.front().getSimilarity();
      if (cn.d_sim >= opts.similarityThreshold) {
        cn.d_res = std::make_shared<RascalResult>(res.front());
      }
    }
  }
  return cn;
}

// As details::buildProximityGraph, on the calling thread.
ProxGraph build_proximity_graph(const std::vector<const RDKit::ROMol *> &mols,
                                const Rascal::RascalClusterOptions &clusOpts) {
  if (mols.size() < 2) {
    return ProxGraph();
  }
  ProxGraph proxGraph(mols.size(),
                      std::vector<ClusNode>(mols.size(), ClusNode()));
  std::vector<std::pair<unsigned int, unsigned int>> toDo;
  Rascal::RascalOptions opts;
  opts.similarityThreshold = clusOpts.similarityCutoff;
  for (unsigned int i = 0; i < mols.size() - 1; ++i) {
    for (unsigned int j = i + 1; j < mols.size(); ++j) {
      toDo.emplace_back(i, j);
    }
  }
  std::vector<ClusNode> molSims(toDo.size());
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  for (std::size_t k = 0; k < toDo.size(); ++k) {
    const auto [i, j] = toDo[k];
    molSims[k] =
        calc_mol_mol_similarity(*mols[i], *mols[j], i, j, opts, clusOpts);
  }
  for (const auto &cn : molSims) {
    proxGraph[cn.d_mol1Num][cn.d_mol2Num] =
        proxGraph[cn.d_mol2Num][cn.d_mol1Num] = cn;
  }
  return proxGraph;
}

// As details::disconnectProximityGraphs.
std::vector<std::vector<unsigned int>> disconnect_proximity_graphs(
    ProxGraph &proxGraph) {
  std::vector<std::vector<unsigned int>> subGraphs;
  std::vector<bool> done(proxGraph.size(), false);
  auto nextStart = std::find(done.begin(), done.end(), false);
  while (nextStart != done.end()) {
    std::list<unsigned int> nodes;
    std::list<unsigned int> toDo(1, std::distance(done.begin(), nextStart));
    while (!toDo.empty()) {
      auto nextNode = toDo.front();
      toDo.pop_front();
      if (!done[nextNode]) {
        nodes.push_back(nextNode);
      }
      done[nextNode] = true;
      for (size_t i = 0; i < proxGraph.size(); ++i) {
        if (!done[i] && proxGraph[nextNode][i].d_res) {
          toDo.push_back(i);
          nodes.push_back(i);
          done[i] = true;
        }
      }
    }
    nodes.sort();
    subGraphs.push_back(std::vector(nodes.begin(), nodes.end()));
    nextStart = std::find(done.begin(), done.end(), false);
  }
  return subGraphs;
}

// As details::g_ij.
double g_ij(const std::shared_ptr<RDKit::ROMol> &mol, double a, double b,
            unsigned int p) {
  auto molFrags = RDKit::MolOps::getMolFrags(*mol, false);
  int numBigFrags = 0;
  for (const auto &mf : molFrags) {
    if (mf->getNumBonds() > p) {
      ++numBigFrags;
    }
  }
  numBigFrags = numBigFrags == 0 ? molFrags.size() : numBigFrags;
  double g = mol->getNumAtoms();
  g += b * (1.0 - a * (numBigFrags - 1)) * mol->getNumBonds();
  return g;
}

// As details::makeSubClusters, whose MCESs match substructures too.
std::vector<std::vector<unsigned int>> make_sub_clusters(
    const std::vector<ClusNode> &nbors,
    const Rascal::RascalClusterOptions &clusOpts) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::vector<std::vector<unsigned int>> subClusters;

  std::vector<const ClusNode *> tmpNbors;
  for (const auto &n : nbors) {
    tmpNbors.push_back(&n);
  }

  while (!tmpNbors.empty()) {
    subClusters.push_back(std::vector<unsigned int>{
        tmpNbors.front()->d_mol1Num, tmpNbors.front()->d_mol2Num});
    auto m1 = tmpNbors.front()->d_res->getMcesMol();
    auto g_12 = g_ij(m1, clusOpts.a, clusOpts.b, clusOpts.minFragSize);
    for (size_t i = 1; i < tmpNbors.size(); ++i) {
      auto m2 = tmpNbors[i]->d_res->getMcesMol();
      auto g_13 = g_ij(m2, clusOpts.a, clusOpts.b, clusOpts.minFragSize);

      auto results = Rascal::rascalMCES(*m1, *m2);
      if (results.empty() || results.front().getBondMatches().empty()) {
        continue;
      }
      auto res = results.front();
      auto g_12_13 =
          g_ij(res.getMcesMol(), clusOpts.a, clusOpts.b, clusOpts.minFragSize);
      double sim = g_12_13 / std::min(g_12, g_13);
      if (sim > clusOpts.minIntraClusterSim) {
        subClusters.back().push_back(tmpNbors[i]->d_mol2Num);
        subClusters.back().push_back(tmpNbors[i]->d_mol1Num);
        tmpNbors[i] = nullptr;
      }
    }
    tmpNbors.front() = nullptr;
    tmpNbors.erase(std::remove(tmpNbors.begin(), tmpNbors.end(), nullptr),
                   tmpNbors.end());
    std::sort(subClusters.back().begin(), subClusters.back().end());
    subClusters.back().erase(
        std::unique(subClusters.back().begin(), subClusters.back().end()),
        subClusters.back().end());
  }
  return subClusters;
}

// As details::formInitialClusters.
std::vector<std::vector<unsigned int>> form_initial_clusters(
    const std::vector<unsigned int> &subGraph, const ProxGraph &proxGraph,
    const Rascal::RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> clusters;
  if (subGraph.size() < 2) {
    return clusters;
  }
  for (auto i : subGraph) {
    std::vector<ClusNode> nbors;
    for (auto j : subGraph) {
      if (proxGraph[i][j].d_res) {
        nbors.push_back(proxGraph[i][j]);
      }
    }
    std::sort(nbors.begin(), nbors.end(),
              [](const ClusNode &c1, const ClusNode &c2) -> bool {
                return c1.d_sim > c2.d_sim;
              });
    if (!nbors.empty()) {
      auto subClusters = make_sub_clusters(nbors, clusOpts);
      clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());
    }
  }
  std::sort(clusters.begin(), clusters.end(),
            [](const std::vector<unsigned int> &c1,
               const std::vector<unsigned int> &c2) -> bool {
              if (c1.size() == c2.size()) {
                return c1.front() < c2.front();
              } else {
                return c1.size() > c2.size();
              }
            });
  clusters.erase(std::unique(clusters.begin(), clusters.end()), clusters.end());
  return clusters;
}

// As details::mergeClusters.
std::vector<std::vector<unsigned int>> merge_clusters(
    const std::vector<std::vector<unsigned int>> &clusters,
    const Rascal::RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> outClusters(clusters);

  if (outClusters.size() < 2) {
    return outClusters;
  }

  for (size_t i = 0; i < outClusters.size() - 1; ++i) {
    for (size_t j = i + 1; j < outClusters.size(); ++j) {
      std::vector<int> inCommon;
      std::set_intersection(outClusters[i].begin(), outClusters[i].end(),
                            outClusters[j].begin(), outClusters[j].end(),
                            std::back_inserter(inCommon));
      double s =
          double(inCommon.size()) / std::min(double(outClusters[i].size()),
                                             double(outClusters[j].size()));
      if (s > clusOpts.clusterMergeSim) {
        outClusters[i].insert(outClusters[i].end(), outClusters[j].begin(),
                              outClusters[j].end());
        outClusters[j].clear();
        std::sort(outClusters[i].begin(), outClusters[i].end());
        outClusters[i].erase(
            std::unique(outClusters[i].begin(), outClusters[i].end()),
            outClusters[i].end());
      }
    }
    outClusters.erase(
        std::remove_if(outClusters.begin(), outClusters.end(),
                       [](const std::vector<unsigned int> &c) -> bool {
                         return c.empty();
                       }),
        outClusters.end());
  }

  return outClusters;
}

// As details::sortClusterMembersByMeanSim.
void sort_cluster_members_by_mean_sim(
    const ProxGraph &proxGraph,
    std::vector<std::vector<unsigned int>> &clusters) {
  for (auto &clus : clusters) {
    std::vector<std::pair<unsigned int, double>> clusSims;
    for (unsigned int i = 0U; i < clus.size(); ++i) {
      double totSim = 0.0;
      for (unsigned int j = 0U; j < clus.size(); ++j) {
        if (i != j) {
          totSim += proxGraph[clus[i]][clus[j]].d_sim;
        }
      }
      clusSims.push_back({clus[i], totSim / (clus.size() - 1)});
    }
    std::sort(clusSims.begin(), clusSims.end(),
              [](const std::pair<unsigned int, double> &p1,
                 const std::pair<unsigned int, double> &p2) -> bool {
                return p1.second > p2.second;
              });
    std::transform(
        clusSims.begin(), clusSims.end(), clus.begin(),
        [](const std::pair<unsigned int, double> &p) -> unsigned int {
          return p.first;
        });
  }
}

// As details::makeClusters.
std::vector<std::vector<unsigned int>> make_clusters(
    const std::vector<std::vector<unsigned int>> &subGraphs,
    const ProxGraph &proxGraph, const Rascal::RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> clusters;
  for (const auto &sg : subGraphs) {
    auto theseClusters = form_initial_clusters(sg, proxGraph, clusOpts);
    auto mergedClusters = merge_clusters(theseClusters, clusOpts);
    clusters.insert(clusters.end(), mergedClusters.begin(),
                    mergedClusters.end());
  }
  std::sort(clusters.begin(), clusters.end(),
            [](const std::vector<unsigned int> &c1,
               const std::vector<unsigned int> &c2) -> bool {
              return c1.size() > c2.size();
            });
  return clusters;
}

// As details::collectSingletons.
std::vector<unsigned int> collect_singletons(const ProxGraph &proxGraph) {
  std::vector<unsigned int> singletons;
  for (size_t i = 0; i < proxGraph.size(); ++i) {
    bool single = true;
    for (const auto &cn : proxGraph[i]) {
      if (cn.d_res) {
        single = false;
        break;
      }
    }
    if (single) {
      singletons.push_back(i);
    }
  }
  return singletons;
}

// As details::buildNborLists in RascalButinaCluster.cpp.
std::vector<std::vector<unsigned int>> build_nbor_lists(
    const ProxGraph &proxGraph) {
  std::vector<std::vector<unsigned int>> nborLists;
  for (size_t i = 0; i < proxGraph.size(); ++i) {
    std::vector<std::pair<unsigned int, double>> tmpList;
    for (const auto &cn : proxGraph[i]) {
      if (cn.d_res) {
        if (i == cn.d_mol1Num) {
          tmpList.push_back({cn.d_mol2Num, cn.d_sim});
        } else {
          tmpList.push_back({cn.d_mol1Num, cn.d_sim});
        }
      }
    }
    std::sort(tmpList.begin(), tmpList.end(),
              [](const std::pair<unsigned int, double> &p1,
                 const std::pair<unsigned int, double> &p2) -> bool {
                return p1.second > p2.second;
              });
    std::vector<unsigned int> nborList(tmpList.size() + 1, 0);
    nborList[0] = i;
    std::transform(
        tmpList.begin(), tmpList.end(), nborList.begin() + 1,
        [](const std::pair<unsigned int, double> &p) -> unsigned int {
          return p.first;
        });
    nborLists.push_back(nborList);
  }
  std::sort(nborLists.begin(), nborLists.end(),
            [](const std::vector<unsigned int> &nl1,
               const std::vector<unsigned int> &nl2) -> bool {
              if (nl1.size() == nl2.size()) {
                return nl1 > nl2;
              } else {
                return nl1.size() > nl2.size();
              }
            });
  return nborLists;
}

// As details::formClusters in RascalButinaCluster.cpp; destroys nborLists.
std::vector<std::vector<unsigned int>> form_clusters(
    std::vector<std::vector<unsigned int>> &nborLists) {
  std::vector<std::vector<unsigned int>> clusters;

  while (!nborLists.empty()) {
    clusters.push_back(nborLists.front());
    std::set<unsigned int> inNborList(nborLists.front().begin(),
                                      nborLists.front().end());
    nborLists.front().clear();
    for (auto &nborList : nborLists) {
      for (auto &n : nborList) {
        if (inNborList.find(n) != inNborList.end()) {
          n = std::numeric_limits<unsigned int>::max();
        }
      }
      nborList.erase(std::remove(nborList.begin(), nborList.end(),
                                 std::numeric_limits<unsigned int>::max()),
                     nborList.end());
    }
    nborLists.erase(
        std::remove_if(nborLists.begin(), nborLists.end(),
                       [](const std::vector<unsigned int> &nl) -> bool {
                         return nl.empty();
                       }),
        nborLists.end());
    std::sort(nborLists.begin(), nborLists.end(),
              [](const std::vector<unsigned int> &nl1,
                 const std::vector<unsigned int> &nl2) -> bool {
                if (nl1.size() == nl2.size()) {
                  return nl1 > nl2;
                } else {
                  return nl1.size() > nl2.size();
                }
              });
  }
  return clusters;
}

std::vector<const RDKit::ROMol *> mol_pointers(BorrowedMols &mols) {
  std::vector<const RDKit::ROMol *> res;
  for (const auto &mol : mols.mols()) {
    res.push_back(mol.get());
  }
  check_rascal_mols(res);
  return res;
}

}  // namespace

// RDKit trusts the options (RascalMCES.cpp):
// - makeModularProduct computes the bond distance matrices only for a
//   maxFragSeparation above -1 but reads them for any value but -1.
// - calcLowerBound casts a product with similarityThreshold to unsigned
//   int; a NaN threshold passes every comparison before it.
// - assignEquivalentAtoms dereferences MolFromSmarts's null for a class of
//   equivalentAtoms that does not parse.
// completeSmallestRings implies completeAromaticRings (RascalOptions.h), and
// buildPairs reads the ring SMILES that only completeAromaticRings makes, so
// it is set here as documented.
Rascal::RascalOptions rascal_options_from(const RascalOptionsData &d) {
  Rascal::RascalOptions o;
  o.similarityThreshold = d.similarity_threshold;
  o.completeAromaticRings =
      d.complete_aromatic_rings || d.complete_smallest_rings;
  o.ringMatchesRingOnly = d.ring_matches_ring_only;
  o.completeSmallestRings = d.complete_smallest_rings;
  o.exactConnectionsMatch = d.exact_connections_match;
  o.singleLargestFrag = d.single_largest_frag;
  o.minFragSize = d.min_frag_size;
  o.maxFragSeparation = d.max_frag_separation;
  o.allBestMCESs = d.all_best_mcess;
  o.maxBestMCESs = d.max_best_mcess;
  o.timeout = d.timeout;
  o.returnEmptyMCES = d.return_empty_mces;
  o.maxBondMatchPairs = d.max_bond_match_pairs;
  o.equivalentAtoms = std::string(d.equivalent_atoms);
  o.ignoreBondOrders = d.ignore_bond_orders;
  o.ignoreAtomAromaticity = d.ignore_atom_aromaticity;
  o.minCliqueSize = d.min_clique_size;

  if (o.maxFragSeparation < -1) {
    throw ValueErrorException("maxFragSeparation must be -1 or more");
  }
  if (std::isnan(o.similarityThreshold)) {
    throw ValueErrorException("similarityThreshold is NaN");
  }
  // As assignEquivalentAtoms splits the classes; it refuses more than 9.
  std::vector<std::string> classes;
  boost::split(classes, o.equivalentAtoms, boost::is_any_of(" "),
               boost::token_compress_on);
  if (!o.equivalentAtoms.empty() && classes.size() <= 9) {
    for (const auto &smarts : classes) {
      if (!smarts.empty() && !RDKit::v2::SmilesParse::MolFromSmarts(smarts)) {
        throw RDKit::SmilesParse::SmilesParseException(
            "could not parse the equivalentAtoms SMARTS '" + smarts + "'");
      }
    }
  }
  return o;
}

// The results remove atoms and bonds from copies of the molecules
// (makeCliqueFrags, extractRings, makeMolFrags), which trusts what
// check_removals checks.
void check_rascal_mols(const std::vector<const RDKit::ROMol *> &mols) {
  int large = 0;
  for (const auto *mol : mols) {
    check_matrix_size(*mol);
    check_removals(*mol);
    if (std::int64_t{mol->getNumAtoms()} + mol->getNumBonds() >
        kRascalMaxSquared) {
      ++large;
    }
  }
  if (large > 1) {
    throw ValueErrorException(
        "RDKit squares a count of matched atoms and bonds in int, which "
        "overflows when two molecules each have more than 46340 atoms and "
        "bonds");
  }
}

RascalOptionsData rascal_options_default() {
  return to_data(Rascal::RascalOptions());
}

RascalClusterOptionsData rascal_cluster_options_default() {
  return to_data(Rascal::RascalClusterOptions());
}

std::unique_ptr<RascalResultList> find_mces(const RWMol &mol1,
                                            const RWMol &mol2,
                                            const RascalOptionsData &opts) {
  auto options = rascal_options_from(opts);
  check_rascal_mols({&mol1, &mol2});
  auto list = std::make_unique<RascalResultList>();
  // The MCES matches substructures (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  list->results = Rascal::rascalMCES(mol1, mol2, options);
  return list;
}

std::size_t rascal_result_list_size(const RascalResultList &list) {
  return list.results.size();
}

std::unique_ptr<RascalResult> rascal_result_list_take(RascalResultList &list,
                                                      std::size_t i) {
  return std::make_unique<RascalResult>(std::move(list.results.at(i)));
}

// The copy constructor copies the result's molecules.
std::unique_ptr<RascalResult> rascal_result_clone(const RascalResult &result) {
  return std::make_unique<RascalResult>(result);
}

std::unique_ptr<std::string> rascal_result_smarts_string(
    const RascalResult &result) {
  return std::make_unique<std::string>(result.getSmarts());
}

// Bond and atom indices, which are never negative.
rust::Vec<IndexPair> rascal_result_bond_matches(const RascalResult &result) {
  return to_pairs(result.getBondMatches());
}

rust::Vec<IndexPair> rascal_result_atom_matches(const RascalResult &result) {
  return to_pairs(result.getAtomMatches());
}

// A result returnEmptyMCES makes has no molecules, and largestFragOnly
// dereferences the fragment molecule it then lacks (makeMolFrags's null).
// With no atom matches the call changes nothing anyway: the fragments are
// none.
void rascal_result_largest_fragment_only(RascalResult &result) {
  if (result.getAtomMatches().empty()) {
    return;
  }
  result.largestFragOnly();
}

double rascal_result_similarity(const RascalResult &result) {
  return result.getSimilarity();
}

int rascal_result_num_fragments(const RascalResult &result) {
  return result.getNumFrags();
}

std::uint32_t rascal_result_largest_fragment_size(const RascalResult &result) {
  return result.getLargestFragSize();
}

double rascal_result_tier1_sim(const RascalResult &result) {
  return result.getTier1Sim();
}

double rascal_result_tier2_sim(const RascalResult &result) {
  return result.getTier2Sim();
}

bool rascal_result_timed_out(const RascalResult &result) {
  return result.getTimedOut();
}

// As rascalCluster in RascalCluster.cpp.
rust::Vec<IndexList> rascal_cluster(BorrowedMols &mols,
                                    const RascalClusterOptionsData &opts) {
  const auto clusOpts = cluster_options_from(opts);
  const auto ptrs = mol_pointers(mols);
  auto proxGraph = build_proximity_graph(ptrs, clusOpts);
  auto subGraphs = disconnect_proximity_graphs(proxGraph);
  auto clusters = make_clusters(subGraphs, proxGraph, clusOpts);
  auto singletons = collect_singletons(proxGraph);
  clusters.push_back(singletons);
  sort_cluster_members_by_mean_sim(proxGraph, clusters);
  return to_lists(clusters);
}

// As rascalButinaCluster in RascalButinaCluster.cpp.
rust::Vec<IndexList> rascal_butina_cluster(
    BorrowedMols &mols, const RascalClusterOptionsData &opts) {
  const auto clusOpts = cluster_options_from(opts);
  const auto ptrs = mol_pointers(mols);
  auto proxGraph = build_proximity_graph(ptrs, clusOpts);
  auto nborLists = build_nbor_lists(proxGraph);
  return to_lists(form_clusters(nborLists));
}

}  // namespace rdkit_cxx
