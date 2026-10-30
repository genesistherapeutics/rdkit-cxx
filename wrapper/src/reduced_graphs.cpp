// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/reduced_graphs.h"

#include <mutex>
#include <string>
#include <vector>

#include <GraphMol/ReducedGraphs/ReducedGraphs.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

// The feature types RDKit's default patterns give (nFeatures in
// ReducedGraphs.cpp), plus one for aromatic ring centroids.
constexpr int num_types = 6;

// generateErGFingerprintForReducedGraph computes each bin from the atoms'
// _ErGAtomTypes in int (overflowing for large types) and from the path
// bounds in unsigned (wrapping for negative ones), and sizes the vector as
// 21 times the number of bins in unsigned.
void check_erg_arguments(const RDKit::ROMol &mol, std::int32_t min_path,
                         std::int32_t max_path) {
  if (min_path < 0 || max_path <= min_path) {
    throw ValueErrorException("ErG paths need 0 <= minPath < maxPath");
  }
  std::uint64_t bins = static_cast<std::uint64_t>(max_path) - min_path + 1;
  if (bins * (num_types * (num_types + 1) / 2) > UINT32_MAX) {
    throw ValueErrorException("maxPath - minPath is too large");
  }
  for (const auto atom : mol.atoms()) {
    std::vector<int> types;
    if (atom->getPropIfPresent("_ErGAtomTypes", types)) {
      for (auto type : types) {
        if (type < 0 || type >= num_types) {
          throw ValueErrorException("atom " + std::to_string(atom->getIdx()) +
                                    " has an _ErGAtomTypes entry of " +
                                    std::to_string(type) +
                                    ", which is not an ErG type");
        }
      }
    }
  }
}

// The reduced graph marks features found by substructure matches
// (control_c_lock.h).
// It removes atoms from a copy of the molecule (check_removals).
std::unique_ptr<RWMol> reduced_graph(const RWMol &mol) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::unique_ptr<RDKit::ROMol> graph(
      RDKit::ReducedGraphs::generateMolExtendedReducedGraph(mol));
  return std::make_unique<RWMol>(*graph);
}

rust::Vec<double> erg_fingerprint(const RDKit::ROMol &graph,
                                  double fuzz_increment, std::int32_t min_path,
                                  std::int32_t max_path) {
  check_erg_arguments(graph, min_path, max_path);
  // The fingerprint reads the reduced graph's distance matrix (U3, U6).
  check_matrix_size(graph);
  MatrixCacheGuard guard(graph);
  std::unique_ptr<RDNumeric::DoubleVector> fp(
      RDKit::ReducedGraphs::generateErGFingerprintForReducedGraph(
          graph, nullptr, fuzz_increment, min_path, max_path));
  rust::Vec<double> result;
  result.reserve(fp->size());
  for (unsigned int i = 0; i < fp->size(); ++i) {
    result.push_back((*fp)[i]);
  }
  return result;
}

}  // namespace

std::unique_ptr<RWMol> generate_mol_extended_reduced_graph(const RWMol &mol) {
  // The copy RDKit edits keeps the input's computed properties, but no
  // molecule Rust holds has a cached matrix (U6).
  auto graph = reduced_graph(mol);
  clear_matrix_caches(*graph);
  return graph;
}

rust::Vec<double> generate_erg_fingerprint_for_reduced_graph(
    const RWMol &mol, double fuzz_increment, std::int32_t min_path,
    std::int32_t max_path) {
  return erg_fingerprint(mol, fuzz_increment, min_path, max_path);
}

rust::Vec<double> get_erg_fingerprint(const RWMol &mol, double fuzz_increment,
                                      std::int32_t min_path,
                                      std::int32_t max_path) {
  // As getErGFingerprint, which makes the reduced graph and fingerprints it;
  // the steps are separate here so that the graph's size is checked.
  auto graph = reduced_graph(mol);
  clear_matrix_caches(*graph);
  return erg_fingerprint(*graph, fuzz_increment, min_path, max_path);
}

}  // namespace rdkit_cxx
