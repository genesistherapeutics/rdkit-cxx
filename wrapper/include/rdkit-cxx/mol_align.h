/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Wrappers for GraphMol/MolAlign (rdkit.Chem.rdMolAlign): RMSD alignment of
// molecules and conformers, and Open3DAlign.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <boost/shared_ptr.hpp>

#include "rdkit-cxx/force_field.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/chirality.rs and
// src/ffi/mol_align.rs).
struct IndexPair;
struct BestAlignmentParamsData;
struct O3AOptionsData;

// One O3A alignment. RDKit's O3A keeps pointers to both molecules, which the
// Rust side keeps borrowed for as long as this lives.
class O3AHolder {
 public:
  explicit O3AHolder(boost::shared_ptr<RDKit::MolAlign::O3A> o3a)
      : d_o3a(std::move(o3a)) {}
  RDKit::MolAlign::O3A &get() const { return *d_o3a; }

 private:
  boost::shared_ptr<RDKit::MolAlign::O3A> d_o3a;
};

class O3AList {
 public:
  std::vector<boost::shared_ptr<RDKit::MolAlign::O3A>> items;
};

BestAlignmentParamsData best_alignment_params_default();

double get_alignment_transform(const RWMol &prb_mol, const RWMol &ref_mol,
                               std::int32_t prb_cid, std::int32_t ref_cid,
                               rust::Slice<const IndexPair> atom_map,
                               rust::Slice<const double> weights, bool reflect,
                               std::uint32_t max_iters,
                               rust::Vec<double> &trans);
double align_mol(RWMol &prb_mol, const RWMol &ref_mol, std::int32_t prb_cid,
                 std::int32_t ref_cid, rust::Slice<const IndexPair> atom_map,
                 rust::Slice<const double> weights, bool reflect,
                 std::uint32_t max_iters);
double get_best_alignment_transform(const RWMol &prb_mol, const RWMol &ref_mol,
                                    const BestAlignmentParamsData &params,
                                    rust::Vec<double> &trans,
                                    rust::Vec<IndexPair> &best_match);
double get_best_rms(RWMol &prb_mol, const RWMol &ref_mol,
                    const BestAlignmentParamsData &params);
rust::Vec<double> get_all_conformer_best_rms(
    const RWMol &mol, const BestAlignmentParamsData &params);
rust::Vec<double> get_all_conformer_best_rms_to_ref(
    const RWMol &prb_mol, const RWMol &ref_mol,
    const BestAlignmentParamsData &params);
double calc_rms(const RWMol &prb_mol, const RWMol &ref_mol,
                const BestAlignmentParamsData &params);
rust::Vec<double> align_mol_conformers(
    RWMol &mol, rust::Slice<const std::uint32_t> atom_ids,
    rust::Slice<const std::uint32_t> conf_ids,
    rust::Slice<const double> weights, bool reflect, std::uint32_t max_iters);
void random_transform(RWMol &mol, std::int32_t cid, std::int32_t seed);

std::unique_ptr<O3AList> get_mmff_o3a(RWMol &prb_mol, const RWMol &ref_mol,
                                      const MMFFMolProperties *prb_props,
                                      const MMFFMolProperties *ref_props,
                                      bool for_probe_confs,
                                      const O3AOptionsData &options);
std::unique_ptr<O3AList> get_crippen_o3a(RWMol &prb_mol, const RWMol &ref_mol,
                                         bool has_prb_contribs,
                                         rust::Slice<const double> prb_contribs,
                                         bool has_ref_contribs,
                                         rust::Slice<const double> ref_contribs,
                                         bool for_probe_confs,
                                         const O3AOptionsData &options);
std::size_t o3a_list_size(const O3AList &list);
std::unique_ptr<O3AHolder> o3a_list_take(O3AList &list, std::size_t i);
double o3a_align(O3AHolder &o3a);
double o3a_trans(const O3AHolder &o3a, rust::Vec<double> &trans);
double o3a_score(const O3AHolder &o3a);
rust::Vec<IndexPair> o3a_matches(const O3AHolder &o3a);
rust::Vec<double> o3a_weights(const O3AHolder &o3a);

}  // namespace rdkit_cxx
