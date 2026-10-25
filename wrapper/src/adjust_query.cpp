/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/adjust_query.h"

#include <GraphMol/GenericGroups/GenericGroups.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/adjust_query.rs.h"

namespace rdkit_cxx {

namespace {

AdjustQueryParametersData to_data(
    const RDKit::MolOps::AdjustQueryParameters &params) {
  return AdjustQueryParametersData{
      params.adjustDegree,
      params.adjustDegreeFlags,
      params.adjustHeavyDegree,
      params.adjustHeavyDegreeFlags,
      params.adjustRingCount,
      params.adjustRingCountFlags,
      params.makeDummiesQueries,
      params.aromatizeIfPossible,
      params.makeBondsGeneric,
      params.makeBondsGenericFlags,
      params.makeAtomsGeneric,
      params.makeAtomsGenericFlags,
      params.adjustRingChain,
      params.adjustRingChainFlags,
      params.useStereoCareForBonds,
      params.adjustConjugatedFiveRings,
      params.setMDLFiveRingAromaticity,
      params.adjustSingleBondsToDegreeOneNeighbors,
      params.adjustSingleBondsBetweenAromaticAtoms};
}

RDKit::MolOps::AdjustQueryParameters from_data(
    const AdjustQueryParametersData &data) {
  RDKit::MolOps::AdjustQueryParameters params;
  params.adjustDegree = data.adjust_degree;
  params.adjustDegreeFlags = data.adjust_degree_flags;
  params.adjustHeavyDegree = data.adjust_heavy_degree;
  params.adjustHeavyDegreeFlags = data.adjust_heavy_degree_flags;
  params.adjustRingCount = data.adjust_ring_count;
  params.adjustRingCountFlags = data.adjust_ring_count_flags;
  params.makeDummiesQueries = data.make_dummies_queries;
  params.aromatizeIfPossible = data.aromatize_if_possible;
  params.makeBondsGeneric = data.make_bonds_generic;
  params.makeBondsGenericFlags = data.make_bonds_generic_flags;
  params.makeAtomsGeneric = data.make_atoms_generic;
  params.makeAtomsGenericFlags = data.make_atoms_generic_flags;
  params.adjustRingChain = data.adjust_ring_chain;
  params.adjustRingChainFlags = data.adjust_ring_chain_flags;
  params.useStereoCareForBonds = data.use_stereo_care_for_bonds;
  params.adjustConjugatedFiveRings = data.adjust_conjugated_five_rings;
  params.setMDLFiveRingAromaticity = data.set_mdl_five_ring_aromaticity;
  params.adjustSingleBondsToDegreeOneNeighbors =
      data.adjust_single_bonds_to_degree_one_neighbors;
  params.adjustSingleBondsBetweenAromaticAtoms =
      data.adjust_single_bonds_between_aromatic_atoms;
  return params;
}

}  // namespace

AdjustQueryParametersData adjust_query_params_to_data(
    const RDKit::MolOps::AdjustQueryParameters &params) {
  return to_data(params);
}

RDKit::MolOps::AdjustQueryParameters adjust_query_params_from_data(
    const AdjustQueryParametersData &data) {
  return from_data(data);
}

AdjustQueryParametersData adjust_query_parameters_default() {
  return to_data(RDKit::MolOps::AdjustQueryParameters());
}

AdjustQueryParametersData adjust_query_parameters_no_adjustments() {
  return to_data(RDKit::MolOps::AdjustQueryParameters::noAdjustments());
}

std::unique_ptr<RWMol> adjust_query_properties(
    const RWMol &mol, const AdjustQueryParametersData &data) {
  auto params = from_data(data);
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::adjustQueryProperties(*res, &params);
  return res;
}

std::unique_ptr<RWMol> adjust_query_properties_with_generic_groups(
    const RWMol &mol, const AdjustQueryParametersData &data) {
  auto params = from_data(data);
  std::unique_ptr<RDKit::ROMol> res(
      RDKit::GenericGroups::adjustQueryPropertiesWithGenericGroups(mol,
                                                                   &params));
  return std::make_unique<RWMol>(*res);
}

void set_generic_queries_from_properties(RWMol &mol, bool use_atom_labels,
                                         bool use_s_groups) {
  RDKit::GenericGroups::setGenericQueriesFromProperties(mol, use_atom_labels,
                                                        use_s_groups);
}

void convert_generic_queries_to_substance_groups(RWMol &mol) {
  RDKit::GenericGroups::convertGenericQueriesToSubstanceGroups(mol);
}

}  // namespace rdkit_cxx
