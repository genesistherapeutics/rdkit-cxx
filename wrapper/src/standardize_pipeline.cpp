// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_pipeline.h"

#include <exception>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

#include <GraphMol/MolOps.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/MolStandardize/Metal.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/standardize_pipeline.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/standardize_normalize.h"
#include "rdkit-cxx/streams.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;

namespace {

std::string text(const rust::Vec<std::uint8_t> &bytes) {
  return to_std(rust::Slice<const std::uint8_t>(bytes.data(), bytes.size()));
}

MS::PipelineOptions to_options(const PipelineOptionsData &data) {
  MS::PipelineOptions options;
  options.strictParsing = data.strict_parsing;
  options.reportAllFailures = data.report_all_failures;
  options.allowEmptyMolecules = data.allow_empty_molecules;
  options.allowEnhancedStereo = data.allow_enhanced_stereo;
  options.allowAromaticBondType = data.allow_aromatic_bond_type;
  options.allowDativeBondType = data.allow_dative_bond_type;
  options.is2DZeroThreshold = data.is_2d_zero_threshold;
  options.atomClashLimit = data.atom_clash_limit;
  options.minMedianBondLength = data.min_median_bond_length;
  options.bondLengthLimit = data.bond_length_limit;
  options.allowLongBondsInRings = data.allow_long_bonds_in_rings;
  options.allowAtomBondClashExemption = data.allow_atom_bond_clash_exemption;
  options.metalNof = text(data.metal_nof);
  options.metalNon = text(data.metal_non);
  options.normalizerData = text(data.normalizer_data);
  options.normalizerMaxRestarts = data.normalizer_max_restarts;
  options.scaledMedianBondLength = data.scaled_median_bond_length;
  options.outputV2000 = data.output_v2000;
  return options;
}

void check_metal_smarts(const std::string &smarts, const char *what) {
  check_smiles_cx(smarts, true, true);
  std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
  if (!query) {
    throw ValueErrorException(std::string(what) +
                              " SMARTS do not parse: " + smarts);
  }
  check_metal_query(*query);
}

// As Operations::parse, and then check_removals on the molecule it reads:
// a mol block sets its substance groups' PARENT chains, and the validation
// and standardization steps remove atoms and bonds from it and its copies
// (metal disconnection, fragment removal, getMolFrags).
RDKit::RWMOL_SPTR checked_parse(const std::string &molblock,
                                MS::PipelineResult &result,
                                const MS::PipelineOptions &options) {
  auto mol = MS::Operations::parse(molblock, result, options);
  if (mol) {
    check_removals(*mol);
  }
  return mol;
}

// As Operations::standardize, which static_casts the ROMol that
// Normalizer::normalize makes to an RWMol, a type it is not; this copies it
// into an RWMol instead. The rest is RDKit's, step for step.
RDKit::RWMOL_SPTR checked_standardize(RDKit::RWMOL_SPTR mol,
                                      MS::PipelineResult &result,
                                      const MS::PipelineOptions &options) {
  auto smiles = RDKit::MolToSmiles(*mol);
  auto reference = smiles;

  // pipeline_new checked that both SMARTS parse.
  try {
    MS::MetalDisconnectorOptions mdOpts;
    MS::MetalDisconnector metalDisconnector(mdOpts);
    std::unique_ptr<RDKit::ROMol> metalNof{
        RDKit::SmartsToMol(options.metalNof)};
    metalDisconnector.setMetalNof(*metalNof);
    std::unique_ptr<RDKit::ROMol> metalNon{
        RDKit::SmartsToMol(options.metalNon)};
    metalDisconnector.setMetalNon(*metalNon);
    metalDisconnector.disconnectInPlace(*mol);
  } catch (...) {
    result.append(
        MS::METAL_STANDARDIZATION_ERROR,
        "An error occurred while processing the bonding of metal species.");
    return mol;
  }

  smiles = RDKit::MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(MS::METALS_DISCONNECTED,
                  "One or more metal atoms were disconnected.");
  }
  reference = smiles;

  try {
    std::unique_ptr<MS::Normalizer> normalizer{};
    if (options.normalizerData.empty()) {
      normalizer.reset(new MS::Normalizer);
    } else {
      std::istringstream sstr(options.normalizerData);
      normalizer.reset(new MS::Normalizer(sstr, options.normalizerMaxRestarts));
    }
    std::unique_ptr<RDKit::ROMol> normalized(normalizer->normalize(*mol));
    mol.reset(new RDKit::RWMol(*normalized));
    mol->updatePropertyCache(false);
  } catch (...) {
    result.append(MS::NORMALIZER_STANDARDIZATION_ERROR,
                  "An error occurred while normalizing the representation of "
                  "some functional groups");
    return mol;
  }

  smiles = RDKit::MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(MS::NORMALIZATION_APPLIED,
                  "The representation of some functional groups was adjusted.");
  }
  reference = smiles;

  try {
    MS::LargestFragmentChooser fragmentChooser;
    fragmentChooser.chooseInPlace(*mol);
  } catch (...) {
    result.append(
        MS::FRAGMENT_STANDARDIZATION_ERROR,
        "An error occurred while removing the disconnected fragments");
    return mol;
  }

  smiles = RDKit::MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(MS::FRAGMENTS_REMOVED,
                  "One or more disconnected fragments (e.g., counterions) were "
                  "removed.");
  }

  constexpr bool cleanIt = true;
  constexpr bool force = true;
  constexpr bool flagPossible = true;
  RDKit::MolOps::assignStereochemistry(*mol, cleanIt, force, flagPossible);

  return mol;
}

// As Operations::makeParent, which sums the formal charges of the molecule
// and of its parent in int (U9); the uncharger it runs first only shrinks
// their magnitudes.
MS::RWMOL_SPTR_PAIR checked_make_parent(RDKit::RWMOL_SPTR mol,
                                        MS::PipelineResult &result,
                                        const MS::PipelineOptions &options) {
  check_formal_charge_sum(*mol);
  return MS::Operations::makeParent(mol, result, options);
}

// As Operations::serialize, whose mol-block writers compute 2D coordinates
// for a molecule without a conformer (prepareMol in MolFileWriter.cpp), so
// it holds the depiction lock then and checks the molecule under it (U4).
// serialize records its writers' exceptions as output errors; so does this
// for the check's. Serializing matches nothing, and the depiction lock comes
// before control_c_mutex (U4), so pipeline_run's is released for it.
void locked_serialize(MS::RWMOL_SPTR_PAIR output, MS::PipelineResult &result,
                      const MS::PipelineOptions &options) {
  ReleaseControlC release;
  std::unique_lock<std::mutex> lock;
  if (!output.first->getNumConformers() || !output.second->getNumConformers()) {
    lock = lock_depiction();
    try {
      for (const auto *mol : {output.first.get(), output.second.get()}) {
        if (!mol->getNumConformers()) {
          check_compute_2d_coords_mol(*mol, false);
        }
      }
    } catch (const std::exception &e) {
      result.append(
          MS::OUTPUT_ERROR,
          "Can't write molecule to output format: " + std::string(e.what()));
      return;
    }
  }
  MS::Operations::serialize(output, result, options);
}

}  // namespace

PipelineOptionsData pipeline_options_default() {
  MS::PipelineOptions options;
  PipelineOptionsData data;
  data.strict_parsing = options.strictParsing;
  data.report_all_failures = options.reportAllFailures;
  data.allow_empty_molecules = options.allowEmptyMolecules;
  data.allow_enhanced_stereo = options.allowEnhancedStereo;
  data.allow_aromatic_bond_type = options.allowAromaticBondType;
  data.allow_dative_bond_type = options.allowDativeBondType;
  data.is_2d_zero_threshold = options.is2DZeroThreshold;
  data.atom_clash_limit = options.atomClashLimit;
  data.min_median_bond_length = options.minMedianBondLength;
  data.bond_length_limit = options.bondLengthLimit;
  data.allow_long_bonds_in_rings = options.allowLongBondsInRings;
  data.allow_atom_bond_clash_exemption = options.allowAtomBondClashExemption;
  data.metal_nof = to_bytes(options.metalNof);
  data.metal_non = to_bytes(options.metalNon);
  data.normalizer_data = to_bytes(options.normalizerData);
  data.normalizer_max_restarts = options.normalizerMaxRestarts;
  data.scaled_median_bond_length = options.scaledMedianBondLength;
  data.output_v2000 = options.outputV2000;
  return data;
}

std::unique_ptr<Pipeline> pipeline_new(const PipelineOptionsData &data) {
  auto options = to_options(data);
  check_metal_smarts(options.metalNof, "metalNof");
  check_metal_smarts(options.metalNon, "metalNon");
  // The normalizer reads these reactions with their CXSMILES extensions
  // (file_checks.h).
  check_transformation_lines(options.normalizerData);
  auto pipeline = std::make_unique<Pipeline>(options);
  auto steps = MS::Operations::standardizationSteps;
  bool replaced = false;
  for (auto &step : steps) {
    if (step.second == &MS::Operations::standardize) {
      step.second = &checked_standardize;
      replaced = true;
    }
  }
  if (!replaced) {
    throw std::logic_error("RDKit's Pipeline has no standardize step");
  }
  pipeline->setStandardizationSteps(steps);
  pipeline->setParse(&checked_parse);
  pipeline->setMakeParent(&checked_make_parent);
  pipeline->setSerialize(&locked_serialize);
  return pipeline;
}

PipelineResultData pipeline_run(const Pipeline &pipeline, rust::Str mol_data) {
  std::string block(mol_data);
  check_input_text(block);
  // Validating and standardizing match substructures (control_c_lock.h).
  auto result = [&] {
    HoldControlC control_c;
    return pipeline.run(block);
  }();
  PipelineResultData data;
  data.status = static_cast<std::uint32_t>(result.status);
  data.stage = result.stage;
  for (const auto &entry : result.log) {
    data.log.push_back(PipelineLogEntryData{
        static_cast<std::uint32_t>(entry.status), to_bytes(entry.detail)});
  }
  data.input_mol_data = to_bytes(result.inputMolData);
  data.output_mol_data = to_bytes(result.outputMolData);
  data.parent_mol_data = to_bytes(result.parentMolData);
  return data;
}

}  // namespace rdkit_cxx
