// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize.h"

#include <algorithm>
#include <climits>
#include <exception>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/MolStandardize/Charge.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/standardize.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/threads.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;

namespace {

std::string text(const rust::Vec<std::uint8_t> &bytes) {
  return to_std(rust::Slice<const std::uint8_t>(bytes.data(), bytes.size()));
}

std::string path(const rust::Vec<std::uint8_t> &bytes) {
  return to_path(rust::Slice<const std::uint8_t>(bytes.data(), bytes.size()));
}

CleanupDataEntry entry(const std::string &name, const std::string &first,
                       const std::string &second = {},
                       const std::string &third = {}) {
  return CleanupDataEntry{to_bytes(name), to_bytes(first), to_bytes(second),
                          to_bytes(third)};
}

bool uses_tautomers(StandardizeOp op) {
  return op == StandardizeOp::CanonicalTautomer ||
         op == StandardizeOp::TautomerParent ||
         op == StandardizeOp::SuperParent;
}

// Whether the op reionizes or uncharges the molecule, which sums its charges.
bool sums_charges(StandardizeOp op, bool skip_standardize) {
  switch (op) {
    case StandardizeOp::Normalize:
    case StandardizeOp::RemoveFragments:
    case StandardizeOp::CanonicalTautomer:
      return false;
    case StandardizeOp::FragmentParent:
    case StandardizeOp::StereoParent:
    case StandardizeOp::IsotopeParent:
      return !skip_standardize;
    default:
      return true;
  }
}

// As canonicalTautomerInPlace in MolStandardize.cpp, with a checked
// enumerator.
void canonical_tautomer_in_place(RWMol &mol,
                                 const MS::TautomerEnumerator &enumerator) {
  check_tautomer_score(mol);
  enumerator.canonicalizeInPlace(mol);
}

// As tautomerParentInPlace in MolStandardize.cpp, which makes its enumerator
// from the parameters (canonicalTautomerInPlace); here it is checked.
void tautomer_parent_in_place(RWMol &mol, const CleanupParameters &params,
                              bool skip_standardize,
                              const MS::TautomerEnumerator &enumerator) {
  if (!skip_standardize) {
    MS::cleanupInPlace(mol, params);
  }
  canonical_tautomer_in_place(mol, enumerator);
  MS::cleanupInPlace(mol, params);
}

// As superParentInPlace in MolStandardize.cpp, with tautomer_parent_in_place.
void super_parent_in_place(RWMol &mol, const CleanupParameters &params,
                           bool skip_standardize,
                           const MS::TautomerEnumerator &enumerator) {
  if (!skip_standardize) {
    MS::cleanupInPlace(mol, params);
  }
  MS::chargeParentInPlace(mol, params, true);
  MS::isotopeParentInPlace(mol, params, true);
  MS::stereoParentInPlace(mol, params, true);
  tautomer_parent_in_place(mol, params, true, enumerator);
  MS::cleanupInPlace(mol, params);
}

// One op with its parameters, and the objects it shares across molecules.
// Normalize, Reionize, and RemoveFragments share one normalizer, reionizer,
// or remover across the threads of a multi-molecule call, as RDKit's own
// multi-molecule forms do, and the tautomer ops one checked enumerator.
class Standardizer {
 public:
  Standardizer(StandardizeOp op, const CleanupParametersData &data,
               bool skip_standardize)
      : d_op(op),
        d_params(to_cleanup_parameters(data)),
        d_skip(skip_standardize) {
    switch (op) {
      case StandardizeOp::Normalize:
        d_normalizer.reset(MS::normalizerFromParams(d_params));
        break;
      case StandardizeOp::Reionize:
        d_reionizer.reset(MS::reionizerFromParams(d_params));
        break;
      case StandardizeOp::RemoveFragments:
        d_remover.reset(MS::fragmentRemoverFromParams(d_params));
        break;
      default:
        break;
    }
    if (uses_tautomers(op)) {
      d_enumerator = checked_enumerator(d_params);
    }
  }

  // Standardizing removes hydrogens and fragments, splits the molecule into
  // its fragments (getMolFrags), and writes tautomers as SMILES, all of
  // which remove atoms from the molecule or a copy (check_removals); in
  // place, this is checked before the first edit.
  void check(const RWMol &mol) const {
    check_removals(mol);
    if (sums_charges(d_op, d_skip)) {
      check_reionizer_atoms(mol);
    }
  }

  void in_place(RWMol &mol) const {
    switch (d_op) {
      case StandardizeOp::Cleanup:
        MS::cleanupInPlace(mol, d_params);
        return;
      case StandardizeOp::Normalize:
        d_normalizer->normalizeInPlace(mol);
        return;
      case StandardizeOp::Reionize:
        d_reionizer->reionizeInPlace(mol);
        return;
      case StandardizeOp::RemoveFragments:
        d_remover->removeInPlace(mol);
        return;
      case StandardizeOp::CanonicalTautomer:
        canonical_tautomer_in_place(mol, *d_enumerator);
        return;
      case StandardizeOp::TautomerParent:
        tautomer_parent_in_place(mol, d_params, d_skip, *d_enumerator);
        return;
      case StandardizeOp::FragmentParent:
        MS::fragmentParentInPlace(mol, d_params, d_skip);
        return;
      case StandardizeOp::StereoParent:
        MS::stereoParentInPlace(mol, d_params, d_skip);
        return;
      case StandardizeOp::IsotopeParent:
        MS::isotopeParentInPlace(mol, d_params, d_skip);
        return;
      case StandardizeOp::ChargeParent:
        MS::chargeParentInPlace(mol, d_params, d_skip);
        return;
      case StandardizeOp::SuperParent:
        super_parent_in_place(mol, d_params, d_skip, *d_enumerator);
        return;
    }
    throw ValueErrorException("unknown standardization");
  }

  // As the functions in MolStandardize.cpp that return a new molecule.
  // normalize and canonicalTautomer static_cast the ROMol their objects
  // return to RWMol, which it is not, so the objects are called here.
  std::unique_ptr<RWMol> copy(const RWMol &mol) const {
    switch (d_op) {
      case StandardizeOp::Normalize:
        return adopt_mol(d_normalizer->normalize(mol));
      case StandardizeOp::Reionize:
        return adopt_mol(d_reionizer->reionize(mol));
      case StandardizeOp::RemoveFragments:
        return adopt_mol(d_remover->remove(mol));
      case StandardizeOp::CanonicalTautomer:
        check_tautomer_score(mol);
        return adopt_mol(d_enumerator->canonicalize(mol));
      default: {
        auto result = std::make_unique<RWMol>(mol);
        in_place(*result);
        return result;
      }
    }
  }

 private:
  StandardizeOp d_op;
  CleanupParameters d_params;
  bool d_skip;
  std::unique_ptr<MS::Normalizer> d_normalizer;
  std::unique_ptr<MS::Reionizer> d_reionizer;
  std::unique_ptr<MS::FragmentRemover> d_remover;
  std::unique_ptr<MS::TautomerEnumerator> d_enumerator;
};

// Runs work(i) for each i below n on up to num_threads threads, as RDKit's
// getNumThreadsToUse counts them, each i on one thread. If one throws, the
// rest of that thread's share is skipped, and the first exception (by
// thread) is rethrown once every thread has finished. getNumThreadsToUse
// itself negates a negative count in int, undefined for INT_MIN, so the
// count is num_threads_to_use's.
template <class Work>
void run_threads(std::size_t n, int num_threads, const Work &work) {
  std::size_t count = std::min<std::size_t>(n, num_threads_to_use(num_threads));
  if (count <= 1) {
    for (std::size_t i = 0; i < n; ++i) {
      work(i);
    }
    return;
  }
  std::vector<std::exception_ptr> errors(count);
  std::vector<std::thread> threads;
  auto join = [&threads] {
    for (auto &thread : threads) {
      thread.join();
    }
  };
  try {
    for (std::size_t t = 0; t < count; ++t) {
      threads.emplace_back([&, t] {
        try {
          for (std::size_t i = t; i < n; i += count) {
            work(i);
          }
        } catch (...) {
          errors[t] = std::current_exception();
        }
      });
    }
  } catch (...) {
    // A thread that could not start: join the others before unwinding,
    // since destroying a joinable thread calls std::terminate.
    join();
    throw;
  }
  join();
  for (const auto &error : errors) {
    if (error) {
      std::rethrow_exception(error);
    }
  }
}

RWMol &as_rwmol(RDKit::ROMol &mol) {
  // Every molecule Rust lends is an RWMol (CONVENTIONS.md W1).
  auto *rwmol = dynamic_cast<RWMol *>(&mol);
  if (!rwmol) {
    throw ValueErrorException("expected an RWMol");
  }
  return *rwmol;
}

}  // namespace

std::unique_ptr<RWMol> adopt_mol(RDKit::ROMol *mol) {
  std::unique_ptr<RDKit::ROMol> owned(mol);
  if (!owned) {
    throw ValueErrorException("RDKit returned no molecule");
  }
  std::unique_ptr<RWMol> result;
  if (auto *rwmol = dynamic_cast<RWMol *>(owned.get())) {
    owned.release();
    result.reset(rwmol);
  } else {
    result = std::make_unique<RWMol>(*owned);
  }
  clear_matrix_caches(*result);
  return result;
}

namespace {

// readTransformations (TransformCatalogUtils.cpp) on a normalizations file's
// text: lines read 511 characters at a time until one fails, an empty line
// or a "//" comment skipped, and each other's first two tab-separated
// fields, spaces erased, are its name and reaction SMARTS. A line without a
// second field RDKit reads past its tokens for, so it is refused here.
std::vector<std::pair<std::string, std::string>> normalization_lines(
    const std::string &text) {
  std::vector<std::pair<std::string, std::string>> data;
  std::istringstream in(text);
  char line[512];
  while (!in.eof() && !in.fail()) {
    in.getline(line, sizeof(line), '\n');
    std::string tmp(line);
    if (tmp.empty() || tmp.compare(0, 2, "//") == 0) {
      continue;
    }
    std::vector<std::string> fields;
    boost::split(fields, tmp, boost::is_any_of("\t"), boost::token_compress_on);
    fields.erase(std::remove(fields.begin(), fields.end(), ""), fields.end());
    if (fields.size() < 2) {
      throw ValueErrorException(
          "a normalization line needs a name and a reaction SMARTS, separated "
          "by a tab: " +
          tmp);
    }
    boost::erase_all(fields[0], " ");
    boost::erase_all(fields[1], " ");
    check_reaction_cx(fields[1], true);
    data.emplace_back(fields[0], fields[1]);
  }
  return data;
}

}  // namespace

// RDKit reads a normalizations file itself, each line's reaction SMARTS with
// its CXSMILES extension (file_checks.h); with `read_normalizations` (for
// the calls that normalize), the wrapper reads it once instead, checks it,
// and hands RDKit its transformations as normalizationData, so RDKit never
// reads a file that changed after the check. (The fragment, acid-base, and
// tautomer files' readers erase the spaces from their SMARTS, which leaves
// no extension.)
CleanupParameters to_cleanup_parameters(const CleanupParametersData &data,
                                        bool read_normalizations) {
  CleanupParameters params;
  params.rdbase = text(data.rdbase);
  params.normalizations = path(data.normalizations);
  params.acidbaseFile = path(data.acidbase_file);
  params.fragmentFile = path(data.fragment_file);
  params.tautomerTransforms = path(data.tautomer_transforms);
  params.maxRestarts = data.max_restarts;
  params.preferOrganic = data.prefer_organic;
  params.doCanonical = data.do_canonical;
  params.maxTautomers = data.max_tautomers;
  params.maxTransforms = data.max_transforms;
  params.tautomerRemoveSp3Stereo = data.tautomer_remove_sp3_stereo;
  params.tautomerRemoveBondStereo = data.tautomer_remove_bond_stereo;
  params.tautomerRemoveIsotopicHs = data.tautomer_remove_isotopic_hs;
  params.tautomerReassignStereo = data.tautomer_reassign_stereo;
  params.largestFragmentChooserUseAtomCount =
      data.largest_fragment_chooser_use_atom_count;
  params.largestFragmentChooserCountHeavyAtomsOnly =
      data.largest_fragment_chooser_count_heavy_atoms_only;
  // RDKit parses these as they are: the normalizations as reactions, the
  // rest as SMARTS (file_checks.h).
  for (const auto &e : data.normalization_data) {
    params.normalizationData.emplace_back(text(e.name), text(e.first));
    check_reaction_cx(params.normalizationData.back().second, true);
  }
  for (const auto &e : data.fragment_data) {
    params.fragmentData.emplace_back(text(e.name), text(e.first));
    check_smiles_cx(params.fragmentData.back().second, true, true);
  }
  for (const auto &e : data.acidbase_data) {
    params.acidbaseData.emplace_back(text(e.name), text(e.first),
                                     text(e.second));
    check_smiles_cx(std::get<1>(params.acidbaseData.back()), true, true);
    check_smiles_cx(std::get<2>(params.acidbaseData.back()), true, true);
  }
  for (const auto &e : data.tautomer_transform_data) {
    params.tautomerTransformData.emplace_back(text(e.name), text(e.first),
                                              text(e.second), text(e.third));
    check_smiles_cx(std::get<1>(params.tautomerTransformData.back()), true,
                    true);
  }
  // normalizerFromParams reads the file only without normalizationData.
  if (read_normalizations && params.normalizationData.empty() &&
      !params.normalizations.empty()) {
    params.normalizationData =
        normalization_lines(read_file(params.normalizations));
    if (params.normalizationData.empty()) {
      // RDKit's empty catalog from an empty file, not its defaults.
      throw ValueErrorException("the normalizations file " +
                                params.normalizations +
                                " has no normalizations");
    }
  }
  return params;
}

CleanupParametersData from_cleanup_parameters(const CleanupParameters &params) {
  CleanupParametersData data;
  data.rdbase = to_bytes(params.rdbase);
  data.normalizations = to_bytes(params.normalizations);
  data.acidbase_file = to_bytes(params.acidbaseFile);
  data.fragment_file = to_bytes(params.fragmentFile);
  data.tautomer_transforms = to_bytes(params.tautomerTransforms);
  data.max_restarts = params.maxRestarts;
  data.prefer_organic = params.preferOrganic;
  data.do_canonical = params.doCanonical;
  data.max_tautomers = params.maxTautomers;
  data.max_transforms = params.maxTransforms;
  data.tautomer_remove_sp3_stereo = params.tautomerRemoveSp3Stereo;
  data.tautomer_remove_bond_stereo = params.tautomerRemoveBondStereo;
  data.tautomer_remove_isotopic_hs = params.tautomerRemoveIsotopicHs;
  data.tautomer_reassign_stereo = params.tautomerReassignStereo;
  data.largest_fragment_chooser_use_atom_count =
      params.largestFragmentChooserUseAtomCount;
  data.largest_fragment_chooser_count_heavy_atoms_only =
      params.largestFragmentChooserCountHeavyAtomsOnly;
  for (const auto &[name, smarts] : params.normalizationData) {
    data.normalization_data.push_back(entry(name, smarts));
  }
  for (const auto &[name, smarts] : params.fragmentData) {
    data.fragment_data.push_back(entry(name, smarts));
  }
  for (const auto &[name, acid, base] : params.acidbaseData) {
    data.acidbase_data.push_back(entry(name, acid, base));
  }
  for (const auto &[name, smarts, bonds, charges] :
       params.tautomerTransformData) {
    data.tautomer_transform_data.push_back(entry(name, smarts, bonds, charges));
  }
  return data;
}

void check_reionizer_atoms(const RDKit::ROMol &mol) {
  if (mol.getNumAtoms() > INT_MAX / 255) {
    throw ValueErrorException(
        "RDKit's reionizer sums formal charges in int, which overflows for a "
        "molecule of more than " +
        std::to_string(INT_MAX / 255) + " atoms");
  }
}

void check_tautomer_score(const RDKit::ROMol &mol, std::uint64_t rings) {
  std::uint64_t atoms = mol.getNumAtoms();
  std::uint64_t bonds = mol.getNumBonds();
  // 47000 is the default terms' |score|s, summed, times 1000 matches; the
  // directly counted patterns add up to 8 per bond and 1 per atom.
  std::uint64_t bound = 250 * rings + 47000 + 8 * bonds + atoms + 382 * atoms;
  if (bound > static_cast<std::uint64_t>(INT_MAX)) {
    throw ValueErrorException(
        "RDKit adds tautomer scores in int, which overflows for a molecule "
        "with this many atoms, bonds, or rings");
  }
}

void check_tautomer_score(const RDKit::ROMol &mol) {
  // The atoms and bonds alone may break the bound, before the rings cost
  // anything to count.
  check_tautomer_score(mol, 0);
  RDKit::ROMol copy(mol);
  RDKit::MolOps::symmetrizeSSSR(copy);
  check_tautomer_score(mol, copy.getRingInfo()->numRings());
}

std::unique_ptr<MS::TautomerEnumerator> checked_enumerator(
    const CleanupParameters &params) {
  // As TautomerEnumerator(const CleanupParameters &) in Tautomer.cpp.
  std::unique_ptr<MS::TautomerCatalogParams> catalog;
  if (params.tautomerTransformData.empty()) {
    catalog =
        std::make_unique<MS::TautomerCatalogParams>(params.tautomerTransforms);
  } else {
    catalog = std::make_unique<MS::TautomerCatalogParams>(
        params.tautomerTransformData);
  }
  for (const auto &transform : catalog->getTransforms()) {
    if (!transform.Mol) {
      throw ValueErrorException("a tautomer transform has no SMARTS");
    }
    std::string name;
    transform.Mol->getPropIfPresent(RDKit::common_properties::_Name, name);
    if (!transform.BondTypes.empty() &&
        transform.BondTypes.size() < transform.Mol->getNumBonds()) {
      throw ValueErrorException(
          "tautomer transform '" + name + "' has " +
          std::to_string(transform.BondTypes.size()) + " bond types for " +
          std::to_string(transform.Mol->getNumBonds()) + " bonds");
    }
    if (!transform.Charges.empty() &&
        transform.Charges.size() < transform.Mol->getNumAtoms()) {
      throw ValueErrorException(
          "tautomer transform '" + name + "' has " +
          std::to_string(transform.Charges.size()) + " charges for " +
          std::to_string(transform.Mol->getNumAtoms()) + " atoms");
    }
  }
  auto enumerator = std::make_unique<MS::TautomerEnumerator>(
      new MS::TautomerCatalog(catalog.get()));
  enumerator->setMaxTautomers(static_cast<unsigned>(params.maxTautomers));
  enumerator->setMaxTransforms(static_cast<unsigned>(params.maxTransforms));
  enumerator->setRemoveSp3Stereo(params.tautomerRemoveSp3Stereo);
  enumerator->setRemoveBondStereo(params.tautomerRemoveBondStereo);
  enumerator->setRemoveIsotopicHs(params.tautomerRemoveIsotopicHs);
  enumerator->setReassignStereo(params.tautomerReassignStereo);
  return enumerator;
}

CleanupParametersData cleanup_parameters_default() {
  return from_cleanup_parameters(CleanupParameters());
}

CleanupParametersData update_params_from_json(const CleanupParametersData &data,
                                              rust::Str json) {
  auto params = to_cleanup_parameters(data, false);
  MS::updateCleanupParamsFromJSON(params, std::string(json));
  return from_cleanup_parameters(params);
}

std::unique_ptr<RWMol> standardize_copy(StandardizeOp op, const RWMol &mol,
                                        const CleanupParametersData &data,
                                        bool skip_standardize) {
  Standardizer standardizer(op, data, skip_standardize);
  standardizer.check(mol);
  // Standardizing matches substructures (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return standardizer.copy(mol);
}

void standardize_in_place(StandardizeOp op, RWMol &mol,
                          const CleanupParametersData &data,
                          bool skip_standardize) {
  Standardizer standardizer(op, data, skip_standardize);
  standardizer.check(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  standardizer.in_place(mol);
}

// As the multi-molecule forms in MolStandardize.cpp
// (standardizeMultipleMolsInPlace), whose threads stride by numThreads
// rather than by the threads they run, so 0 loops forever and a negative
// count has two threads standardize one molecule at once, and which let an
// exception escape a thread (std::terminate). RDKit's check for a molecule
// listed twice is not needed: Rust lends each one mutably once. Each thread
// holds control_c_mutex while it standardizes a molecule, so the threads
// take turns at RDKit's matching (control_c_lock.h).
void standardize_in_place_mols(StandardizeOp op, BorrowedMols &mols,
                               std::int32_t num_threads,
                               const CleanupParametersData &data,
                               bool skip_standardize) {
  std::vector<RWMol *> list;
  for (const auto &mol : mols.mols()) {
    list.push_back(&as_rwmol(*mol));
  }
  Standardizer standardizer(op, data, skip_standardize);
  for (const auto *mol : list) {
    standardizer.check(*mol);
  }
  run_threads(list.size(), num_threads, [&](std::size_t i) {
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    standardizer.in_place(*list[i]);
  });
}

// As standardizeSmiles in MolStandardize.cpp, with the size check.
std::unique_ptr<std::string> standardize_smiles(rust::Str smiles) {
  std::string input(smiles);
  check_smiles_cx(input, true, true);
  std::unique_ptr<RWMol> mol{RDKit::SmilesToMol(input, 0, false)};
  if (!mol) {
    throw ValueErrorException("SMILES Parse Error: syntax error for input: " +
                              input);
  }
  check_reionizer_atoms(*mol);
  // CXSMILES can give the substance groups PARENT chains (Standardizer).
  check_removals(*mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  MS::cleanupInPlace(*mol);
  return std::make_unique<std::string>(RDKit::MolToSmiles(*mol));
}

}  // namespace rdkit_cxx
