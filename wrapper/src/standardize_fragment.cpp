// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_fragment.h"

#include <mutex>
#include <sstream>
#include <string>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/standardize.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;

namespace {

MS::LargestFragmentChooser chooser(bool prefer_organic, bool use_atom_count,
                                   bool count_heavy_atoms_only) {
  MS::CleanupParameters params;
  params.preferOrganic = prefer_organic;
  params.largestFragmentChooserUseAtomCount = use_atom_count;
  params.largestFragmentChooserCountHeavyAtomsOnly = count_heavy_atoms_only;
  return MS::LargestFragmentChooser(params);
}

// The remover's methods do not change it, but are not const.
FragmentRemover &mut(const FragmentRemover &remover) {
  return const_cast<FragmentRemover &>(remover);
}

}  // namespace

std::unique_ptr<FragmentRemover> fragment_remover_default() {
  return std::make_unique<FragmentRemover>();
}

std::unique_ptr<FragmentRemover> fragment_remover_new(
    rust::Slice<const std::uint8_t> fragment_filename, bool leave_last,
    bool skip_if_all_match) {
  // RDKit's fragment reader parses no CXSMILES extension: it erases the
  // spaces from the SMARTS it splits at tabs (read_transformation_lines in
  // file_checks.h), so the file needs no check.
  return std::make_unique<FragmentRemover>(to_path(fragment_filename),
                                           leave_last, skip_if_all_match);
}

std::unique_ptr<FragmentRemover> fragment_remover_from_data(
    rust::Str fragment_data, bool leave_last, bool skip_if_all_match) {
  // As fragment_remover_new, the data needs no check.
  std::istringstream in{std::string(fragment_data)};
  return std::make_unique<FragmentRemover>(in, leave_last, skip_if_all_match);
}

// The remover matches its fragments (control_c_lock.h). Both remove atoms
// (check_removals), in place after this check.
std::unique_ptr<RWMol> fragment_remover_remove(const FragmentRemover &remover,
                                               const RWMol &mol) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return adopt_mol(mut(remover).remove(mol));
}

void fragment_remover_remove_in_place(const FragmentRemover &remover,
                                      RWMol &mol) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  mut(remover).removeInPlace(mol);
}

std::unique_ptr<RWMol> largest_fragment_chooser_choose(
    const RWMol &mol, bool prefer_organic, bool use_atom_count,
    bool count_heavy_atoms_only) {
  check_removals(mol);
  return adopt_mol(
      chooser(prefer_organic, use_atom_count, count_heavy_atoms_only)
          .choose(mol));
}

void largest_fragment_chooser_choose_in_place(RWMol &mol, bool prefer_organic,
                                              bool use_atom_count,
                                              bool count_heavy_atoms_only) {
  check_removals(mol);
  chooser(prefer_organic, use_atom_count, count_heavy_atoms_only)
      .chooseInPlace(mol);
}

}  // namespace rdkit_cxx
