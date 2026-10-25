// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_fragment.h"

#include <sstream>
#include <string>

#include "rdkit-cxx/errors.h"
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
  return std::make_unique<FragmentRemover>(to_path(fragment_filename),
                                           leave_last, skip_if_all_match);
}

std::unique_ptr<FragmentRemover> fragment_remover_from_data(
    rust::Str fragment_data, bool leave_last, bool skip_if_all_match) {
  std::istringstream in{std::string(fragment_data)};
  return std::make_unique<FragmentRemover>(in, leave_last, skip_if_all_match);
}

std::unique_ptr<RWMol> fragment_remover_remove(const FragmentRemover &remover,
                                               const RWMol &mol) {
  return adopt_mol(mut(remover).remove(mol));
}

void fragment_remover_remove_in_place(const FragmentRemover &remover,
                                      RWMol &mol) {
  mut(remover).removeInPlace(mol);
}

std::unique_ptr<RWMol> largest_fragment_chooser_choose(
    const RWMol &mol, bool prefer_organic, bool use_atom_count,
    bool count_heavy_atoms_only) {
  return adopt_mol(
      chooser(prefer_organic, use_atom_count, count_heavy_atoms_only)
          .choose(mol));
}

void largest_fragment_chooser_choose_in_place(RWMol &mol, bool prefer_organic,
                                              bool use_atom_count,
                                              bool count_heavy_atoms_only) {
  chooser(prefer_organic, use_atom_count, count_heavy_atoms_only)
      .chooseInPlace(mol);
}

}  // namespace rdkit_cxx
