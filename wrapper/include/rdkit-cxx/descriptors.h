/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolDescriptors's 2D descriptors. The wrappers whose RDKit calls read or
// write a descriptor cache clear the molecule's caches before and after the
// call (CONVENTIONS.md U8), so Python's force arguments are not bound.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/Descriptors/Lipinski.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using NumRotatableBondsOptions = RDKit::Descriptors::NumRotatableBondsOptions;

// Shared structs, defined by the bridge (src/ffi/descriptors.rs).
struct AtomPairsConstants;

// Crippen's (logP, MR).
rust::Vec<double> calc_crippen_descriptors(const RWMol &mol, bool include_hs);
double calc_labute_asa(const RWMol &mol, bool include_hs);
double calc_tpsa(const RWMol &mol, bool include_s_and_p);
double calc_exact_mol_wt(const RWMol &mol, bool only_heavy);
double calc_amw(const RWMol &mol, bool only_heavy);
rust::Vec<double> calc_labute_asa_contribs(const RWMol &mol, bool include_hs);
rust::Vec<double> calc_crippen_contribs(const RWMol &mol);
std::unique_ptr<std::string> calc_mol_formula(const RWMol &mol,
                                              bool separate_isotopes,
                                              bool abbreviate_h_isotopes);

std::uint32_t calc_num_lipinski_hbd(const RWMol &mol);
std::uint32_t calc_num_lipinski_hba(const RWMol &mol);
std::uint32_t calc_num_hbd(const RWMol &mol);
std::uint32_t calc_num_hba(const RWMol &mol);
std::uint32_t calc_num_rotatable_bonds(const RWMol &mol,
                                       NumRotatableBondsOptions strict);
std::uint32_t calc_num_rings(const RWMol &mol);
std::uint32_t calc_num_aromatic_rings(const RWMol &mol);
std::uint32_t calc_num_saturated_rings(const RWMol &mol);
std::uint32_t calc_num_heterocycles(const RWMol &mol);
std::uint32_t calc_num_aromatic_heterocycles(const RWMol &mol);
std::uint32_t calc_num_aromatic_carbocycles(const RWMol &mol);
std::uint32_t calc_num_saturated_heterocycles(const RWMol &mol);
std::uint32_t calc_num_saturated_carbocycles(const RWMol &mol);
std::uint32_t calc_num_aliphatic_rings(const RWMol &mol);
std::uint32_t calc_num_aliphatic_heterocycles(const RWMol &mol);
std::uint32_t calc_num_aliphatic_carbocycles(const RWMol &mol);
std::uint32_t calc_num_heavy_atoms(const RWMol &mol);
std::uint32_t calc_num_atoms(const RWMol &mol);
std::uint32_t calc_num_heteroatoms(const RWMol &mol);
std::uint32_t calc_num_amide_bonds(const RWMol &mol);
double calc_fraction_csp3(const RWMol &mol);
std::uint32_t calc_num_atom_stereo_centers(const RWMol &mol);
std::uint32_t calc_num_unspecified_atom_stereo_centers(const RWMol &mol);
// The atoms go to `atoms`, as Python's atoms list.
std::uint32_t calc_num_spiro_atoms(const RWMol &mol,
                                   rust::Vec<std::uint32_t> &atoms);
std::uint32_t calc_num_bridgehead_atoms(const RWMol &mol,
                                        rust::Vec<std::uint32_t> &atoms);

// Empty bins are RDKit's defaults, as Python's empty lists are.
rust::Vec<double> slogp_vsa(const RWMol &mol, rust::Slice<const double> bins);
rust::Vec<double> smr_vsa(const RWMol &mol, rust::Slice<const double> bins);
rust::Vec<double> peoe_vsa(const RWMol &mol, rust::Slice<const double> bins);
rust::Vec<double> custom_prop_vsa(const RWMol &mol, rust::Str custom_prop_name,
                                  rust::Slice<const double> bins);
rust::Vec<std::uint32_t> mqns(const RWMol &mol);

double calc_chi0v(const RWMol &mol);
double calc_chi1v(const RWMol &mol);
double calc_chi2v(const RWMol &mol);
double calc_chi3v(const RWMol &mol);
double calc_chi4v(const RWMol &mol);
double calc_chi_nv(const RWMol &mol, std::uint32_t n);
double calc_chi0n(const RWMol &mol);
double calc_chi1n(const RWMol &mol);
double calc_chi2n(const RWMol &mol);
double calc_chi3n(const RWMol &mol);
double calc_chi4n(const RWMol &mol);
double calc_chi_nn(const RWMol &mol, std::uint32_t n);
// Each atom's contribution goes to `atom_contribs`, as Python's list.
double calc_hall_kier_alpha(const RWMol &mol, rust::Vec<double> &atom_contribs);
double calc_kappa1(const RWMol &mol);
double calc_kappa2(const RWMol &mol);
double calc_kappa3(const RWMol &mol);
double calc_phi(const RWMol &mol);

rust::Vec<double> get_atom_features(const RWMol &mol, std::uint32_t atomid,
                                    bool addchiral);
void calc_oxidation_numbers(RWMol &mol);
rust::Vec<double> calc_autocorr_2d(const RWMol &mol,
                                   rust::Str custom_atom_property);
rust::Vec<double> bcut_2d(const RWMol &mol);
// The (high, low) eigenvalues.
rust::Vec<double> bcut_2d_atom_props(const RWMol &mol,
                                     rust::Slice<const double> atom_props);
rust::Vec<double> bcut_2d_atom_propname(const RWMol &mol,
                                        rust::Str atom_propname);

// The fingerprints and atom codes rdMolDescriptors binds that RDKit has not
// deprecated.
std::unique_ptr<BitVect> get_maccs_keys_fingerprint(const RWMol &mol);
std::uint32_t get_atom_pair_atom_code(const Atom &atom,
                                      std::uint32_t branch_subtract,
                                      bool include_chirality);
std::uint32_t get_atom_pair_code(std::uint32_t atom1_code,
                                 std::uint32_t atom2_code,
                                 std::uint32_t distance,
                                 bool include_chirality);
rust::Vec<std::uint32_t> get_connectivity_invariants(
    const RWMol &mol, bool include_ring_membership);
rust::Vec<std::uint32_t> get_feature_invariants(const RWMol &mol);
// AtomPairs.h's constants, as AtomPairsParameters's properties.
AtomPairsConstants atom_pairs_constants();
std::unique_ptr<std::string> atom_pairs_version();

// For the other descriptor wrappers. MolData3Ddescriptors indexes its
// property tables by atomic number - 1, and they cover 1 to 110: a dummy
// atom reads before the start, and an atom past darmstadtium past the end.
// Throws ValueErrorException for any other atom.
void check_property_table_atoms(const RWMol &mol);

}  // namespace rdkit_cxx
