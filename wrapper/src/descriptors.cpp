// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/descriptors.h"

#include <cmath>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/Descriptors/AtomFeat.h>
#include <GraphMol/Descriptors/BCUT.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/OxidationNumbers.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/descriptors.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace Desc = RDKit::Descriptors;

// The VSA descriptors bin each atom with std::upper_bound, which needs the
// bins sorted: for others, and for a NaN, its result is undefined.
std::optional<std::vector<double>> checked_bins(
    rust::Slice<const double> bins) {
  if (bins.empty()) {
    return std::nullopt;
  }
  std::vector<double> result(bins.begin(), bins.end());
  for (std::size_t i = 0; i < result.size(); ++i) {
    if (std::isnan(result[i]) || (i > 0 && result[i] < result[i - 1])) {
      throw ValueErrorException("bins must be sorted, without NaN");
    }
  }
  return result;
}

rust::Vec<double> pair_vec(const std::pair<double, double> &pair) {
  rust::Vec<double> result;
  result.push_back(pair.first);
  result.push_back(pair.second);
  return result;
}

rust::Vec<std::uint32_t> index_vec(const std::vector<unsigned int> &indices) {
  rust::Vec<std::uint32_t> result;
  result.reserve(indices.size());
  for (auto idx : indices) {
    result.push_back(idx);
  }
  return result;
}

// A path of n + 1 atoms visits at least n atoms (its last may close a ring),
// so ChiNv and ChiNn have no paths for n above the number of atoms. RDKit's
// path search runs to the length it is given regardless, and n + 1 wraps
// for the largest n, so those return 0 here without it.
bool chi_has_paths(const RWMol &mol, std::uint32_t n) {
  return n <= mol.getNumAtoms();
}

}  // namespace

void check_property_table_atoms(const RWMol &mol) {
  for (const auto *atom : mol.atoms()) {
    auto num = atom->getAtomicNum();
    if (num < 1 || num > 110) {
      throw ValueErrorException(
          "these descriptors' property tables cover atomic numbers 1 to 110; "
          "atom " +
          std::to_string(atom->getIdx()) + " has atomic number " +
          std::to_string(num));
    }
  }
}

// Crippen's contributions type the atoms by substructure matches, as do the
// descriptors below that read them (SlogP_VSA, SMR_VSA, BCUT2D), the counts
// Lipinski.cpp makes from SMARTS, and the MQNs, which include those counts
// (control_c_lock.h).
rust::Vec<double> calc_crippen_descriptors(const RWMol &mol, bool include_hs) {
  DescriptorCacheGuard guard(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  double logp = 0.0;
  double mr = 0.0;
  Desc::calcCrippenDescriptors(mol, logp, mr, include_hs);
  return pair_vec({logp, mr});
}

double calc_labute_asa(const RWMol &mol, bool include_hs) {
  DescriptorCacheGuard guard(mol);
  return Desc::calcLabuteASA(mol, include_hs);
}

double calc_tpsa(const RWMol &mol, bool include_s_and_p) {
  DescriptorCacheGuard guard(mol);
  return Desc::calcTPSA(mol, false, include_s_and_p);
}

double calc_exact_mol_wt(const RWMol &mol, bool only_heavy) {
  return Desc::calcExactMW(mol, only_heavy);
}

double calc_amw(const RWMol &mol, bool only_heavy) {
  return Desc::calcAMW(mol, only_heavy);
}

// Python's _CalcLabuteASAContribs: the implicit hydrogens' contribution,
// then each atom's.
rust::Vec<double> calc_labute_asa_contribs(const RWMol &mol, bool include_hs) {
  DescriptorCacheGuard guard(mol);
  std::vector<double> contribs(mol.getNumAtoms());
  double h_contrib = 0.0;
  Desc::getLabuteAtomContribs(mol, contribs, h_contrib, include_hs, true);
  rust::Vec<double> result;
  result.reserve(contribs.size() + 1);
  result.push_back(h_contrib);
  for (auto contrib : contribs) {
    result.push_back(contrib);
  }
  return result;
}

// Python's _CalcCrippenContribs: each atom's logP and MR contributions, in
// pairs. The atoms are typed by substructure matches, as in
// calc_crippen_descriptors.
rust::Vec<double> calc_crippen_contribs(const RWMol &mol) {
  DescriptorCacheGuard guard(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  std::vector<double> logp(mol.getNumAtoms());
  std::vector<double> mr(mol.getNumAtoms());
  Desc::getCrippenAtomContribs(mol, logp, mr, true);
  rust::Vec<double> result;
  result.reserve(2 * logp.size());
  for (std::size_t i = 0; i < logp.size(); ++i) {
    result.push_back(logp[i]);
    result.push_back(mr[i]);
  }
  return result;
}

// getMolFormula (MolProps.cpp) sums the formal charges in int and negates a
// negative sum (U9).
std::unique_ptr<std::string> calc_mol_formula(const RWMol &mol,
                                              bool separate_isotopes,
                                              bool abbreviate_h_isotopes) {
  check_formal_charge_sum(mol);
  return std::make_unique<std::string>(
      Desc::calcMolFormula(mol, separate_isotopes, abbreviate_h_isotopes));
}

#define RDKIT_CXX_COUNT(name, function) \
  std::uint32_t name(const RWMol &mol) { return function(mol); }
#define RDKIT_CXX_COUNT_MATCHES(name, function)               \
  std::uint32_t name(const RWMol &mol) {                      \
    std::lock_guard<std::mutex> control_c(control_c_mutex()); \
    return function(mol);                                     \
  }

RDKIT_CXX_COUNT(calc_num_lipinski_hbd, Desc::calcLipinskiHBD)
RDKIT_CXX_COUNT(calc_num_lipinski_hba, Desc::calcLipinskiHBA)
RDKIT_CXX_COUNT_MATCHES(calc_num_hbd, Desc::calcNumHBD)
RDKIT_CXX_COUNT_MATCHES(calc_num_hba, Desc::calcNumHBA)
RDKIT_CXX_COUNT(calc_num_rings, Desc::calcNumRings)
RDKIT_CXX_COUNT(calc_num_aromatic_rings, Desc::calcNumAromaticRings)
RDKIT_CXX_COUNT(calc_num_saturated_rings, Desc::calcNumSaturatedRings)
RDKIT_CXX_COUNT(calc_num_heterocycles, Desc::calcNumHeterocycles)
RDKIT_CXX_COUNT(calc_num_aromatic_heterocycles,
                Desc::calcNumAromaticHeterocycles)
RDKIT_CXX_COUNT(calc_num_aromatic_carbocycles, Desc::calcNumAromaticCarbocycles)
RDKIT_CXX_COUNT(calc_num_saturated_heterocycles,
                Desc::calcNumSaturatedHeterocycles)
RDKIT_CXX_COUNT(calc_num_saturated_carbocycles,
                Desc::calcNumSaturatedCarbocycles)
RDKIT_CXX_COUNT(calc_num_aliphatic_rings, Desc::calcNumAliphaticRings)
RDKIT_CXX_COUNT(calc_num_aliphatic_heterocycles,
                Desc::calcNumAliphaticHeterocycles)
RDKIT_CXX_COUNT(calc_num_aliphatic_carbocycles,
                Desc::calcNumAliphaticCarbocycles)
RDKIT_CXX_COUNT(calc_num_heavy_atoms, Desc::calcNumHeavyAtoms)
RDKIT_CXX_COUNT(calc_num_atoms, Desc::calcNumAtoms)
RDKIT_CXX_COUNT_MATCHES(calc_num_heteroatoms, Desc::calcNumHeteroatoms)
RDKIT_CXX_COUNT_MATCHES(calc_num_amide_bonds, Desc::calcNumAmideBonds)
RDKIT_CXX_COUNT(calc_num_atom_stereo_centers, Desc::numAtomStereoCenters)
RDKIT_CXX_COUNT(calc_num_unspecified_atom_stereo_centers,
                Desc::numUnspecifiedAtomStereoCenters)

#undef RDKIT_CXX_COUNT
#undef RDKIT_CXX_COUNT_MATCHES

std::uint32_t calc_num_rotatable_bonds(const RWMol &mol,
                                       NumRotatableBondsOptions strict) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return Desc::calcNumRotatableBonds(mol, strict);
}

double calc_fraction_csp3(const RWMol &mol) {
  return Desc::calcFractionCSP3(mol);
}

std::uint32_t calc_num_spiro_atoms(const RWMol &mol,
                                   rust::Vec<std::uint32_t> &atoms) {
  std::vector<unsigned int> found;
  auto count = Desc::calcNumSpiroAtoms(mol, &found);
  atoms = index_vec(found);
  return count;
}

std::uint32_t calc_num_bridgehead_atoms(const RWMol &mol,
                                        rust::Vec<std::uint32_t> &atoms) {
  std::vector<unsigned int> found;
  auto count = Desc::calcNumBridgeheadAtoms(mol, &found);
  atoms = index_vec(found);
  return count;
}

// The VSA descriptors read the Labute and Crippen contributions through
// their caches (U8).
rust::Vec<double> slogp_vsa(const RWMol &mol, rust::Slice<const double> bins) {
  auto lbins = checked_bins(bins);
  DescriptorCacheGuard guard(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return to_vec(Desc::calcSlogP_VSA(mol, lbins ? &*lbins : nullptr));
}

rust::Vec<double> smr_vsa(const RWMol &mol, rust::Slice<const double> bins) {
  auto lbins = checked_bins(bins);
  DescriptorCacheGuard guard(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return to_vec(Desc::calcSMR_VSA(mol, lbins ? &*lbins : nullptr));
}

rust::Vec<double> peoe_vsa(const RWMol &mol, rust::Slice<const double> bins) {
  auto lbins = checked_bins(bins);
  DescriptorCacheGuard guard(mol);
  return to_vec(Desc::calcPEOE_VSA(mol, lbins ? &*lbins : nullptr));
}

rust::Vec<double> custom_prop_vsa(const RWMol &mol, rust::Str custom_prop_name,
                                  rust::Slice<const double> bins) {
  auto lbins = checked_bins(bins);
  DescriptorCacheGuard guard(mol);
  return to_vec(
      Desc::calcCustomProp_VSA(mol, std::string(custom_prop_name),
                               lbins ? *lbins : std::vector<double>()));
}

rust::Vec<std::uint32_t> mqns(const RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return to_vec(Desc::calcMQNs(mol));
}

// The Chi descriptors read their per-atom deltas through the cache, by atom,
// without checking the cached vector's length (U8).
#define RDKIT_CXX_CHI(name, function) \
  double name(const RWMol &mol) {     \
    DescriptorCacheGuard guard(mol);  \
    return function(mol);             \
  }

RDKIT_CXX_CHI(calc_chi0v, Desc::calcChi0v)
RDKIT_CXX_CHI(calc_chi1v, Desc::calcChi1v)
RDKIT_CXX_CHI(calc_chi2v, Desc::calcChi2v)
RDKIT_CXX_CHI(calc_chi3v, Desc::calcChi3v)
RDKIT_CXX_CHI(calc_chi4v, Desc::calcChi4v)
RDKIT_CXX_CHI(calc_chi0n, Desc::calcChi0n)
RDKIT_CXX_CHI(calc_chi1n, Desc::calcChi1n)
RDKIT_CXX_CHI(calc_chi2n, Desc::calcChi2n)
RDKIT_CXX_CHI(calc_chi3n, Desc::calcChi3n)
RDKIT_CXX_CHI(calc_chi4n, Desc::calcChi4n)

#undef RDKIT_CXX_CHI

double calc_chi_nv(const RWMol &mol, std::uint32_t n) {
  if (!chi_has_paths(mol, n)) {
    return 0.0;
  }
  DescriptorCacheGuard guard(mol);
  return Desc::calcChiNv(mol, n);
}

double calc_chi_nn(const RWMol &mol, std::uint32_t n) {
  if (!chi_has_paths(mol, n)) {
    return 0.0;
  }
  DescriptorCacheGuard guard(mol);
  return Desc::calcChiNn(mol, n);
}

double calc_hall_kier_alpha(const RWMol &mol,
                            rust::Vec<double> &atom_contribs) {
  std::vector<double> contribs(mol.getNumAtoms());
  auto alpha = Desc::calcHallKierAlpha(mol, &contribs);
  atom_contribs = to_vec(contribs);
  return alpha;
}

double calc_kappa1(const RWMol &mol) { return Desc::calcKappa1(mol); }
double calc_kappa2(const RWMol &mol) { return Desc::calcKappa2(mol); }
double calc_kappa3(const RWMol &mol) { return Desc::calcKappa3(mol); }
double calc_phi(const RWMol &mol) { return Desc::calcPhi(mol); }

rust::Vec<double> get_atom_features(const RWMol &mol, std::uint32_t atomid,
                                    bool addchiral) {
  check_index(atomid, mol.getNumAtoms());
  std::vector<double> result;
  Desc::AtomFeatVect(mol, result, static_cast<int>(atomid), addchiral);
  return to_vec(result);
}

void calc_oxidation_numbers(RWMol &mol) { Desc::calcOxidationNumbers(mol); }

// AUTOCORR2D reads the topological distance matrix (U6) and, without a
// custom property, the property tables.
rust::Vec<double> calc_autocorr_2d(const RWMol &mol,
                                   rust::Str custom_atom_property) {
  std::string property(custom_atom_property);
  check_matrix_size(mol);
  if (property.empty()) {
    check_property_table_atoms(mol);
  }
  MatrixCacheGuard guard(mol);
  std::vector<double> result;
  Desc::AUTOCORR2D(mol, result, property);
  return to_vec(result);
}

// BCUT2D(mol) reads Crippen contributions through their cache, on a copy
// that keeps the molecule's computed properties (U8), and makes that copy
// with removeAllHs (check_removals).
rust::Vec<double> bcut_2d(const RWMol &mol) {
  check_removals(mol);
  DescriptorCacheGuard guard(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return to_vec(Desc::BCUT2D(mol));
}

rust::Vec<double> bcut_2d_atom_props(const RWMol &mol,
                                     rust::Slice<const double> atom_props) {
  std::vector<double> props(atom_props.begin(), atom_props.end());
  return pair_vec(Desc::BCUT2D(mol, props));
}

rust::Vec<double> bcut_2d_atom_propname(const RWMol &mol,
                                        rust::Str atom_propname) {
  return pair_vec(Desc::BCUT2D(mol, std::string(atom_propname)));
}

// MACCS.cpp makes its patterns on first use and assigns them to a global
// flyweight without a lock, so two first calls race (U4). The first call is
// made once here, on an empty molecule, which makes the patterns and returns;
// later calls only read them. Every call matches the patterns
// (control_c_lock.h).
std::unique_ptr<BitVect> get_maccs_keys_fingerprint(const RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  static std::once_flag patterns_made;
  std::call_once(patterns_made, [] {
    RDKit::ROMol empty;
    std::unique_ptr<ExplicitBitVect> fp(
        RDKit::MACCSFingerprints::getFingerprintAsBitVect(empty));
  });
  return wrap_bit_vect(std::unique_ptr<ExplicitBitVect>(
      RDKit::MACCSFingerprints::getFingerprintAsBitVect(mol)));
}

std::uint32_t get_atom_pair_atom_code(const Atom &atom,
                                      std::uint32_t branch_subtract,
                                      bool include_chirality) {
  // getAtomCode may run the CIP labeler on the atom's molecule.
  std::unique_lock<std::mutex> control_c;
  if (include_chirality) {
    control_c = lock_control_c_for_cip();
  }
  return RDKit::AtomPairs::getAtomCode(&atom, branch_subtract,
                                       include_chirality);
}

// The codes shift by at most 16 bits, so any 32-bit codes are defined.
std::uint32_t get_atom_pair_code(std::uint32_t atom1_code,
                                 std::uint32_t atom2_code,
                                 std::uint32_t distance,
                                 bool include_chirality) {
  return RDKit::AtomPairs::getAtomPairCode(atom1_code, atom2_code, distance,
                                           include_chirality);
}

rust::Vec<std::uint32_t> get_connectivity_invariants(
    const RWMol &mol, bool include_ring_membership) {
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  RDKit::MorganFingerprints::getConnectivityInvariants(mol, invars,
                                                       include_ring_membership);
  return to_vec(invars);
}

// The feature invariants come from substructure matches (control_c_lock.h).
rust::Vec<std::uint32_t> get_feature_invariants(const RWMol &mol) {
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  RDKit::MorganFingerprints::getFeatureInvariants(mol, invars);
  return to_vec(invars);
}

// As atomPairTypes in Code/GraphMol/Descriptors/nbWrap/rdMolDescriptors.cpp:
// the whole table, whose last entry is unused.
AtomPairsConstants atom_pairs_constants() {
  namespace AP = RDKit::AtomPairs;
  AtomPairsConstants constants;
  constants.num_type_bits = AP::numTypeBits;
  constants.num_pi_bits = AP::numPiBits;
  constants.num_branch_bits = AP::numBranchBits;
  constants.num_chiral_bits = AP::numChiralBits;
  constants.code_size = AP::codeSize;
  constants.num_path_bits = AP::numPathBits;
  constants.num_atom_pair_fingerprint_bits = AP::numAtomPairFingerprintBits;
  for (auto type : AP::atomNumberTypes) {
    constants.atom_types.push_back(type);
  }
  return constants;
}

std::unique_ptr<std::string> atom_pairs_version() {
  return std::make_unique<std::string>(RDKit::AtomPairs::atomPairsVersion);
}

}  // namespace rdkit_cxx
