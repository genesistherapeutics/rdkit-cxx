// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/enumerate.h"

#include <climits>
#include <limits>
#include <mutex>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <GraphMol/ChemReactions/Enumerate/CartesianProduct.h>
#include <GraphMol/ChemReactions/Enumerate/EnumerationPickler.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSample.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSampleAllBBs.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chem_reactions.rs.h"
#include "rdkit-cxx/ffi/enumerate.rs.h"
#include "rdkit-cxx/mol_checks.h"
// RDKit builds EvenSamplePairsStrategy into its library but does not install
// its header; this is RDKit 2025.03's (CMake checks the release series).
#include "rdkit-cxx/vendor/EvenSamplePairs.h"

namespace rdkit_cxx {

namespace {

using RDKit::MOL_SPTR_VECT;
using RDKit::ROMOL_SPTR;
using RDKit::EnumerationTypes::BBS;
using RDKit::EnumerationTypes::RGROUPS;
using StrategyPtr = boost::shared_ptr<RDKit::EnumerationStrategyBase>;

constexpr std::uint64_t kOverflow =
    RDKit::EnumerationStrategyBase::EnumerationOverflow;

// The strategies' and the library's protected members, reached through
// pointers to members a subclass may form.
struct StrategyAccess : RDKit::EnumerationStrategyBase {
  static const RGROUPS &permutation(const RDKit::EnumerationStrategyBase &s) {
    return s.*(&StrategyAccess::m_permutation);
  }
  static const RGROUPS &sizes(const RDKit::EnumerationStrategyBase &s) {
    return s.*(&StrategyAccess::m_permutationSizes);
  }
};

struct LibraryAccess : RDKit::EnumerateLibraryBase {
  static StrategyPtr &enumerator(RDKit::EnumerateLibraryBase &lib) {
    return lib.*(&LibraryAccess::m_enumerator);
  }
  static const StrategyPtr &enumerator(const RDKit::EnumerateLibraryBase &lib) {
    return lib.*(&LibraryAccess::m_enumerator);
  }
  static const StrategyPtr &initial(const RDKit::EnumerateLibraryBase &lib) {
    return lib.*(&LibraryAccess::m_initialEnumerator);
  }
};

// EvenSamplePairsStrategy's state is private. An explicit instantiation may
// name private members ([temp.spec.general]), so Steal<Tag, &C::member>
// defines a friend that returns the pointer to the member, found through
// its tag.
template <class Tag, typename Tag::type Member>
struct Steal {
  friend typename Tag::type member(Tag) { return Member; }
};

#define RDKIT_CXX_EVEN_MEMBER(Tag, Type, name)           \
  struct Tag {                                           \
    using type = Type RDKit::EvenSamplePairsStrategy::*; \
    friend type member(Tag);                             \
  };                                                     \
  template struct Steal<Tag, &RDKit::EvenSamplePairsStrategy::name>;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
RDKIT_CXX_EVEN_MEMBER(UsedCount, std::vector<boost::int64_t>, used_count)
RDKIT_CXX_EVEN_MEMBER(VarUsed, std::vector<std::vector<boost::uint64_t>>,
                      var_used)
RDKIT_CXX_EVEN_MEMBER(PairUsed, std::vector<std::vector<boost::uint64_t>>,
                      pair_used)
RDKIT_CXX_EVEN_MEMBER(PairCounts, std::vector<std::vector<boost::uint64_t>>,
                      pair_counts)
RDKIT_CXX_EVEN_MEMBER(Modulus, boost::uint64_t, M)
#pragma GCC diagnostic pop
#undef RDKIT_CXX_EVEN_MEMBER

template <class Tag>
const auto &even(const RDKit::EvenSamplePairsStrategy &s) {
  return s.*member(Tag());
}

// The random strategies draw from uniform_int_distribution<int>(0, size -
// 1): an empty group is an empty range, which boost does not allow (and
// RandomSampleAllBBsStrategy takes it modulo), and a group over 2^31
// overflows int.
void check_random_sizes(const RGROUPS &sizes, const char *name) {
  for (auto size : sizes) {
    if (size == 0 || size - 1 > static_cast<std::uint64_t>(INT_MAX)) {
      throw ValueErrorException(
          std::string(name) +
          " needs between 1 and 2^31 building blocks in every group");
    }
  }
}

// What a strategy needs of the building blocks it is initialized with.
void check_initialize(const RDKit::EnumerationStrategyBase &strategy,
                      const RGROUPS &sizes) {
  if (dynamic_cast<const RDKit::RandomSampleAllBBsStrategy *>(&strategy)) {
    // Its initializeStrategy dereferences max_element of the sizes.
    if (sizes.empty()) {
      throw ValueErrorException(
          "RandomSampleAllBBsStrategy needs a group of building blocks");
    }
    check_random_sizes(sizes, "RandomSampleAllBBsStrategy");
  } else if (dynamic_cast<const RDKit::RandomSampleStrategy *>(&strategy)) {
    check_random_sizes(sizes, "RandomSampleStrategy");
  } else if (dynamic_cast<const RDKit::EvenSamplePairsStrategy *>(&strategy)) {
    // Its initializeStrategy doubles a uint64 modulus until it reaches the
    // number of permutations, which never ends past 2^63.
    auto n = RDKit::computeNumProducts(sizes);
    if (n != kOverflow && n > (std::uint64_t{1} << 63)) {
      throw ValueErrorException(
          "EvenSamplePairsStrategy takes at most 2^63 permutations");
    }
  }
}

void bad_state(const std::string &what) {
  throw ValueErrorException("bad enumeration state: " + what);
}

void check_shape(const std::vector<std::vector<boost::uint64_t>> &rows,
                 const std::vector<std::uint64_t> &lengths, const char *name) {
  if (rows.size() != lengths.size()) {
    bad_state(name);
  }
  for (std::size_t i = 0; i < rows.size(); ++i) {
    if (rows[i].size() != lengths[i]) {
      bad_state(name);
    }
  }
}

// A strategy read from a pickle must be one RDKit could have reached by
// enumerating the library's building blocks, `sizes`: the strategies index
// their state, and the building blocks, with it unchecked. The pickles are
// Boost.Serialization archives, whose framing Boost trusts, so the Rust
// callers take only archives RDKit wrote (`unsafe fn`, CONVENTIONS.md U2);
// these checks catch an archive written for another library.
void check_state(const RDKit::EnumerationStrategyBase &strategy,
                 const RGROUPS &sizes) {
  const auto &position = StrategyAccess::permutation(strategy);
  if (StrategyAccess::sizes(strategy) != sizes) {
    bad_state("its building block counts are not the library's");
  }
  if (position.size() != sizes.size()) {
    bad_state("its position has the wrong length");
  }
  auto n = strategy.getNumPermutations();
  if (n != RDKit::computeNumProducts(sizes)) {
    bad_state("its permutation count is wrong");
  }
  if (const auto *cartesian =
          dynamic_cast<const RDKit::CartesianProductStrategy *>(&strategy)) {
    // CartesianProductStrategy::next carries from each group to the next,
    // past the last one, unless its position is the one its count of steps
    // reaches.
    auto processed = cartesian->getPermutationIdx();
    if (n == kOverflow || processed < n) {
      std::uint64_t step = processed ? processed - 1 : 0;
      for (std::size_t i = 0; i < sizes.size(); ++i) {
        if (position[i] != step % sizes[i]) {
          bad_state("its position does not match its count");
        }
        step /= sizes[i];
      }
    }
  } else if (dynamic_cast<const RDKit::RandomSampleAllBBsStrategy *>(
                 &strategy)) {
    check_random_sizes(sizes, "RandomSampleAllBBsStrategy");
  } else if (dynamic_cast<const RDKit::RandomSampleStrategy *>(&strategy)) {
    check_random_sizes(sizes, "RandomSampleStrategy");
  } else if (const auto *pairs =
                 dynamic_cast<const RDKit::EvenSamplePairsStrategy *>(
                     &strategy)) {
    // try_add indexes these by the decoded building blocks, and next takes
    // seeds modulo M.
    if (n == kOverflow) {
      bad_state("too many permutations");
    }
    std::uint64_t monomers = 0;
    for (auto size : sizes) {
      monomers += size;
    }
    if (even<UsedCount>(*pairs).size() != sizes.size()) {
      bad_state("used_count");
    }
    check_shape(even<VarUsed>(*pairs), sizes, "var_used");
    check_shape(even<PairUsed>(*pairs),
                std::vector<std::uint64_t>(monomers, monomers), "pair_used");
    check_shape(even<PairCounts>(*pairs),
                std::vector<std::uint64_t>(sizes.size(), sizes.size()),
                "pair_counts");
    // try_add casts a pair's count over a double to int.
    for (const auto &row : even<PairUsed>(*pairs)) {
      for (auto count : row) {
        if (count > static_cast<std::uint64_t>(INT_MAX)) {
          bad_state("pair_used");
        }
      }
    }
    if (even<Modulus>(*pairs) == 0) {
      bad_state("a zero modulus");
    }
  } else {
    bad_state("an unknown strategy");
  }
}

// Groups the lent molecules: copies the library owns, or the molecules
// themselves, lent for one call.
BBS to_bbs(BorrowedMols &mols, rust::Slice<const std::uint32_t> group_sizes,
           bool copy) {
  BBS bbs;
  std::size_t next = 0;
  for (auto size : group_sizes) {
    MOL_SPTR_VECT group;
    for (std::uint32_t j = 0; j < size; ++j) {
      const auto &mol = mols.mols().at(next++);
      if (copy) {
        auto *owned = new RWMol(*mol);
        clear_matrix_caches(*owned);
        group.emplace_back(owned);
      } else {
        group.push_back(mol);
      }
    }
    bbs.push_back(std::move(group));
  }
  if (next != mols.mols().size()) {
    throw std::logic_error("the building block groups do not cover the list");
  }
  return bbs;
}

RDKit::EnumerationParams params_from(const EnumerationParamsData &data) {
  RDKit::EnumerationParams params;
  params.reagentMaxMatchCount = data.reagent_max_match_count;
  params.sanePartialProducts = data.sane_partial_products;
  return params;
}

rust::Vec<std::uint64_t> to_vec(const RGROUPS &values) {
  rust::Vec<std::uint64_t> res;
  for (auto value : values) {
    res.push_back(value);
  }
  return res;
}

const RDKit::EnumerationStrategyBase &strategy_of(const EnumerateLibrary &lib) {
  const auto &strategy = LibraryAccess::enumerator(lib);
  // RDKit's getPosition dereferences it unchecked (a segfault in Python for
  // an EnumerateLibrary()).
  if (!strategy) {
    throw ValueErrorException("the library has no enumeration strategy");
  }
  return *strategy;
}

}  // namespace

std::unique_ptr<EnumerationStrategyBase> enumeration_strategy_new(
    std::uint8_t kind) {
  switch (kind) {
    case 0:
      return std::make_unique<RDKit::CartesianProductStrategy>();
    case 1:
      return std::make_unique<RDKit::RandomSampleStrategy>();
    case 2:
      return std::make_unique<RDKit::RandomSampleAllBBsStrategy>();
    case 3:
      return std::make_unique<RDKit::EvenSamplePairsStrategy>();
    default:
      throw std::logic_error("unknown enumeration strategy");
  }
}

std::unique_ptr<EnumerationStrategyBase> enumeration_strategy_clone(
    const EnumerationStrategyBase &strategy) {
  return std::unique_ptr<EnumerationStrategyBase>(strategy.copy());
}

void enumeration_strategy_initialize(
    EnumerationStrategyBase &strategy, const ChemicalReaction &rxn,
    BorrowedMols &mols, rust::Slice<const std::uint32_t> group_sizes) {
  // The strategies keep only the groups' sizes.
  auto bbs = to_bbs(mols, group_sizes, false);
  check_initialize(strategy, RDKit::getSizesFromBBs(bbs));
  strategy.initialize(rxn, bbs);
}

rust::Vec<std::uint64_t> enumeration_strategy_get_position(
    const EnumerationStrategyBase &strategy) {
  return to_vec(strategy.getPosition());
}

std::uint64_t enumeration_strategy_get_num_permutations(
    const EnumerationStrategyBase &strategy) {
  return strategy.getNumPermutations();
}

std::unique_ptr<std::string> enumeration_strategy_type(
    const EnumerationStrategyBase &strategy) {
  return std::make_unique<std::string>(strategy.type());
}

bool enumeration_strategy_has_next(const EnumerationStrategyBase &strategy) {
  return static_cast<bool>(strategy);
}

rust::Vec<std::uint64_t> enumeration_strategy_next(
    EnumerationStrategyBase &strategy) {
  return to_vec(strategy.next());
}

bool enumeration_strategy_skip(EnumerationStrategyBase &strategy,
                               std::uint64_t skip_count) {
  return strategy.skip(skip_count);
}

std::unique_ptr<std::string> even_sample_pairs_strategy_stats(
    const EnumerationStrategyBase &strategy) {
  const auto *even =
      dynamic_cast<const RDKit::EvenSamplePairsStrategy *>(&strategy);
  if (!even) {
    throw std::logic_error("not an EvenSamplePairsStrategy");
  }
  return std::make_unique<std::string>(even->stats());
}

EnumerationParamsData enumeration_params_default() {
  RDKit::EnumerationParams params;
  return EnumerationParamsData{params.reagentMaxMatchCount,
                               params.sanePartialProducts};
}

bool enumerate_library_can_serialize() {
  return RDKit::EnumerateLibraryCanSerialize();
}

std::unique_ptr<EnumerateLibrary> enumerate_library_new() {
  return std::make_unique<EnumerateLibrary>();
}

std::unique_ptr<EnumerateLibrary> enumerate_library_new_with_params(
    const ChemicalReaction &rxn, BorrowedMols &mols,
    rust::Slice<const std::uint32_t> group_sizes,
    const EnumerationParamsData &data) {
  // Python's library shares Python's molecules without keeping them alive;
  // this one owns copies. Construction drops the building blocks the
  // reaction does not match, and running it matches too (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return std::make_unique<EnumerateLibrary>(
      rxn, to_bbs(mols, group_sizes, true), params_from(data));
}

std::unique_ptr<EnumerateLibrary> enumerate_library_new_with_enumerator(
    const ChemicalReaction &rxn, BorrowedMols &mols,
    rust::Slice<const std::uint32_t> group_sizes,
    const EnumerationStrategyBase &strategy,
    const EnumerationParamsData &data) {
  auto params = params_from(data);
  auto bbs = to_bbs(mols, group_sizes, true);
  // The constructor initializes the strategy with the building blocks its
  // initialized reaction matches; check those first.
  ChemicalReaction initialized(rxn);
  initialized.initReactantMatchers(true);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  check_initialize(strategy,
                   RDKit::getSizesFromBBs(RDKit::removeNonmatchingReagents(
                       initialized, bbs, params)));
  return std::make_unique<EnumerateLibrary>(rxn, bbs, strategy, params);
}

std::unique_ptr<ProductSets> enumerate_library_get_reagents(
    const EnumerateLibrary &lib) {
  return to_product_sets(lib.getReagents());
}

std::unique_ptr<ChemicalReaction> enumerate_library_get_reaction(
    const EnumerateLibrary &lib) {
  auto rxn = std::make_unique<ChemicalReaction>(lib.getReaction());
  normalize_templates(*rxn);
  return rxn;
}

std::unique_ptr<EnumerationStrategyBase> enumerate_library_get_enumerator(
    const EnumerateLibrary &lib) {
  return std::unique_ptr<EnumerationStrategyBase>(strategy_of(lib).copy());
}

rust::Vec<std::uint64_t> enumerate_library_get_position(
    const EnumerateLibrary &lib) {
  return to_vec(strategy_of(lib).getPosition());
}

bool enumerate_library_has_next(const EnumerateLibrary &lib) {
  return static_cast<bool>(strategy_of(lib));
}

std::unique_ptr<ProductSets> enumerate_library_next(EnumerateLibrary &lib) {
  auto &strategy = LibraryAccess::enumerator(lib);
  if (!strategy) {
    throw ValueErrorException("the library has no enumeration strategy");
  }
  if (!static_cast<bool>(*strategy)) {
    // As EnumerateLibrary::next's PRECONDITION.
    throw std::runtime_error("No more enumerations");
  }
  const auto &indices = strategy->next();
  const auto &bbs = lib.getReagents();
  // RDKit's next indexes the building blocks with these unchecked.
  if (indices.size() != bbs.size()) {
    bad_state("the position does not fit the building blocks");
  }
  MOL_SPTR_VECT reactants;
  for (std::size_t i = 0; i < bbs.size(); ++i) {
    if (indices[i] >= bbs[i].size()) {
      bad_state("the position does not fit the building blocks");
    }
    reactants.push_back(bbs[i][indices[i]]);
  }
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_product_sets(lib.getReaction().runReactants(reactants));
}

std::unique_ptr<std::string> enumerate_library_get_state(
    const EnumerateLibrary &lib) {
  strategy_of(lib);
  return std::make_unique<std::string>(lib.getState());
}

void enumerate_library_set_state(EnumerateLibrary &lib, rust::Str state) {
  auto strategy =
      RDKit::EnumerationStrategyPickler::fromPickle(std::string(state));
  if (!strategy) {
    throw ValueErrorException("the state holds no enumeration strategy");
  }
  check_state(*strategy, RDKit::getSizesFromBBs(lib.getReagents()));
  LibraryAccess::enumerator(lib) = strategy;
}

void enumerate_library_reset_state(EnumerateLibrary &lib) { lib.resetState(); }

// The archive pickles the reaction and the building blocks
// (check_picklable).
rust::Vec<std::uint8_t> enumerate_library_serialize(
    const EnumerateLibrary &lib) {
  check_picklable_templates(lib.getReaction());
  for (const auto &bbs : lib.getReagents()) {
    for (const auto &bb : bbs) {
      check_picklable(*bb);
    }
  }
  auto text = lib.Serialize();
  rust::Vec<std::uint8_t> res;
  res.reserve(text.size());
  for (char c : text) {
    res.push_back(static_cast<std::uint8_t>(c));
  }
  return res;
}

bool enumerate_library_is_empty(const EnumerateLibrary &lib) {
  // As EnumerateLibraryBase::load's check.
  const auto &rxn = lib.getReaction();
  return !rxn.getNumReactantTemplates() && !rxn.getNumProductTemplates() &&
         !rxn.getNumAgentTemplates();
}

std::unique_ptr<EnumerateLibrary> enumerate_library_from_string(
    rust::Slice<const std::uint8_t> data) {
  // Loaded into a new library, so a load that throws partway leaves nothing
  // behind. The reaction and building blocks are molecule pickles (U2); the
  // strategies are checked against the building blocks.
  auto lib = std::make_unique<EnumerateLibrary>();
  lib->initFromString(
      std::string(reinterpret_cast<const char *>(data.data()), data.size()));
  const EnumerateLibrary &loaded = *lib;
  // The reaction and the building blocks are molecule pickles (U10).
  check_reaction_monomer_types(loaded.getReaction());
  for (const auto &group : loaded.getReagents()) {
    for (const auto &mol : group) {
      if (mol) {
        check_pickled_mol(*mol);
      }
    }
  }
  auto sizes = RDKit::getSizesFromBBs(loaded.getReagents());
  for (const auto *strategy :
       {&LibraryAccess::enumerator(loaded), &LibraryAccess::initial(loaded)}) {
    if (*strategy) {
      check_state(**strategy, sizes);
    }
  }
  return lib;
}

}  // namespace rdkit_cxx
