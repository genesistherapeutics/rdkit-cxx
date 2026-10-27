// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/synthon_space_search.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/StreamOps.h>
#include <boost/algorithm/string.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/rascal_mces.rs.h"
#include "rdkit-cxx/ffi/substruct.rs.h"
#include "rdkit-cxx/ffi/synthon_space_search.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/rascal_mces.h"
#include "rdkit-cxx/streams.h"
#include "rdkit-cxx/substruct.h"
#include "rdkit-cxx/symmetrize_lock.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace SSS = RDKit::SynthonSpaceSearch;

constexpr double kTwoTo32 = 4294967296.0;

// ---------------------------------------------------------------------------
// The parameters.

SynthonSpaceSearchParamsData to_data(const SSS::SynthonSpaceSearchParams &p) {
  SynthonSpaceSearchParamsData d{};
  d.max_hits = p.maxHits;
  d.max_num_frags = p.maxNumFragSets;
  d.hit_start = p.hitStart;
  d.random_sample = p.randomSample;
  d.random_seed = p.randomSeed;
  d.build_hits = p.buildHits;
  d.similarity_cutoff = p.similarityCutoff;
  d.frag_similarity_adjuster = p.fragSimilarityAdjuster;
  d.approx_similarity_adjuster = p.approxSimilarityAdjuster;
  d.min_hit_heavy_atoms = p.minHitHeavyAtoms;
  d.max_hit_heavy_atoms = p.maxHitHeavyAtoms;
  d.min_hit_mol_wt = p.minHitMolWt;
  d.max_hit_mol_wt = p.maxHitMolWt;
  d.min_hit_chiral_atoms = p.minHitChiralAtoms;
  d.max_hit_chiral_atoms = p.maxHitChiralAtoms;
  d.time_out = p.timeOut;
  d.to_try_chunk_size = p.toTryChunkSize;
  d.num_threads = p.numThreads;
  return d;
}

// SynthonSpaceSearcher (SynthonSpaceSearcher.cpp) ends its search at
// Clock::now() + std::chrono::seconds(timeOut), in nanoseconds of int64,
// which overflows (undefined) for a timeOut that reaches past the clock's
// range.
void check_time_out(std::uint64_t seconds) {
  using Clock = std::chrono::steady_clock;
  auto room = std::chrono::duration_cast<std::chrono::seconds>(
                  Clock::duration::max() - Clock::now().time_since_epoch())
                  .count();
  if (seconds > static_cast<std::uint64_t>(room)) {
    throw ValueErrorException(
        "timeOut is too large: RDKit adds it to the clock in nanoseconds of "
        "64 bits");
  }
}

// The parameters RDKit runs a search with, checked for what it trusts:
// - buildAllHits erases the results from begin() + maxHits + hitStart and
//   clears those before begin() + hitStart, so a maxHits below -1 or a
//   negative hitStart reach before the vector's beginning (a segfault in
//   Python), and it adds the two in int64 (signed overflow);
// - the time limit (check_time_out).
// RDKit's own threads (numThreads) would each construct the ControlCHandlers
// of their substructure matches at once, which control_c_mutex cannot order
// (control_c_lock.h), and call std::terminate on an exception, so the search
// runs on the calling thread: the hits do not depend on the thread count.
SSS::SynthonSpaceSearchParams params_from(
    const SynthonSpaceSearchParamsData &d) {
  SSS::SynthonSpaceSearchParams p;
  p.maxHits = d.max_hits;
  p.maxNumFragSets = d.max_num_frags;
  p.hitStart = d.hit_start;
  p.randomSample = d.random_sample;
  p.randomSeed = d.random_seed;
  p.buildHits = d.build_hits;
  p.similarityCutoff = d.similarity_cutoff;
  p.fragSimilarityAdjuster = d.frag_similarity_adjuster;
  p.approxSimilarityAdjuster = d.approx_similarity_adjuster;
  p.minHitHeavyAtoms = d.min_hit_heavy_atoms;
  p.maxHitHeavyAtoms = d.max_hit_heavy_atoms;
  p.minHitMolWt = d.min_hit_mol_wt;
  p.maxHitMolWt = d.max_hit_mol_wt;
  p.minHitChiralAtoms = d.min_hit_chiral_atoms;
  p.maxHitChiralAtoms = d.max_hit_chiral_atoms;
  p.timeOut = d.time_out;
  p.toTryChunkSize = d.to_try_chunk_size;
  p.numThreads = 1;

  if (p.maxHits < -1) {
    throw ValueErrorException(
        "maxHits must be -1, for no limit, or more: RDKit erases its results "
        "from before their beginning for a smaller one");
  }
  if (p.hitStart < 0) {
    throw ValueErrorException(
        "hitStart must not be negative: RDKit clears results from before "
        "their beginning for a negative one");
  }
  if (p.maxHits >= 0 &&
      p.hitStart > std::numeric_limits<std::int64_t>::max() - p.maxHits) {
    throw ValueErrorException(
        "maxHits + hitStart must fit in 64 bits: RDKit adds them in int64");
  }
  check_time_out(p.timeOut);
  return p;
}

// ---------------------------------------------------------------------------
// What the space's reactions hold.

// Every reaction RDKit keeps, with its synthon sets.
template <class F>
void for_each_reaction(SSS::SynthonSpace &space, F &&f) {
  for (const auto &name : space.getReactionNames()) {
    f(name, *space.getReaction(name));
  }
}

// RWMol's removals of atoms and bonds follow each substance group's PARENT
// to the group with that "index" (check_substance_group_parents), and the
// reader and the searches remove atoms and bonds from
// products that combineMols builds from a synthon of each of a reaction's
// sets, keeping every group's properties, so a chain can run from one
// synthon's groups to another's. Every chain of a reaction's products runs
// along the links from each group of its synthons to each one with the index
// its PARENT names, which must have no cycle and no chain longer than
// MAX_PARENT_CHAIN. This can refuse a reaction none of whose products has
// such a chain, as a product's group with an index hides an earlier one's.
class ParentChains {
 public:
  void add(const RDKit::ROMol &mol) {
    for (const auto &sg : RDKit::getSubstanceGroups(mol)) {
      Group group;
      unsigned value;
      if (sg.getPropIfPresent("index", value)) {
        group.index = value;
      }
      if (sg.getPropIfPresent("PARENT", value)) {
        group.parent = value;
      }
      if (group.index || group.parent) {
        d_groups.push_back(group);
      }
    }
  }

  void add(const ParentChains &other) {
    d_groups.insert(d_groups.end(), other.d_groups.begin(),
                    other.d_groups.end());
  }

  // The longest chain from each group, by a depth-first search of its own
  // stack (the chains can be long).
  void check(const std::string &reaction) const {
    std::unordered_map<unsigned, std::vector<std::size_t>> by_index;
    bool parents = false;
    for (std::size_t i = 0; i < d_groups.size(); ++i) {
      if (d_groups[i].index) {
        by_index[*d_groups[i].index].push_back(i);
      }
      parents = parents || d_groups[i].parent;
    }
    if (!parents) {
      return;
    }
    static const std::vector<std::size_t> none;
    auto links = [&](std::size_t i) -> const std::vector<std::size_t> & {
      if (!d_groups[i].parent) {
        return none;
      }
      auto it = by_index.find(*d_groups[i].parent);
      return it == by_index.end() ? none : it->second;
    };
    enum State : char { kNew, kOpen, kDone };
    std::vector<State> state(d_groups.size(), kNew);
    std::vector<std::uint64_t> longest(d_groups.size(), 0);
    std::vector<std::pair<std::size_t, std::size_t>> stack;
    for (std::size_t start = 0; start < d_groups.size(); ++start) {
      if (state[start] != kNew) {
        continue;
      }
      state[start] = kOpen;
      stack.emplace_back(start, 0);
      while (!stack.empty()) {
        auto &[group, next] = stack.back();
        const auto &to = links(group);
        if (next < to.size()) {
          auto link = to[next++];
          if (state[link] == kOpen) {
            throw ValueErrorException(
                "reaction " + reaction +
                " has synthons whose substance groups' PARENT chains can make "
                "a cycle in a product, which RDKit follows forever when it "
                "removes an atom or bond");
          }
          if (state[link] == kNew) {
            state[link] = kOpen;
            stack.emplace_back(link, 0);
          } else {
            longest[group] = std::max(longest[group], longest[link] + 1);
          }
          continue;
        }
        state[group] = kDone;
        auto done = group;
        stack.pop_back();
        if (longest[done] > MAX_PARENT_CHAIN) {
          throw ValueErrorException(
              "reaction " + reaction +
              " has synthons whose substance groups' PARENT chains can be "
              "longer than " +
              std::to_string(MAX_PARENT_CHAIN) +
              " in a product, which RDKit follows recursively");
        }
        if (!stack.empty()) {
          auto parent = stack.back().first;
          longest[parent] = std::max(longest[parent], longest[done] + 1);
        }
      }
    }
  }

 private:
  struct Group {
    std::optional<unsigned> index;
    std::optional<unsigned> parent;
  };
  std::vector<Group> d_groups;
};

// What RDKit's code after a read trusts about the space it holds:
// - summarize counts reactions by their number of synthon sets in an array
//   of MAX_CONNECTOR_NUM + 1 (a write past it for more sets), and
//   buildAddAndSubtractFPs converts a product of up to 79 samples per set,
//   times 0.75, to int (undefined past INT_MAX, from 5 sets);
// - the searches, enumeration, and fingerprint building step through each
//   reaction's sets from the first, and read each set's first synthon
//   (front() of an empty vector);
// - every synthon has the molecules and pattern fingerprint the readers
//   make, which the searches and the writer dereference;
// - the searches remove atoms and bonds from the
//   synthons' molecules and their products (check_removals, ParentChains).
// A text file passes these by check_synthon_text; a database file, whatever
// RDKit wrote, is checked too.
void check_space(SSS::SynthonSpace &space) {
  for_each_reaction(
      space, [](const std::string &name, const SSS::SynthonSet &reaction) {
        const auto &sets = reaction.getSynthons();
        ParentChains chains;
        if (sets.empty() || sets.size() > SSS::MAX_CONNECTOR_NUM) {
          throw ValueErrorException("reaction " + name + " has " +
                                    std::to_string(sets.size()) +
                                    " synthon sets; RDKit takes 1 to " +
                                    std::to_string(SSS::MAX_CONNECTOR_NUM));
        }
        for (const auto &set : sets) {
          if (set.empty()) {
            throw ValueErrorException("reaction " + name +
                                      " has an empty synthon set");
          }
          for (const auto &[id, synthon] : set) {
            if (!synthon || !synthon->getOrigMol() ||
                !synthon->getSearchMol() || !synthon->getPattFP()) {
              throw ValueErrorException("reaction " + name + " has synthon " +
                                        id + " without its molecules");
            }
            for (const auto *mol :
                 {synthon->getOrigMol().get(), synthon->getSearchMol().get()}) {
              check_removals(*mol);
              chains.add(*mol);
            }
          }
        }
        chains.check(name);
      });
}

// writeDBFile pickles each synthon's molecules and connector regions and
// each reaction's connector regions (Synthon.cpp, SynthonSet.cpp), so a
// space it writes holds only molecules a pickle holds (check_picklable).
void check_space_picklable(SSS::SynthonSpace &space) {
  for_each_reaction(
      space, [](const std::string &, const SSS::SynthonSet &reaction) {
        for (const auto &region : reaction.getConnectorRegions()) {
          check_picklable(*region);
        }
        for (const auto &set : reaction.getSynthons()) {
          for (const auto &[id, synthon] : set) {
            check_picklable(*synthon->getOrigMol());
            check_picklable(*synthon->getSearchMol());
            for (const auto &region : synthon->getConnRegions()) {
              check_picklable(*region);
            }
          }
        }
      });
}

// The most atoms a product of the space can have: each reaction's largest
// synthon from each set, joined.
std::uint64_t largest_product(SSS::SynthonSpace &space) {
  std::uint64_t largest = 0;
  for_each_reaction(
      space, [&](const std::string &, const SSS::SynthonSet &reaction) {
        std::uint64_t atoms = 0;
        for (const auto &set : reaction.getSynthons()) {
          unsigned most = 0;
          for (const auto &[id, synthon] : set) {
            most = std::max({most, synthon->getOrigMol()->getNumAtoms(),
                             synthon->getSearchMol()->getNumAtoms()});
          }
          atoms += most;
        }
        largest = std::max(largest, atoms);
      });
  return largest;
}

// For RDKit code that computes each molecule's n * n matrices in int: the
// space's synthons and every product, which are at most largest_product
// atoms (check_matrix_size's limit, U3).
void check_space_matrix_sizes(SSS::SynthonSpace &space) {
  if (largest_product(space) > 46340) {
    throw ValueErrorException(
        "the space's products can have more than 46340 atoms, whose n * n "
        "matrices RDKit computes in int");
  }
}

void clear_space_matrix_caches(SSS::SynthonSpace &space) {
  for_each_reaction(space,
                    [](const std::string &, const SSS::SynthonSet &reaction) {
                      for (const auto &set : reaction.getSynthons()) {
                        for (const auto &[id, synthon] : set) {
                          clear_matrix_caches(*synthon->getOrigMol());
                          clear_matrix_caches(*synthon->getSearchMol());
                        }
                      }
                    });
}

// The fingerprint and RASCAL searches compute distance matrices of the
// space's synthons, which RDKit caches in their properties (U6); they are
// cleared before and after, as a MatrixCacheGuard clears a molecule's.
class SpaceMatrixCacheGuard {
 public:
  explicit SpaceMatrixCacheGuard(SSS::SynthonSpace &space) : d_space(space) {
    clear_space_matrix_caches(d_space);
  }
  ~SpaceMatrixCacheGuard() { clear_space_matrix_caches(d_space); }
  SpaceMatrixCacheGuard(const SpaceMatrixCacheGuard &) = delete;
  SpaceMatrixCacheGuard &operator=(const SpaceMatrixCacheGuard &) = delete;

 private:
  SSS::SynthonSpace &d_space;
};

// ---------------------------------------------------------------------------
// The text reader, readStream in SynthonSpace.cpp, which this follows line by
// line. What it trusts:
// - SynthonSet::addSynthon indexes the reaction's synthon sets by the line's
//   synthon number, resizing the sets to the number plus 1 in int: a negative
//   number writes before the vector (a segfault in Python), INT_MAX
//   overflows, and a large one sizes a vector of that many sets for each
//   reaction before the empty ones are dropped;
// - at most MAX_CONNECTOR_NUM synthon sets a reaction (check_space);
// - makeSynthonSearchMols zips a sample product from each synthon and then
//   finds the synthon's fragment of it, indexing at -1 when it finds none,
//   which a synthon of only dummy atoms (or none, an empty SMILES) gives (a
//   segfault in Python);
// - it counts lines in int (U12);
// - it splits lines with std::regex, whose matcher recurses once per
//   character of a run of separators (U12's limit);
// - it removes atoms and bonds from each synthon and from the sample
//   products (check_removals, ParentChains).

constexpr int kMaxSynthonNumber = 15;

std::vector<std::string> split_line(const std::string &str,
                                    const std::regex &re) {
  return {std::sregex_token_iterator(str.begin(), str.end(), re, -1),
          std::sregex_token_iterator()};
}

// As deduceFormat in SynthonSpace.cpp, once a line's separators are checked.
int deduce_format(const std::string &line) {
  const std::vector<std::vector<std::string>> first_line_opts{
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id"},
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id",
                               "release"},
      std::vector<std::string>{"SMILES", "synton_id", "synton_role",
                               "reaction_id"}};
  const std::regex regext("\\t");
  auto parts = split_line(line, regext);
  for (std::size_t i = 0; i < first_line_opts.size(); ++i) {
    if (parts == first_line_opts[i]) {
      return static_cast<int>(i + 6);
    }
  }
  const std::regex regexws("\\s+");
  parts = split_line(line, regexws);
  for (std::size_t i = 0; i < first_line_opts.size(); ++i) {
    if (parts == first_line_opts[i]) {
      return static_cast<int>(i);
    }
  }
  const std::regex regexc(",+");
  parts = split_line(line, regexc);
  for (std::size_t i = 0; i < first_line_opts.size(); ++i) {
    if (parts == first_line_opts[i]) {
      return static_cast<int>(i + first_line_opts.size());
    }
  }
  return -1;
}

// As getSynthonNum in SynthonSpace.cpp; throws as it does.
int synthon_number(int format, const std::string &field) {
  int num = std::numeric_limits<int>::max();
  if (format == 0 || format == 1 || format == 3 || format == 4 || format == 6 ||
      format == 7) {
    num = std::stoi(field);
  } else if (format == 2 || format == 5 || format == 8) {
    num = std::stoi(field.substr(7));
  }
  return num;
}

// As fixConnectors in SynthonSpace.cpp, whose patterns are these literals.
void fix_connectors(std::string &smiles) {
  for (unsigned i = 0; i < SSS::MAX_CONNECTOR_NUM; ++i) {
    const auto &from = SSS::CONNECTOR_SYMBOLS[i];
    const std::string to = "[" + std::to_string(i + 1) + "*]";
    std::string out;
    std::size_t pos = 0;
    for (auto at = smiles.find(from); at != std::string::npos;
         at = smiles.find(from, pos)) {
      out.append(smiles, pos, at - pos);
      out += to;
      pos = at + from.size();
    }
    out.append(smiles, pos, std::string::npos);
    smiles = std::move(out);
  }
}

void check_separator_runs(const std::string &line, std::uint64_t line_num) {
  std::size_t run = 0;
  for (char c : line) {
    if (c == ' ' || c == '\t' || c == '\v' || c == '\f' || c == '\r' ||
        c == ',') {
      if (++run > kMaxStereoGroupLine) {
        throw ValueErrorException(
            "line " + std::to_string(line_num) +
            " of the synthon file has a run of more than " +
            std::to_string(kMaxStereoGroupLine) +
            " separators, which RDKit's std::regex matches recursively");
      }
    } else {
      run = 0;
    }
  }
}

// The synthon's substance groups, parsed as RDKit's Synthon constructor
// parses it; none when the constructor will throw for the SMILES.
std::optional<ParentChains> check_synthon_smiles(const std::string &smiles,
                                                 std::uint64_t line_num) {
  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  check_cx_atom_props(smiles);
  auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles, params);
  if (!mol) {
    return std::nullopt;
  }
  auto where =
      "the synthon on line " + std::to_string(line_num) + " (" + smiles + ")";
  bool heavy = false;
  for (const auto *atom : mol->atoms()) {
    heavy = heavy || atom->getAtomicNum() != 0;
  }
  if (!heavy) {
    throw ValueErrorException(
        where +
        " has no atom but dummies: RDKit indexes its fragment of a product "
        "at -1");
  }
  try {
    check_removals(*mol);
  } catch (const ValueErrorException &e) {
    throw ValueErrorException(where + ": " + e.what());
  }
  ParentChains chains;
  chains.add(*mol);
  return chains;
}

// Lines as std::getline gives them to readStream: split at '\n', with the
// empty line after a final '\n', and as getLine strips a '\r' at the end.
void check_synthon_text(std::string_view text) {
  // Compiled per call, not as function statics: their thread-safe
  // initialization synchronizes through atomics and a futex that Valgrind's
  // DRD cannot see, so the races job would report every read of them.
  const std::regex regext("\\t+");
  const std::regex regexws("\\s+");
  const std::regex regexc(",+");
  int format = -1;
  std::uint64_t num_lines = 0;
  std::map<std::string, std::set<int>> sets;
  std::map<std::string, ParentChains> chains;
  std::unordered_map<std::string, ParentChains> checked;
  std::size_t pos = 0;
  while (true) {
    auto end = text.find('\n', pos);
    std::string line(text.substr(pos, end == std::string_view::npos
                                          ? std::string_view::npos
                                          : end - pos));
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (++num_lines > kMaxInputLines) {
      throw ValueErrorException("the synthon file has more than " +
                                std::to_string(kMaxInputLines) +
                                " lines, which RDKit counts in int");
    }
    check_separator_runs(line, num_lines);
    if (!line.empty() && line[0] != '#') {
      if (format == -1) {
        format = deduce_format(line);
        if (format == -1) {
          return;  // RDKit refuses the file here.
        }
      } else {
        std::vector<std::string> fields;
        if (format < 3) {
          fields = split_line(line, regexws);
        } else if (format > 3 && format < 6) {
          fields = split_line(line, regexc);
        } else if (format > 5) {
          fields = split_line(line, regext);
        }
        if (fields.size() < 4) {
          return;  // RDKit refuses the line.
        }
        int num;
        try {
          num = synthon_number(format, fields[2]);
        } catch (const std::exception &) {
          return;  // RDKit throws the same here.
        }
        if (num < 0 || num > kMaxSynthonNumber) {
          throw ValueErrorException(
              "the synthon number on line " + std::to_string(num_lines) +
              " is " + std::to_string(num) + "; RDKit indexes synthon sets " +
              "by it, which must be from 0 to " +
              std::to_string(kMaxSynthonNumber));
        }
        auto &numbers = sets[fields[3]];
        numbers.insert(num);
        if (numbers.size() > SSS::MAX_CONNECTOR_NUM) {
          throw ValueErrorException(
              "reaction " + fields[3] + " has more than " +
              std::to_string(SSS::MAX_CONNECTOR_NUM) +
              " synthon sets; RDKit counts a reaction's sets in an array of " +
              std::to_string(SSS::MAX_CONNECTOR_NUM + 1));
        }
        auto smiles = fields[0];
        fix_connectors(smiles);
        auto it = checked.find(smiles);
        if (it == checked.end()) {
          auto groups = check_synthon_smiles(smiles, num_lines);
          if (!groups) {
            return;  // RDKit refuses the SMILES.
          }
          it = checked.emplace(smiles, std::move(*groups)).first;
        }
        chains[fields[3]].add(it->second);
      }
    }
    if (end == std::string_view::npos) {
      break;
    }
    pos = end + 1;
  }
  for (const auto &[reaction, reaction_chains] : chains) {
    reaction_chains.check(reaction);
  }
}

void read_text(SynthonSpace &space, const std::string &name,
               const char *cancelled_what) {
  auto text = read_file(name);
  check_synthon_text(text);
  auto fresh = std::make_unique<SSS::SynthonSpace>();
  bool cancelled = false;
  {
    std::istringstream in(std::move(text));
    // readStream constructs a ControlCHandler, and the synthons' pattern
    // fingerprints match patterns (control_c_lock.h).
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    fresh->readStream(in, cancelled);
  }
  if (cancelled) {
    throw std::runtime_error(cancelled_what);
  }
  check_space(*fresh);
  space.space = std::move(fresh);
}

// ---------------------------------------------------------------------------
// The database files (writeDBFile and readDBFile in SynthonSpace.cpp). The
// reader trusts the file RDKit wrote: it reads the counts and offsets of the
// synthon and reaction records from the header, maps the file, and reads
// record i from offset i to offset i + 1 (the end of the file for the last),
// with no check against the file's size, and it reads the first reaction's
// offset even when there are no reactions (the file of an empty space; a
// segfault in Python). The records hold molecule and bit-vector pickles,
// counts, and bools RDKit trusts, which is why the reader is unsafe in Rust.
// This reads the header and offsets as RDKit does and refuses a file whose
// offsets do not ascend through it, which also tells a complete file from one
// the writer, which ignores a failed stream, left short: it writes zeros for
// the offsets first and the real ones last.

constexpr std::int32_t kEndianId = 0xa100f;
constexpr std::int32_t kVersionMajor = 3;
constexpr std::int32_t kVersionMinor = 1;

struct DbCounts {
  std::uint64_t synthons = 0;
  std::uint64_t reactions = 0;
};

class DbReader {
 public:
  DbReader(std::ifstream &in, std::uint64_t size, const std::string &name)
      : d_in(in), d_size(size), d_name(name) {}

  template <class T>
  T read() {
    unsigned char bytes[sizeof(T)];
    need(sizeof(T));
    d_in.read(reinterpret_cast<char *>(bytes), sizeof(T));
    if (!d_in) {
      truncated();
    }
    d_pos += sizeof(T);
    // Little-endian, as streamWrite writes.
    std::uint64_t value = 0;
    for (std::size_t i = 0; i < sizeof(T); ++i) {
      value |= std::uint64_t{bytes[i]} << (8 * i);
    }
    return static_cast<T>(value);
  }

  void skip(std::uint64_t n) {
    need(n);
    d_in.seekg(static_cast<std::streamoff>(n), std::ios::cur);
    d_pos += n;
  }

  std::uint64_t pos() const { return d_pos; }
  std::uint64_t size() const { return d_size; }

  [[noreturn]] void fail(const std::string &why) const {
    throw ValueErrorException("the synthon database " + d_name + " " + why);
  }

 private:
  void need(std::uint64_t n) const {
    if (n > d_size - d_pos) {
      truncated();
    }
  }
  [[noreturn]] void truncated() const { fail("is cut short"); }

  std::ifstream &d_in;
  std::uint64_t d_size;
  std::uint64_t d_pos = 0;
  const std::string &d_name;
};

DbCounts check_db_file(const std::string &name, bool reading) {
  std::ifstream in(name, std::ios::binary);
  if (!in || in.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  in.seekg(0, std::ios::end);
  auto end = in.tellg();
  if (end < 0) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  in.seekg(0, std::ios::beg);
  DbReader r(in, static_cast<std::uint64_t>(end), name);
  if (r.read<std::int32_t>() != kEndianId) {
    r.fail("is not a synthon database of this endianness");
  }
  auto major = r.read<std::int32_t>();
  auto minor = r.read<std::int32_t>();
  // RDKit refuses versions below 3.0 and reads later ones as its own after
  // a warning; a later RDKit's records may differ (a state from another
  // library, U2).
  if (major != kVersionMajor || minor < 0 || minor > kVersionMinor) {
    r.fail("has format version " + std::to_string(major) + "." +
           std::to_string(minor) + "; this RDKit reads 3.0 to 3.1");
  }
  // A bool RDKit reads byte for byte: any value but 0 or 1 is undefined.
  auto has_fps = r.read<std::uint8_t>();
  if (has_fps > 1) {
    r.fail("has a bool that is neither 0 nor 1");
  }
  if (has_fps) {
    r.skip(r.read<std::uint32_t>());
  }
  DbCounts counts;
  counts.synthons = r.read<std::uint64_t>();
  counts.reactions = r.read<std::uint64_t>();
  r.read<std::uint64_t>();
  if (reading && counts.reactions == 0) {
    r.fail(
        "has no reactions: RDKit reads the first reaction's offset all the "
        "same");
  }
  auto remaining = (r.size() - r.pos()) / 8;
  if (counts.synthons > remaining ||
      counts.reactions > remaining - counts.synthons) {
    r.fail("is cut short");
  }
  auto table_end = r.pos() + 8 * (counts.synthons + counts.reactions);
  std::uint64_t previous = 0;
  for (std::uint64_t i = 0; i < counts.synthons + counts.reactions; ++i) {
    auto offset = r.read<std::uint64_t>();
    bool ok = i == 0 ? offset == table_end : offset > previous;
    if (!ok || offset >= r.size()) {
      r.fail("has record offsets that do not ascend through the file");
    }
    previous = offset;
  }
  return counts;
}

// ---------------------------------------------------------------------------
// Files RDKit writes by name. It ignores a stream that fails (U12).

// writeDBFile opens the file itself, writes zeros for the offsets, the
// records, and then the offsets: the file is truncated first, so a write that
// fails leaves a file check_db_file refuses, not an older one.
// The synthons writeDBFile writes, one record for each in the space's pool,
// which is keyed on their SMILES and holds those of every reaction's sets.
// RDKit 2025.09 has no getNumSynthons.
std::uint64_t count_synthons(SSS::SynthonSpace &space) {
  std::set<std::string> smiles;
  for_each_reaction(space,
                    [&](const std::string &, const SSS::SynthonSet &reaction) {
                      for (const auto &set : reaction.getSynthons()) {
                        for (const auto &[id, synthon] : set) {
                          smiles.insert(synthon->getSmiles());
                        }
                      }
                    });
  return smiles.size();
}

void write_db(SSS::SynthonSpace &space, const std::string &name) {
  {
    std::ofstream out(name, std::ios::binary | std::ios::trunc);
    if (!out) {
      throw RDKit::BadFileException("Bad output file " + name);
    }
  }
  space.writeDBFile(name);
  try {
    auto counts = check_db_file(name, false);
    if (counts.synthons != count_synthons(space) ||
        counts.reactions != space.getNumReactions()) {
      throw ValueErrorException("does not hold the space");
    }
  } catch (const ValueErrorException &e) {
    throw RDKit::BadFileException("could not write " + name + ": " + e.what());
  }
}

// ---------------------------------------------------------------------------
// The queries, which the searches split into fragments, removing bonds and
// atoms (check_removals). An extended query's search splits each of its
// molecules, or each tautomer of its tautomer queries (extendedSearch in
// SynthonSpace.cpp).

void check_query(const RDKit::ROMol &query) { check_removals(query); }

void check_query(const ExtendedQueryMol &query) {
  using XQ = ExtendedQueryMol;
  auto check_tautomers = [](const RDKit::TautomerQuery &tautomers) {
    check_removals(tautomers.getTemplateMolecule());
    for (const auto &tautomer : tautomers.getTautomers()) {
      check_removals(*tautomer);
    }
  };
  if (const auto *mol = std::get_if<XQ::RWMol_T>(&query.xqmol)) {
    if (*mol) {
      check_removals(**mol);
    }
  } else if (const auto *bundle = std::get_if<XQ::MolBundle_T>(&query.xqmol)) {
    if (*bundle) {
      for (const auto &mol : (*bundle)->getMols()) {
        check_removals(*mol);
      }
    }
  } else if (const auto *tq = std::get_if<XQ::TautomerQuery_T>(&query.xqmol)) {
    if (*tq) {
      check_tautomers(**tq);
    }
  } else if (const auto *tqs =
                 std::get_if<XQ::TautomerBundle_T>(&query.xqmol)) {
    if (*tqs) {
      for (const auto &tq : **tqs) {
        check_tautomers(*tq);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Results and callbacks.

// Python raises KeyboardInterrupt for a search RDKit's SIGINT handler
// canceled.
std::unique_ptr<SearchResults> finish(SSS::SearchResults results,
                                      const char *what) {
  if (results.getCancelled()) {
    throw std::runtime_error(std::string(what) + " canceled");
  }
  return std::make_unique<SearchResults>(std::move(results));
}

// Hit molecules for Rust, which owns an RWMol (W1); a copy carries RDKit's
// cached matrices, which are dropped (U6).
std::unique_ptr<RWMol> copy_hit(const RDKit::ROMol &mol) {
  auto copy = std::make_unique<RWMol>(mol);
  clear_matrix_caches(*copy);
  return copy;
}

bool call_rust(SearchCallback &callback,
               std::vector<std::unique_ptr<RDKit::ROMol>> &hits) {
  auto list = std::make_unique<MolList>();
  for (const auto &hit : hits) {
    list->push(copy_hit(*hit));
  }
  hits.clear();
  return call_search_callback(callback, std::move(list));
}

// Marks a closure RDKit calls with a ControlCHandler alive, as fmcs.cpp's
// does for its final match checks: the lock stays held, so a wrapper the
// closure calls that would take it fails (control_c_lock.h).
class InControlCCallback {
 public:
  InControlCCallback() : d_outer(in_control_c_callback) {
    in_control_c_callback = true;
  }
  ~InControlCCallback() { in_control_c_callback = d_outer; }
  InControlCCallback(const InControlCCallback &) = delete;
  InControlCCallback &operator=(const InControlCCallback &) = delete;

 private:
  bool d_outer;
};

// Python's CallbackAdapter, which hands each batch of hits to the closure;
// the closure's error (or panic) is the rust::Error it throws, which RDKit,
// on this thread, lets through. SynthonSpace's substructure and fingerprint
// searches construct a ControlCHandler for the whole search, so the closure
// runs with control_c_mutex held.
SSS::SearchResultCallback held_callback(SearchCallback &callback) {
  return [&callback](std::vector<std::unique_ptr<RDKit::ROMol>> &hits) {
    InControlCCallback in_callback;
    return call_rust(callback, hits);
  };
}

// The RASCAL search constructs no ControlCHandler of its own, and calls the
// closure between its matches, so the closure runs with the lock released.
SSS::SearchResultCallback released_callback(SearchCallback &callback) {
  return [&callback](std::vector<std::unique_ptr<RDKit::ROMol>> &hits) {
    ReleaseControlC release;
    return call_rust(callback, hits);
  };
}

// ---------------------------------------------------------------------------
// The fingerprint and RASCAL searches' thresholds.

// findSynthonSearchStartF and getHitSynthons
// (SynthonSpaceFingerprintSearcher.cpp) convert the fragment threshold,
// similarityCutoff - fragSimilarityAdjuster, times a fragment's on-bits, and
// 1 + those bits over it, to unsigned int: undefined outside its range, for
// a negative or NaN threshold, or one so small (but not 0, which it handles)
// or large that the bits it bounds pass 2^32.
void check_fingerprint_threshold(const SSS::SynthonSpaceSearchParams &params,
                                 const FpGenerator &gen) {
  double t = params.similarityCutoff - params.fragSimilarityAdjuster;
  double bits = gen.getOptions()->d_fpSize;
  bool ok =
      t == 0.0 || (t > 0.0 && bits / t + 1.0 < kTwoTo32 && t * bits < kTwoTo32);
  if (!ok) {
    throw ValueErrorException(
        "similarityCutoff - fragSimilarityAdjuster must be 0, or above 0 and "
        "small enough that fpSize times it, and fpSize over it, stay below "
        "2^32: RDKit converts them to unsigned int");
  }
}

// What the fingerprint searches fingerprint: the query (whose matrix caches
// the generator fills, U6), its fragments (check_query), which have its
// atoms and at most 8 dummies, and the synthons and products, which atom
// pairs need small enough for their distance matrices.
void check_fingerprint_search(SSS::SynthonSpace *space, const RWMol &query,
                              const FpGenerator &gen,
                              const SSS::SynthonSpaceSearchParams &params) {
  check_query(query);
  check_fingerprint_mol(gen, query);
  if (fingerprint_reads_distance_matrix(gen)) {
    if (query.getNumAtoms() > 46340 - 8) {
      throw ValueErrorException(
          "the query's fragments can have more than 46340 atoms, whose n * n "
          "matrices RDKit computes in int");
    }
    if (space) {
      check_space_matrix_sizes(*space);
    }
  }
  check_fingerprint_threshold(params, gen);
}

// getHitSynthons and findSynthonSearchStartR (SynthonSpaceRascalSearcher.cpp)
// convert the fragment threshold, the options' similarityThreshold -
// fragSimilarityAdjuster, times a fragment's atoms and bonds, and 1 + those
// over it, to unsigned int, with no case for a threshold of 0 or less, and
// RASCAL then runs with that threshold. Fragments (check_query) have the
// query's atoms and bonds and at most 8 dummies with their bonds.
void check_rascal_search(SSS::SynthonSpace *space, const RWMol &query,
                         const RDKit::RascalMCES::RascalOptions &options,
                         const SSS::SynthonSpaceSearchParams &params) {
  check_query(query);
  double n =
      static_cast<double>(query.getNumAtoms()) + query.getNumBonds() + 16.0;
  if (n > static_cast<double>(kRascalMaxSquared)) {
    throw ValueErrorException(
        "the query's fragments can have more than 46340 atoms and bonds, "
        "whose matched counts RASCAL squares in int");
  }
  if (space) {
    check_space_matrix_sizes(*space);
  }
  double t = options.similarityThreshold - params.fragSimilarityAdjuster;
  if (!(t > 0.0 && n / t + 1.0 < kTwoTo32 && t * n < kTwoTo32)) {
    throw ValueErrorException(
        "similarityThreshold - fragSimilarityAdjuster must be above 0, and "
        "small enough that the query's atoms and bonds times it, and over it, "
        "stay below 2^32: RDKit converts them to unsigned int");
  }
}

}  // namespace

SynthonSpace::SynthonSpace()
    : space(std::make_unique<RDKit::SynthonSpaceSearch::SynthonSpace>()) {}

SynthonSpaceSearchParamsData synthon_space_search_params_default() {
  return to_data(SSS::SynthonSpaceSearchParams());
}

std::unique_ptr<SynthonSpace> synthon_space_new() {
  return std::make_unique<SynthonSpace>();
}

void synthon_space_read_text_file(SynthonSpace &space,
                                  rust::Slice<const std::uint8_t> path) {
  read_text(space, to_path(path), "Database read canceled");
}

void synthon_space_read_db_file(SynthonSpace &space,
                                rust::Slice<const std::uint8_t> path) {
  auto name = to_path(path);
  check_db_file(name, true);
  auto fresh = std::make_unique<SSS::SynthonSpace>();
  // RDKit's reader threads (numThreads) would call std::terminate on an
  // exception; the space does not depend on how many read it.
  fresh->readDBFile(name, 1);
  check_space(*fresh);
  space.space = std::move(fresh);
}

void synthon_space_write_db_file(const SynthonSpace &space,
                                 rust::Slice<const std::uint8_t> path) {
  check_space_picklable(*space.space);
  write_db(*space.space, to_path(path));
}

// As writeEnumeratedFile, which opens the file itself; the wrapper opens it,
// has RDKit write to the stream, and checks the stream (U12).
void synthon_space_write_enumerated_file(const SynthonSpace &space,
                                         rust::Slice<const std::uint8_t> path) {
  auto name = to_path(path);
  std::ofstream out(name);
  if (!out) {
    throw RDKit::BadFileException("Bad output file " + name);
  }
  BOOST_LOG(rdWarningLog) << "Writing the enumerated file may take quite"
                          << " some time and result in a large file."
                          << std::endl;
  space.space->enumerateToStream(out);
  out.close();
  if (out.fail()) {
    throw RDKit::BadFileException("could not write " + name);
  }
}

std::uint64_t synthon_space_get_num_reactions(const SynthonSpace &space) {
  return space.space->getNumReactions();
}

std::uint64_t synthon_space_get_num_products(const SynthonSpace &space) {
  return space.space->getNumProducts();
}

void synthon_space_summarise(const SynthonSpace &space) {
  space.space->summarise(std::cout);
}

std::unique_ptr<std::string> synthon_space_get_synthon_fingerprint_type(
    const SynthonSpace &space) {
  return std::make_unique<std::string>(
      space.space->getSynthonFingerprintType());
}

// buildSynthonFingerprints reads the Ctrl-C flag, and the generator may
// match patterns or run the CIP labeler (control_c_lock.h).
void synthon_space_build_synthon_fingerprints(SynthonSpace &space,
                                              const FpGenerator &gen) {
  if (fingerprint_reads_distance_matrix(gen)) {
    check_space_matrix_sizes(*space.space);
  }
  SpaceMatrixCacheGuard guard(*space.space);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  space.space->buildSynthonFingerprints(gen);
}

// The substructure searches construct a ControlCHandler, and match.

std::unique_ptr<SearchResults> synthon_space_substructure_search(
    SynthonSpace &space, const RWMol &query,
    const SubstructMatchParametersData &match_params,
    const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  auto match = substruct_params_from_data(match_params);
  check_query(query);
  SSS::SearchResults results;
  {
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    results = space.space->substructureSearch(query, match, params);
  }
  return finish(std::move(results), "SubstructureSearch");
}

std::unique_ptr<SearchResults> synthon_space_substructure_search_extended(
    SynthonSpace &space, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &match_params,
    const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  auto match = substruct_params_from_data(match_params);
  check_query(query);
  SSS::SearchResults results;
  {
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    results = space.space->substructureSearch(query, match, params);
  }
  return finish(std::move(results), "SubstructureSearch");
}

void synthon_space_substructure_search_incremental(
    SynthonSpace &space, const RWMol &query, SearchCallback &callback,
    const SubstructMatchParametersData &match_params,
    const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  auto match = substruct_params_from_data(match_params);
  check_query(query);
  {
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    space.space->substructureSearch(query, held_callback(callback), match,
                                    params);
  }
}

// The fingerprint searches construct a ControlCHandler, read the Ctrl-C
// flag, and fingerprint with a generator that may match patterns or run the
// CIP labeler.

std::unique_ptr<SearchResults> synthon_space_fingerprint_search(
    SynthonSpace &space, const RWMol &query, const FpGenerator &gen,
    const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  check_fingerprint_search(space.space.get(), query, gen, params);
  SSS::SearchResults results;
  {
    MatrixCacheGuard query_guard(query);
    SpaceMatrixCacheGuard space_guard(*space.space);
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    results = space.space->fingerprintSearch(query, gen, params);
  }
  return finish(std::move(results), "FingerprintSearch");
}

void synthon_space_fingerprint_search_incremental(
    SynthonSpace &space, const RWMol &query, const FpGenerator &gen,
    SearchCallback &callback, const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  check_fingerprint_search(space.space.get(), query, gen, params);
  {
    MatrixCacheGuard query_guard(query);
    SpaceMatrixCacheGuard space_guard(*space.space);
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    space.space->fingerprintSearch(query, gen, held_callback(callback), params);
  }
}

// The RASCAL searches match substructures (deltaYExchangePossible in
// RascalMCES.cpp), as do their hits' fragments' pattern fingerprints.

std::unique_ptr<SearchResults> synthon_space_rascal_search(
    SynthonSpace &space, const RWMol &query,
    const RascalOptionsData &rascal_options,
    const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  auto options = rascal_options_from(rascal_options);
  check_rascal_search(space.space.get(), query, options, params);
  SSS::SearchResults results;
  {
    SpaceMatrixCacheGuard space_guard(*space.space);
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    results = space.space->rascalSearch(query, options, params);
  }
  return finish(std::move(results), "RascalSearch");
}

void synthon_space_rascal_search_incremental(
    SynthonSpace &space, const RWMol &query,
    const RascalOptionsData &rascal_options, SearchCallback &callback,
    const SynthonSpaceSearchParamsData &data) {
  auto params = params_from(data);
  auto options = rascal_options_from(rascal_options);
  check_rascal_search(space.space.get(), query, options, params);
  {
    SpaceMatrixCacheGuard space_guard(*space.space);
    HoldControlC control_c;
    space.space->rascalSearch(query, options, released_callback(callback),
                              params);
  }
}

std::unique_ptr<MolList> substructure_result_get_hit_molecules(
    const SearchResults &results) {
  auto list = std::make_unique<MolList>();
  for (const auto &hit : results.getHitMolecules()) {
    list->push(copy_hit(*hit));
  }
  return list;
}

std::uint64_t substructure_result_get_max_num_results(
    const SearchResults &results) {
  return results.getMaxNumResults();
}

bool substructure_result_get_timed_out(const SearchResults &results) {
  return results.getTimedOut();
}

bool substructure_result_get_cancelled(const SearchResults &results) {
  return results.getCancelled();
}

std::unique_ptr<std::string> formatted_integer_string(std::int64_t value) {
  return std::make_unique<std::string>(
      RDKit::SynthonSpaceSearch::formattedIntegerString(value));
}

}  // namespace rdkit_cxx
