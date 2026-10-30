// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/filter_catalog.h"

#include <algorithm>
#include <climits>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/filter_catalog.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/threads.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

using RDKit::ExclusionList;
using RDKit::FilterCatalogEntry;
using RDKit::FilterHierarchyMatcher;
using RDKit::FilterMatch;
using RDKit::FilterMatcherBase;
using RDKit::SmartsMatcher;
using RDKit::FilterMatchOps::And;
using RDKit::FilterMatchOps::Not;
using RDKit::FilterMatchOps::Or;
using MatcherPtr = boost::shared_ptr<FilterMatcherBase>;
using NodePtr = boost::shared_ptr<FilterHierarchyMatcher>;

// The composite matchers keep their parts private. An explicit
// instantiation may name private members ([temp.spec.general]), so
// Steal<Tag, &C::member> defines a friend that returns the pointer to the
// member, found through its tag.
template <class Tag, typename Tag::type Member>
struct Steal {
  friend typename Tag::type member(Tag) { return Member; }
};

#define RDKIT_CXX_MATCHER_MEMBER(Tag, Class, Type, name) \
  struct Tag {                                           \
    using type = Type Class::*;                          \
    friend type member(Tag);                             \
  };                                                     \
  template struct Steal<Tag, &Class::name>;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif
RDKIT_CXX_MATCHER_MEMBER(AndArg1, And, MatcherPtr, arg1)
RDKIT_CXX_MATCHER_MEMBER(AndArg2, And, MatcherPtr, arg2)
RDKIT_CXX_MATCHER_MEMBER(OrArg1, Or, MatcherPtr, arg1)
RDKIT_CXX_MATCHER_MEMBER(OrArg2, Or, MatcherPtr, arg2)
RDKIT_CXX_MATCHER_MEMBER(NotArg1, Not, MatcherPtr, arg1)
RDKIT_CXX_MATCHER_MEMBER(OffPatterns, ExclusionList, std::vector<MatcherPtr>,
                         d_offPatterns)
RDKIT_CXX_MATCHER_MEMBER(NodeChildren, FilterHierarchyMatcher,
                         std::vector<NodePtr>, d_children)
RDKIT_CXX_MATCHER_MEMBER(NodeMatcher, FilterHierarchyMatcher, MatcherPtr,
                         d_matcher)
RDKIT_CXX_MATCHER_MEMBER(EntryMatcher, FilterCatalogEntry, MatcherPtr,
                         d_matcher)
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#undef RDKIT_CXX_MATCHER_MEMBER

template <class Tag, class Class>
const auto &part(const Class &matcher) {
  return matcher.*member(Tag());
}

// A FilterHierarchyMatcher() has no matcher of its own, which its isValid,
// hasMatch, and getMatches dereference (its documentation calls it a root
// that passes matches to its children).
bool is_empty_node(const FilterMatcherBase &matcher) {
  const auto *node = dynamic_cast<const FilterHierarchyMatcher *>(&matcher);
  return node && !part<NodeMatcher>(*node);
}

// Every substructure match constructs a ControlCHandler, so the wrappers
// whose RDKit call can match hold control_c_mutex for the whole call
// through HoldControlC (control_c_lock.h, CONVENTIONS.md U4). A Rust
// matcher's methods may call such wrappers themselves, and RDKit calls them
// only between its matches, with no ControlCHandler alive (one lives only
// inside a SubstructMatch, which calls no matcher), so each runs inside a
// ReleaseControlC, which releases the lock and takes it back after, even
// when they throw. A callback can run while it is already released: one
// that calls a wrapper that does not match (filter_matcher_is_valid) on
// another Rust matcher.

// A Rust callback gets its own copy of the molecule, as a Rust molecule,
// without cached matrices or descriptors (U6, U8).
std::unique_ptr<RWMol> callback_mol(const RDKit::ROMol &mol) {
  auto copy = std::make_unique<RWMol>(mol);
  clear_matrix_caches(*copy);
  clear_descriptor_caches(*copy);
  return copy;
}

// A matcher whose methods call a Rust callback, which is Send + Sync since
// RDKit's catalogs share their matchers between threads. Copies share it. A
// panic in it becomes an exception in Rust's half of the bridge, which
// RDKit's matching passes on.
class RustFilterMatcher : public FilterMatcherBase {
 public:
  RustFilterMatcher(std::string name,
                    std::shared_ptr<rust::Box<FilterMatcherCallback>> callback)
      : FilterMatcherBase(std::move(name)), d_callback(std::move(callback)) {}

  bool isValid() const override {
    ReleaseControlC released;
    return filter_matcher_callback_is_valid(**d_callback);
  }

  bool hasMatch(const RDKit::ROMol &mol) const override {
    auto copy = callback_mol(mol);
    ReleaseControlC released;
    return filter_matcher_callback_has_match(**d_callback, std::move(copy));
  }

  bool getMatches(const RDKit::ROMol &mol,
                  std::vector<FilterMatch> &matches) const override {
    FilterMatchList list;
    auto copy = callback_mol(mol);
    bool matched;
    {
      ReleaseControlC released;
      matched = filter_matcher_callback_get_matches(**d_callback,
                                                    std::move(copy), list);
    }
    // As the Python trampoline, matches count only when the matcher matched.
    if (matched) {
      matches.insert(matches.end(), list.matches.begin(), list.matches.end());
    }
    return matched;
  }

  boost::shared_ptr<FilterMatcherBase> copy() const override {
    return boost::make_shared<RustFilterMatcher>(*this);
  }

 private:
  std::shared_ptr<rust::Box<FilterMatcherCallback>> d_callback;
};

// What a walk of a matcher's parts checks.
enum class Walk {
  // A matcher that another will hold, whose methods RDKit then calls: no
  // part may be an empty hierarchy node.
  Use,
  // A matcher to pickle: SmartsMatcher's save dereferences its pattern, and
  // pickles it (check_picklable).
  Pickle,
  // A matcher read from an archive, which can hold what Rust cannot make:
  // empty hierarchy nodes or children, empty exclusion patterns, cycles,
  // and patterns whose monomer infos RDKit's pickler would misread (U10).
  Load,
};

// The parts a matcher holds, and checks of the matcher itself.
std::vector<const FilterMatcherBase *> parts(const FilterMatcherBase &matcher,
                                             Walk walk) {
  std::vector<const FilterMatcherBase *> result;
  auto add = [&](const MatcherPtr &p) {
    if (p) {
      result.push_back(p.get());
    }
  };
  if (walk == Walk::Pickle &&
      dynamic_cast<const RustFilterMatcher *>(&matcher)) {
    // Boost knows only RDKit's matcher classes.
    throw ValueErrorException("a matcher written in Rust cannot be pickled");
  }
  if (const auto *m = dynamic_cast<const SmartsMatcher *>(&matcher)) {
    const auto &pattern = m->getPattern();
    if (walk == Walk::Pickle && !pattern) {
      throw ValueErrorException(
          "a SmartsMatcher without a pattern cannot be pickled");
    }
    if (walk == Walk::Pickle) {
      check_picklable(*pattern);
    }
    if (walk == Walk::Load && pattern) {
      check_pickled_mol(*pattern);
    }
  } else if (const auto *m = dynamic_cast<const And *>(&matcher)) {
    // A missing argument makes the matcher invalid, which RDKit checks.
    add(part<AndArg1>(*m));
    add(part<AndArg2>(*m));
  } else if (const auto *m = dynamic_cast<const Or *>(&matcher)) {
    add(part<OrArg1>(*m));
    add(part<OrArg2>(*m));
  } else if (const auto *m = dynamic_cast<const Not *>(&matcher)) {
    add(part<NotArg1>(*m));
  } else if (const auto *m = dynamic_cast<const ExclusionList *>(&matcher)) {
    for (const auto &p : part<OffPatterns>(*m)) {
      // isValid dereferences each pattern.
      if (!p) {
        throw ValueErrorException("an ExclusionList has an empty pattern");
      }
      add(p);
    }
  } else if (const auto *m =
                 dynamic_cast<const FilterHierarchyMatcher *>(&matcher)) {
    if (walk != Walk::Pickle && !part<NodeMatcher>(*m)) {
      throw ValueErrorException(
          "a FilterHierarchyMatcher without a pattern cannot be used inside "
          "another matcher or an entry");
    }
    add(part<NodeMatcher>(*m));
    for (const auto &child : part<NodeChildren>(*m)) {
      // getMatches dereferences each child.
      if (!child) {
        throw ValueErrorException(
            "a FilterHierarchyMatcher has an empty child");
      }
      result.push_back(child.get());
    }
  }
  return result;
}

// Walks a matcher and its parts depth first, with its own stack, visiting
// each once; a part that holds a matcher still being walked is a cycle,
// which RDKit's recursive methods would follow forever.
void check_matcher(const FilterMatcherBase &root, Walk walk) {
  struct Frame {
    std::vector<const FilterMatcherBase *> parts;
    std::size_t next = 0;
    const FilterMatcherBase *matcher;
  };
  std::unordered_map<const FilterMatcherBase *, bool> done;
  std::vector<Frame> stack;
  auto enter = [&](const FilterMatcherBase *matcher) {
    done[matcher] = false;
    stack.push_back(Frame{parts(*matcher, walk), 0, matcher});
  };
  enter(&root);
  while (!stack.empty()) {
    auto &frame = stack.back();
    if (frame.next == frame.parts.size()) {
      done[frame.matcher] = true;
      stack.pop_back();
      continue;
    }
    const auto *next = frame.parts[frame.next++];
    auto seen = done.find(next);
    if (seen == done.end()) {
      enter(next);
    } else if (!seen->second) {
      throw ValueErrorException("a matcher contains itself");
    }
  }
}

void check_use(const FilterMatcherHandle &matcher) {
  check_matcher(*matcher.ptr, Walk::Use);
}

// SubstructMatch records a recursive query's matches in the query, holding
// its mutex for the whole match (RDK_BUILD_THREADSAFE_SSS), while
// RecursiveStructureQuery::copy reads them without it. The crate's own
// matches are serialized by control_c_mutex, which this holds too, but the
// query locks also order it against any match RDKit makes without it. A pattern
// shared with matchers other threads run (RDKit's copies share their patterns,
// as the functional-group hierarchy shares its patterns with its flattened map)
// is copied with each recursive query locked, in the order SubstructMatch locks
// them: each atom's query tree in atom order, a recursive query before its
// molecule's atoms' queries, and those before its children.
std::unique_ptr<RWMol> copy_pattern(const RDKit::ROMol &pattern) {
  using AtomQuery = RDKit::Atom::QUERYATOM_QUERY;
  // The crate's matches hold control_c_mutex and then, in SubstructMatch,
  // the query locks; taking them in the same order cannot deadlock.
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::vector<std::unique_lock<std::mutex>> locks;
  std::unordered_set<const void *> seen;
  std::vector<const AtomQuery *> todo;
  auto push_atoms = [&](const RDKit::ROMol &mol) {
    std::vector<const AtomQuery *> queries;
    for (const auto *atom : mol.atoms()) {
      if (atom->hasQuery()) {
        queries.push_back(atom->getQuery());
      }
    }
    todo.insert(todo.end(), queries.rbegin(), queries.rend());
  };
  push_atoms(pattern);
  while (!todo.empty()) {
    const auto *query = todo.back();
    todo.pop_back();
    if (!seen.insert(query).second) {
      continue;
    }
    std::vector<const AtomQuery *> children;
    for (auto child = query->beginChildren(); child != query->endChildren();
         ++child) {
      children.push_back(child->get());
    }
    todo.insert(todo.end(), children.rbegin(), children.rend());
    const auto *recursive =
        query->getDescription() == "RecursiveStructure"
            ? dynamic_cast<const RDKit::RecursiveStructureQuery *>(query)
            : nullptr;
    if (recursive) {
      locks.emplace_back(
          const_cast<RDKit::RecursiveStructureQuery *>(recursive)->d_mutex);
      const auto *inner = recursive->getQueryMol();
      if (inner && seen.insert(inner).second) {
        push_atoms(*inner);
      }
    }
  }
  auto copy = std::make_unique<RWMol>(pattern);
  clear_matrix_caches(*copy);
  return copy;
}

std::unique_ptr<FilterMatcherHandle> handle(MatcherPtr ptr) {
  auto result = std::make_unique<FilterMatcherHandle>();
  result->ptr = std::move(ptr);
  return result;
}

template <class T>
T &as(FilterMatcherHandle &matcher) {
  auto *result = dynamic_cast<T *>(matcher.ptr.get());
  if (!result) {
    throw std::logic_error("a matcher handle holds another class");
  }
  return *result;
}

template <class T>
const T &as(const FilterMatcherHandle &matcher) {
  return as<T>(const_cast<FilterMatcherHandle &>(matcher));
}

std::unique_ptr<FilterMatchList> match_list(std::vector<FilterMatch> matches) {
  auto list = std::make_unique<FilterMatchList>();
  list->matches = std::move(matches);
  return list;
}

std::unique_ptr<FilterCatalogEntryHandle> entry_handle(
    boost::shared_ptr<const FilterCatalogEntry> entry) {
  auto result = std::make_unique<FilterCatalogEntryHandle>();
  result->ptr = boost::const_pointer_cast<FilterCatalogEntry>(std::move(entry));
  return result;
}

// The entry, copied first unless this handle is its only owner: catalogs,
// match results, and other handles share entries, which other threads may
// read. A count of one means no other owner exists to copy the pointer
// (FilterCatalogEntry is not enable_shared_from_this).
FilterCatalogEntry &entry_mut(FilterCatalogEntryHandle &entry) {
  if (entry.ptr.use_count() != 1) {
    entry.ptr = boost::make_shared<FilterCatalogEntry>(*entry.ptr);
  }
  return *entry.ptr;
}

// An entry without a matcher is invalid, which RDKit checks.
void check_entry(const FilterCatalogEntry &entry, Walk walk) {
  if (const auto &matcher = part<EntryMatcher>(entry)) {
    check_matcher(*matcher, walk);
  }
}

// Boost.Serialization writes the strings' lengths as size_t; nothing else
// the entries pickle is sized. FilterCatalogCanSerialize is true in builds
// with Boost.Serialization; without, RDKit's precondition reports it.
void check_can_serialize() {
  if (!RDKit::FilterCatalogCanSerialize()) {
    throw std::runtime_error(
        "Pickling of FilterCatalog instances is not enabled");
  }
}

std::unique_ptr<std::string> copy_string(const std::string &text) {
  return std::make_unique<std::string>(text);
}

// The catalogs' sets are bits of FilterCatalogs, which RDKit casts to the
// enum: a value outside the enum's range (its bits, 0 to 2^14 - 1) is
// undefined, and addCatalog ignores the bits outside ALL anyway.
FilterCatalogParams::FilterCatalogs catalogs_of(std::uint32_t bits) {
  return static_cast<FilterCatalogParams::FilterCatalogs>(
      bits & FilterCatalogParams::ALL);
}

// RDKit's RunFilterCatalog returns this entry for SMILES that do not parse.
const boost::shared_ptr<FilterCatalogEntry> &bad_smiles_entry() {
  static const boost::shared_ptr<FilterCatalogEntry> entry(
      new FilterCatalogEntry("no valid RDKit molecule", MatcherPtr()));
  return entry;
}

}  // namespace

std::unique_ptr<FilterMatcherHandle> filter_matcher_new(
    rust::Str name, rust::Box<FilterMatcherCallback> callback) {
  return handle(boost::make_shared<RustFilterMatcher>(
      std::string(name),
      std::make_shared<rust::Box<FilterMatcherCallback>>(std::move(callback))));
}

std::unique_ptr<FilterMatcherHandle> filter_matcher_copy(
    const FilterMatcherHandle &matcher) {
  return handle(matcher.ptr->copy());
}

std::unique_ptr<std::string> filter_matcher_get_name(
    const FilterMatcherHandle &matcher) {
  return copy_string(matcher.ptr->getName());
}

bool filter_matcher_is_valid(const FilterMatcherHandle &matcher) {
  return !is_empty_node(*matcher.ptr) && matcher.ptr->isValid();
}

bool filter_matcher_has_match(const FilterMatcherHandle &matcher,
                              const RWMol &mol) {
  if (is_empty_node(*matcher.ptr)) {
    throw ValueErrorException("a FilterHierarchyMatcher without a pattern");
  }
  HoldControlC control_c;
  return matcher.ptr->hasMatch(mol);
}

// As FilterMatcherBaseGetMatches in Wrap/FilterCatalog.cpp: the matches
// only when it matches.
std::unique_ptr<FilterMatchList> filter_matcher_get_matches(
    const FilterMatcherHandle &matcher, const RWMol &mol) {
  if (is_empty_node(*matcher.ptr)) {
    throw ValueErrorException("a FilterHierarchyMatcher without a pattern");
  }
  std::vector<FilterMatch> matches;
  HoldControlC control_c;
  if (!matcher.ptr->getMatches(mol, matches)) {
    matches.clear();
  }
  return match_list(std::move(matches));
}

std::unique_ptr<FilterMatcherHandle> smarts_matcher_new(rust::Str name) {
  return handle(boost::make_shared<SmartsMatcher>(std::string(name)));
}

std::unique_ptr<FilterMatcherHandle> smarts_matcher_from_mol(const RWMol &mol) {
  return handle(boost::make_shared<SmartsMatcher>(mol));
}

std::unique_ptr<FilterMatcherHandle> smarts_matcher_new_mol(
    rust::Str name, const RWMol &mol, std::uint32_t min_count,
    std::uint32_t max_count) {
  return handle(boost::make_shared<SmartsMatcher>(std::string(name), mol,
                                                  min_count, max_count));
}

std::unique_ptr<FilterMatcherHandle> smarts_matcher_new_smarts(
    rust::Str name, rust::Str smarts, std::uint32_t min_count,
    std::uint32_t max_count) {
  return handle(boost::make_shared<SmartsMatcher>(
      std::string(name), std::string(smarts), min_count, max_count));
}

std::unique_ptr<RWMol> smarts_matcher_get_pattern(
    const FilterMatcherHandle &matcher) {
  const auto &pattern = as<SmartsMatcher>(matcher).getPattern();
  if (!pattern) {
    return nullptr;
  }
  return copy_pattern(*pattern);
}

void smarts_matcher_set_pattern_mol(FilterMatcherHandle &matcher,
                                    const RWMol &mol) {
  as<SmartsMatcher>(matcher).setPattern(static_cast<const RDKit::ROMol &>(mol));
}

void smarts_matcher_set_pattern_smarts(FilterMatcherHandle &matcher,
                                       rust::Str smarts) {
  as<SmartsMatcher>(matcher).setPattern(std::string(smarts));
}

std::uint32_t smarts_matcher_get_min_count(const FilterMatcherHandle &matcher) {
  return as<SmartsMatcher>(matcher).getMinCount();
}

void smarts_matcher_set_min_count(FilterMatcherHandle &matcher,
                                  std::uint32_t count) {
  as<SmartsMatcher>(matcher).setMinCount(count);
}

std::uint32_t smarts_matcher_get_max_count(const FilterMatcherHandle &matcher) {
  return as<SmartsMatcher>(matcher).getMaxCount();
}

void smarts_matcher_set_max_count(FilterMatcherHandle &matcher,
                                  std::uint32_t count) {
  as<SmartsMatcher>(matcher).setMaxCount(count);
}

std::unique_ptr<FilterMatcherHandle> exclusion_list_new() {
  return handle(boost::make_shared<ExclusionList>());
}

// RDKit checks that the pattern is valid, which an empty hierarchy node
// cannot say.
void exclusion_list_add_pattern(FilterMatcherHandle &list,
                                const FilterMatcherHandle &base) {
  check_use(base);
  as<ExclusionList>(list).addPattern(*base.ptr);
}

std::unique_ptr<FilterMatcherList> filter_matcher_list_new() {
  return std::make_unique<FilterMatcherList>();
}

// As SetOffPatterns in Wrap/FilterCatalog.cpp, which copies each.
void filter_matcher_list_push(FilterMatcherList &list,
                              const FilterMatcherHandle &matcher) {
  check_use(matcher);
  list.matchers.push_back(matcher.ptr->copy());
}

void exclusion_list_set_exclusion_patterns(FilterMatcherHandle &list,
                                           const FilterMatcherList &patterns) {
  as<ExclusionList>(list).setExclusionPatterns(patterns.matchers);
}

std::unique_ptr<FilterMatcherHandle> filter_hierarchy_matcher_new() {
  return handle(boost::make_shared<FilterHierarchyMatcher>());
}

std::unique_ptr<FilterMatcherHandle> filter_hierarchy_matcher_from_matcher(
    const FilterMatcherHandle &matcher) {
  check_use(matcher);
  return handle(boost::make_shared<FilterHierarchyMatcher>(*matcher.ptr));
}

void filter_hierarchy_matcher_set_pattern(FilterMatcherHandle &node,
                                          const FilterMatcherHandle &matcher) {
  check_use(matcher);
  as<FilterHierarchyMatcher>(node).setPattern(*matcher.ptr);
}

std::unique_ptr<FilterMatcherHandle> filter_hierarchy_matcher_add_child(
    FilterMatcherHandle &node, const FilterMatcherHandle &hierarchy) {
  // RDKit checks the child has a valid pattern.
  return handle(as<FilterHierarchyMatcher>(node).addChild(
      as<FilterHierarchyMatcher>(hierarchy)));
}

std::unique_ptr<FilterMatcherHandle> filter_match_ops_and(
    const FilterMatcherHandle &arg1, const FilterMatcherHandle &arg2) {
  check_use(arg1);
  check_use(arg2);
  return handle(boost::make_shared<And>(*arg1.ptr, *arg2.ptr));
}

std::unique_ptr<FilterMatcherHandle> filter_match_ops_or(
    const FilterMatcherHandle &arg1, const FilterMatcherHandle &arg2) {
  check_use(arg1);
  check_use(arg2);
  return handle(boost::make_shared<Or>(*arg1.ptr, *arg2.ptr));
}

std::unique_ptr<FilterMatcherHandle> filter_match_ops_not(
    const FilterMatcherHandle &arg1) {
  check_use(arg1);
  return handle(boost::make_shared<Not>(*arg1.ptr));
}

std::size_t filter_match_list_size(const FilterMatchList &list) {
  return list.matches.size();
}

std::unique_ptr<FilterMatcherHandle> filter_match_list_matcher(
    const FilterMatchList &list, std::size_t i) {
  return handle(list.matches.at(i).filterMatch);
}

rust::Vec<IndexPair> filter_match_list_atom_pairs(const FilterMatchList &list,
                                                  std::size_t i) {
  rust::Vec<IndexPair> pairs;
  for (const auto &pair : list.matches.at(i).atomPairs) {
    pairs.push_back(IndexPair{static_cast<std::uint32_t>(pair.first),
                              static_cast<std::uint32_t>(pair.second)});
  }
  return pairs;
}

void filter_match_list_push(FilterMatchList &list,
                            const FilterMatcherHandle &matcher,
                            rust::Slice<const IndexPair> atom_pairs) {
  RDKit::MatchVectType pairs;
  for (const auto &pair : atom_pairs) {
    pairs.emplace_back(static_cast<int>(pair.first),
                       static_cast<int>(pair.second));
  }
  list.matches.emplace_back(matcher.ptr, std::move(pairs));
}

std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_new() {
  auto result = std::make_unique<FilterCatalogEntryHandle>();
  result->ptr = boost::make_shared<FilterCatalogEntry>();
  return result;
}

std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_new_matcher(
    rust::Str name, const FilterMatcherHandle &matcher) {
  check_use(matcher);
  auto result = std::make_unique<FilterCatalogEntryHandle>();
  result->ptr =
      boost::make_shared<FilterCatalogEntry>(std::string(name), *matcher.ptr);
  return result;
}

std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_share(
    const FilterCatalogEntryHandle &entry) {
  auto result = std::make_unique<FilterCatalogEntryHandle>();
  result->ptr = entry.ptr;
  return result;
}

bool filter_catalog_entry_same(const FilterCatalogEntryHandle &a,
                               const FilterCatalogEntryHandle &b) {
  return a.ptr == b.ptr;
}

bool filter_catalog_entry_is_valid(const FilterCatalogEntryHandle &entry) {
  return entry.ptr->isValid();
}

std::unique_ptr<std::string> filter_catalog_entry_get_description(
    const FilterCatalogEntryHandle &entry) {
  return copy_string(entry.ptr->getDescription());
}

void filter_catalog_entry_set_description(FilterCatalogEntryHandle &entry,
                                          rust::Str description) {
  entry_mut(entry).setDescription(std::string(description));
}

// As FilterCatalogEntryGetMatches in Wrap/FilterCatalog.cpp.
std::unique_ptr<FilterMatchList> filter_catalog_entry_get_filter_matches(
    const FilterCatalogEntryHandle &entry, const RWMol &mol) {
  std::vector<FilterMatch> matches;
  HoldControlC control_c;
  if (!entry.ptr->getFilterMatches(mol, matches)) {
    matches.clear();
  }
  return match_list(std::move(matches));
}

bool filter_catalog_entry_has_filter_match(
    const FilterCatalogEntryHandle &entry, const RWMol &mol) {
  HoldControlC control_c;
  return entry.ptr->hasFilterMatch(mol);
}

std::unique_ptr<std::vector<std::string>> filter_catalog_entry_get_prop_list(
    const FilterCatalogEntryHandle &entry) {
  return std::make_unique<std::vector<std::string>>(entry.ptr->getPropList());
}

void filter_catalog_entry_set_prop(FilterCatalogEntryHandle &entry,
                                   rust::Str key, rust::Str val) {
  entry_mut(entry).setProp<std::string>(std::string(key), std::string(val));
}

std::unique_ptr<std::string> filter_catalog_entry_get_prop(
    const FilterCatalogEntryHandle &entry, rust::Str key) {
  return copy_string(entry.ptr->getProp<std::string>(std::string(key)));
}

void filter_catalog_entry_clear_prop(FilterCatalogEntryHandle &entry,
                                     rust::Str key) {
  entry_mut(entry).clearProp(std::string(key));
}

std::size_t filter_catalog_entry_list_size(const FilterCatalogEntryList &list) {
  return list.entries.size();
}

std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_entry_list_get(
    const FilterCatalogEntryList &list, std::size_t i) {
  return entry_handle(list.entries.at(i));
}

rust::Vec<std::size_t> filter_catalog_entry_list_row_ends(
    const FilterCatalogEntryList &list) {
  return to_vec(list.row_ends);
}

std::unique_ptr<FilterCatalogParams> filter_catalog_params_new() {
  return std::make_unique<FilterCatalogParams>();
}

std::unique_ptr<FilterCatalogParams> filter_catalog_params_from_catalogs(
    std::uint32_t catalogs) {
  return std::make_unique<FilterCatalogParams>(catalogs_of(catalogs));
}

bool filter_catalog_params_add_catalog(FilterCatalogParams &params,
                                       std::uint32_t catalogs) {
  return params.addCatalog(catalogs_of(catalogs));
}

rust::Vec<std::uint32_t> filter_catalogs_values() {
  using P = FilterCatalogParams;
  return to_vec(std::vector<std::uint32_t>{
      P::PAINS_A, P::PAINS_B, P::PAINS_C, P::PAINS, P::BRENK, P::NIH, P::ZINC,
      P::CHEMBL_Glaxo, P::CHEMBL_Dundee, P::CHEMBL_BMS, P::CHEMBL_SureChEMBL,
      P::CHEMBL_MLSMR, P::CHEMBL_Inpharmatica, P::CHEMBL_LINT, P::CHEMBL,
      P::ALL});
}

std::unique_ptr<FilterCatalog> filter_catalog_new() {
  return std::make_unique<FilterCatalog>();
}

std::unique_ptr<FilterCatalog> filter_catalog_from_params(
    const FilterCatalogParams &params) {
  return std::make_unique<FilterCatalog>(params);
}

std::unique_ptr<FilterCatalog> filter_catalog_from_catalogs(
    std::uint32_t catalogs) {
  return std::make_unique<FilterCatalog>(catalogs_of(catalogs));
}

// The pickle is a Boost.Serialization archive, which Boost reads trusting
// its class and object ids, so the Rust caller promises bytes RDKit wrote
// (CONVENTIONS.md U2). Those can still hold matchers RDKit's methods
// dereference or recurse through forever, and patterns with monomer infos
// its pickler misreads, which a Load walk refuses.
std::unique_ptr<FilterCatalog> filter_catalog_from_pickle(
    rust::Slice<const std::uint8_t> pickle) {
  check_can_serialize();
  auto catalog = std::make_unique<FilterCatalog>(to_std(pickle));
  for (unsigned i = 0; i < catalog->getNumEntries(); ++i) {
    check_entry(*catalog->getEntryWithIdx(i), Walk::Load);
  }
  return catalog;
}

rust::Vec<std::uint8_t> filter_catalog_serialize(const FilterCatalog &catalog) {
  check_can_serialize();
  for (unsigned i = 0; i < catalog.getNumEntries(); ++i) {
    check_entry(*catalog.getEntryWithIdx(i), Walk::Pickle);
  }
  return to_bytes(catalog.Serialize());
}

rust::Vec<std::uint8_t> filter_catalog_entry_serialize(
    const FilterCatalogEntryHandle &entry) {
  check_can_serialize();
  check_entry(*entry.ptr, Walk::Pickle);
  return to_bytes(entry.ptr->Serialize());
}

// As filter_catalog_add_entry in Wrap/FilterCatalog.cpp, the catalog gets
// a copy.
void filter_catalog_add_entry(FilterCatalog &catalog,
                              const FilterCatalogEntryHandle &entry) {
  catalog.addEntry(boost::make_shared<FilterCatalogEntry>(*entry.ptr));
}

bool filter_catalog_remove_entry(FilterCatalog &catalog, std::uint32_t idx) {
  return catalog.removeEntry(idx);
}

// As FilterCatalogRemoveEntry: by identity, so only an entry the catalog
// returned (and not changed since, which copies it) is found.
bool filter_catalog_remove_entry_entry(FilterCatalog &catalog,
                                       const FilterCatalogEntryHandle &entry) {
  return catalog.removeEntry(catalog.getIdxForEntry(entry.ptr.get()));
}

std::uint32_t filter_catalog_get_num_entries(const FilterCatalog &catalog) {
  return catalog.getNumEntries();
}

std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_get_entry(
    const FilterCatalog &catalog, std::uint32_t idx) {
  check_index(idx, catalog.getNumEntries());
  return entry_handle(catalog.getEntry(idx));
}

bool filter_catalog_has_match(const FilterCatalog &catalog, const RWMol &mol) {
  HoldControlC control_c;
  return catalog.hasMatch(mol);
}

std::unique_ptr<FilterCatalogEntryHandle> filter_catalog_get_first_match(
    const FilterCatalog &catalog, const RWMol &mol) {
  boost::shared_ptr<const FilterCatalogEntry> entry;
  {
    HoldControlC control_c;
    entry = catalog.getFirstMatch(mol);
  }
  if (!entry) {
    return nullptr;
  }
  return entry_handle(std::move(entry));
}

std::unique_ptr<FilterCatalogEntryList> filter_catalog_get_matches(
    const FilterCatalog &catalog, const RWMol &mol) {
  auto list = std::make_unique<FilterCatalogEntryList>();
  std::vector<FilterCatalog::CONST_SENTRY> entries;
  {
    HoldControlC control_c;
    entries = catalog.getMatches(mol);
  }
  for (auto &entry : entries) {
    list->entries.push_back(
        boost::const_pointer_cast<FilterCatalogEntry>(entry));
  }
  list->row_ends.push_back(list->entries.size());
  return list;
}

std::unique_ptr<FilterMatchList> filter_catalog_get_filter_matches(
    const FilterCatalog &catalog, const RWMol &mol) {
  HoldControlC control_c;
  return match_list(catalog.getFilterMatches(mol));
}

bool filter_catalog_can_serialize() {
  return RDKit::FilterCatalogCanSerialize();
}

// RDKit's RunFilterCatalog, on the crate's threads (threads.h). RDKit's own
// negates a nonpositive thread count in int (undefined for INT_MIN) and
// steps each thread's index in unsigned int, which wraps for 2^32 SMILES, so
// threads process the same SMILES at once, and its threads' matches would
// race on ControlCHandler's statics, which control_c_mutex cannot order.
// Here each thread parses its SMILES and matches each molecule holding
// control_c_mutex, so only the parsing runs in parallel. As RDKit's, a
// SMILES that does not parse gets the "no valid RDKit molecule" entry, and
// a sanitization error is the call's error.
std::unique_ptr<FilterCatalogEntryList> run_filter_catalog(
    const FilterCatalog &catalog, rust::Slice<const rust::Str> smiles,
    std::int32_t num_threads) {
  std::vector<std::string> texts;
  texts.reserve(smiles.size());
  for (auto text : smiles) {
    texts.emplace_back(text);
  }
  std::vector<std::vector<FilterCatalog::CONST_SENTRY>> results(texts.size());
  // More threads than SMILES would have nothing to do.
  auto n = static_cast<std::size_t>(num_threads_to_use(num_threads));
  n = std::max<std::size_t>(1, std::min(n, texts.size()));
  run_on_threads(static_cast<unsigned>(n), [&](unsigned t) {
    for (std::size_t i = t; i < texts.size(); i += n) {
      std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(texts[i]));
      if (mol) {
        HoldControlC control_c;
        results[i] = catalog.getMatches(*mol);
      } else {
        results[i].push_back(bad_smiles_entry());
      }
    }
  });
  auto list = std::make_unique<FilterCatalogEntryList>();
  for (auto &row : results) {
    for (auto &entry : row) {
      list->entries.push_back(
          boost::const_pointer_cast<FilterCatalogEntry>(entry));
    }
    list->row_ends.push_back(list->entries.size());
  }
  return list;
}

// RDKit returns its global hierarchy, which Python then changes in place;
// this is a catalog sharing its entries, which are not changed (entry_mut
// copies them).
std::unique_ptr<FilterCatalog> get_functional_group_hierarchy() {
  const auto &hierarchy = RDKit::GetFunctionalGroupHierarchy();
  auto catalog = std::make_unique<FilterCatalog>();
  for (unsigned i = 0; i < hierarchy.getNumEntries(); ++i) {
    catalog->addEntry(
        boost::const_pointer_cast<FilterCatalogEntry>(hierarchy.getEntry(i)));
  }
  return catalog;
}

// The patterns are the hierarchy's own, which its matchers (and their
// copies, in other threads' catalogs) match with, so they are copied
// locked.
std::unique_ptr<MolList> get_flattened_functional_group_hierarchy(
    bool normalized) {
  auto list = std::make_unique<MolList>();
  for (const auto &[name, pattern] :
       RDKit::GetFlattenedFunctionalGroupHierarchy(normalized)) {
    list->push_named(name, copy_pattern(*pattern));
  }
  return list;
}

}  // namespace rdkit_cxx
