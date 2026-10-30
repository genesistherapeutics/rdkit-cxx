// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/inchi.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/inchi.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/inchi.rs.h"

// The InChI library keeps no state between calls: inchi_dll.c works on the
// caller's input and its own locals, so threads may call it at once (its
// mol2inchi demo runs it on pthreads). inchi.cpp has no mutable statics
// either; the static its findAlternatingBonds comment mentions is now the
// caller's _visited set, and its recursion is at most seven bonds deep.

namespace rdkit_cxx {

namespace {

// inchi_api.h's limits, for the fields MolToInchi narrows into.
// AT_NUM (num_atoms, neighbor[], num_stereo0D) is a signed short, and InChI
// takes at most MAX_ATOMS (32766, with LargeMolecules) atoms.
constexpr unsigned kMaxInchiAtoms = 32766;
constexpr std::size_t kMaxInchiStereo0D = 32767;
// num_iso_H is a signed char.
constexpr unsigned kMaxInchiHs = 127;
// MAXVAL: the entries of inchi_Atom's neighbor, bond_type, and bond_stereo.
constexpr unsigned kMaxInchiBonds = 20;
// ISOTOPIC_SHIFT_MAX: isotopic_mass is ISOTOPIC_SHIFT_FLAG (10000) plus the
// isotope's difference from the average mass, and InChI reads a value
// outside 10000 +/- 100 as an absolute mass instead (inchi_dll.c), whose
// difference it then casts to a signed char.
constexpr std::int64_t kMaxIsotopeShift = 100;
// See check_coordinates.
constexpr double kMaxInchiCoordinate = 1e18;
// The longest timeouts accepted, in ms; see check_timeout.
constexpr std::uint64_t kMaxTimeoutMsec = 1000000000000ULL;

// The structure-perception options, which InChI's parser
// (set_common_options_by_parg and ReadCommandLineParms in ichiparm.c) matches
// whole and ignoring case. Other options steer InChI's file input and output,
// and some of those corrupt memory through the API RDKit calls: "/OutputSDF
// /Key" makes GetINCHI free() a pointer it never allocated.
constexpr std::array<std::string_view, 22> kInchiOptions = {
    "SNon",      "SRel",         "SRac",
    "SUCF",      "ChiralFlagON", "ChiralFlagOFF",
    "NEWPSOFF",  "DoNotAddH",    "SUU",
    "SLUUD",     "FixedH",       "RecMet",
    "KET",       "15T",          "LargeMolecules",
    "Polymers",  "Polymers105",  "NPZz",
    "SAtZZ",     "AuxNone",      "WarnOnEmptyStructure",
    "NoWarnings"};

bool equals_ignoring_case(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    auto lower = [](char c) {
      return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
    };
    if (lower(a[i]) != lower(b[i])) {
      return false;
    }
  }
  return true;
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// Reads digits from text at pos, as a value capped just past limit.
std::uint64_t read_digits(std::string_view text, std::size_t &pos,
                          std::uint64_t limit) {
  std::uint64_t value = 0;
  while (pos < text.size() && is_digit(text[pos])) {
    value = std::min(value * 10 + static_cast<unsigned>(text[pos] - '0'),
                     limit + 1);
    ++pos;
  }
  return value;
}

// A timeout: W<seconds>, digits with an optional fraction, or
// WM<milliseconds>, digits. InChI reads seconds with strtod and converts
// t * 1000 to long, undefined for NaN ("Wnan") and for 2^63 exactly, which
// its range check lets through; it then converts the milliseconds to clock
// ticks (INCHI_CLOCK_T in ichicano.c, a double cast to clock_t) and adds
// them to the clock, which overflows for timeouts near LONG_MAX. Timeouts
// are therefore plain decimals of at most 10^12 ms (about 31 years).
bool is_timeout(std::string_view name) {
  if (name.size() < 2 || (name[0] != 'W' && name[0] != 'w')) {
    return false;
  }
  bool msec = name[1] == 'M' || name[1] == 'm';
  std::size_t pos = msec ? 2 : 1;
  std::size_t start = pos;
  std::uint64_t limit = msec ? kMaxTimeoutMsec : kMaxTimeoutMsec / 1000;
  std::uint64_t value = read_digits(name, pos, limit);
  if (pos == start || value > limit) {
    return false;
  }
  if (!msec && pos < name.size() && name[pos] == '.') {
    ++pos;
    std::size_t fraction = pos;
    read_digits(name, pos, 0);
    if (pos == fraction) {
      return false;
    }
  }
  return pos == name.size();
}

// InChI parses the options itself (parse_options_string in inchi_dll.c,
// which splits at spaces and tabs), so only options known to be safe pass:
// each token is '/' or '-' and an option above, or a timeout.
void check_options(std::string_view options) {
  std::size_t pos = 0;
  while (pos < options.size()) {
    if (options[pos] == ' ' || options[pos] == '\t') {
      ++pos;
      continue;
    }
    std::size_t end = options.find_first_of(" \t", pos);
    if (end == std::string_view::npos) {
      end = options.size();
    }
    std::string_view token = options.substr(pos, end - pos);
    pos = end;
    if (token[0] == '/' || token[0] == '-') {
      std::string_view name = token.substr(1);
      if (is_timeout(name)) {
        continue;
      }
      bool known = false;
      for (auto option : kInchiOptions) {
        known = known || equals_ignoring_case(name, option);
      }
      if (known) {
        continue;
      }
    }
    throw ValueErrorException(
        "unsupported InChI option '" + std::string(token) +
        "': options are the structure-perception options SNon, SRel, SRac, "
        "SUCF, ChiralFlagON, ChiralFlagOFF, NEWPSOFF, DoNotAddH, SUU, SLUUD, "
        "FixedH, RecMet, KET, 15T, LargeMolecules, Polymers, Polymers105, "
        "NPZz, SAtZZ, AuxNone, WarnOnEmptyStructure, NoWarnings, and the "
        "timeouts W<seconds> and WM<milliseconds> (at most 10^12 ms), each "
        "prefixed with / or -");
  }
}

bool is_tetrahedral(const RDKit::Atom &atom) {
  return atom.getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW ||
         atom.getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW;
}

// InChI computes with the coordinates' differences, products, and norms in
// double and converts angles and directions made from them to integers
// (GetMinDistDistribution in strutil.c, which indexes an array with the
// result, and the z_dir values in ichister.c), which is undefined for NaN.
// A non-finite coordinate makes NaN, and so does one large enough for those
// products to overflow (inf - inf, 0 * inf); below 1e18, products of up to
// sixteen coordinates stay finite.
void check_coordinates(const RDKit::ROMol &mol) {
  if (!mol.getNumConformers()) {
    return;
  }
  const RDKit::Conformer &conf = **mol.beginConformers();
  for (unsigned i = 0; i < mol.getNumAtoms(); ++i) {
    const RDGeom::Point3D &p = conf.getAtomPos(i);
    for (double c : {p.x, p.y, p.z}) {
      if (!std::isfinite(c) || std::fabs(c) > kMaxInchiCoordinate) {
        throw ValueErrorException(
            "atom " + std::to_string(i) +
            " has a coordinate InChI cannot use: coordinates must be finite "
            "and at most 1e18 in magnitude");
      }
    }
  }
}

// An atom's implicit hydrogens as MolToInchi counts them on its kekulized
// copy (calcImplicitValence, then getNumImplicitHs); 0 when the copy could
// not be kekulized, as MolToInchi then fails first, or the count fails.
unsigned implicit_hs(std::optional<RDKit::RWMol> &kekulized, unsigned idx) {
  if (!kekulized) {
    return 0;
  }
  RDKit::Atom *atom = kekulized->getAtomWithIdx(idx);
  try {
    atom->calcImplicitValence();
  } catch (const std::exception &) {
    return 0;  // MolToInchi throws there too
  }
  return atom->getNumImplicitHs();
}

// What MolToInchi narrows into InChI's fields (inchi_Atom and inchi_Input in
// inchi_api.h) without checking it, read as MolToInchi reads it.
void check_inchi_input(const RWMol &mol) {
  // num_atoms and every neighbor index are signed shorts, which more atoms
  // would wrap.
  if (mol.getNumAtoms() > kMaxInchiAtoms) {
    throw ValueErrorException(
        "InChI takes at most 32766 atoms; the molecule has " +
        std::to_string(mol.getNumAtoms()));
  }
  check_coordinates(mol);

  // MolToInchi reads hydrogen counts from its copy, after updating the
  // property cache when the molecule needs it; so does this one.
  RDKit::RWMol copy(mol);
  if (mol.needsUpdatePropertyCache()) {
    copy.updatePropertyCache(false);
  }
  // MolToInchi kekulizes its copy (failing, before any stereo, if it cannot)
  // and counts implicit hydrogens on it.
  std::optional<RDKit::RWMol> kekulized(copy);
  try {
    RDKit::MolOps::Kekulize(*kekulized, false);
  } catch (const RDKit::MolSanitizeException &) {
    kekulized.reset();
  }
  const RDKit::PeriodicTable *table = RDKit::PeriodicTable::getTable();
  std::size_t num_stereo = 0;
  for (const RDKit::Atom *atom : copy.atoms()) {
    unsigned idx = atom->getIdx();
    // isotopic_mass is ISOTOPIC_SHIFT_FLAG + isotope - round(average mass),
    // computed in int (which a large isotope overflows) and stored in a
    // short. A shift beyond ISOTOPIC_SHIFT_MAX is read as an absolute mass,
    // so isotope 113 for carbon gives /i1+105 and 70000 gives /i1+104.
    if (unsigned isotope = atom->getIsotope()) {
      std::int64_t average =
          static_cast<int>(table->getAtomicWeight(atom->getAtomicNum()) + 0.5);
      std::int64_t shift = static_cast<std::int64_t>(isotope) - average;
      if (shift > kMaxIsotopeShift || shift < -kMaxIsotopeShift) {
        throw ValueErrorException(
            "atom " + std::to_string(idx) + "'s isotope " +
            std::to_string(isotope) +
            " is more than 100 from its element's average mass, which InChI "
            "cannot represent");
      }
    }
    // num_iso_H[0] is a signed char. MolToInchi lets InChI count the
    // hydrogens of C, N, O, F, Cl, Br, and I (-1), and passes every other
    // atom's, and every radical's, total, which 128 to 255 make negative
    // (Python raises) and larger counts wrap around again.
    switch (atom->getAtomicNum()) {
      case 6:
      case 7:
      case 8:
      case 9:
      case 17:
      case 35:
      case 53:
        if (!atom->getNumRadicalElectrons()) {
          break;
        }
        [[fallthrough]];
      default:
        if (atom->getTotalNumHs() > kMaxInchiHs) {
          throw ValueErrorException(
              "atom " + std::to_string(idx) + " has " +
              std::to_string(atom->getTotalNumHs()) +
              " hydrogens; InChI takes at most 127 per atom");
        }
    }
    // MolToInchi (RDKit 2024.03, External/INCHI-API/inchi.cpp) makes a
    // stereo0D for every atom with a chiral tag other than CHI_UNSPECIFIED,
    // or with a "molParity" property, and with at least three neighbors and
    // implicit hydrogens, writing each neighbor into its four-entry neighbor
    // array: five or more neighbors write past it (2026.03 ignores such
    // atoms).
    if ((atom->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED ||
         atom->hasProp("molParity")) &&
        atom->getDegree() > 4) {
      throw ValueErrorException(
          "atom " + std::to_string(idx) +
          " has a chiral tag or molParity and " +
          std::to_string(atom->getDegree()) +
          " neighbors; RDKit 2024.03's InChI writer takes at most 4");
    }
    if (is_tetrahedral(*atom)) {
      ++num_stereo;
      // With fewer than three neighbor atoms it fills in only those (and the
      // atom itself for the fourth when there are three), and passes a CW
      // or CCW center whose implicit hydrogens (counted on its kekulized
      // copy, after calcImplicitValence) bring the total to three or more,
      // of any size, with uninitialized neighbor indices (2026.03 ignores
      // totals above four). Counting all hydrogens, as the total degree
      // does, refuses those and centers of explicit hydrogens too.
      if (atom->getDegree() < 3 &&
          (atom->getTotalDegree() >= 3 ||
           atom->getDegree() + implicit_hs(kekulized, idx) >= 3)) {
        throw ValueErrorException(
            "atom " + std::to_string(idx) +
            " is a tetrahedral stereo center with " +
            std::to_string(atom->getDegree()) +
            " neighbor atoms and implicit hydrogens; InChI needs at least 3 "
            "neighbor atoms");
      }
    }
  }
  // MolToInchi stores each bond in the neighbor arrays of its lower-index
  // atom, which hold MAXVAL (20) entries, without a bound (2026.03 fails
  // instead): more overflow the heap array of inchi_Atoms.
  std::vector<unsigned> lower_bonds(copy.getNumAtoms(), 0);
  for (const RDKit::Bond *bond : copy.bonds()) {
    unsigned lower = std::min(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    if (++lower_bonds[lower] > kMaxInchiBonds) {
      throw ValueErrorException(
          "atom " + std::to_string(lower) +
          " is the lower-index atom of more than 20 bonds; RDKit 2024.03's "
          "InChI writer stores at most 20 per atom");
    }
  }
  for (const RDKit::Bond *bond : copy.bonds()) {
    if (bond->getStereo() > RDKit::Bond::STEREOANY &&
        bond->getStereoAtoms().size() >= 2) {
      // MolToInchi passes the first two stereo atoms as stereo0D neighbors
      // and checks only the first against the molecule.
      for (std::size_t i = 0; i < 2; ++i) {
        int stereo_atom = bond->getStereoAtoms()[i];
        if (stereo_atom < 0 ||
            static_cast<unsigned>(stereo_atom) >= copy.getNumAtoms()) {
          throw ValueErrorException(
              "bond " + std::to_string(bond->getIdx()) + "'s stereo atom " +
              std::to_string(stereo_atom) + " is not in the molecule");
        }
      }
      ++num_stereo;
    } else if (bond->getStereo() == RDKit::Bond::STEREOANY) {
      ++num_stereo;
    }
  }
  // num_stereo0D is a signed short too. The count is an upper bound on the
  // stereo0D elements MolToInchi makes: one per tetrahedral center and one
  // per double bond with stereo.
  if (num_stereo > kMaxInchiStereo0D) {
    throw ValueErrorException(
        "InChI takes at most 32767 stereo elements; the molecule has up to " +
        std::to_string(num_stereo));
  }
}

}  // namespace

namespace {

// InChI keeps each atom's hydrogen counts, isotopic ones too, in signed
// chars (inchi_Atom::num_iso_H), and InchiToMol copies the isotopic ones
// into unsigned ints (ISOTOPES_t in inchi.cpp): a count of 128 or more
// reads as about 4.3 billion, and InchiToMol adds that many hydrogen atoms.
// The counts follow H, D, or T in the hydrogen layers (/h) and the isotopic
// ones (/i and the sublayers after it); the formula, with its total count,
// is the second layer, and the fixed-H (/f) and reconnected (/r) layers
// open with formulas too.
void check_inchi_h_counts(std::string_view inchi) {
  std::size_t layer = 0;
  std::size_t start = 0;
  while (start <= inchi.size()) {
    auto end = inchi.find('/', start);
    if (end == std::string_view::npos) {
      end = inchi.size();
    }
    auto text = inchi.substr(start, end - start);
    if (layer >= 2 && !text.empty() && (text[0] == 'h' || text[0] == 'i')) {
      for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] != 'H' && text[i] != 'D' && text[i] != 'T') {
          continue;
        }
        unsigned count = 0;
        std::size_t j = i + 1;
        for (; j < text.size() && text[j] >= '0' && text[j] <= '9'; ++j) {
          count = std::min(count * 10 + (text[j] - '0'), 1000u);
        }
        if (count >= 128) {
          throw ValueErrorException(
              "InChI keeps an atom's hydrogen counts in a signed char, so it "
              "cannot hold " +
              std::string(text.substr(i, j - i)));
        }
        i = j - 1;
      }
    }
    ++layer;
    start = end + 1;
  }
}

}  // namespace

InchiToMolData inchi_to_mol(rust::Str inchi, bool sanitize, bool remove_hs) {
  check_inchi_h_counts(std::string_view(inchi.data(), inchi.size()));
  RDKit::ExtraInchiReturnValues rv;
  std::unique_ptr<RWMol> mol;
  {
    // InchiToMol's cleanup matches patterns, and each match constructs a
    // ControlCHandler (control_c_lock.h, CONVENTIONS.md U4).
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    mol.reset(RDKit::InchiToMol(std::string(inchi), rv, sanitize, remove_hs));
  }
  return InchiToMolData{std::move(mol), rv.returnCode,
                        std::make_unique<std::string>(rv.messagePtr),
                        std::make_unique<std::string>(rv.logPtr)};
}

MolToInchiData mol_to_inchi(const RWMol &mol, rust::Str options) {
  std::string opts(options);
  check_options(opts);
  check_inchi_input(mol);
  RDKit::ExtraInchiReturnValues rv;
  std::string inchi;
  {
    // MolToInchi's reverse cleanup matches patterns (U4).
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    inchi = RDKit::MolToInchi(mol, rv, opts.empty() ? nullptr : opts.c_str());
  }
  return MolToInchiData{std::make_unique<std::string>(std::move(inchi)),
                        rv.returnCode,
                        std::make_unique<std::string>(rv.messagePtr),
                        std::make_unique<std::string>(rv.logPtr),
                        std::make_unique<std::string>(rv.auxInfoPtr)};
}

std::unique_ptr<std::string> mol_to_inchi_key(const RWMol &mol,
                                              rust::Str options) {
  std::string opts(options);
  check_options(opts);
  check_inchi_input(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return std::make_unique<std::string>(
      RDKit::MolToInchiKey(mol, opts.empty() ? nullptr : opts.c_str()));
}

std::unique_ptr<std::string> inchi_to_inchi_key(rust::Str inchi) {
  return std::make_unique<std::string>(
      RDKit::InchiToInchiKey(std::string(inchi)));
}

}  // namespace rdkit_cxx
