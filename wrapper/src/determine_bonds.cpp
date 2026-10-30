// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/determine_bonds.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/DetermineBonds/DetermineBonds.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDConfig.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"

#ifdef RDK_BUILD_YAEHMOP_SUPPORT
extern "C" {
// YAeHMOP's built-in parameter table (yaehmop/tightbind/eht_parms.h), one
// line per orbital, ending with "END".
extern const char *defaultParms[];
}
#endif

namespace rdkit_cxx {

namespace {

// RDKit sums formal charges in int (MolOps::getFormalCharge, U9) after
// setAtomicCharges, or the Hückel option, sets some, each stored in 8 bits;
// headroom for every atom's covers them. setAtomicRadicals takes std::abs of
// `charge`, undefined for INT_MIN.
void check_charges(const RWMol &mol, std::int32_t charge) {
  if (charge == INT_MIN) {
    throw ValueErrorException("charge must be above INT_MIN");
  }
  check_formal_charge_sum(mol, 128 * std::int64_t{mol.getNumAtoms()});
}

// ConnectTheDots_Large (FileParsers/ProximityBonds.cpp) hashes each atom's
// coordinates as int(x / 5.45) * 571 + ..., undefined for a coordinate that
// is not finite or whose hash overflows int; within 1e7 it fits.
void check_proximity_coordinates(const RWMol &mol) {
  if (!mol.getNumConformers()) {
    return;  // ConnectTheDots then adds no bonds.
  }
  for (const auto &pos : mol.getConformer().getPositions()) {
    for (double v : {pos.x, pos.y, pos.z}) {
      if (!(std::abs(v) <= 1e7)) {
        throw ValueErrorException(
            "connect-the-dots needs finite coordinates within 1e7");
      }
    }
  }
}

// ConnectTheDots_Large, pruning a hydrogen bonded to several atoms, reads each
// one's monomer info as an AtomPDBResidueInfo without checking even its type
// (CONVENTIONS.md U10), past the end of any other kind. Those are set aside for
// the call, so that RDKit sees no info there, and put back after it, even when
// it throws.
class OtherMonomerInfoSetAside {
 public:
  explicit OtherMonomerInfoSetAside(RWMol &mol) {
    for (auto atom : mol.atoms()) {
      auto *info = atom->getMonomerInfo();
      if (info && !dynamic_cast<RDKit::AtomPDBResidueInfo *>(info)) {
        d_saved.emplace_back(
            atom, std::unique_ptr<RDKit::AtomMonomerInfo>(info->copy()));
        atom->setMonomerInfo(nullptr);
      }
    }
  }
  ~OtherMonomerInfoSetAside() {
    for (auto &[atom, info] : d_saved) {
      atom->setMonomerInfo(info.release());
    }
  }
  OtherMonomerInfoSetAside(const OtherMonomerInfoSetAside &) = delete;
  OtherMonomerInfoSetAside &operator=(const OtherMonomerInfoSetAside &) =
      delete;

 private:
  std::vector<std::pair<RDKit::Atom *, std::unique_ptr<RDKit::AtomMonomerInfo>>>
      d_saved;
};

#ifdef RDK_BUILD_YAEHMOP_SUPPORT

// The Hückel option runs YAeHMOP (EHTTools::runMol), whose C code exits the
// process, or reads past its buffers, for input it does not expect. What it
// reads is checked here first, as fill_atomic_parms (tightbind/fileio.c)
// and check_for_errors (tightbind/genutil.c) read it.

std::string upper(std::string text) {
  for (auto &c : text) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  }
  return text;
}

// The parameter lines, without comments, as runMol and fill_atomic_parms
// choose them: from the file $BIND_PARM_FILE names, else from
// $RDBASE/Data/eht_parms.dat when it opens; from the built-in table when the
// file does not open. Without either variable, fill_atomic_parms opens a
// path fixed when YAeHMOP was built, and when that does not open, as
// conda-forge's does not, it frees its copy of the name twice; that case is
// refused. The second value says whether the lines end with the table's
// "END", which YAeHMOP needs after an element's lines when it reads a file.
std::pair<std::vector<std::string>, bool> eht_parameter_lines() {
  std::string name;
  if (const char *bind = std::getenv("BIND_PARM_FILE")) {
    name = bind;
  } else if (const char *rdbase = std::getenv("RDBASE")) {
    std::string candidate = std::string(rdbase) + "/Data/eht_parms.dat";
    if (std::ifstream(candidate).good()) {
      name = candidate;
    }
  }
  if (name.empty()) {
    throw ValueErrorException(
        "the Hückel method needs RDBASE set to RDKit's share/RDKit directory, "
        "or BIND_PARM_FILE to an EHT parameter file: without them, YAeHMOP "
        "frees memory twice when the parameter file it was built with is "
        "missing");
  }
  std::vector<std::string> lines;
  if (std::ifstream in(name); in) {
    std::string line;
    while (std::getline(in, line)) {
      // fgets reads up to 2047 characters into a 240-byte buffer.
      if (line.size() > 238) {
        throw ValueErrorException(
            "the EHT parameter file has a line longer than YAeHMOP's "
            "239-character buffer");
      }
      // skipcomments skips lines that start, after spaces, with ';' or
      // end there.
      auto start = line.find_first_not_of(' ');
      if (start == std::string::npos || line[start] == ';') {
        continue;
      }
      std::istringstream fields(line);
      std::string label;
      if (!(fields >> label)) {
        continue;
      }
      lines.push_back(line);
    }
    return {lines, false};
  }
  for (const char **line = defaultParms; std::strcmp(*line, "END") != 0;
       ++line) {
    lines.emplace_back(*line);
  }
  return {lines, true};
}

struct EhtElement {
  int valence = 0;
  int orbitals = 0;
};

// An element's valence electrons and orbitals, from the first run of lines
// whose label starts as its symbol does (both upper case, compared over two
// characters, as strncmp(symb, tstring, 2)).
EhtElement eht_element(const std::vector<std::string> &lines, bool ends,
                       const std::string &symbol) {
  auto key = upper(symbol).substr(0, 2);
  bool found = false;
  bool has_orbital[4] = {false, false, false, false};
  int at_number = 0;
  EhtElement element;
  for (std::size_t i = 0; i < lines.size(); ++i) {
    std::istringstream fields(lines[i]);
    std::string label;
    fields >> label;
    // sscanf("%s") reads the label into ten bytes.
    if (label.size() > 9) {
      throw ValueErrorException(
          "the EHT parameters have a label longer than YAeHMOP's buffer");
    }
    bool matches = upper(label).substr(0, 2) == key;
    if (!matches) {
      if (found) {
        break;
      }
      continue;
    }
    found = true;
    int nzeta = 0;
    int nquant = 0;
    std::string ang;
    double hii = 0.0;
    double rest[4];
    if (!(fields >> at_number >> element.valence >> nzeta >> nquant >> ang >>
          hii >> rest[0] >> rest[1] >> rest[2] >> rest[3]) ||
        ang.size() != 1) {
      // sscanf would leave values unset, and reads the orbital into two
      // bytes.
      throw ValueErrorException("the EHT parameters for " + symbol +
                                " are malformed");
    }
    auto orbital = std::string("spdf").find(
        static_cast<char>(std::tolower(static_cast<unsigned char>(ang[0]))));
    if (orbital != std::string::npos) {
      has_orbital[orbital] = hii != 0.0 && nquant != 0;
    }
    if (i + 1 == lines.size() && !ends) {
      // At the end of the file, fill_atomic_parms keeps the last label and
      // reads the same line forever.
      throw ValueErrorException("the EHT parameter file ends with " + symbol +
                                "'s parameters, which YAeHMOP reads forever");
    }
  }
  if (!found) {
    // fill_atomic_parms calls fatal, which exits the process.
    throw ValueErrorException("YAeHMOP has no parameters for " + symbol);
  }
  if (at_number >= 0) {
    const int sizes[4] = {1, 3, 5, 7};
    for (int i = 0; i < 4; ++i) {
      element.orbitals += has_orbital[i] ? sizes[i] : 0;
    }
  }
  return element;
}

// runMol copies each atom's symbol and reads its parameters; a symbol
// starting with '*', RDKit's dummy atom, makes fill_atomic_parms read
// special parameters from its input file, which is null. check_for_errors
// exits the process unless there are some electrons and at most two per
// orbital, and YAeHMOP sizes its matrices as int(orbitals * orbitals).
void check_hueckel(const RWMol &mol, std::int32_t charge) {
  check_charges(mol, charge);
  check_matrix_size(mol);
  if (!mol.getNumAtoms()) {
    return;  // connectivityHueckel throws on its own.
  }
  auto [lines, ends] = eht_parameter_lines();
  std::int64_t valence = 0;
  std::int64_t orbitals = 0;
  // connectivityHueckel sets the first atom's formal charge to `charge`,
  // which Atom keeps in 8 bits.
  std::int64_t total_charge = static_cast<std::int8_t>(charge);
  bool first = true;
  for (const auto atom : mol.atoms()) {
    if (atom->getAtomicNum() == 0) {
      throw ValueErrorException(
          "YAeHMOP reads a dummy atom's parameters from a null file");
    }
    auto element = eht_element(lines, ends, atom->getSymbol());
    valence += element.valence;
    orbitals += element.orbitals;
    if (!first) {
      total_charge += atom->getFormalCharge();
    }
    first = false;
  }
  if (orbitals > 46340) {
    throw ValueErrorException(
        "YAeHMOP sizes its matrices in int, which overflows above 46340 "
        "orbitals");
  }
  auto electrons = valence - total_charge;
  if (electrons <= 0 || electrons > 2 * orbitals) {
    throw ValueErrorException(
        "YAeHMOP needs at least one electron and at most two per orbital; "
        "the charge leaves " +
        std::to_string(electrons) + " for " + std::to_string(orbitals) +
        " orbitals");
  }
}

#endif

// What the connectivity step trusts, for the method it uses.
void check_connectivity(const RWMol &mol, bool use_hueckel, std::int32_t charge,
                        bool use_vdw) {
  if (use_hueckel) {
#ifdef RDK_BUILD_YAEHMOP_SUPPORT
    check_hueckel(mol, charge);
#endif
    // Without YAeHMOP, RDKit's own ValueError.
  } else if (use_vdw) {
    // connectivityVdW reads the n * n 3D distance matrix.
    check_matrix_size(mol);
  } else {
    check_proximity_coordinates(mol);
    // determineConnectivity removes every bond first.
    check_proximity_hydrogens(mol, true);
  }
}

}  // namespace

// determineConnectivity first removes the molecule's bonds, which trusts
// what check_removals checks, checked before the first edit.
void determine_connectivity(RWMol &mol, bool use_hueckel, std::int32_t charge,
                            double cov_factor, bool use_vdw) {
  check_connectivity(mol, use_hueckel, charge, use_vdw);
  check_removals(mol);
  MatrixCacheGuard guard(mol);
  OtherMonomerInfoSetAside set_aside(mol);
  RDKit::determineConnectivity(mol, use_hueckel, charge, cov_factor, use_vdw);
}

// determineBondOrders constructs a ControlCHandler (control_c_mutex), and
// sanitizes the molecule when it embeds chirality (U6).
void determine_bond_orders(RWMol &mol, std::int32_t charge,
                           bool allow_charged_fragments, bool embed_chiral,
                           bool use_atom_map) {
  check_charges(mol, charge);
  MatrixCacheGuard guard(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  RDKit::determineBondOrders(mol, charge, allow_charged_fragments, embed_chiral,
                             use_atom_map);
}

// As determineBonds, which leaves a molecule of one atom or none as it is.
void determine_bonds(RWMol &mol, bool use_hueckel, std::int32_t charge,
                     double cov_factor, bool allow_charged_fragments,
                     bool embed_chiral, bool use_atom_map, bool use_vdw) {
  if (mol.getNumAtoms() <= 1) {
    return;
  }
  check_connectivity(mol, use_hueckel, charge, use_vdw);
  check_removals(mol);
  // The bond orders start from the charge the Hückel option set.
  check_charges(mol, charge);
  MatrixCacheGuard guard(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  OtherMonomerInfoSetAside set_aside(mol);
  RDKit::determineBonds(mol, use_hueckel, charge, cov_factor,
                        allow_charged_fragments, embed_chiral, use_atom_map,
                        use_vdw);
}

bool hueckel_enabled() {
#ifdef RDK_BUILD_YAEHMOP_SUPPORT
  return true;
#else
  return false;
#endif
}

}  // namespace rdkit_cxx
