// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/errors.h"

#include <stdexcept>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

namespace rdkit_cxx {

namespace {

// Python's rdBase initializes RDKit's logs on import; this does the same when
// the library loads. Every bridge links this file (E5), so it always runs.
const bool logs_initialized = [] {
  RDLog::InitLogs();
  return true;
}();

std::string encode(const std::string &variant, const std::string &payload,
                   const std::string &what) {
  return "\x1Frdkit\x1F" + variant + "\x1F" + payload + "\x1F" + what;
}

}  // namespace

IndexError::IndexError(unsigned index, unsigned len)
    : d_index(index),
      d_len(len),
      d_msg("index " + std::to_string(index) + " is out of range for length " +
            std::to_string(len)) {}

void check_index(unsigned index, unsigned len) {
  if (index >= len) {
    throw IndexError(index, len);
  }
}

// The most derived types come first (E4).
std::string encode_current_exception() {
  try {
    throw;
  } catch (const RDKit::AtomValenceException &e) {
    return encode("AtomValence", std::to_string(e.getAtomIdx()), e.what());
  } catch (const RDKit::AtomKekulizeException &e) {
    return encode("AtomKekulize", std::to_string(e.getAtomIdx()), e.what());
  } catch (const RDKit::KekulizeException &e) {
    std::string indices;
    for (auto idx : e.getAtomIndices()) {
      if (!indices.empty()) {
        indices += " ";
      }
      indices += std::to_string(idx);
    }
    return encode("Kekulize", indices, e.what());
  } catch (const RDKit::MolSanitizeException &e) {
    return encode("MolSanitize", "", e.what());
  } catch (const IndexError &e) {
    return encode("Index",
                  std::to_string(e.index()) + " " + std::to_string(e.len()),
                  e.what());
  } catch (const KeyErrorException &e) {
    return encode("Key", e.key(), e.what());
  } catch (const ValueErrorException &e) {
    return encode("Value", "", e.what());
  } catch (const RDKit::ConformerException &e) {
    return encode("Value", "", e.what());
  } catch (const RDKit::ChemicalReactionException &e) {
    // Python's ChemicalReactionException subclasses ValueError.
    return encode("Value", "", e.what());
  } catch (const RDKit::BadFileException &e) {
    return encode("Os", "", e.what());
  } catch (const RDKit::SmilesParse::SmilesParseException &e) {
    return encode("Parse", "", e.what());
  } catch (const RDKit::FileParseException &e) {
    return encode("Parse", "", e.what());
  } catch (const RDKit::SLNParseException &e) {
    return encode("Parse", "", e.what());
  } catch (const RDKit::ChemicalReactionParserException &e) {
    return encode("Parse", "", e.what());
  } catch (const Invar::Invariant &e) {
    // As nanobind's translator in Code/RDBoost/nbWrap/RDBase.cpp.
    return encode("Runtime", "", e.toUserString());
  } catch (const std::exception &e) {
    return encode("Runtime", "", e.what());
  } catch (...) {
    return encode("Runtime", "", "unknown C++ exception");
  }
}

void raise_for_tests(rust::Str kind_str) {
  std::string kind(kind_str);
  std::string message = "test " + kind;
  if (kind == "atom_valence") {
    throw RDKit::AtomValenceException(message, 3);
  } else if (kind == "atom_kekulize") {
    throw RDKit::AtomKekulizeException(message, 2);
  } else if (kind == "kekulize") {
    throw RDKit::KekulizeException(message,
                                   std::vector<unsigned>{0, 1, 2, 3, 4});
  } else if (kind == "mol_sanitize") {
    throw RDKit::MolSanitizeException(message);
  } else if (kind == "index") {
    throw IndexError(9, 3);
  } else if (kind == "key") {
    throw KeyErrorException(std::string("nope"));
  } else if (kind == "value") {
    throw ValueErrorException(message);
  } else if (kind == "conformer") {
    throw RDKit::ConformerException(message);
  } else if (kind == "bad_file") {
    throw RDKit::BadFileException(message);
  } else if (kind == "smiles_parse") {
    throw RDKit::SmilesParse::SmilesParseException(message);
  } else if (kind == "sln_parse") {
    throw RDKit::SLNParseException(message);
  } else if (kind == "reaction") {
    throw RDKit::ChemicalReactionException(message);
  } else if (kind == "reaction_parse") {
    throw RDKit::ChemicalReactionParserException(message);
  } else if (kind == "invariant") {
    throw Invar::Invariant("Invariant Violation", message, "false", __FILE__,
                           __LINE__);
  } else if (kind == "std") {
    throw std::runtime_error(message);
  }
  throw std::logic_error("raise_for_tests: unknown kind " + kind);
}

}  // namespace rdkit_cxx
