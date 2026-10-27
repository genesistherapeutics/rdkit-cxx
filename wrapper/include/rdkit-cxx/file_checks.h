/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks shared by the file-format wrappers (rdmolfiles). Each guards an
// input RDKit's readers or writers trust; the definitions say what RDKit
// would do without it (CONVENTIONS.md U1).

#pragma once

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rdkit_cxx {

// The bytes of a file, read once, so RDKit never rereads a file that changed
// after a check or after it indexed the records (U2). Throws
// RDKit::BadFileException("Bad input file <name>") when it cannot be opened,
// as RDKit's readers do.
std::string read_file(const std::string &name);

// As read_file, and also refuses an empty file with
// RDKit::BadFileException("Invalid input file <name>"), as RDKit's
// openAndCheckStream (MolSupplier.h) does.
std::string read_supplier_file(const std::string &name);

// The most digits the writers format a coordinate with. libstdc++ formats a
// double in a stack buffer as long as its output (alloca in
// locale_facets.tcc), so a larger precision overflows the stack.
constexpr unsigned kMaxPrecision = 1000;
void check_precision(unsigned precision);

// RDKit keeps the SMILES and SMARTS parsers' debug level in globals
// (yysmiles_debug, yysmarts_debug) that every parse reads without a lock;
// MolFromSmiles and MolFromSmarts write a nonzero debugParse there, racing
// with parses on other threads. The parsers therefore take only 0.
void check_debug_parse(int debug_parse);

// RDKit's Marvin reader (MarvinMolBase::parseAtomsAndBonds in
// MarvinDefs.cpp), which reads molecules and reactions alike, copies the
// second atom reference of each bond before it checks that there are two,
// reading past the end of its token vector for a bond with one or none (a
// segfault in Python). This parses the XML as RDKit does, which throws for
// malformed XML, and checks every bond.
void check_mrv(const std::string &text);

// RDKit 2025.03's CXSMILES parser (CXSmilesOps.cpp) reads past the end of
// an extension's text for an "atomProp:" list without a "|" or "," after
// it, and, when strictCXSMILES is off, copies past it for a data or polymer
// substance group whose fields run to the end; it also keeps a stereo group
// that lists an atom twice (check_cx_stereo_groups). check_cx_extensions
// throws ValueErrorException for the extension `cx` (from its "|") if RDKit
// would. RDKit's SMILES and TDT suppliers parse no CXSMILES.
void check_cx_extensions(std::string_view cx, bool strict);

// The extension MolFromSmiles and MolFromSmarts parse (preprocessSmiles in
// SmilesParse.cpp): after the text's first space or tab, trimmed, if it
// starts with "|" and allowCXSMILES is set.
void check_smiles_cx(std::string_view text, bool allow_cx = true,
                     bool strict = true);

// The extension the reaction SMARTS and SMILES parsers parse (parseReaction
// in DaylightParser.cpp): from the text's first "|" after its start,
// trimmed, parsed strictly.
void check_reaction_cx(std::string_view text, bool allow_cx = true);

// The reactions MolStandardize's transformation readers parse from text
// (readTransformations in TransformCatalogUtils.cpp, which the Normalizer's
// stream constructor uses): the second tab-separated field of each line, of
// at most 511 characters, without spaces. (Its fragment, acid-base, and
// tautomer readers erase the spaces from their SMARTS too, and a SMARTS
// without spaces or tabs has no extension.)
void check_transformation_lines(const std::string &text);

// The (name, reaction SMARTS) pairs readTransformations would read from
// `text`, for RDKit's data constructors. A line with a name and no reaction,
// whose second field RDKit reads from a finished tokenizer, is refused
// (ValueErrorException).
std::vector<std::pair<std::string, std::string>> transformation_data(
    const std::string &text);

}  // namespace rdkit_cxx
