/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
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

namespace RDKit {
class ROMol;
}  // namespace RDKit

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

// RDKit 2025.09's CXSMILES parser (parse_it in CXSmilesOps.cpp) reads its
// text, the extension from its opening "|", up to end() and in places past
// it:
// - an "atomProp:" block (parse_atom_props) reads the byte after the text's
//   terminating null when the text ends inside it or its last value runs to
//   the end (rdkit/rdkit#9361), whatever strictCXSMILES says;
// - a data substance group's ("SgD:") field without a ":" after it, its
//   coordinates without a ")", and a polymer group's ("Sg:") type without a
//   ":" step past the end, and the parse fails there; without
//   strictCXSMILES, MolFromSmiles and MolFromSmarts then copy the text up to
//   where it stopped into _CXSMILES_Data, past the terminating null;
// - an enhanced stereo group ("a:", "o<n>:", "&<n>:") may list an atom twice
//   (parse_enhanced_stereo merges sections by type and number), and the
//   parser's own removal of hydrogens then leaves the duplicate dangling
//   (cleanupStereoGroups reads the freed atom).
// Each such block is read here from wherever it appears in the text, as the
// parser would read it there, which covers wherever the parser does; a
// refusal throws ValueErrorException. (As rdkit-2024.09's and
// rdkit-2025.03's checks.)
void check_cx_extensions(std::string_view cx, bool strict);

// The extension MolFromSmiles and MolFromSmarts parse (preprocessSmiles in
// SmilesParse.cpp): after the text's first space or tab, trimmed, if it
// starts with "|" and allowCXSMILES is set.
void check_smiles_cx(const std::string &text, bool allow_cx, bool strict);

// The extension the reaction SMARTS and SMILES parsers parse (parseReaction
// in DaylightParser.cpp): from the text's first "|" after its start,
// trimmed; RDKit throws when it fails.
void check_reaction_cx(const std::string &text, bool allow_cx);

// The reactions MolStandardize's transformation readers parse from text
// (readTransformations in TransformCatalogUtils.cpp, which Normalizer's
// stream constructor uses): the second tab-separated field of each line, of
// at most 511 characters, without spaces. (Its fragment, acid-base, and
// tautomer readers erase the spaces from their SMARTS too, and a SMARTS
// without spaces or tabs has no extension.)
void check_transformation_lines(const std::string &text);

// addRecursiveQueries parses an atom's `prop` value as SMARTS when it names
// no query (ChemTransforms.cpp).
void check_query_props(const RDKit::ROMol &mol, const std::string &prop);

}  // namespace rdkit_cxx
