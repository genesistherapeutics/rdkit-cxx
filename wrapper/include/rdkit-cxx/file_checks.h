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

// RDKit 2024.03's CXSMILES extension parser (parse_it and the parsers it
// calls in SmilesParse/CXSmilesOps.cpp) moves an iterator `first` through
// the extension text up to `last`, its end(), and reads `*first` while
// `first <= last`, so it reads the string's terminator, and it steps one
// past it where a field runs to the end. Three fields go further:
//
// - a data S-group (`SgD:`) reads five fields, each with
//   parse_data_sgroup_attr, which reads `*first` while `first != last` and
//   steps past the ':' that ends the field whether or not it found one, so
//   a field that runs to the end leaves `first` at end() + 2 and the next
//   reads past the string (2026.03 adds a precondition);
// - an atom property (`atomProp:`) whose value runs to the end leaves
//   `first` at end() + 1, which parse_atom_props then dereferences (2026.03
//   returns first);
// - a polymer S-group (`Sg:`) whose type runs to the end steps `first` to
//   end() + 2, which MolFromSmiles then copies the text up to, as
//   "_CXSMILES_Data", when strictCXSMILES is off (handleCXPartAndName in
//   SmilesParse.cpp).
//
// cx_extension_problem follows each of those fields wherever it could start
// in the text, as RDKit would read it there, and says why it is refused, or
// returns an empty string. Other fields stop at end() + 1 at most, and
// strings keep their terminator, so they need no check.
std::string cx_extension_problem(std::string_view cx);

// The extension text MolFromSmiles and MolFromSmarts parse
// (preprocessSmiles in SmilesParse.cpp): with allowCXSMILES, the text after
// the first space or tab, trimmed, when it starts with "|". Throws
// ValueErrorException for one cx_extension_problem refuses.
void check_smiles_cx(const std::string &text, bool allow_cxsmiles = true);

// As check_smiles_cx, for the extension text the reaction SMARTS and SMILES
// parsers parse (parseReaction in ChemReactions/DaylightParser.cpp): from
// the first "|", when it is not the first character, trimmed.
void check_reaction_cx(const std::string &text, bool allow_cxsmiles = true);

// RDKit's Marvin reader (MarvinMolBase::parseAtomsAndBonds in
// MarvinDefs.cpp), which reads molecules and reactions alike, copies the
// second atom reference of each bond before it checks that there are two,
// reading past the end of its token vector for a bond with one or none (a
// segfault in Python). This parses the XML as RDKit does, which throws for
// malformed XML, and checks every bond.
void check_mrv(const std::string &text);

}  // namespace rdkit_cxx
