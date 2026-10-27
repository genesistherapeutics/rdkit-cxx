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

// RDKit's Marvin reader (MarvinMolBase::parseAtomsAndBonds in
// MarvinDefs.cpp), which reads molecules and reactions alike, copies the
// second atom reference of each bond before it checks that there are two,
// reading past the end of its token vector for a bond with one or none (a
// segfault in Python). This parses the XML as RDKit does, which throws for
// malformed XML, and checks every bond.
void check_mrv(const std::string &text);

// RDKit 2026.03's CXSMILES parser (parse_it in CXSmilesOps.cpp) reads a
// data substance group's ("SgD:") field without a ":" after it, its
// coordinates without a ")", and a polymer group's ("Sg:") type without a
// ":" (read_text_to) through the text's terminating null, and steps past
// it; the parse fails there. Without strictCXSMILES, MolFromSmiles and
// MolFromSmarts then copy the text up to where it stopped into
// _CXSMILES_Data (handleCXPartAndName in SmilesParse.cpp), up to 2 bytes
// past its end: for a CX part too long for std::string's inline buffer, 1
// byte past its allocation. Each such block is read here from wherever it
// appears in the text, as the parser would read it there, which covers
// wherever the parser does; a refusal throws ValueErrorException. With
// strictCXSMILES RDKit rethrows before it copies, so nothing is checked.
void check_cx_extensions(std::string_view cx, bool strict);

// The extension MolFromSmiles and MolFromSmarts parse (preprocessSmiles in
// SmilesParse.cpp): after the text's first space or tab, trimmed, if it
// starts with "|" and allowCXSMILES is set.
void check_smiles_cx(const std::string &text, bool allow_cx, bool strict);

}  // namespace rdkit_cxx
