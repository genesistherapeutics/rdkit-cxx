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

#include <GraphMol/ROMol.h>

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

// RDKit's PDB, sequence, and Maestro writers cast an atom's monomer info to
// AtomPDBResidueInfo when it is tagged PDBRESIDUE (MaeWriter.cpp: whenever it
// is present), but a pickle or an AtomMonomerInfo constructed with that type
// carries the tag on the base class, which the cast would read past.
// `every_info` checks every present info rather than only tagged ones.
void check_pdb_residue_infos(const RDKit::ROMol &mol, bool every_info);

}  // namespace rdkit_cxx
