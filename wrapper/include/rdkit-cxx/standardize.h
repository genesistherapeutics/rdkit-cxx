/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's CleanupParameters and the module functions that take
// them (MolStandardize.h): cleanup, normalization, reionization, fragment
// removal, the canonical tautomer, and the parents.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Tautomer.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared types, defined by the bridge (src/ffi/standardize.rs).
struct CleanupDataEntry;
struct CleanupParametersData;
enum class StandardizeOp : std::uint8_t;

using CleanupParameters = RDKit::MolStandardize::CleanupParameters;

// Conversions between the bridge's struct and RDKit's. File names holding a
// NUL byte are refused (util.h's to_path).
CleanupParameters to_cleanup_parameters(const CleanupParametersData &data);
CleanupParametersData from_cleanup_parameters(const CleanupParameters &params);
// Reads the normalizations file into the data RDKit would read from it, for
// a call that normalizes (U1: RDKit would otherwise open it unchecked).
void read_normalizations(CleanupParameters &params);

// Takes ownership of a molecule an RDKit object returned, which may be an
// ROMol or an RWMol, as an RWMol without cached matrices (CONVENTIONS.md U6).
std::unique_ptr<RWMol> adopt_mol(RDKit::ROMol *mol);

// RDKit's reionizer sums a molecule's formal charges in int before and after
// it applies its charge corrections, which set charges of up to 128 in
// magnitude, and subtracts the sums (Reionizer::reionizeInPlace in
// Charge.cpp), so the difference is up to 255 per atom; the uncharger sums
// the charges again. Throws ValueErrorException for a molecule of more than
// INT_MAX / 255 atoms, which could overflow them (CONVENTIONS.md U9).
void check_reionizer_atoms(const RDKit::ROMol &mol);

// RDKit's default tautomer score adds, in int, 100 or 250 for each aromatic
// ring of the symmetrized SSSR (scoreRings), whose relevant cycles can be
// exponentially many; each substructure term's score times its match count
// (at most 1000 per pattern, or one per bond or atom for the patterns
// scoreSubstructsFiltered counts directly); and minus the hydrogens on P, S,
// Se, and Te atoms (at most 382 per atom) (Tautomer.cpp). Throws
// ValueErrorException unless those bounds sum to at most INT_MAX for `mol`,
// which has `rings` rings.
void check_tautomer_score(const RDKit::ROMol &mol, std::uint64_t rings);
// As above, counting the rings of a copy of `mol` that has the symmetrized
// SSSR, so `mol` itself is not changed. Tautomers share their molecule's
// graph, and so its rings.
void check_tautomer_score(const RDKit::ROMol &mol);

// A tautomer enumerator for the parameters, as TautomerEnumerator(params)
// makes one, whose transforms have been checked. RDKit indexes a transform's
// bond types and charges by the bonds and atoms of its SMARTS without
// checking their lengths (TautomerEnumerator::enumerate in Tautomer.cpp), so
// shorter lists, from tautomerTransformData or the transforms file, are read
// past their ends: they are refused with ValueErrorException. The transforms
// are checked as RDKit parsed them, from one read of the file.
std::unique_ptr<RDKit::MolStandardize::TautomerEnumerator> checked_enumerator(
    const CleanupParameters &params);

CleanupParametersData cleanup_parameters_default();
CleanupParametersData update_params_from_json(const CleanupParametersData &data,
                                              rust::Str json);

// The module functions, by op. skip_standardize applies to the parents.
std::unique_ptr<RWMol> standardize_copy(StandardizeOp op, const RWMol &mol,
                                        const CleanupParametersData &data,
                                        bool skip_standardize);
void standardize_in_place(StandardizeOp op, RWMol &mol,
                          const CleanupParametersData &data,
                          bool skip_standardize);
void standardize_in_place_mols(StandardizeOp op, BorrowedMols &mols,
                               std::int32_t num_threads,
                               const CleanupParametersData &data,
                               bool skip_standardize);

std::unique_ptr<std::string> standardize_smiles(rust::Str smiles);

}  // namespace rdkit_cxx
