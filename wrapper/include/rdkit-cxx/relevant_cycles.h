/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// MolOps::symmetrizeSSSR's default algorithm, RDL, enumerates the relevant
// cycles (RingInfo::atomRelevantCycles), first reserving RDL_getNofRC's
// double, which it converts to size_t unchecked: undefined for a count of 2^64
// or more and for RDL_INVALID_RC_COUNT (DBL_MAX), which RingDecomposerLib
// returns for a count too large to compute. The count grows exponentially (a
// ring of n diamonds, 3n atoms, has 2^n + n). Every relevant cycle then
// becomes a ring, and the ring info stores ring indices as int, so the
// wrappers refuse more than INT_MAX relevant cycles (U1): a count that fits
// size_t but not int could not be enumerated in memory anyway.
//
// symmetrizeSSSR is reached from far more than the wrappers that call it:
// sanitization's SANITIZE_SYMMRINGS step, aromaticity and stereo perception
// on a molecule whose ring info is not a symmetrized SSSR, the parsers, the
// depictor, and much else. Every wrapper whose RDKit call can reach it calls
// one of the checks here first, on the graph that call symmetrizes, or on a
// graph whose relevant cycles include those of every graph it symmetrizes.
//
// The same check has been proposed to RDKit, for atomRelevantCycles and
// numRelevantCycles to throw ValueErrorException for such a count. Once a
// release the crate binds has it, this file and its calls can go: RDKit's own
// exception then maps to the same Error::Value.

#pragma once

#include <climits>
#include <cstdint>
#include <istream>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace rdkit_cxx {

// Relevant cycles are distinct nonzero elements of the cycle space, so a
// graph with c independent cycles (bonds - atoms + fragments) has at most
// 2^c - 1 of them, and the checks count them only for a graph with more than
// MAX_INDEPENDENT_CYCLES. Adding a bond adds at most one independent cycle,
// so a bound on them with every bond bounds them for any subset of the bonds,
// or for any graph built from these bonds and others, plus one for each of
// the others. Relevant cycles do not shrink with the graph, though: removing
// a bond can make longer cycles relevant. So a check above the bound counts
// the cycles of the graph RDKit will symmetrize, not of a larger one.
constexpr unsigned MAX_INDEPENDENT_CYCLES = 31;

// Whether the graph RDKit gives RingDecomposerLib (every bond but zero-order
// ones, and dative and hydrogen bonds only if included), and every graph of
// some of its bonds, has at most INT_MAX relevant cycles, without counting
// them: it has at most MAX_INDEPENDENT_CYCLES independent cycles, or its
// biconnected blocks, over which relevant cycles add up, have few enough
// (a block of c independent cycles has at most 2^c - 1 relevant ones; the
// blocks of a graph of some of its bonds lie within its blocks, with no more
// independent cycles among them).
bool few_relevant_cycles(const RDKit::ROMol &mol, bool include_dative_bonds,
                         bool include_hydrogen_bonds, bool include_zero_bonds);

// The independent cycles of the molecule's graph with every bond.
std::uint64_t independent_cycles(const RDKit::ROMol &mol);

// Counts the independent cycles of a graph given as its bonds, one at a time,
// with a union-find over its atoms' labels.
class CycleCounter {
 public:
  void add_bond(std::uint64_t begin, std::uint64_t end);
  std::uint64_t cycles() const { return d_cycles; }

 private:
  std::uint64_t root(std::uint64_t node);
  std::unordered_map<std::uint64_t, std::uint64_t> d_parent;
  std::uint64_t d_cycles = 0;
};

// Whether each product of the reaction is its reactants less some bonds and
// atoms, with trees added that each join them by one bond, so it has no more
// independent cycles: the product templates bond mapped atoms only as the
// reactant templates do, and once the mapped atoms are taken as one, their
// other bonds form no cycle. `exactly` also asks that the products keep
// every atom and bond of the reactant templates, as bonds RDKit's ring
// graph takes, so they keep the reactants' biconnected blocks.
bool reaction_keeps_cycles(const RDKit::ChemicalReaction &rxn, bool exactly);

// What a set of reactions applied in turn, as the standardizer's
// normalizations are, does to the cycles of what they apply to.
struct TransformCycles {
  bool keeps = true;          // every one reaction_keeps_cycles
  bool keeps_exactly = true;  // every one reaction_keeps_cycles exactly
};
TransformCycles transform_cycles(
    const std::vector<std::shared_ptr<RDKit::ChemicalReaction>> &rxns);

// Throws ValueErrorException if the standardizer (MolStandardize), on this
// molecule, could symmetrize more than INT_MAX relevant cycles. Its steps
// sanitize and perceive what each one leaves: the molecule less hydrogens
// of one bond, its fragments, less bonds to metals, with the normalizations
// applied, and its tautomers. From a molecule of at most
// MAX_INDEPENDENT_CYCLES independent cycles (with every bond), under
// normalizations that keep cycles, none has more. Past that, a molecule with
// no zero-order, dative or hydrogen bond, and no bond to an atom the metal
// disconnection or organometallic cleanup could cut or make dative, keeps
// its biconnected blocks through every step, under normalizations that keep
// them exactly, so its own relevant cycles bound every step's. Anything
// else is refused: normalizations that bond atoms anew, which apply again
// and again, and, past the bound, the rest (`custom_metals`: metal patterns
// the caller chose, which can cut any bond).
void check_standardize_relevant_cycles(const RDKit::ROMol &mol,
                                       TransformCycles normalizations,
                                       bool custom_metals = false);

// Throws ValueErrorException if symmetrizeSSSR, recalculating, would
// enumerate more than INT_MAX relevant cycles of the molecule's graph.
void check_graph_relevant_cycles(const RDKit::ROMol &mol,
                                 bool include_dative_bonds = false,
                                 bool include_hydrogen_bonds = false);

// Throws ValueErrorException if symmetrizeSSSR, with these options, would
// enumerate more than INT_MAX relevant cycles: those of the molecule's ring
// families when it keeps them (they exist and recalc_sssr is false, or the
// ring info is not initialized), otherwise those of its graph. RDKit's own
// errors (RingDecomposerLib failing) are left to RDKit.
void check_relevant_cycles(const RDKit::ROMol &mol,
                           bool include_dative_bonds = false,
                           bool include_hydrogen_bonds = false,
                           bool recalc_sssr = true);

// As check_relevant_cycles, for the RDKit calls that symmetrize a molecule's
// rings only when they are not already a symmetrized SSSR (`if
// (!ringInfo->isSymmSssr()) symmetrizeSSSR(mol)`: stereo perception,
// detectBondStereochemistry, drawing, the JSON writer, and others), or, with
// unless_initialized, only when they are not found at all (setAromaticity).
// O(1) for a sanitized molecule.
void check_perceive_relevant_cycles(const RDKit::ROMol &mol,
                                    bool unless_initialized = false);

// As check_relevant_cycles, for MolOps::sanitizeMol with these operations,
// which finds the rings of the graph its cleanup steps leave (they can make
// bonds dative), in SANITIZE_SYMMRINGS or in the later steps that perceive
// rings.
void check_sanitize_relevant_cycles(const RDKit::ROMol &mol,
                                    unsigned sanitize_ops);

// As check_sanitize_relevant_cycles, for MolOps::removeHs with these
// parameters, which sanitizes what it leaves when `sanitize` is set
// (removeHigherDegrees removes hydrogens from rings).
void check_remove_hs_relevant_cycles(
    const RDKit::ROMol &mol, const RDKit::MolOps::RemoveHsParameters &ps,
    bool sanitize);

// As check_graph_relevant_cycles, for the RDKit calls that perceive
// stereochemistry (assignStereochemistry, findPotentialStereo, and the
// writers that call them) on the molecule, a copy of it, or each of its
// fragments: they symmetrize the rings when the ring info is not a
// symmetrized SSSR, and a fragment cut from a molecule of several has none.
// The fragments' relevant cycles are among the molecule's.
void check_stereo_relevant_cycles(const RDKit::ROMol &mol);

// As check_stereo_relevant_cycles, for RDKit's SMARTS writers (MolToSmarts,
// and GetAtomSmarts, which the mol-block writers call for query atoms),
// which write each recursive query's molecule with MolToSmarts, perceiving
// its stereochemistry as they do the molecule's. With `include_mol` false,
// only the recursive queries' molecules are checked (for writers that write
// the molecule otherwise).
void check_smarts_relevant_cycles(const RDKit::ROMol &mol,
                                  bool include_mol = true);
// For GetAtomSmarts on one atom: its recursive queries' molecules.
void check_atom_smarts_relevant_cycles(const RDKit::Atom &atom);

// As check_sanitize_relevant_cycles, for a reader that builds a molecule
// without perceiving anything and then, as asked, removes its hydrogens
// (MolOps::removeHs with its default parameters, sanitizing) or sanitizes
// it; `raw` is what it builds, read without either.
void check_parsed_relevant_cycles(const RDKit::ROMol &raw, bool sanitize,
                                  bool remove_hs);

// As check_sanitize_relevant_cycles, for MolOps::getMolFrags with
// sanitizeFrags, which sanitizes each fragment it cuts from the molecule.
// Its callers check what the cutting trusts (check_removals).
void check_frags_sanitize_relevant_cycles(const RDKit::ROMol &mol);

// For ChemicalReaction::runReactants (or, with reactant_idx, runReactant)
// with these reactants. A run sets its products' bond directions when a
// reactant has a directed single bond (generateOneProductSet), with
// setDoubleBondNeighborDirections, which symmetrizes each product's rings;
// with `sanitize` the caller (or RDKit, as removeNonmatchingReagents does)
// sanitizes the products. A product is its template's atoms and bonds, with
// what is left of each reactant glued to them at the atoms its template
// matched, so its independent cycles are at most the product templates'
// bonds, plus each reactant's independent cycles and the atoms of its
// template. Above that bound the reaction is run first on copies of the
// reactants without bond directions, which perceives nothing and makes the
// same products, and they are checked. The reaction is initialized; the
// caller holds control_c_mutex().
void check_products_relevant_cycles(const RDKit::ChemicalReaction &rxn,
                                    const RDKit::MOL_SPTR_VECT &reactants,
                                    unsigned max_products, bool sanitize,
                                    unsigned reactant_idx = UINT_MAX);

// For MolOps::removeAllHs, which removes every hydrogen (rings' too) and
// sanitizes what is left when `sanitize` is set.
void check_remove_all_hs_relevant_cycles(const RDKit::ROMol &mol,
                                         bool sanitize);

// Half the characters that can write a ring closure (digits and '%') in
// `text`: a bound on the independent cycles of the molecule of any SMILES in
// it, whose bonds but its ring closures form a spanning forest.
std::size_t smiles_cycle_bound(std::string_view text);
// The same for one SMILES, not SMARTS (whose recursive queries close rings
// inside brackets), counting only labels outside brackets.
std::size_t smiles_closure_bound(std::string_view smiles);

// As check_sanitize_relevant_cycles, for SmilesToMol, which sanitizes what it
// parses, or removes hydrogens and perceives stereochemistry. Text whose
// SMILES part bounds the cycles is parsed once, by RDKit; other text is
// parsed first without sanitizing, and that molecule checked.
void check_smiles_relevant_cycles(const std::string &text,
                                  const RDKit::SmilesParserParams &params);
void check_smiles_relevant_cycles(
    const std::string &text,
    const RDKit::v2::SmilesParse::SmilesParserParams &params);

// Mol blocks. RDKit's reader (MolFromMolDataStream) reads a block's CTAB
// without perceiving anything, then (finishMolProcessing) symmetrizes the
// rings of the graph it read: detectBondStereochemistry does, for every
// block with a conformer (all of them), after only the cleanup steps
// (SANITIZE_CLEANUP, which changes no bond's ring family) when it sanitizes
// and removes hydrogens, and without sanitizing; sanitization symmetrizes
// the graph its organometallic cleanup leaves, after hydrogens are removed
// if they are. The independent cycles of the bonds a block lists bound them
// all.

// Throws ValueErrorException if RDKit, reading `block` (the text
// MolFromMolDataStream reads: header, counts line, CTAB through "M  END")
// with these parameters, would symmetrize more than INT_MAX relevant cycles.
// It reads the CTAB as RDKit does, with RDKit's own CTAB reader, and checks
// each graph RDKit would symmetrize; a block RDKit's reader fails on before
// it perceives anything passes, since RDKit fails the same way.
void check_mol_block_relevant_cycles(
    std::string_view block,
    const RDKit::v2::FileParsers::MolFileParserParams &params);

// The independent cycles of the bonds of the Mol2 block RDKit reads from
// `text` (MolFromMol2DataStream: the first molecule's counts line and bond
// lines), or UINT64_MAX when its bond lines cannot be told apart as RDKit
// finds them; 0 for a block RDKit fails on before it builds the molecule.
std::uint64_t mol2_cycle_bound(const std::string &text);

// Throws ValueErrorException if RDKit's SCSR reader (MolFromSCSRDataStream),
// reading `text` with these parameters, would symmetrize more than INT_MAX
// relevant cycles: it finishes the main molecule as the mol-block reader
// does, finishes each template without sanitizing, and finishes the molecule
// it expands them into with the parameters. The expansion is counted on a
// copy RDKit expands with every bond made dative (or, for its hydrogen
// bonds, left as is), so it perceives no rings; the copy's dative bonds are
// the expansion's. A block with dative bonds of its own, or whose expansion
// RDKit's organometallic cleanup could change, is refused when the copy has
// more than MAX_INDEPENDENT_CYCLES independent cycles: its rings cannot be
// counted without RDKit perceiving them.
void check_scsr_relevant_cycles(
    const std::string &text,
    const RDKit::v2::FileParsers::MolFileParserParams &params,
    const RDKit::v2::FileParsers::MolFromSCSRParams &scsr);

// Throws ValueErrorException if InchiToMol, reading `inchi`, would
// symmetrize more than INT_MAX relevant cycles. It perceives stereochemistry
// even unsanitized, so there is no reading the structure without perceiving
// it; but its graph, but for hydrogens, is the one the connection layers
// ("/c", for the main structure and after "/r" for the reconnected one)
// write. Each component's layer names its atoms by number along chains, so
// it has as many independent cycles as numbers past the first, less distinct
// atoms past the first; above MAX_INDEPENDENT_CYCLES, the graph is built
// from the layer and counted. Sanitization could make a bond to a metal
// dative, so a structure above the bound with one is refused.
void check_inchi_relevant_cycles(std::string_view inchi, bool sanitize);

// Molecule pickles. MolPickler::molFromPickle assigns stereochemistry to a
// pickle older than format 4.0 as it reads it (assignStereochemistry), which
// symmetrizes its rings before any wrapper sees its bonds, so such pickles
// are refused. Every molecule pickle starts with endianId, the VERSION tag,
// and its major, minor, and patch versions, and RDKit reads each one it finds
// in a reaction, bundle, tautomer or extended query, catalog, library, or
// synthon space, and each recursive query's, from a run of the bytes it is
// given (once decoded, in extended-query JSON's base64 and PNG metadata's
// compressed chunks). So these throw ValueErrorException for any run of the
// bytes that reads as an older pickle's header; bytes that only look like
// one, in a property's value say, are refused too.
void check_mol_pickle_versions(std::string_view bytes);
// The same for the bytes of a stream, read to its end.
void check_mol_pickle_versions(std::istream &in);

// Throws ValueErrorException if RDKit's CDXML reader (the ChemDraw
// library's, under CDXMLToChemicalReactions), reading `text` with these
// parameters, would symmetrize more than INT_MAX relevant cycles of one of
// its fragments: it gives every fragment a conformer, so
// detectBondStereochemistry symmetrizes each, sanitized or not, besides
// sanitization, and it fuses abbreviations with molzip. The fragments are
// counted on a copy whose bonds are all hydrogen bonds (or, for its dative
// and hydrogen bonds, dative ones), which RDKit reads without perceiving
// any rings; the copy's hydrogen bonds are the fragments' graphs. Binary
// CDX input is refused, as is a document whose fragments have more than
// MAX_INDEPENDENT_CYCLES independent cycles and bonds to metals that
// sanitization could make dative.
void check_cdxml_relevant_cycles(const std::string &text, bool sanitize,
                                 bool remove_hs);

// Throws ValueErrorException if RDKit's Marvin reader, reading the molecule
// of `text` with these parameters, would symmetrize more than INT_MAX
// relevant cycles. It symmetrizes them whenever the molecule has
// coordinates (detectBondStereochemistry, as the mol-block reader does),
// sanitized or not, and in sanitization. The bonds of the document's
// bondArrays bound them, unless RDKit copies some (a MultipleSgroup); a
// document over the bound is read again with its coordinates removed, which
// RDKit reads without perceiving anything, and that molecule checked.
void check_mrv_relevant_cycles(const std::string &text, bool sanitize,
                               bool remove_hs);
// For a Marvin reaction (ReactionFromMrvDataStream), whose templates RDKit
// reads as it reads molecules.
void check_mrv_reaction_relevant_cycles(const std::string &text, bool sanitize,
                                        bool remove_hs);

// Throws ValueErrorException if RDKit's SVG reader (RDKitSVGToMol), reading
// the molecules in the metadata of `svg` with these parameters, would
// symmetrize more than INT_MAX relevant cycles of one: it removes their
// hydrogens or sanitizes each, from its rdkit:bond elements.
void check_svg_relevant_cycles(const std::string &svg, bool sanitize,
                               bool remove_hs);

// As check_mol_block_relevant_cycles, for each template of the rxn block
// `text` (ReactionFromRxnBlock), which RDKit reads as it reads mol blocks:
// a V2000 block's from "$MOL" lines, a V3000 block's as CTABs, each finished
// (agents of a V2000 block without sanitizing).
void check_rxn_block_relevant_cycles(
    const std::string &text,
    const RDKit::v2::FileParsers::MolFileParserParams &params);

// Follows text line by line as RDKit's mol-block readers read it, to bound
// the independent cycles of each block's bonds: the header and counts line,
// the atom lines, and the bond lines, whose atoms it counts cycles over
// (V2000's columns and V3000's tokens, read as RDKit reads them); then the
// properties up to the line that ends the CTAB: V2000's "M  END", which the
// lines that an "A", "G", or "S  SKP" line takes do not end, as in RDKit's
// reader; V3000's first line that does not start "M  V30 ", which RDKit
// accepts only as "M  END".
//
// In SD text, RDKit starts a record only after a "$$$$" line (or at the
// start), but not after every one: a record's header lines may be "$$$$",
// and a record RDKit fails on is skipped up to the next. So a record is
// followed from after each "$$$$" line, alongside any still being read.
//
// A record's text is kept until its bonds bound its cycles. A record over the
// bound is kept, and its lines after the bonds held, through the line that
// ends its CTAB, when it is checked (check_mol_block_relevant_cycles): RDKit
// symmetrizes nothing before it has read that line.
class MolRecordScanner {
 public:
  // With `sd`, records follow "$$$$" lines; otherwise the text is one block.
  MolRecordScanner(const RDKit::v2::FileParsers::MolFileParserParams &params,
                   bool sd);

  // A physical line, without its "\n", which starts at offset `at`. Returns
  // false once a record is refused (error() says why).
  bool line(std::string_view text, std::uint64_t at);
  // Where held lines start: those of a record over the bound, after its
  // bonds, until it is checked; UINT64_MAX when nothing is held. A refusal
  // takes effect there.
  std::uint64_t held_from() const;
  const std::string &error() const { return d_error; }
  std::uint64_t refuse_at() const { return d_refuse_at; }

 private:
  enum class State : std::uint8_t {
    kHeader,
    kCounts,
    kAtoms,
    kBonds,
    kProperties,
    kV3000,
  };
  // A record being followed: where RDKit's reader would stand in it, the
  // cycles of its bonds, and its text while kept.
  struct Record {
    State state = State::kHeader;
    std::uint64_t lines = 0;  // lines read in the current state
    std::uint64_t atoms = 0;
    std::uint64_t bonds = 0;
    CycleCounter cycles;
    // V2000 properties: lines an "A", "G", or "S  SKP" line takes, and
    // whether the first property line (which may be blank) is still to come.
    std::uint64_t skip = 0;
    bool first_property = true;
    // V3000: logical lines read, and the one being joined across "-"s.
    std::uint64_t logical = 0;
    std::string joined;
    bool bonds_done = false;
    std::string text;
    bool keep = true;
    std::uint64_t held_from = UINT64_MAX;
  };
  // Feeds a line to a record; false once the record is done (its CTAB read,
  // or failed as RDKit's reader fails).
  bool feed(Record &record, std::string_view raw, std::uint64_t at);
  void end_bonds(Record &record, std::uint64_t next_line);
  bool end_ctab(Record &record);

  RDKit::v2::FileParsers::MolFileParserParams d_params;
  bool d_sd;
  std::vector<Record> d_records;
  std::string d_error;
  std::uint64_t d_refuse_at = UINT64_MAX;
};

// Follows the lines of a SMILES file as RDKit's SMILES suppliers read them
// (SmilesMolSupplier, MultithreadedSmilesMolSupplier): blank lines and lines
// that start with "#" are skipped, the first other line is a title line when
// the parameters say so, and every other line is a record, whose SMILES is a
// field of it. A line with at most 2 * MAX_INDEPENDENT_CYCLES characters that
// can write ring closures (smiles_cycle_bound) is not checked. A longer one
// is held from the character past that many, until it ends and its SMILES is
// checked (check_smiles_relevant_cycles): RDKit parses a line only once it
// has read it whole, and the line before that character, should a refusal
// cut the input there, is one RDKit can parse.
class SmilesLineScanner {
 public:
  explicit SmilesLineScanner(
      const RDKit::v2::FileParsers::SmilesMolSupplierParams &params)
      : d_params(params), d_title(params.titleLine) {}
  // For a TDT file read with these parameters (TDTMolSupplier), whose
  // records' SMILES are the lines that start with "$SMI<", between the
  // first "<" and the last ">".
  explicit SmilesLineScanner(
      const RDKit::v2::FileParsers::TDTMolSupplierParams &params)
      : d_tdt(true), d_tdt_params(params.parseParameters), d_title(false) {}

  // A character of the line being read, at offset `at`.
  void byte(char c, std::uint64_t at) {
    if (d_tdt && d_column < 5) {
      d_tdt_line = d_tdt_line && c == "$SMI<"[d_column];
      ++d_column;
      return;
    }
    if ((!d_tdt || d_tdt_line) && (c == '%' || (c >= '0' && c <= '9')) &&
        ++d_closures == 2 * MAX_INDEPENDENT_CYCLES + 1) {
      d_held_from = at;
    }
  }
  // The end of a line, whose text is `text` (without its "\n"). Returns
  // false once a line is refused (error() says why).
  bool line(std::string_view text);
  std::uint64_t held_from() const { return d_held_from; }
  const std::string &error() const { return d_error; }
  std::uint64_t refuse_at() const { return d_refuse_at; }

 private:
  RDKit::v2::FileParsers::SmilesMolSupplierParams d_params;
  bool d_tdt = false;
  RDKit::v2::SmilesParse::SmilesParserParams d_tdt_params;
  unsigned d_column = 0;
  bool d_tdt_line = true;
  bool d_title;
  std::uint64_t d_closures = 0;
  std::uint64_t d_held_from = UINT64_MAX;
  std::string d_error;
  std::uint64_t d_refuse_at = UINT64_MAX;
};

// As check_smiles_relevant_cycles, for a record line of a SMILES file read
// with these parameters (SmilesLineScanner).
void check_smiles_line_relevant_cycles(
    std::string_view line,
    const RDKit::v2::FileParsers::SmilesMolSupplierParams &params);

// Throws ValueErrorException if RDKit, reading the mol blocks of `text` (one
// block, or SD records) with these parameters, would symmetrize more than
// INT_MAX relevant cycles (MolRecordScanner).
void check_mol_text_relevant_cycles(
    std::string_view text,
    const RDKit::v2::FileParsers::MolFileParserParams &params);

}  // namespace rdkit_cxx
