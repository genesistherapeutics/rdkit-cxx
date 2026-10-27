// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! CXSMILES extensions RDKit 2025.09 reads past or leaves dangling
//! (wrapper/include/rdkit-cxx/file_checks.h), at every entry point that hands
//! RDKit SMILES or SMARTS from the caller, and the Marvin atom ids that do the
//! same for stereo groups.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::{
    CleanupParameters, Normalizer, Pipeline, PipelineOptions, SubstructTerm, normalizer_from_data,
    normalizer_from_params, score_substructs_with_terms,
};
use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
use rdkit_cxx::chem::rd_deprotect::DeprotectData;
use rdkit_cxx::chem::rd_substruct_library::CachedSmilesMolHolder;
use rdkit_cxx::chem::rdfiltercatalog::SmartsMatcher;
use rdkit_cxx::chem::rdmolops::{mol_add_recursive_queries, parse_mol_query_def_file_text};
use rdkit_cxx::chem::scaffolds::rd_scaffold_network::ScaffoldNetworkParams;
use rdkit_cxx::chem::*;

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

const OVERRUN: &str = "|atomProp:0.p.v";

#[test]
fn smiles_parsers() {
    let lax = SmilesParserParams {
        strict_cxsmiles: false,
        ..Default::default()
    };
    // Python's values: RDKit 2025.09 parses these.
    let mol = mol_from_smiles_with_params("C |SgD:0:field:value::::|", &lax).unwrap();
    assert_eq!(mol.get_num_atoms(), 1);
    assert_eq!(
        mol_from_smiles("C[C@H](F)Cl |&1:1|")
            .unwrap()
            .get_stereo_groups()
            .len(),
        1
    );
    // A data group's field without its ":" and a polymer group's type
    // without one, which RDKit 2025.09 copies past the end without
    // strictCXSMILES (and refuses itself with it).
    for text in ["C |SgD:0:field", "C |Sg:n"] {
        assert!(is_value(mol_from_smiles_with_params(text, &lax)), "{text}");
        assert!(mol_from_smiles(text).is_err(), "{text}");
    }
    // A stereo group that lists an atom twice, which dangles once the
    // parser removes the hydrogen it names (valgrind: an invalid read).
    for text in [
        "[H][C@](C)(F)Cl |&1:0,0,1|",
        "C[C@H](F)Cl |&1:1,1|",
        "C[C@H](F)Cl |o1:1,o1:1|",
    ] {
        assert!(is_value(mol_from_smiles(text)), "{text}");
        assert!(is_value(mol_from_smarts(text)), "{text}");
    }
    assert!(is_value(reaction_from_smarts("[C:1]>>[C:1] |&1:0,0|")));
}

#[test]
fn every_entry_point() {
    let text = format!("[C:1]>>[C:1] {OVERRUN}");
    assert!(is_value(reaction_from_smarts(&text)));
    assert!(is_value(DeprotectData::new("x", &text, "a", "b")));
    assert!(is_value(ScaffoldNetworkParams::new(&[&text])));

    let smarts = format!("[#6] {OVERRUN}");
    let terms = [SubstructTerm::new("x", &smarts, 1)];
    assert!(is_value(score_substructs_with_terms(
        &mol_from_smiles("C").unwrap(),
        &terms
    )));
    let mut holder = CachedSmilesMolHolder::default();
    assert!(is_value(holder.add_smiles(&format!("C {OVERRUN}"))));
    assert_eq!(holder.add_smiles("C").unwrap(), 0);

    // SmartsMatcher's refused SMARTS leaves it without a pattern, as one that
    // does not parse does.
    assert!(!SmartsMatcher::from_smarts("x", &smarts).is_valid());
    assert!(SmartsMatcher::from_smarts("x", "[#6]").is_valid());
    let mut matcher = SmartsMatcher::from_smarts("x", "[#6]");
    matcher.set_pattern_str(&smarts);
    assert!(!matcher.is_valid());

    // Query definitions, from text and from a file.
    let defs = format!("AcidChloride\t{smarts}\n");
    assert!(is_value(parse_mol_query_def_file_text(&defs)));
    let path =
        std::env::temp_dir().join(format!("rdkit-cxx-query-defs-{}.txt", std::process::id()));
    std::fs::write(&path, &defs).unwrap();
    assert!(is_value(parse_mol_query_def_file(&path)));
    std::fs::write(&path, "AcidChloride\t[#6]\n").unwrap();
    assert_eq!(parse_mol_query_def_file(&path).unwrap().len(), 1);
    std::fs::remove_file(&path).unwrap();

    // An atom's query property that names no query is parsed as SMARTS.
    let mut mol = mol_from_smiles("CC").unwrap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("query", &smarts);
    assert!(is_value(mol_add_recursive_queries(
        &mut mol,
        &BTreeMap::new(),
        "query"
    )));
}

#[test]
fn normalizations() {
    let line = format!("x\t[C:1]>>[C:1] {OVERRUN}");
    let params = CleanupParameters::default();
    assert!(is_value(normalizer_from_data(&line, &params)));
    let path = std::env::temp_dir().join(format!(
        "rdkit-cxx-normalizations-{}.txt",
        std::process::id()
    ));
    std::fs::write(&path, &line).unwrap();
    assert!(is_value(Normalizer::new(&path, 10)));
    // A cleanup parameters' normalizations file is read and checked by the
    // wrapper, which hands RDKit its lines.
    let from_file = CleanupParameters {
        normalizations_file: path.clone(),
        ..Default::default()
    };
    assert!(is_value(normalizer_from_params(&from_file)));
    std::fs::write(&path, "Nitro to N+(O-)=O\t[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]\n").unwrap();
    let normalizer = normalizer_from_params(&from_file).unwrap();
    let nitro = normalizer
        .normalize(&mol_from_smiles("CN(=O)=O").unwrap())
        .unwrap();
    assert_eq!(mol_to_smiles(&nitro).unwrap(), "C[N+](=O)[O-]");
    std::fs::remove_file(&path).unwrap();
    let data = CleanupParameters {
        normalization_data: vec![("x".to_owned(), format!("[C:1]>>[C:1] {OVERRUN}"))],
        ..Default::default()
    };
    assert!(is_value(normalizer_from_params(&data)));
    let pipeline = PipelineOptions {
        normalizer_data: line,
        ..Default::default()
    };
    assert!(is_value(Pipeline::new(&pipeline)));
}

/// A Marvin block for `C[C@H](F)Cl` in an AND stereo group, as RDKit 2025.09.6
/// writes it, with the fluorine's id given as `id`.
fn mrv(id: &str) -> String {
    format!(
        r#"<cml><MDocument><MChemicalStruct><molecule molID="m1"><atomArray><atom id="a1" elementType="C" x2="0.00000" y2="0.00000"/><atom id="a2" elementType="C" x2="1.29904" y2="0.75" mrvStereoGroup="and1"/><atom id="{id}" elementType="F" x2="2.59808" y2="0.00000" mrvStereoGroup="and1"/><atom id="a4" elementType="Cl" x2="1.29904" y2="2.25"/></atomArray><bondArray><bond id="b1" atomRefs2="a2 a1" order="1"><bondStereo>W</bondStereo></bond><bond id="b2" atomRefs2="a2 {id}" order="1"/><bond id="b3" atomRefs2="a2 a4" order="1"/></bondArray></molecule></MChemicalStruct></MDocument></cml>"#
    )
}

#[test]
fn marvin_atom_ids() {
    assert_eq!(mol_from_mrv_block(&mrv("a3")).unwrap().get_num_atoms(), 4);
    // Two atoms with one id put the first in its stereo group twice, which
    // dangles once removed (a segfault in RDKit 2025.09.6's Python).
    assert!(mol_from_mrv_block(&mrv("a2")).is_err());
}
