// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Inputs RDKit 2024.09 mishandles and later releases refuse or handle: each
// is refused with Error::Value or Error::Parse before RDKit reads it, and
// the nearest well-formed input still works.

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::rd_fmcs::*;
use rdkit_cxx::chem::rd_rascal_mces::*;
use rdkit_cxx::chem::rd_substruct_library::*;
use rdkit_cxx::chem::scaffolds::rd_scaffold_network::*;
use rdkit_cxx::chem::*;

mod common;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn is_value<T: std::fmt::Debug>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

fn is_parse<T: std::fmt::Debug>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Parse { .. }))
}

#[test]
fn cml_parity_needs_four_neighbors() {
    // The writer reads four neighbors for a tetrahedral atom's parity.
    let mol = smiles("C[C@H](O)N");
    assert!(is_value(mol_to_cml_block(&mol)));
    let hs = add_hs(&mol).unwrap();
    assert!(mol_to_cml_block(&hs).unwrap().contains("atomParity"));
}

#[test]
fn fragmenting_refuses_repeated_bonds_and_counts_past_the_bonds() {
    let mol = smiles("CCC");
    assert!(is_value(fragment_on_bonds(&mol, &[0, 0])));
    assert!(is_value(fragment_on_some_bonds(&mol, &[1, 0, 1])));
    // More bonds to break than given, up to 64 and beyond, make nothing.
    for num_to_break in [3, 64, 100] {
        let params = FragmentOnSomeBondsParams {
            num_to_break,
            ..Default::default()
        };
        assert!(
            fragment_on_some_bonds_with_params(&mol, &[0, 1], &params)
                .unwrap()
                .is_empty()
        );
    }
}

fn mol2(bond: &str) -> String {
    format!(
        "@<TRIPOS>MOLECULE\nethane\n 2 1 0 0 0\nSMALL\nNO_CHARGES\n\n\
         @<TRIPOS>ATOM\n\
         \x20     1 C1         0.0000    0.0000    0.0000 C.3   1  UNL1  0.0000\n\
         \x20     2 C2         1.5000    0.0000    0.0000 C.3   1  UNL1  0.0000\n\
         @<TRIPOS>BOND\n     1 {bond} 1\n"
    )
}

#[test]
fn mol2_bonds_must_join_atoms() {
    assert_eq!(
        mol_from_mol2_block(&mol2("1 2")).unwrap().get_num_bonds(),
        1
    );
    // RDKit 2024.09 requires only one end to be an atom, and 0 wraps.
    for bond in ["1 3", "3 1", "0 2", "1 4294967295"] {
        assert!(is_parse(mol_from_mol2_block(&mol2(bond))), "{bond}");
    }
    let path = common::temp_path("bad.mol2");
    std::fs::write(&path, mol2("1 3")).unwrap();
    assert!(is_parse(mol_from_mol2_file(&path)));
    std::fs::write(&path, mol2("1 2")).unwrap();
    assert_eq!(mol_from_mol2_file(&path).unwrap().get_num_atoms(), 2);
}

fn v3000(sgroups: usize, lines: &[&str]) -> String {
    let mut text = format!(
        "\n  RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
         M  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 {sgroups} 0 0\nM  V30 BEGIN ATOM\n\
         M  V30 1 C 0 0 0 0\nM  V30 2 C 1.3 0 0 0\nM  V30 3 C 2.6 0 0 0\n\
         M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\n\
         M  V30 END BOND\n"
    );
    if sgroups > 0 || !lines.is_empty() {
        text += "M  V30 BEGIN SGROUP\n";
        for line in lines {
            text += &format!("M  V30 {line}\n");
        }
        text += "M  V30 END SGROUP\n";
    }
    text + "M  V30 END CTAB\nM  END\n"
}

#[test]
fn v3000_substance_group_fields_need_their_values() {
    let fine = v3000(
        1,
        &["1 SUP 0 ATOMS=(1 3) XBONDS=(1 2) LABEL=Me SAP=(3 3 2 1)"],
    );
    let mol = mol_from_mol_block(&fine).unwrap();
    assert_eq!(get_mol_substance_groups(&mol).len(), 1);
    let defaults = v3000(1, &["DEFAULT CLASS=x", "1 DAT 0 ATOMS=(1 1) FIELDNAME=a"]);
    mol_from_mol_block(&defaults).unwrap();

    for bad in [
        // RDKit 2024.09 pops the last character of an empty attachment id.
        "1 SUP 0 ATOMS=(1 3) SAP=(3 3 2)",
        "1 SUP 0 ATOMS=(1 3) SAP=(",
        // It leaves the external id, a CSTATE's count, a COMPNO, and an
        // array's first value unset when there is nothing to read.
        "1 DAT",
        "1 SUP 0 CSTATE=(",
        "1 DAT 0 CSTATE=(1",
        "1 DAT 0 COMPNO=",
        "1 DAT 0 ATOMS=(1",
    ] {
        assert!(is_value(mol_from_mol_block(&v3000(1, &[bad]))), "{bad}");
    }
    // The END SGROUP line is read as a group's when the count says two.
    assert!(is_value(mol_from_mol_block(&v3000(
        2,
        &["1 DAT 0 ATOMS=(1 1)"]
    ))));
}

#[test]
fn v3000_refusals_reach_streams_cut_within_a_group_line() {
    // A group line that crosses the scanner's 64 KiB chunks, with a short
    // attachment point after the boundary: RDKit has read the line's start
    // when the stream is refused, and must not parse it cut off.
    let pad = "x".repeat(70_000);
    let line = format!("1 SUP 0 ATOMS=(1 3) LABEL={pad} SAP=(3 3 2)");
    let text = v3000(1, &[&line]) + "$$$$\n";
    let path = common::temp_path("long_sgroup.sdf");
    std::fs::write(&path, &text).unwrap();
    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert!(suppl.mols().all(|m| m.is_err()));
    let mut forward = ForwardSDMolSupplier::new(&path).unwrap();
    assert!(forward.mols().all(|m| m.is_err()));
    assert!(is_value(mol_from_mol_file(&path)));
}

const BRIDGE_PDB: &str = "\
HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
HETATM    2  H1  UNL     1       1.000   0.000   0.000  1.00  0.00           H  \n\
HETATM    3  C2  UNL     1       2.050   0.000   0.000  1.00  0.00           C  \n\
END\n";

#[test]
fn proximity_bonding_refuses_shared_hydrogens() {
    // Proximity bonding would bond the hydrogen to both carbons, and RDKit
    // 2024.09 removes one bond while it iterates over them.
    assert!(is_value(mol_from_pdb_block(BRIDGE_PDB)));
    let no_proximity = MolFromPDBBlockParams {
        proximity_bonding: false,
        remove_hs: false,
        ..Default::default()
    };
    let mol = mol_from_pdb_block_with_params(BRIDGE_PDB, &no_proximity).unwrap();
    assert_eq!((mol.get_num_atoms(), mol.get_num_bonds()), (3, 0));
}

#[test]
fn mcs_verbose_output_with_store_all_is_refused() {
    let mols = [smiles("CCO"), smiles("CCN")];
    let both = MCSParameters {
        verbose: true,
        store_all: true,
        ..Default::default()
    };
    assert!(is_value(find_mcs_with_params(&mols, &both)));
    let verbose = MCSParameters {
        verbose: true,
        ..Default::default()
    };
    assert_eq!(
        find_mcs_with_params(&mols, &verbose).unwrap().num_atoms(),
        2
    );
    // A seed of one atom ends the search before verbose output reads its
    // bonds.
    let seeded = MCSParameters {
        verbose: true,
        initial_seed: "[#6]".to_string(),
        ..Default::default()
    };
    find_mcs_with_params(&mols, &seeded).unwrap();
    let one_atom = [smiles("C"), smiles("CO")];
    assert_eq!(
        find_mcs_with_params(&one_atom, &verbose)
            .unwrap()
            .num_atoms(),
        1
    );
    assert_eq!(find_mcs(&mols).unwrap().num_atoms(), 2);
}

#[test]
fn zipping_with_coordinates_refuses_what_it_would_not_zip() {
    let coords = MolzipParams {
        generate_coordinates: true,
        ..Default::default()
    };
    // RDKit zips one of two R groups of the same SMILES, and reads the
    // other's atoms in the zipped molecule.
    let mut repeated = [
        smiles("c1ccc([*:1])cc1[*:2]"),
        smiles("C[*:1]"),
        smiles("C[*:1]"),
    ];
    assert!(is_value(molzip_fragments_with_params(
        &mut repeated,
        &coords
    )));
    let mut fine = [
        smiles("c1ccc([*:1])cc1[*:2]"),
        smiles("C[*:1]"),
        smiles("N[*:2]"),
    ];
    let zipped = molzip_fragments_with_params(&mut fine, &coords)
        .unwrap()
        .unwrap();
    assert_eq!(mol_to_smiles(&zipped).unwrap(), "Cc1cccc(N)c1");
    assert_eq!(zipped.get_num_conformers(), 1);
    let fragment_on_bonds = MolzipParams {
        label: MolzipLabel::FragmentOnBonds,
        ..coords
    };
    let mut two = [smiles("C[1*].[2*]C"), smiles("C[*:1]")];
    assert!(is_value(molzip_fragments_with_params(
        &mut two,
        &fragment_on_bonds
    )));
}

#[test]
fn replacing_a_core_needs_the_neighbors_of_shared_atoms() {
    // Core atoms 0 and 2 share molecule atom 0, and the match lacks their
    // neighbor.
    let (mol, core) = (smiles("CC"), smiles("CCC"));
    assert!(is_value(replace_core(&mol, &core, &[(0, 0), (2, 0)])));
}

#[test]
fn dropping_a_multithreaded_supplier_mid_read_is_safe() {
    let block = mol_to_mol_block(&smiles("CCO")).unwrap() + "$$$$\n";
    let path = common::temp_path("many.sdf");
    std::fs::write(&path, block.repeat(2000)).unwrap();
    let params = MultithreadedSDMolSupplierParams {
        num_writer_threads: 4,
        size_input_queue: 2,
        size_output_queue: 2,
        ..Default::default()
    };
    for _ in 0..20 {
        let mut suppl = MultithreadedSDMolSupplier::new_with_params(&path, &params).unwrap();
        suppl.mols().next().unwrap().unwrap();
        drop(suppl);
    }
}

#[test]
fn stereo_groups_listing_an_atom_twice_are_refused() {
    // RDKit 2024.09 keeps the duplicate, which dangles once the atom is
    // removed: removing hydrogens, Python crashes on the second.
    for text in [
        "C[C@H](F)Cl |&1:1,1|",
        "[H][C@](C)(F)Cl |&1:0,0,1|",
        "C[C@H](F)C[C@H](F)Cl |o1:1,o1:4,o1:1|",
        "C[C@H](F)Cl |a:1,a:1|",
    ] {
        assert!(is_value(mol_from_smiles(text)), "{text}");
    }
    for text in ["C[C@H](F)Cl |&1:1,&2:1|", "C[C@H](F)Cl |a:1,o1:1|"] {
        assert!(mol_from_smiles(text).is_ok(), "{text}");
    }

    let block = mol_to_v3k_mol_block(&smiles("C[C@H](F)Cl |&1:1|")).unwrap();
    assert!(mol_from_mol_block(&block).is_ok());
    // A list shorter than its count repeats its last atom, or reads one
    // unset (still so in RDKit 2026.09).
    for atoms in ["(2 2 2)", "(2 2)", "(1 )", "(2 )"] {
        let bad = block.replace("ATOMS=(1 2)", &format!("ATOMS={atoms}"));
        assert_ne!(bad, block);
        assert!(is_value(mol_from_mol_block(&bad)), "{atoms}");
    }
    let two = mol_to_v3k_mol_block(&smiles("C[C@H](F)[C@H](F)Cl |&1:1,3|")).unwrap();
    assert!(two.contains("ATOMS=(2 2 4)"), "{two}");
    assert!(mol_from_mol_block(&two).is_ok());
    let repeated = two.replace("ATOMS=(2 2 4)", "ATOMS=(3 2 2 4)");
    assert!(is_value(mol_from_mol_block(&repeated)));
    // Only a collection block's lines are read as stereo groups.
    let elsewhere = block.replace(
        "M  V30 END CTAB",
        "M  V30 BEGIN OBJ3D\nM  V30 MDLV30/STERAC1 ATOMS=(1 )\nM  V30 END OBJ3D\nM  V30 END CTAB",
    );
    assert_ne!(elsewhere, block);
    assert!(mol_from_mol_block(&elsewhere).is_ok());

    let mrv = mol_to_mrv_block(&smiles("C[C@H](F)Cl |&1:1|")).unwrap();
    let twice = mrv.replace(
        "</atomArray>",
        "<atom id=\"a2\" elementType=\"Br\" x2=\"3\" y2=\"3\" mrvStereoGroup=\"and1\"/></atomArray>",
    );
    assert_ne!(twice, mrv);
    assert!(is_parse(mol_from_mrv_block(&twice)));

    let json = mol_to_json(&smiles("C[C@H](F)Cl |&1:1|")).unwrap();
    let doubled = json.replace("\"atoms\":[1]", "\"atoms\":[1,1]");
    assert_ne!(doubled, json);
    assert!(is_value(json_to_mols(&doubled)));
    assert_eq!(json_to_mols(&json).unwrap().len(), 1);
}

// RDKit 2024.09's CXSMILES parser reads past the end of the text for an
// atomProp block whose last value ends it (#9361); every entry point that
// parses caller text as SMILES or SMARTS with CX extensions checks it.
// ("|atomProp:0.a", a name without a value, RDKit refuses.)
const OVERRUN: &str = "|atomProp:0.p.v";

#[test]
fn cx_extensions_are_checked_at_every_entry_point() {
    let rascal = RascalOptions {
        equivalent_atoms: format!("[#6]\t{OVERRUN}"),
        ..Default::default()
    };
    assert!(is_value(find_mces_with_opts(
        &smiles("CCO"),
        &smiles("CCN"),
        &rascal
    )));
    let fine = RascalOptions {
        equivalent_atoms: "[#6] [#7]".to_string(),
        ..Default::default()
    };
    find_mces_with_opts(&smiles("CCO"), &smiles("CCN"), &fine).unwrap();

    let breaker = "[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]";
    assert!(is_value(ScaffoldNetworkParams::new(&[&format!(
        "{breaker}{OVERRUN}"
    )])));
    ScaffoldNetworkParams::new(&[breaker]).unwrap();

    let mut holder = CachedSmilesMolHolder::default();
    assert!(is_value(holder.add_smiles(&format!("CC {OVERRUN}"))));
    holder.add_smiles("CC |$_AV:a;b$|").unwrap();

    assert!(is_value(parse_mol_query_def_file_text(&format!(
        "Alkyl\tC {OVERRUN}\n"
    ))));
    assert_eq!(
        parse_mol_query_def_file_text("Alkyl\t[CX4]\n")
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn mol_block_smarts_are_checked() {
    let block = mol_to_mol_block(&smiles("CC")).unwrap();
    let with = |line: &str| block.replace("M  END", &format!("{line}\nM  END"));
    // An "M  MRV SMA" line's SMARTS, from its 16th byte.
    assert!(is_value(mol_from_mol_block(&with(&format!(
        "M  MRV SMA   1 C {OVERRUN}"
    )))));
    let sma = mol_from_mol_block(&with("M  MRV SMA   1 [#6]")).unwrap();
    assert_eq!(sma.get_num_atoms(), 2);
    // A data group's field, which a SMARTSQ query parses, also when it runs
    // across an SCD line into the SED's.
    let data = |fields: &[&str]| {
        let mut lines = vec![
            "M  STY  1   1 DAT".to_string(),
            "M  SAL   1  1   1".to_string(),
            "M  SDT   1 MRV_QUERY                 F".to_string(),
        ];
        lines.extend(fields.iter().map(|f| f.to_string()));
        with(&lines.join("\n"))
    };
    assert!(is_value(mol_from_mol_block(&data(&[&format!(
        "M  SED   1 C {OVERRUN}"
    )]))));
    assert!(is_value(mol_from_mol_block(&data(&[
        &format!("M  SCD   1 {} |atomP", "C".repeat(62)),
        "M  SED   1 rop:0.p.v",
    ]))));
    mol_from_mol_block(&data(&["M  SED   1 [#6]"])).unwrap();
    // Lines of another group within a field would be skipped by RDKit.
    assert!(is_value(mol_from_mol_block(&data(&[
        "M  SCD   1 C",
        "M  SED   2 C",
    ]))));
    // A V3000 FIELDDATA value.
    let v3k = v3000(
        1,
        &[&format!(
            "1 DAT 0 ATOMS=(1 1) FIELDNAME=MRV_QUERY QUERYTYPE=SMARTSQ FIELDDATA=\"C {OVERRUN}\""
        )],
    );
    assert!(is_value(mol_from_mol_block(&v3k)));
    let v3k = v3000(
        1,
        &["1 DAT 0 ATOMS=(1 1) FIELDNAME=MRV_QUERY QUERYTYPE=SMARTSQ FIELDDATA=\"[#6]\""],
    );
    mol_from_mol_block(&v3k).unwrap();
}

#[test]
fn normalization_files_are_read_and_checked_once() {
    let path = common::temp_path("normalizations.txt");
    std::fs::write(&path, format!("Bad\t[C:1]>>[N:1]{OVERRUN}\n")).unwrap();
    assert!(is_value(Normalizer::new(&path, 200)));
    let params = CleanupParameters {
        normalizations_file: path.clone(),
        ..Default::default()
    };
    assert!(is_value(normalize_with_params(&smiles("CCO"), &params)));
    assert!(is_value(cleanup_with_params(&smiles("CCO"), &params)));
    assert!(is_value(normalizer_from_params(&params)));
    // The ops that do not normalize do not read it.
    reionize_with_params(&smiles("CC(=O)O"), &params).unwrap();

    std::fs::write(&path, "// a comment\nAmine\t[O:1]>>[N:1]\n").unwrap();
    let normalizer = Normalizer::new(&path, 200).unwrap();
    assert_eq!(
        mol_to_smiles(&normalizer.normalize(&smiles("CCO")).unwrap()).unwrap(),
        "CCN"
    );
    let params = CleanupParameters {
        normalizations_file: path.clone(),
        ..Default::default()
    };
    let normalized = normalize_with_params(&smiles("CCO"), &params).unwrap();
    assert_eq!(mol_to_smiles(&normalized).unwrap(), "CCN");
    std::fs::write(&path, "// nothing\n").unwrap();
    assert!(is_value(normalize_with_params(&smiles("CCO"), &params)));
}
