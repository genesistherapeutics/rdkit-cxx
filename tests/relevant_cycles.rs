// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Every path to RDKit's relevant-cycle enumeration refuses a molecule with
//! more relevant cycles than RDKit can enumerate (relevant_cycles.h).
//!
//! The molecules are rings of "diamonds": `n` hub atoms, each joined to the
//! next by two paths of two bonds, which have 2^n + n relevant cycles. 64
//! diamonds give 2^64 + 64, whose conversion to size_t in
//! RingInfo::atomRelevantCycles is undefined; on x86-64 RDKit then
//! enumerates cycles until memory runs out, so no test may reach it
//! unguarded. 10 diamonds give 1034, which every path reads.

mod common;

use std::io::Cursor;

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::standardize_smiles;
use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};
use rdkit_cxx::chem::rd_mol_interchange::mol_to_json;
use rdkit_cxx::chem::rdinchi::inchi_to_mol;
use rdkit_cxx::chem::*;

/// The atoms of a ring of `n` diamonds: hub `i` is atom `3i`, and its paths
/// to the next hub go through atoms `3i + 1` and `3i + 2`.
fn diamond_bonds(n: u32) -> Vec<(u32, u32)> {
    let mut bonds = Vec::new();
    for i in 0..n {
        let next = 3 * ((i + 1) % n);
        for middle in [3 * i + 1, 3 * i + 2] {
            bonds.push((3 * i, middle));
            bonds.push((middle, next));
        }
    }
    bonds
}

/// A ring of `n` diamonds, with its rings found by `fast_find_rings`, so
/// not a symmetrized SSSR.
fn diamonds(n: u32) -> RWMol {
    let mut rw = RWMol::default();
    for _ in 0..3 * n {
        rw.add_atom(&Atom::from_num(6));
    }
    for (a, b) in diamond_bonds(n) {
        rw.add_bond_with_order(a, b, BondType::Single).unwrap();
    }
    rw.update_property_cache().unwrap();
    fast_find_rings(&mut rw).unwrap();
    rw
}

/// SMILES for a ring of `n` diamonds.
fn diamonds_smiles(n: u32) -> String {
    // The hubs' second paths close on the next hub through %12 and %13 in
    // turn, and the last hub's paths on the first through %10 and %11.
    let mut smiles = String::from("C%10%11(C%12)C");
    for i in 1..n {
        let (close, open) = if i % 2 == 1 { (12, 13) } else { (13, 12) };
        if i == n - 1 {
            smiles += &format!("C%{close}(C%10)C%11");
        } else {
            smiles += &format!("C%{close}(C%{open})C");
        }
    }
    smiles
}

/// A V2000 mol block for a ring of `n` diamonds.
fn diamonds_v2000(n: u32) -> String {
    let bonds = diamond_bonds(n);
    let mut block = format!(
        "diamonds\n     RDKit          2D\n\n{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
        3 * n,
        bonds.len()
    );
    for i in 0..3 * n {
        let angle = f64::from(i) * std::f64::consts::TAU / f64::from(3 * n);
        block += &format!(
            "{:10.4}{:10.4}{:10.4} C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            20.0 * angle.cos(),
            20.0 * angle.sin(),
            0.0
        );
    }
    for (a, b) in bonds {
        block += &format!("{:3}{:3}  1  0\n", a + 1, b + 1);
    }
    block + "M  END\n"
}

/// A V3000 mol block for a ring of `n` diamonds, with each bond line
/// continued across two physical lines.
fn diamonds_v3000(n: u32) -> String {
    let bonds = diamond_bonds(n);
    let mut block = String::from(
        "diamonds\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n",
    );
    block += "M  V30 BEGIN CTAB\n";
    block += &format!("M  V30 COUNTS {} {} 0 0 0\n", 3 * n, bonds.len());
    block += "M  V30 BEGIN ATOM\n";
    for i in 0..3 * n {
        block += &format!("M  V30 {} C {}.0 0.0 0.0 0\n", i + 1, i);
    }
    block += "M  V30 END ATOM\nM  V30 BEGIN BOND\n";
    for (i, (a, b)) in bonds.into_iter().enumerate() {
        block += &format!("M  V30 {} 1 {} -\nM  V30 {}\n", i + 1, a + 1, b + 1);
    }
    block + "M  V30 END BOND\nM  V30 END CTAB\nM  END\n"
}

/// A chain of `n` diamonds, the last hub not joined to the first: `n`
/// independent cycles, but only `n` relevant ones, which every path reads
/// however many independent cycles there are.
fn diamond_chain(n: u32) -> RWMol {
    let mut rw = RWMol::default();
    for _ in 0..3 * n + 1 {
        rw.add_atom(&Atom::from_num(6));
    }
    for i in 0..n {
        for middle in [3 * i + 1, 3 * i + 2] {
            rw.add_bond_with_order(3 * i, middle, BondType::Single)
                .unwrap();
            rw.add_bond_with_order(middle, 3 * i + 3, BondType::Single)
                .unwrap();
        }
    }
    rw.update_property_cache().unwrap();
    rw
}

/// A ring of `n` diamonds without ring info.
fn diamonds_raw(n: u32) -> RWMol {
    let mut rw = RWMol::default();
    for _ in 0..3 * n {
        rw.add_atom(&Atom::from_num(6));
    }
    for (a, b) in diamond_bonds(n) {
        rw.add_bond_with_order(a, b, BondType::Single).unwrap();
    }
    rw.update_property_cache().unwrap();
    rw
}

/// Atom positions on a circle wide enough that no two are within bonding
/// distance (the PDB reader's proximity bonds).
fn position(i: u32, atoms: u32) -> (f64, f64) {
    let radius = f64::from(atoms) * 3.0 / std::f64::consts::TAU;
    let angle = f64::from(i) * std::f64::consts::TAU / f64::from(atoms);
    (radius * angle.cos(), radius * angle.sin())
}

/// A Mol2 block for a ring of `n` diamonds.
fn diamonds_mol2(n: u32) -> String {
    let bonds = diamond_bonds(n);
    let mut block = format!(
        "@<TRIPOS>MOLECULE\ndiamonds\n{} {} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n",
        3 * n,
        bonds.len()
    );
    for i in 0..3 * n {
        let (x, y) = position(i, 3 * n);
        block += &format!(
            "{:>7} C{:<6}{x:>10.4}{y:>10.4}{:>10.4} C.3 1 RES1 0.0000\n",
            i + 1,
            i + 1,
            0.0
        );
    }
    block += "@<TRIPOS>BOND\n";
    for (i, (a, b)) in bonds.into_iter().enumerate() {
        block += &format!("{:>6}{:>6}{:>6} 1\n", i + 1, a + 1, b + 1);
    }
    block
}

/// A PDB block for a ring of `n` diamonds, its bonds in CONECT records.
fn diamonds_pdb(n: u32) -> String {
    let mut block = String::new();
    for i in 0..3 * n {
        let (x, y) = position(i, 3 * n);
        block += &format!(
            "HETATM{:>5}  C   UNL A   1    {x:>8.3}{y:>8.3}{:>8.3}  1.00  0.00           C\n",
            i + 1,
            0.0
        );
    }
    let mut neighbors = vec![Vec::new(); 3 * n as usize];
    for (a, b) in diamond_bonds(n) {
        neighbors[a as usize].push(b + 1);
        neighbors[b as usize].push(a + 1);
    }
    for (i, list) in neighbors.iter().enumerate() {
        for chunk in list.chunks(4) {
            block += &format!("CONECT{:>5}", i + 1);
            for j in chunk {
                block += &format!("{j:>5}");
            }
            block += "\n";
        }
    }
    block + "END\n"
}

const INCHI_10: &str = "InChI=1S/C30H40/c1-21-2-22(1)5-24(6-22)9-26(10-24)13-28(14-26)17-30(18-28)19-29(20-30)15-27(16-29)11-25(12-27)7-23(3-21,4-21)8-25/h1-20H2";
const INCHI_64: &str = concat!(
    "InChI=1S/C192H256/c1-129-2-130(1)5-132(6-130)9-134(10-132)13-136(14-134)",
    "17-138(18-136)21-140(22-138)25-142(26-140)29-144(30-142)33-146(34-144)37",
    "-148(38-146)41-150(42-148)45-152(46-150)49-154(50-152)53-156(54-154)57-1",
    "58(58-156)61-160(62-158)65-162(66-160)69-164(70-162)73-166(74-164)77-168",
    "(78-166)81-170(82-168)85-172(86-170)89-174(90-172)93-176(94-174)97-178(9",
    "8-176)101-180(102-178)105-182(106-180)109-184(110-182)113-186(114-184)11",
    "7-188(118-186)121-190(122-188)125-192(126-190)127-191(128-192)123-189(12",
    "4-191)119-187(120-189)115-185(116-187)111-183(112-185)107-181(108-183)10",
    "3-179(104-181)99-177(100-179)95-175(96-177)91-173(92-175)87-171(88-173)8",
    "3-169(84-171)79-167(80-169)75-165(76-167)71-163(72-165)67-161(68-163)63-",
    "159(64-161)59-157(60-159)55-155(56-157)51-153(52-155)47-151(48-153)43-14",
    "9(44-151)39-147(40-149)35-145(36-147)31-143(32-145)27-141(28-143)23-139(",
    "24-141)19-137(20-139)15-135(16-137)11-133(12-135)7-131(3-129,4-129)8-133",
    "/h1-128H2",
);

fn refused<T: std::fmt::Debug>(result: Result<T, Error>) -> bool {
    match result {
        Err(Error::Value { message }) => message.contains("relevant cycles"),
        other => panic!("expected a refusal, got {other:?}"),
    }
}

fn is_refusal<T>(result: &Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { message }) if message.contains("relevant cycles"))
}

fn mol_params(sanitize: bool, remove_hs: bool) -> MolFromMolBlockParams {
    MolFromMolBlockParams {
        sanitize,
        remove_hs,
        ..Default::default()
    }
}

const PARAMS: [(bool, bool); 3] = [(true, true), (true, false), (false, false)];

#[test]
fn mol_blocks() {
    for block in [diamonds_v2000(10), diamonds_v3000(10)] {
        for (sanitize, remove_hs) in PARAMS {
            let mol =
                mol_from_mol_block_with_params(&block, &mol_params(sanitize, remove_hs)).unwrap();
            assert_eq!(mol.get_num_atoms(), 30);
            assert_eq!(mol.get_ring_info().num_rings().unwrap(), 1034);
        }
    }
    for block in [diamonds_v2000(64), diamonds_v3000(64)] {
        for (sanitize, remove_hs) in PARAMS {
            assert!(refused(mol_from_mol_block_with_params(
                &block,
                &mol_params(sanitize, remove_hs)
            )));
        }
        let path = common::temp_path("diamonds.mol");
        std::fs::write(&path, &block).unwrap();
        assert!(refused(mol_from_mol_file(&path)));
    }
    // Past the bound, the exact count lets a block with few relevant cycles
    // through.
    let chain = mol_to_mol_block(&diamond_chain(40).into()).unwrap();
    let mol = mol_from_mol_block(&chain).unwrap();
    assert_eq!(mol.get_ring_info().num_rings().unwrap(), 40);
}

/// SD text whose `bad` record follows enough small records to straddle the
/// checked stream's 64 KiB chunks, then a last small record.
fn sd_text(bad: &str) -> String {
    let small = diamonds_v2000(2) + "$$$$\n";
    let mut sd = String::new();
    while sd.len() + small.len() < (1 << 16) - bad.len() / 2 {
        sd += &small;
    }
    sd + bad + "$$$$\n" + &small
}

#[test]
fn sd_suppliers() {
    for block in [diamonds_v2000(64), diamonds_v3000(64)] {
        let sd = sd_text(&block);
        assert!(sd.find("diamonds\n").is_some());
        let mut suppl = SDMolSupplier::default();
        suppl.set_data(&sd).unwrap();
        assert!(suppl.mols().any(|mol| is_refusal(&mol)));
        let mut forward = ForwardSDMolSupplier::from_fileobj(Cursor::new(sd.clone())).unwrap();
        assert!(forward.mols().any(|mol| is_refusal(&mol)));
        let path = common::temp_path("diamonds.sdf");
        std::fs::write(&path, &sd).unwrap();
        let mut forward = ForwardSDMolSupplier::new(&path).unwrap();
        assert!(forward.mols().any(|mol| is_refusal(&mol)));
        let mut threaded = MultithreadedSDMolSupplier::new(&path).unwrap();
        assert!(threaded.mols().any(|mol| is_refusal(&mol)));
    }
    for block in [diamonds_v2000(10), diamonds_v3000(10)] {
        let sd = sd_text(&block);
        let mut suppl = SDMolSupplier::default();
        suppl.set_data(&sd).unwrap();
        let rings: Vec<u32> = suppl
            .mols()
            .map(|mol| mol.unwrap().get_ring_info().num_rings().unwrap())
            .collect();
        assert!(rings.contains(&1034));
        let path = common::temp_path("diamonds10.sdf");
        std::fs::write(&path, &sd).unwrap();
        let mut threaded = MultithreadedSDMolSupplier::new(&path).unwrap();
        assert!(threaded.mols().all(|mol| mol.is_ok()));
    }
    // A record held for its exact count, across the chunk, comes through.
    let chain = mol_to_mol_block(&diamond_chain(40).into()).unwrap();
    let sd = sd_text(&chain);
    let rings = |mols: Vec<Result<Mol, Error>>| {
        let rings: Vec<u32> = mols
            .into_iter()
            .map(|mol| mol.unwrap().get_ring_info().num_rings().unwrap())
            .collect();
        assert!(rings.contains(&40));
    };
    let mut suppl = SDMolSupplier::default();
    suppl.set_data(&sd).unwrap();
    rings(suppl.mols().collect());
    let mut forward = ForwardSDMolSupplier::from_fileobj(Cursor::new(sd.clone())).unwrap();
    rings(forward.mols().collect());
    let path = common::temp_path("chain.sdf");
    std::fs::write(&path, &sd).unwrap();
    let mut threaded = MultithreadedSDMolSupplier::new(&path).unwrap();
    rings(threaded.mols().collect());
}

#[test]
fn smiles() {
    assert!(refused(mol_from_smiles(&diamonds_smiles(64))));
    assert_eq!(
        mol_from_smiles(&diamonds_smiles(10))
            .unwrap()
            .get_ring_info()
            .num_rings()
            .unwrap(),
        1034
    );
    // Many independent cycles but few relevant ones parse whatever the
    // bound.
    let chain = mol_to_smiles(&diamond_chain(40).into()).unwrap();
    let mol = mol_from_smiles(&chain).unwrap();
    assert_eq!(mol.get_ring_info().num_rings().unwrap(), 40);
}

/// Checks `$e`: refused for a ring of 64 diamonds, otherwise accepted.
macro_rules! check {
    ($big:expr, $e:expr) => {
        if $big {
            assert!(refused($e));
        } else {
            $e.unwrap();
        }
    };
}

#[test]
fn perception() {
    for n in [10, 64] {
        let big = n == 64;
        let mol: Mol = diamonds(n).into();
        check!(big, mol_to_smiles(&mol));
        check!(big, assign_stereochemistry(&mut mol.clone()));
        check!(big, find_potential_stereo(&mut mol.clone()));
        check!(big, mol_to_json(&mol));
        check!(big, mol_hash(&mol, HashFunction::CanonicalSmiles));
        check!(big, mol_hash(&mol, HashFunction::AnonymousGraph));
        let morgan = get_morgan_generator_with_params(&GetMorganGeneratorParams {
            include_chirality: true,
            ..GetMorganGeneratorParams::default()
        })
        .unwrap();
        check!(big, morgan.get_fingerprint_as_num_py(&mol));
        let mut raw: Mol = diamonds_raw(n).into();
        check!(big, set_aromaticity(&mut raw));
        check!(big, standardize_smiles(&diamonds_smiles(n)));
    }
}

#[test]
fn formats() {
    for n in [10, 64] {
        let big = n == 64;
        let rings = |mol: Mol| assert_eq!(mol.get_ring_info().num_rings().unwrap(), 1034);
        if big {
            assert!(refused(mol_from_mol2_block(&diamonds_mol2(n))));
            assert!(refused(mol_from_pdb_block(&diamonds_pdb(n))));
            assert!(refused(inchi_to_mol(INCHI_64)));
        } else {
            rings(mol_from_mol2_block(&diamonds_mol2(n)).unwrap());
            rings(mol_from_pdb_block(&diamonds_pdb(n)).unwrap());
            rings(inchi_to_mol(INCHI_10).unwrap().mol.unwrap());
        }
    }
}

#[test]
fn reaction_products() {
    let mut rxn = reaction_from_smarts("[C:1]>>[C:1]").unwrap();
    for n in [10, 64] {
        // A directed single bond has the products' double bonds' neighbors'
        // directions set, which symmetrizes their rings.
        let mut mol: Mol = diamonds(n).into();
        mol.get_bond_with_idx_mut(0)
            .unwrap()
            .set_bond_dir(BondDir::Endupright);
        check!(n == 64, rxn.run_reactants(&[mol]));
    }
}

#[test]
fn molzip() {
    for n in [10, 64] {
        let mut frag = diamonds(n);
        let dummy = frag.add_atom(&Atom::from_num(0));
        frag.get_atom_with_idx_mut(dummy)
            .unwrap()
            .set_atom_map_num(1);
        frag.add_bond_with_order(1, dummy, BondType::Single)
            .unwrap();
        let mut parts = vec![mol_from_smiles("C[*:1]").unwrap(), frag.into()];
        check!(n == 64, molzip_fragments(&mut parts));
    }
}

#[test]
fn smiles_and_tdt_suppliers() {
    for n in [10, 64] {
        let big = n == 64;
        let text = format!(
            "smiles name\nC first\n{} ring\nCC last\n",
            diamonds_smiles(n)
        );
        let mut suppl = SmilesMolSupplier::default();
        suppl.set_data(&text).unwrap();
        let results: Vec<_> = suppl.mols().collect();
        assert_eq!(results.iter().any(is_refusal), big);
        let path = common::temp_path(&format!("diamonds{n}.smi"));
        std::fs::write(&path, &text).unwrap();
        let mut threaded = MultithreadedSmilesMolSupplier::new(&path).unwrap();
        let results: Vec<_> = threaded.mols().collect();
        assert_eq!(results.iter().any(is_refusal), big);
        let tdt = format!("$SMI<C>\n|\n$SMI<{}>\n|\n$SMI<CC>\n|\n", diamonds_smiles(n));
        let mut suppl = TDTMolSupplier::default();
        suppl.set_data(&tdt).unwrap();
        let results: Vec<_> = suppl.mols().collect();
        assert_eq!(results.iter().any(is_refusal), big);
    }
}
