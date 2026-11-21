// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Manipulating Molecules" chapter.

use rdkit_cxx::Result;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::{LargestFragmentChooser, Uncharger};
use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
use rdkit_cxx::chem::{
    Mol, RWMol, fragment_on_bonds, fragment_on_some_bonds, get_mol_frags_as_mols, mol_from_smarts,
    mol_from_smiles, mol_to_smiles, replace_substructs,
};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[
    ("RDKitCB_7", "Create Fragments", create_fragments),
    ("RDKitCB_31", "Largest Fragment", largest_fragment),
    ("RDKitCB_29", "Sidechain-Core Enumeration", sidechain_core_enumeration),
    ("RDKitCB_33", "Neutralizing Molecules", neutralizing_molecules),
];

fn smiles_of(mols: &[Mol]) -> Result<Vec<String>> {
    mols.iter().map(mol_to_smiles).collect()
}

/// Cuts a molecule on bonds, three ways. By Paolo Tosco.
fn create_fragments() -> Result<()> {
    // Explicit bonds in the SMILES make the bond indices easy to follow.
    let mol = mol_from_smiles("O-C-C-C-C-N")?;

    // FragmentOnBonds makes every cut at once, in one molecule, and fills the
    // empty valences with dummy atoms.
    let fragmented = fragment_on_bonds(&mol, &[0, 2, 4])?;
    println!("{}", mol_to_smiles(&fragmented)?);
    // GetMolFrags splits it into one molecule per fragment.
    let fragments = smiles_of(&get_mol_frags_as_mols(&fragmented)?)?;
    println!("{fragments:?}");
    assert_eq!(fragments.len(), 4);

    // FragmentOnSomeBonds makes one cut at a time, returning one molecule per
    // choice of bond.
    let one_cut = smiles_of(&fragment_on_some_bonds(&mol, &[0, 2, 4])?)?;
    println!("{one_cut:?}");
    assert_eq!(one_cut.len(), 3);

    // Removing the bonds leaves the empty valences to implicit hydrogens.
    // A batch edit removes them together, so the bond indices hold until
    // the commit.
    let mut rwmol = RWMol::from(mol.clone());
    rwmol.begin_batch_edit()?;
    for idx in [0, 2, 4] {
        let bond = rwmol.get_bond_with_idx(idx)?;
        let (begin, end) = (bond.get_begin_atom_idx(), bond.get_end_atom_idx());
        rwmol.remove_bond(begin, end)?;
    }
    rwmol.commit_batch_edit();
    let pieces = smiles_of(&get_mol_frags_as_mols(&rwmol)?)?;
    println!("{pieces:?}");
    assert_eq!(pieces, ["O", "CC", "CC", "N"]);
    Ok(())
}

/// Keeps the fragment with the most atoms. By Andrew Dalke and Susan Leung.
fn largest_fragment() -> Result<()> {
    let mol = mol_from_smiles(
        "CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl.CO.\
         C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
    )?;
    let expected = "Nc1nc2ncc(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)nc2c(=O)[nH]1";

    let largest = get_mol_frags_as_mols(&mol)?
        .into_iter()
        .max_by_key(|frag| frag.get_num_atoms())
        .unwrap_or(mol.clone());
    println!("{}", mol_to_smiles(&largest)?);
    assert_eq!(mol_to_smiles(&largest)?, expected);

    // MolStandardize does the same.
    let largest = LargestFragmentChooser::default().choose(&mol)?;
    println!("{}", mol_to_smiles(&largest)?);
    assert_eq!(mol_to_smiles(&largest)?, expected);
    Ok(())
}

/// Attaches sidechains to a core, with a reaction and with
/// ReplaceSubstructs. By Chris Earnshaw, Stephen Roughley, Greg Landrum, and
/// Vincent Scalfani.
fn sidechain_core_enumeration() -> Result<()> {
    // The core is *c1c(C)cccc1(O) and the chain CN*.
    let mut rxn = reaction_from_smarts("[c:1][#0].[#0][*:2]>>[c:1]-[*:2]")?;
    let reactants = [mol_from_smiles("*c1c(C)cccc1(O)")?, mol_from_smiles("CN*")?];
    let products = rxn.run_reactants(&reactants)?;
    println!("{} {}", products.len(), products[0].len());
    assert_eq!((products.len(), products[0].len()), (1, 1));
    let product = mol_to_smiles(&products[0][0])?;
    println!("{product}");
    assert_eq!(product, "CNc1c(C)cccc1O");

    // A sidechain whose attachment atom is atom 0 can replace the dummy atom
    // directly.
    let dummy = mol_from_smarts("[#0]")?;
    let core = mol_from_smiles("*c1c(C)cccc1(O)")?;
    let chain = mol_from_smiles("NC")?;
    let products = replace_substructs(&core, &dummy, &chain)?;
    let product = mol_to_smiles(&products[0])?;
    println!("{product}");
    assert_eq!(product, "CNc1c(C)cccc1O");

    // Alkyl chains on an imidazolium core.
    let core = mol_from_smiles("*[n+]1cc[nH]c1")?;
    let mut product_smiles = Vec::new();
    for chain in ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"] {
        let products = replace_substructs(&core, &dummy, &mol_from_smiles(chain)?)?;
        product_smiles.push(mol_to_smiles(&products[0])?);
    }
    println!("{product_smiles:?}");
    assert_eq!(
        product_smiles,
        [
            "C[n+]1cc[nH]c1",
            "CC[n+]1cc[nH]c1",
            "CCC[n+]1cc[nH]c1",
            "CCCC[n+]1cc[nH]c1",
            "CCCCC[n+]1cc[nH]c1",
            "CCCCCC[n+]1cc[nH]c1",
        ]
    );
    Ok(())
}

/// Neutralizes +1 and -1 atoms by adding or removing hydrogens, leaving
/// charge-separated groups such as nitro alone. By Noel O'Boyle.
fn neutralize_atoms(mol: &mut Mol) -> Result<()> {
    let pattern = mol_from_smarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")?;
    for matched in mol.get_substruct_matches(&pattern)? {
        let mut atom = mol.get_atom_with_idx_mut(matched[0])?;
        let charge = atom.get_formal_charge();
        let num_hs = atom.get_total_num_hs()?;
        atom.set_formal_charge(0);
        let num_hs = num_hs
            .checked_add_signed(-charge)
            .expect("the pattern matches +1 atoms only when they have a hydrogen");
        atom.set_num_explicit_hs(num_hs);
        atom.update_property_cache()?;
    }
    Ok(())
}

/// Neutralizes charged atoms, and compares MolStandardize's Uncharger, which
/// keeps a negative charge that balances a positive one it cannot remove.
/// By Noel O'Boyle, adapted for RDKit by Vincent Scalfani.
fn neutralizing_molecules() -> Result<()> {
    let smiles = [
        "CC(CNC[O-])[N+]([O-])=O",
        "C[N+](C)(C)CCC([O-])=O",
        "[O-]C1=CC=[N+]([O-])C=C1",
        "[O-]CCCN=[N+]=[N-]",
        "C[NH+](C)CC[S-]",
        "CP([O-])(=O)OC[NH3+]",
    ];

    let mut neutralized = Vec::new();
    for smi in smiles {
        let mut mol = mol_from_smiles(smi)?;
        neutralize_atoms(&mut mol)?;
        println!("{}", mol_to_smiles(&mol)?);
        neutralized.push(mol_to_smiles(&mol)?);
    }
    assert_eq!(
        neutralized,
        [
            "CC(CNCO)[N+](=O)[O-]",
            "C[N+](C)(C)CCC(=O)O",
            "[O-][n+]1ccc(O)cc1",
            "[N-]=[N+]=NCCCO",
            "CN(C)CCS",
            "CP(=O)(O)OCN",
        ]
    );

    let uncharger = Uncharger::default();
    let mut uncharged = Vec::new();
    for smi in smiles {
        let mol = uncharger.uncharge(&mol_from_smiles(smi)?)?;
        println!("{}", mol_to_smiles(&mol)?);
        uncharged.push(mol_to_smiles(&mol)?);
    }
    assert_eq!(
        uncharged,
        [
            "CC(CNC[O-])[N+](=O)[O-]",
            "C[N+](C)(C)CCC(=O)[O-]",
            "[O-]c1cc[n+]([O-])cc1",
            "[N-]=[N+]=NCCC[O-]",
            "C[NH+](C)CC[S-]",
            "CP(=O)([O-])OC[NH3+]",
        ]
    );
    Ok(())
}
