use rdkit_cxx::Error;
use rdkit_cxx::chem::{
    Atom, BondType, ChiralType, HybridizationType, Mol, RWMol, ReplaceAtomParams, mol_from_smiles,
    mol_to_smiles,
};

fn smiles(mol: &Mol) -> String {
    mol_to_smiles(mol).unwrap()
}

#[test]
fn build_atom_by_atom() {
    let mut rw = RWMol::default();
    assert_eq!(rw.add_atom(&Atom::from_num(6)), 0);
    assert_eq!(rw.add_atom(&Atom::new("O").unwrap()), 1);
    assert_eq!(rw.add_bond_with_order(0, 1, BondType::Single).unwrap(), 1);
    assert_eq!(smiles(&rw), "CO");
    assert_eq!(
        rw.add_bond(0, 5).unwrap_err(),
        Error::Index { index: 5, len: 2 }
    );
    assert!(matches!(rw.add_bond(0, 1), Err(Error::Runtime { .. })));
    assert_eq!(
        rw.get_bond_with_idx(0).unwrap().get_bond_type(),
        BondType::Single
    );
}

#[test]
fn property_cache_is_required_and_reported() {
    let mut rw = RWMol::default();
    rw.add_atom(&Atom::from_num(6));
    rw.add_atom(&Atom::from_num(8));
    rw.add_bond_with_order(0, 1, BondType::Single).unwrap();
    assert!(
        matches!(rw.get_atom_with_idx(0).unwrap().get_total_num_hs(),
                     Err(Error::Runtime { message }) if message.contains("calcImplicitValence"))
    );
    rw.update_property_cache().unwrap();
    assert_eq!(
        rw.get_atom_with_idx(0).unwrap().get_total_num_hs().unwrap(),
        3
    );
}

#[test]
fn strict_property_cache_reports_bad_valence() {
    let mut mol = mol_from_smiles("CC").unwrap();
    mol.get_atom_with_idx_mut(0).unwrap().set_formal_charge(5);
    mol.get_atom_with_idx_mut(0).unwrap().set_num_explicit_hs(4);
    mol.get_atom_with_idx_mut(0).unwrap().set_atomic_num(7);
    mol.get_atom_with_idx_mut(0).unwrap().set_formal_charge(0);
    assert!(matches!(
        mol.update_property_cache(),
        Err(Error::AtomValence { atom_idx: 0, .. })
    ));
    mol.update_property_cache_with_strict(false).unwrap();
    let mut atom = mol.get_atom_with_idx_mut(0).unwrap();
    assert!(matches!(
        atom.update_property_cache(),
        Err(Error::AtomValence { atom_idx: 0, .. })
    ));
    atom.update_property_cache_with_strict(false).unwrap();
}

#[test]
fn edit_through_views() {
    let mut acid = mol_from_smiles("CC(=O)O").unwrap();
    acid.get_atom_with_idx_mut(3).unwrap().set_formal_charge(-1);
    acid.update_property_cache().unwrap();
    assert_eq!(smiles(&acid), "CC(=O)[O-]");

    let mut ethane = mol_from_smiles("CC").unwrap();
    for i in 0..ethane.get_num_atoms() {
        ethane.get_atom_with_idx_mut(i).unwrap().set_isotope(13);
    }
    assert_eq!(smiles(&ethane), "[13CH3][13CH3]");

    let mut ethene = mol_from_smiles("CC").unwrap();
    ethene
        .get_bond_with_idx_mut(0)
        .unwrap()
        .set_bond_type(BondType::Double);
    ethene.update_property_cache().unwrap();
    assert_eq!(smiles(&ethene), "C=C");

    assert_eq!(
        ethene.get_atom_with_idx_mut(9).unwrap_err(),
        Error::Index { index: 9, len: 2 }
    );
    assert_eq!(
        ethene.get_bond_with_idx_mut(9).unwrap_err(),
        Error::Index { index: 9, len: 1 }
    );
}

#[test]
fn mutable_views_read_and_write() {
    let mut mol = mol_from_smiles("CO").unwrap();
    let mut atom = mol.get_atom_with_idx_mut(1).unwrap();
    assert_eq!(
        (atom.get_symbol().unwrap(), atom.get_degree()),
        ("O".to_string(), 1)
    );
    assert_eq!(atom.get_neighbors()[0].get_idx(), 0);
    atom.set_hybridization(HybridizationType::Sp);
    atom.set_chiral_tag(ChiralType::ChiOther);
    atom.set_is_aromatic(true);
    atom.set_no_implicit(true);
    assert_eq!(atom.get_hybridization(), HybridizationType::Sp);
    assert_eq!(atom.get_chiral_tag(), ChiralType::ChiOther);
    assert!(atom.get_is_aromatic());
    assert_eq!(atom.get_owning_mol().get_num_atoms(), 2);

    let mut bond = mol.get_bond_with_idx_mut(0).unwrap();
    bond.set_is_aromatic(true);
    bond.set_is_conjugated(true);
    assert!(bond.get_is_aromatic() && bond.get_is_conjugated());
    assert_eq!(
        (
            bond.get_begin_atom().get_idx(),
            bond.get_end_atom().get_idx()
        ),
        (0, 1)
    );
    bond.set_bond_type(BondType::Other);
    assert!(matches!(
        bond.get_bond_type_as_double(),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn rwmol_edits() {
    let mut rw = RWMol::from(mol_from_smiles("OCC").unwrap());
    rw.remove_atom(0).unwrap();
    assert_eq!(smiles(&rw), "CC");
    rw.replace_atom(0, &Atom::from_num(7)).unwrap();
    assert_eq!(smiles(&rw), "CN");
    let copy = rw.get_mol();
    rw.remove_bond(0, 1).unwrap();
    assert_eq!(
        (smiles(&rw), smiles(&copy)),
        ("C.N".to_string(), "CN".to_string())
    );
    assert_eq!(
        rw.remove_atom(7).unwrap_err(),
        Error::Index { index: 7, len: 2 }
    );
    assert_eq!(
        rw.remove_bond(0, 7).unwrap_err(),
        Error::Index { index: 7, len: 2 }
    );
    let params = ReplaceAtomParams {
        update_label: true,
        preserve_props: true,
    };
    rw.replace_atom_with_params(1, &Atom::from_num(8), &params)
        .unwrap();
    assert_eq!(smiles(&rw), "N.O");
    assert_eq!(
        ReplaceAtomParams::default(),
        ReplaceAtomParams {
            update_label: false,
            preserve_props: false
        }
    );
    let back: Mol = rw.into();
    assert_eq!(back.get_num_atoms(), 2);
}

#[test]
fn atom_constructors() {
    assert_eq!(Atom::new("C").unwrap().get_atomic_num(), 6);
    assert!(matches!(Atom::new("Xx"), Err(Error::Runtime { message }) if message.contains("Xx")));
    let mut atom = Atom::from_num(6);
    atom.set_isotope(13);
    let copy = atom.clone();
    atom.set_isotope(14);
    assert_eq!((atom.get_isotope(), copy.get_isotope()), (14, 13));
    assert!(matches!(
        Atom::from_num(500).get_symbol(),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn rwmol_derefs_to_mol_and_is_send() {
    fn count(mol: &Mol) -> u32 {
        mol.get_num_atoms()
    }
    assert_eq!(count(&RWMol::from(mol_from_smiles("CC").unwrap())), 2);
    fn assert_send<T: Send>() {}
    assert_send::<RWMol>();
    assert_send::<Atom>();
}
