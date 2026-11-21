use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles, mol_to_smiles};

fn main() -> Result<(), rdkit_cxx::Error> {
    let mol = mol_from_smiles("Oc1ccccc1O")?;
    let hydroxyl = mol_from_smarts("[OX2H]")?;
    println!("{} {:?}", mol_to_smiles(&mol)?, mol.get_substruct_matches(&hydroxyl)?);
    println!("qed {:.4}", rdkit_cxx::chem::qed::qed(&mol)?);
    Ok(())
}
