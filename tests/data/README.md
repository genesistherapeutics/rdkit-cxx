# Test data

From RDKit's source tree, under RDKit's BSD license:

- `zim.head100.fpb`, `zinc_random200.1.patt.fpb`, `zinc_random200.2.patt.fpb`:
  `Code/DataStructs/testData`
- `boronic1.rxn`: `rdkit/Chem/SimpleEnum/test_data`
- `reaction1.no_metadata.png`: `Code/GraphMol/ChemReactions/testData`
- `rxn1.cdxml`: `Code/GraphMol/test_data/CDXML`
- `AlexRxn.mrv`: `Code/GraphMol/MarvinParse/test_data`
- `RiboseFullname.mol`: `Code/GraphMol/FileParsers/test_data/macromols`, with
  trailing whitespace trimmed
- `amide_space.txt`, `amide_space_error.txt`, `triazole_space.txt`,
  `extended_query.csv`, `doebner_miller_space.txt`:
  `Code/GraphMol/SynthonSpaceSearch/data`
- `Pharm2D_BaseFeatures.fdef`, `Pharm2D_orderBug.sdf`: `BaseFeatures.fdef` and
  `orderBug.sdf` in `rdkit/Chem/Pharm2D/test_data`, with a final newline added
  and trailing whitespace trimmed
- `Pharm3D_BaseFeatures.fdef`, `Pharm3D_Issue268.fdef`,
  `Pharm3D_Issue268_Mol1.mol`, `Pharm3D_Issue268_Mol2.mol`, `Pharm3D_mol-r.mol`,
  `Pharm3D_mol-s.mol`, `Pharm3D_mol-r-3.mol`, `Pharm3D_mol-s-3.mol`: the files
  of the same names in `rdkit/Chem/Pharm3D/test_data`, with trailing whitespace
  and blank lines trimmed
- `Subshape_5ht3ligs.sdf`: `5ht3ligs.sdf` in `rdkit/Chem/Subshape/test_data`,
  with trailing whitespace trimmed
- `RP-6306_atrop1.sdf`: `Code/GraphMol/FileParsers/test_data/atropisomers`,
  with trailing whitespace trimmed

Made from RDKit's source tree with its Python package (RDKit 2026.09.1):

- `1DWD_ligand.mol`, `1PPC_ligand.mol`, `1DWD_ligand_reordered.mol`: the PDB
  files of the same names in `rdkit/Chem/test_data`, with bond orders assigned
  from `UnitTestMol3D.py`'s template by `AllChem.AssignBondOrdersFromTemplate`,
  written by `Chem.MolToMolFile`
- `feat_dirs.mol`: an ETKDG conformer (seed 42) of
  `O=C(N)c1ccc(cc1)N(C)CCOC(=O)CS(=O)(=O)N1CCOCC1` with hydrogens, from
  `AllChem.EmbedMolecule`, written by `Chem.MolToMolFile`
- `cdk2-syn-clip100.smi`: the molecules (with their hydrogens) and names
  pickled in `rdkit/Chem/Pharm3D/test_data/cdk2-syn-clip100.pkl.gz`, written by
  `Chem.MolToSmiles`; each SMILES reads back to the same SMILES and atom count
