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
- `RP-6306_atrop1.sdf`: `Code/GraphMol/FileParsers/test_data/atropisomers`,
  with trailing whitespace trimmed

Made from RDKit's source tree with its Python package (RDKit 2026.09.1):

- `1DWD_ligand.mol`, `1PPC_ligand.mol`, `1DWD_ligand_reordered.mol`: the PDB
  files of the same names in `rdkit/Chem/test_data`, with bond orders assigned
  from `UnitTestMol3D.py`'s template by `AllChem.AssignBondOrdersFromTemplate`,
  written by `Chem.MolToMolFile`
