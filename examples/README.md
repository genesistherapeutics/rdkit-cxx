# Examples

The [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)'s recipes,
ported line by line. Each prints what the recipe prints and asserts the
cookbook's output, so a run checks it:

```sh
cargo run --example cookbook_hybridization
```

The drawing recipes are left out, since the crate does not bind
`rdMolDraw2D`; recipes that draw their results along the way print SMILES
instead. Where a recipe calls a pure-Python RDKit function, the example
calls its port in the crate: `rdkit_cxx::chem::find_mol_chiral_centers` and
`rdkit_cxx::chem::lipinski`.

| Recipe | ID | Example |
| --- | --- | --- |
| Hybridization Type and Count | RDKitCB_26 | [`cookbook_hybridization`](cookbook_hybridization.rs) |
| Count Ring Systems | RDKitCB_3 | [`cookbook_ring_systems`](cookbook_ring_systems.rs) |
| Identify Aromatic Rings | RDKitCB_8 | [`cookbook_aromatic_rings`](cookbook_aromatic_rings.rs) |
| Identify Aromatic Atoms | RDKitCB_9 | [`cookbook_aromatic_atoms`](cookbook_aromatic_atoms.rs) |
| Identifying Stereochemistry | RDKitCB_30 | [`cookbook_stereochemistry`](cookbook_stereochemistry.rs) |
| Create Fragments | RDKitCB_7 | [`cookbook_fragments`](cookbook_fragments.rs) |
| Largest Fragment | RDKitCB_31 | [`cookbook_largest_fragment`](cookbook_largest_fragment.rs) |
| Sidechain-Core Enumeration | RDKitCB_29 | [`cookbook_sidechain_core_enumeration`](cookbook_sidechain_core_enumeration.rs) |
| Neutralizing Molecules | RDKitCB_33 | [`cookbook_neutralize`](cookbook_neutralize.rs) |
| Functional Group with SMARTS queries | RDKitCB_10 | [`cookbook_functional_groups`](cookbook_functional_groups.rs) |
| Macrocycles with SMARTS queries | RDKitCB_13 | [`cookbook_macrocycles`](cookbook_macrocycles.rs) |
| Returning Substructure Matches as SMILES | RDKitCB_18 | [`cookbook_substructure_smiles`](cookbook_substructure_smiles.rs) |
| Within the Same Fragment | RDKitCB_20 | [`cookbook_same_fragment_match`](cookbook_same_fragment_match.rs) |
| Molecule Hash Strings | RDKitCB_21 | [`cookbook_mol_hash`](cookbook_mol_hash.rs) |
| Contiguous Rotatable Bonds | RDKitCB_22 | [`cookbook_contiguous_rotatable_bonds`](cookbook_contiguous_rotatable_bonds.rs) |
| Kekule SMILES | RDKitCB_4 | [`cookbook_kekule_smiles`](cookbook_kekule_smiles.rs) |
| Isomeric SMILES without isotopes | RDKitCB_5 | [`cookbook_smiles_without_isotopes`](cookbook_smiles_without_isotopes.rs) |
| Reversing Reactions | RDKitCB_6 | [`cookbook_reverse_reaction`](cookbook_reverse_reaction.rs) |
| Reaction Fingerprints and Similarity | RDKitCB_27 | [`cookbook_reaction_fingerprints`](cookbook_reaction_fingerprints.rs) |
| Explicit Valence Error - Partial Sanitization | RDKitCB_15 | [`cookbook_partial_sanitization`](cookbook_partial_sanitization.rs) |
| Detect Chemistry Problems | RDKitCB_14 | [`cookbook_chemistry_problems`](cookbook_chemistry_problems.rs) |
| Explicit Valence and Number of Hydrogens | RDKitCB_11 | [`cookbook_explicit_valence`](cookbook_explicit_valence.rs) |
| Wiener Index | RDKitCB_12 | [`cookbook_wiener_index`](cookbook_wiener_index.rs) |
| Organometallics with Dative Bonds | RDKitCB_19 | [`cookbook_dative_bonds`](cookbook_dative_bonds.rs) |
| Enumerate SMILES | RDKitCB_24 | [`cookbook_enumerate_smiles`](cookbook_enumerate_smiles.rs) |
| Reorder Atoms | RDKitCB_28 | [`cookbook_reorder_atoms`](cookbook_reorder_atoms.rs) |
| Conformer Generation with ETKDG | RDKitCB_25 | [`cookbook_etkdg`](cookbook_etkdg.rs) |

Two recipes differ from the cookbook's text on purpose:

- Neutralizing Molecules: the cookbook's `Uncharger` comparison prints its
  inputs, because `uncharge` returns a new molecule rather than changing its
  argument. The example checks what `uncharge` returns.
- Create Fragments: Python's `with Chem.RWMol(mol)` batches the bond
  removals, so the example brackets them with `begin_batch_edit` and
  `commit_batch_edit`.
