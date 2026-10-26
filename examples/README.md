# Examples

The [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)'s recipes,
ported line by line, all 40 of them. Each prints what the recipe prints and
asserts the cookbook's output, so a run checks it:

```sh
cargo run --example cookbook_hybridization
```

The drawing recipes write their images to `target/cookbook/`. Where a recipe
shows a molecule in Jupyter, the example draws it as Jupyter does
(`IPythonConsole`'s `molSize` and `drawOptions`), as an SVG, and checks what
the drawing holds: the indices, notes, CIP labels, highlights, and colors
the recipe adds. Where the recipe draws with `MolDraw2DCairo`, the example
writes the PNG and checks its size.

Where a recipe calls code RDKit offers only in Python, the example calls its port in
the crate: `rdkit_cxx::chem::find_mol_chiral_centers`,
`rdkit_cxx::chem::lipinski`, and `rdkit_cxx::chem::draw`'s `mol_to_image`,
`mol_to_svg`, `mols_to_grid_image`, and `draw_mol_with_matches`.

| Recipe | ID | Example |
| --- | --- | --- |
| Include an Atom Index | RDKitCB_0 | [`cookbook_atom_index`](cookbook_atom_index.rs) |
| Include a Bond Index | RDKitCB_40 | [`cookbook_bond_index`](cookbook_bond_index.rs) |
| Include a Calculation | RDKitCB_23 | [`cookbook_calculation`](cookbook_calculation.rs) |
| Include Stereo Annotations | RDKitCB_32 | [`cookbook_stereo_annotations`](cookbook_stereo_annotations.rs) |
| Black and White Molecules | RDKitCB_1 | [`cookbook_black_and_white`](cookbook_black_and_white.rs) |
| Highlight a Substructure in a Molecule | RDKitCB_2 | [`cookbook_highlight_substructure`](cookbook_highlight_substructure.rs) |
| Highlight Molecule Differences | RDKitCB_36 | [`cookbook_highlight_differences`](cookbook_highlight_differences.rs) |
| Highlight Entire Molecule | RDKitCB_38 | [`cookbook_highlight_molecule`](cookbook_highlight_molecule.rs) |
| Highlight Molecule with Multiple Colors | RDKitCB_39 | [`cookbook_highlight_multiple_colors`](cookbook_highlight_multiple_colors.rs) |
| Without Implicit Hydrogens | RDKitCB_17 | [`cookbook_without_implicit_hs`](cookbook_without_implicit_hs.rs) |
| With Abbreviations | RDKitCB_34 | [`cookbook_abbreviations`](cookbook_abbreviations.rs) |
| Using CoordGen Library | RDKitCB_37 | [`cookbook_coordgen`](cookbook_coordgen.rs) |
| On a Plot | RDKitCB_35 | [`cookbook_on_a_plot`](cookbook_on_a_plot.rs) |
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

Some recipes differ from the cookbook's text on purpose:

- Neutralizing Molecules: the cookbook's `Uncharger` comparison prints its
  inputs, because `uncharge` returns a new molecule rather than changing its
  argument. The example checks what `uncharge` returns.
- Create Fragments: Python's `with Chem.RWMol(mol)` batches the bond
  removals, so the example brackets them with `begin_batch_edit` and
  `commit_batch_edit`.
- Highlight a Substructure in a Molecule: after a substructure search,
  Jupyter keeps the query, and when the recipe then sets the atoms to
  highlight by hand, drawing them fails in RDKit 2026.09, since they do not
  fit the query. The example shows the failure, then draws the atoms without
  the query, as the recipe intends.
- On a Plot: the recipe puts the molecule on a matplotlib plot, which has no
  counterpart here; the example draws the image.
