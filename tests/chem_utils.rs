// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The ChemUtils scripts' functions; values from Python's modules in RDKit
// 2024.09.6, run on the same input.

use rdkit_cxx::chem::chem_utils::align_depict::*;
use rdkit_cxx::chem::chem_utils::sdf_to_csv::*;
use rdkit_cxx::chem::chem_utils::template_expand::*;
use rdkit_cxx::chem::rd_depictor::compute_2d_coords;
use rdkit_cxx::chem::{
    Mol, SDMolSupplier, mol_from_mol_block, mol_from_smarts, mol_from_smiles, mol_to_smiles,
};
use rdkit_cxx::{Error, Result};

const SALTS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/witch-salts.sdf");

const CORE: &str = r#"
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END
"#;

const ALIGNED: &str = r#"
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  7  2  1  0
M  END
"#;

const ALIGNED_PATT: &str = r#"
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  7  2  1  0
M  END
"#;

const ACCEPT: &str = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"#;

const CSV: &str = r#"SMILES,Code,Name
Cl,AA,Hydrochloride
[Na+],AB,Sodium
Br,AC,Hydrobromide
I,AD,Hydroiodide
"#;

const CSV_KEY: &str = r#"Code,SMILES,Name
AA,Cl,Hydrochloride
AB,[Na+],Sodium
AC,Br,Hydrobromide
"#;

const SDF: &str = r#"TemplateEnum: Mol_1
     RDKit          2D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  7  1  1  0
  9 10  1  0
 10 11  2  0
  9  4  1  0
M  END
>  <seq_num>
1

>  <reagent_indices>
1_1

>  <building_block_1>
a0

>  <building_block_1_src>
x

>  <building_block_2>
1

$$$$
TemplateEnum: Mol_2
     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  7  1  1  0
  9 10  1  0
  9  4  1  0
M  END
>  <seq_num>
2

>  <reagent_indices>
1_2

>  <building_block_1>
a0

>  <building_block_1_src>
x

>  <building_block_2>
2

$$$$
TemplateEnum: Mol_3
     RDKit          2D

 12 12  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2500   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  8  9  1  0
  7  1  1  0
 10 11  1  0
 11 12  2  0
 10  4  1  0
M  END
>  <seq_num>
3

>  <reagent_indices>
2_1

>  <building_block_1>
a1

>  <building_block_1_src>
x

>  <building_block_2>
1

$$$$
TemplateEnum: Mol_4
     RDKit          2D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2500   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  8  9  1  0
  7  1  1  0
 10 11  1  0
 10  4  1  0
M  END
>  <seq_num>
4

>  <reagent_indices>
2_2

>  <building_block_1>
a1

>  <building_block_1_src>
x

>  <building_block_2>
2

$$$$
TemplateEnum: Mol_5
     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  1  1  0
  8  9  1  0
  9 10  2  0
  8  4  1  0
M  END
>  <seq_num>
5

>  <reagent_indices>
3_1

>  <building_block_1>
a2

>  <building_block_1_src>
x

>  <building_block_2>
1

$$$$
TemplateEnum: Mol_6
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  1  1  0
  8  9  1  0
  8  4  1  0
M  END
>  <seq_num>
6

>  <reagent_indices>
3_2

>  <building_block_1>
a2

>  <building_block_1_src>
x

>  <building_block_2>
2

$$$$
"#;

const SDF_NAMES: &str = r#"tmpl_a0_1
     RDKit          2D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  7  1  1  0
  9 10  1  0
 10 11  2  0
  9  4  1  0
M  END
>  <seq_num>
1

>  <reagent_indices>
1_1

>  <building_block_1>
a0

>  <building_block_1_src>
x

>  <building_block_2>
1

$$$$
tmpl_a0_2
     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  7  1  1  0
  9 10  1  0
  9  4  1  0
M  END
>  <seq_num>
2

>  <reagent_indices>
1_2

>  <building_block_1>
a0

>  <building_block_1_src>
x

>  <building_block_2>
2

$$$$
tmpl_a1_1
     RDKit          2D

 12 12  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2500   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  8  9  1  0
  7  1  1  0
 10 11  1  0
 11 12  2  0
 10  4  1  0
M  END
>  <seq_num>
3

>  <reagent_indices>
2_1

>  <building_block_1>
a1

>  <building_block_1_src>
x

>  <building_block_2>
1

$$$$
tmpl_a1_2
     RDKit          2D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2500   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  8  1  0
  8  9  1  0
  7  1  1  0
 10 11  1  0
 10  4  1  0
M  END
>  <seq_num>
4

>  <reagent_indices>
2_2

>  <building_block_1>
a1

>  <building_block_1_src>
x

>  <building_block_2>
2

$$$$
tmpl_a2_1
     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  1  1  0
  8  9  1  0
  9 10  2  0
  8  4  1  0
M  END
>  <seq_num>
5

>  <reagent_indices>
3_1

>  <building_block_1>
a2

>  <building_block_1_src>
x

>  <building_block_2>
1

$$$$
tmpl_a2_2
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
    1.5000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  7  1  1  0
  8  9  1  0
  8  4  1  0
M  END
>  <seq_num>
6

>  <reagent_indices>
3_2

>  <building_block_1>
a2

>  <building_block_1_src>
x

>  <building_block_2>
2

$$$$
"#;

/// Whether the molecule's depiction is the block's, within the block's
/// rounding: the core here is read back from a rounded block.
fn same_depiction(mol: &Mol, block: &str) {
    let expected = mol_from_mol_block(block).unwrap();
    assert_eq!(
        mol_to_smiles(mol).unwrap(),
        mol_to_smiles(&expected).unwrap()
    );
    let (a, b) = (
        mol.get_conformer().unwrap(),
        expected.get_conformer().unwrap(),
    );
    for i in 0..mol.get_num_atoms() {
        let d = a
            .get_atom_position(i)
            .unwrap()
            .distance(&b.get_atom_position(i).unwrap());
        assert!(d < 1e-3, "atom {i} is {d} away");
    }
}

#[test]
fn aligned_depictions() -> Result<()> {
    let core = mol_from_mol_block(CORE)?;
    let mut mol = mol_from_smiles("Oc1ccccc1C")?;
    align_depict(&mut mol, &core)?;
    same_depiction(&mol, ALIGNED);
    let mut mol = mol_from_smiles("Oc1ccccc1C")?;
    let params = AlignDepictParams {
        core_pattern: Some(mol_from_smarts("c1ccccc1")?),
        ..Default::default()
    };
    align_depict_with_params(&mut mol, &core, &params)?;
    same_depiction(&mol, ALIGNED_PATT);
    let mut mol = mol_from_smiles("CCO")?;
    let params = AlignDepictParams {
        accept_failure: true,
        ..Default::default()
    };
    align_depict_with_params(&mut mol, &core, &params)?;
    same_depiction(&mol, ACCEPT);
    match align_depict(&mut mol_from_smiles("CCO")?, &core) {
        Err(Error::Value { message }) => {
            assert_eq!(message, "Substructure match with core not found.")
        }
        other => panic!("{other:?}"),
    }
    let params = AlignDepictParams {
        core_pattern: Some(mol_from_smarts("cc")?),
        ..Default::default()
    };
    assert!(matches!(
        align_depict_with_params(&mut mol_from_smiles("c1ccccc1")?, &core, &params),
        Err(Error::Value { .. })
    ));
    Ok(())
}

#[test]
fn csv() -> Result<()> {
    let mut out = Vec::new();
    let params = ConvertParams {
        stop_after: 4,
        ..Default::default()
    };
    convert_with_params(&mut SDMolSupplier::new(SALTS)?, &mut out, &params)?;
    assert_eq!(String::from_utf8(out).unwrap(), CSV.replace('\n', "\r\n"));
    let mut out = Vec::new();
    let params = ConvertParams {
        key_col: Some("Code".into()),
        stop_after: 3,
        include_chirality: true,
        ..Default::default()
    };
    convert_with_params(&mut SDMolSupplier::new(SALTS)?, &mut out, &params)?;
    assert_eq!(
        String::from_utf8(out).unwrap(),
        CSV_KEY.replace('\n', "\r\n")
    );
    let mut out = Vec::new();
    convert(&mut SDMolSupplier::new(SALTS)?, &mut out)?;
    assert_eq!(String::from_utf8(out).unwrap().matches("\r\n").count(), 241);
    Ok(())
}

fn named(smiles: &[&str]) -> Vec<Option<Mol>> {
    smiles
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut m = mol_from_smiles(s).unwrap();
            m.set_prop("_Name", &format!("a{i}"));
            m.set_prop("src", "x");
            Some(m)
        })
        .collect()
}

fn chains(sc: &[Sidechain]) -> Vec<(usize, String)> {
    sc.iter()
        .map(|(i, m)| (*i, mol_to_smiles(m).unwrap()))
        .collect()
}

#[test]
fn template_expansion() -> Result<()> {
    let mut tmpl = mol_from_smiles("[1*]c1ccc([2*])cc1")?;
    tmpl.set_prop("_Name", "tmpl");
    compute_2d_coords(&mut tmpl)?;
    let c1 = construct_sidechains(named(&["CO", "CCN", "Cl"]))?.unwrap();
    let s2 = ["C(=O)O", "OC"]
        .iter()
        .map(|s| mol_from_smiles(s).ok())
        .collect::<Vec<_>>();
    let params = ConstructSidechainsParams {
        sma: Some("[OH]".into()),
        replace: false,
        ..Default::default()
    };
    let c2 = construct_sidechains_with_params(s2, &params)?.unwrap();
    assert_eq!(
        chains(&c2),
        [(1, "O=CO"), (2, "CO")].map(|(i, s): (usize, &str)| (i, s.to_owned()))
    );
    let s3 = ["ClCCO", "BrCC"]
        .iter()
        .map(|s| mol_from_smiles(s).ok())
        .collect::<Vec<_>>();
    let params = ConstructSidechainsParams {
        sma: Some("[Cl,Br]".into()),
        ..Default::default()
    };
    assert_eq!(
        chains(&construct_sidechains_with_params(s3, &params)?.unwrap()),
        [(1, "CCO"), (2, "CC")].map(|(i, s): (usize, &str)| (i, s.to_owned()))
    );
    let bad = ConstructSidechainsParams {
        sma: Some("[C".into()),
        ..Default::default()
    };
    assert!(construct_sidechains_with_params(named(&["C"]), &bad)?.is_none());
    let moved: Vec<String> =
        move_dummy_neighbors_to_beginning_with_use_all(&mol_from_smiles("OCC(*)C*")?, true)?
            .iter()
            .map(|m| mol_to_smiles(m).unwrap())
            .collect();
    assert_eq!(moved, ["CCCO", "CCCO"]);
    let sidechains = [c1, c2];
    let mut out = Vec::new();
    explode(&tmpl, &sidechains, &mut out)?;
    assert_eq!(String::from_utf8(out).unwrap(), SDF);
    let mut out = Vec::new();
    explode_with_params(
        &tmpl,
        &sidechains,
        &mut out,
        &ExplodeParams {
            auto_names: false,
            ..Default::default()
        },
    )?;
    assert_eq!(String::from_utf8(out).unwrap(), SDF_NAMES);
    Ok(())
}
