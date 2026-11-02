// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::{mol_from_mol_block, mol_to_smiles};

// Expected values come from Python's rdMolStandardize.

const BENZOATE: &str = "
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    2.5833    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8333    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5833   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3333   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4167    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9167    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6667    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9167   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4167   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.2990    0.0000 Na  0  0  0  0  0 15  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9  4  1  0
M  CHG  2   3  -1  10   1
M  END
";

const NITRO: &str = "
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  2  0
M  END
";

const SODIUM_ACETATE: &str = "
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -1.8187    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5196    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5196   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7794    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0785   -0.0000    0.0000 Na  0  0  0  0  0  1  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  4  5  1  0
M  END
";

const SULFOXIDE: &str = "
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
M  CHG  2   2   1   4  -1
M  END
";

const EMPTY: &str = "
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
";

fn smiles_of(block: &str) -> String {
    mol_to_smiles(&mol_from_mol_block(block).unwrap()).unwrap()
}

#[test]
fn options_defaults() {
    let options = PipelineOptions::default();
    assert!(!options.strict_parsing);
    assert!(options.report_all_failures);
    assert!(!options.allow_empty_molecules);
    assert!(!options.allow_enhanced_stereo);
    assert!(!options.allow_aromatic_bond_type);
    assert!(!options.allow_dative_bond_type);
    assert!((options.is_2d_zero_threshold - 1e-3).abs() < 1e-12);
    assert!((options.atom_clash_limit - 0.03).abs() < 1e-12);
    assert!((options.min_median_bond_length - 1e-3).abs() < 1e-12);
    assert!((options.bond_length_limit - 100.0).abs() < 1e-12);
    assert!(options.allow_long_bonds_in_rings);
    assert!(options.allow_atom_bond_clash_exemption);
    assert_eq!(options.metal_nof, "[Li,Na,K,Rb,Cs,Fr]~[#7,#8,F]");
    assert_eq!(options.metal_non, "");
    assert_eq!(options.normalizer_data.len(), 2560);
    assert!(options.normalizer_data.starts_with("// Name\tSMIRKS\n"));
    assert_eq!(options.normalizer_max_restarts, 200);
    assert!((options.scaled_median_bond_length - 1.0).abs() < 1e-12);
    assert!(!options.output_v2000);
}

#[test]
fn status_values() {
    assert_eq!(PipelineStatus::NO_EVENT.bits(), 0);
    assert_eq!(PipelineStatus::PIPELINE_ERROR.bits(), 8191);
    assert_eq!(PipelineStatus::STRUCTURE_MODIFICATION.bits(), 125829120);
    assert_eq!(PipelineStatus::VALIDATION_ERROR.bits(), 124);
    assert_eq!(PipelineStatus::STANDARDIZATION_ERROR.bits(), 3840);
    assert_eq!(
        PipelineStatus::FRAGMENTS_REMOVED | PipelineStatus::PROTONATION_CHANGED,
        PipelineStatus::from_bits_retain(100663296)
    );
    assert_eq!(
        format!("{:?}", PipelineStatus::OUTPUT_ERROR),
        "PipelineStatus::OUTPUT_ERROR"
    );
}

#[test]
fn running_the_pipeline() {
    let pipeline = Pipeline::default();
    let result = pipeline.run(BENZOATE).unwrap();
    assert_eq!(result.status.bits(), 100663296);
    assert_eq!(result.stage, PipelineStage::Completed);
    let log: Vec<_> = result
        .log
        .iter()
        .map(|entry| (entry.status, entry.detail.as_str()))
        .collect();
    assert_eq!(
        log,
        [
            (
                PipelineStatus::FRAGMENTS_REMOVED,
                "One or more disconnected fragments (e.g., counterions) were removed."
            ),
            (
                PipelineStatus::PROTONATION_CHANGED,
                "The protonation state was adjusted."
            ),
        ]
    );
    assert_eq!(result.input_mol_data, BENZOATE);
    assert!(result.output_mol_data.contains("V3000"));
    assert_eq!(smiles_of(&result.output_mol_data), "O=C(O)c1ccccc1");
    assert_eq!(smiles_of(&result.parent_mol_data), "O=C(O)c1ccccc1");

    let v2000 = Pipeline::new(&PipelineOptions {
        output_v2000: true,
        ..Default::default()
    })
    .unwrap();
    let result = v2000.run(BENZOATE).unwrap();
    assert!(result.output_mol_data.contains("V2000"));
    assert_eq!(smiles_of(&result.output_mol_data), "O=C(O)c1ccccc1");

    let result = pipeline.run(NITRO).unwrap();
    assert_eq!(result.status, PipelineStatus::NO_EVENT);
    assert!(result.log.is_empty());
    assert_eq!(smiles_of(&result.output_mol_data), "C[N+](=O)[O-]");

    let result = pipeline.run("not a mol block").unwrap();
    assert_eq!(result.status, PipelineStatus::INPUT_ERROR);
    assert_eq!(result.stage, PipelineStage::ParsingInput);
    assert_eq!(
        result.log[0].detail,
        "Could not instantiate a valid molecule from input"
    );
    assert!(result.output_mol_data.is_empty());

    let result = pipeline.run(EMPTY).unwrap();
    assert_eq!(result.status, PipelineStatus::BASIC_VALIDATION_ERROR);
    assert_eq!(result.stage, PipelineStage::Completed);
    assert_eq!(result.log[0].detail, "Molecule has no atoms");
    let stop = Pipeline::new(&PipelineOptions {
        report_all_failures: false,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(stop.run(EMPTY).unwrap().stage, PipelineStage::Validation);
    let allow = Pipeline::new(&PipelineOptions {
        allow_empty_molecules: true,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(allow.run(EMPTY).unwrap().status, PipelineStatus::NO_EVENT);
}

#[test]
fn standardization_steps() {
    // The wrapper runs its own copy of RDKit's standardization step, which
    // casts the normalizer's ROMol to an RWMol it is not.
    let statuses = |result: &PipelineResult| {
        result
            .log
            .iter()
            .map(|entry| entry.status.bits())
            .collect::<Vec<_>>()
    };
    let result = Pipeline::default().run(SODIUM_ACETATE).unwrap();
    assert_eq!(result.status.bits(), 109051904);
    assert_eq!(result.stage, PipelineStage::Completed);
    assert_eq!(statuses(&result), [8388608, 33554432, 67108864]);
    assert_eq!(
        result.log[0].detail,
        "One or more metal atoms were disconnected."
    );
    assert_eq!(smiles_of(&result.output_mol_data), "CC(=O)O");
    assert_eq!(smiles_of(&result.parent_mol_data), "CC(=O)O");

    let result = Pipeline::default().run(SULFOXIDE).unwrap();
    assert_eq!(result.status, PipelineStatus::NO_EVENT);
    assert_eq!(smiles_of(&result.output_mol_data), "C[S+](C)[O-]");
    let sulfoxide = Pipeline::new(&PipelineOptions {
        normalizer_data: "// Name\tSMIRKS\nSulfoxide to S=O\t[S+:1][O-:2]>>[S+0:1]=[O-0:2]\n"
            .to_string(),
        ..Default::default()
    })
    .unwrap();
    let result = sulfoxide.run(SULFOXIDE).unwrap();
    assert_eq!(result.status, PipelineStatus::NORMALIZATION_APPLIED);
    assert_eq!(
        result.log[0].detail,
        "The representation of some functional groups was adjusted."
    );
    assert_eq!(smiles_of(&result.output_mol_data), "CS(C)=O");

    let bad = Pipeline::new(&PipelineOptions {
        normalizer_data: "// Name\tSMIRKS\nbad\tC>>\n".to_string(),
        ..Default::default()
    })
    .unwrap();
    let result = bad.run(SODIUM_ACETATE).unwrap();
    assert_eq!(result.status.bits(), 8389120);
    assert_eq!(result.stage, PipelineStage::Completed);
    assert_eq!(statuses(&result), [8388608, 512]);
    assert_eq!(
        result.log[1].detail,
        "An error occurred while normalizing the representation of some functional groups"
    );
}

#[test]
fn pipelines_on_threads() {
    let pipeline = Pipeline::default();
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for _ in 0..5 {
                    let result = pipeline.run(BENZOATE).unwrap();
                    assert_eq!(result.stage, PipelineStage::Completed);
                }
            });
        }
    });
}

#[test]
fn metal_smarts_are_checked() {
    // RDKit dereferences the null that SMARTS that do not parse give, and
    // reads past matches of queries that are not two bonded atoms.
    for (nof, non) in [
        ("[Na", ""),
        ("[Na]", ""),
        ("[Na]~[O]~[C]", ""),
        ("[Li,Na,K,Rb,Cs,Fr]~[#7,#8,F]", "[Hg].[C]"),
    ] {
        let options = PipelineOptions {
            metal_nof: nof.to_string(),
            metal_non: non.to_string(),
            ..Default::default()
        };
        assert!(
            matches!(Pipeline::new(&options), Err(Error::Value { .. })),
            "{nof} {non}"
        );
    }
    let options = PipelineOptions {
        metal_non: "[Hg]~[C]".to_string(),
        ..Default::default()
    };
    assert!(Pipeline::new(&options).is_ok());
}

#[test]
fn long_stereo_group_lines_are_refused() {
    // The mol-block parser matches enhanced-stereo lines with a regex that
    // recurses once per character (U12).
    let block = |padding: usize| {
        format!(
            "\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n\
             M  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\n\
             M  V30 2 O 1.5 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\n\
             M  V30 END BOND\nM  V30 BEGIN COLLECTION\nM  V30 MDLV30/STERAC1 ATOMS=(1 1{})\n\
             M  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n",
            " ".repeat(padding)
        )
    };
    let pipeline = Pipeline::default();
    assert!(pipeline.run(&block(0)).is_ok());
    assert!(matches!(
        pipeline.run(&block(10_000)),
        Err(Error::Value { .. })
    ));
}
