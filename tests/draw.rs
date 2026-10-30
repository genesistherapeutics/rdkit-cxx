// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The ported rdkit.Chem.Draw functions against Python's in RDKit 2026.03.6,
// whose SVGs they reproduce byte for byte (by SHA-1). PNGs come from Cairo
// and zlib, whose bytes vary by build, so only their sizes are checked.

use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{DrawColour, MolDraw2DSVG, MolDrawOptions};
use rdkit_cxx::chem::draw::*;
use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
use rdkit_cxx::chem::rd_fingerprint_generator::{
    AdditionalOutput, GetFingerprintParams, GetMorganGeneratorParams,
    get_morgan_generator_with_params,
};
use rdkit_cxx::chem::{
    Mol, RDKFingerprintParams, mol_from_png_string, mol_from_smarts, mol_from_smiles,
    rdk_fingerprint_with_outputs,
};

fn sha1(text: &str) -> String {
    let mut h = sha1_smol::Sha1::new();
    h.update(text.as_bytes());
    h.digest().to_string()
}

/// A PNG's width and height, from its IHDR chunk.
fn png_size(png: &[u8]) -> (u32, u32) {
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    let word = |at: usize| u32::from_be_bytes(png[at..at + 4].try_into().unwrap());
    (word(16), word(20))
}

fn aspirin() -> Mol {
    mol_from_smiles("c1cc(C(=O)O)c(OC(=O)C)cc1").unwrap()
}

fn atom_indices() -> MolDrawOptions {
    MolDrawOptions {
        add_atom_indices: true,
        ..Default::default()
    }
}

fn grid_mols() -> Vec<Mol> {
    ["CCO", "c1ccccc1O", "C[C@H](N)C(=O)O", "CC(=O)Nc1ccc(O)cc1"]
        .iter()
        .map(|smiles| mol_from_smiles(smiles).unwrap())
        .collect()
}

#[test]
fn mol_to_svg_matches_python() {
    let m = aspirin();
    assert_eq!(
        sha1(&mol_to_svg(&m).unwrap()),
        "91e470f6aa690aeb25fe9ae4b9ad737d50e9cd22"
    );
    let highlighted = MolToSVGParams {
        size: (250, 200),
        highlight_atoms: vec![3, 4, 5],
        legend: "aspirin".to_owned(),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mol_to_svg_with_params(&m, &highlighted).unwrap()),
        "552df835bbac8dd19b158b5d3e7a7eae2250b41c"
    );
    let with_options = MolToSVGParams {
        draw_options: Some(atom_indices()),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mol_to_svg_with_params(&m, &with_options).unwrap()),
        "cd33dc765de68029fa82cafb3836e1ec2d07cdc3"
    );
}

#[test]
fn mols_to_grid_image_matches_python() {
    let ms = grid_mols();
    assert_eq!(
        sha1(&mols_to_grid_image_use_svg(&ms).unwrap()),
        "545d04dce8891ad8de975efb5a6a30e624942157"
    );
    // A legend too many is cut, as in Python.
    let params = MolsToGridImageParams {
        mols_per_row: 2,
        sub_img_size: (150, 150),
        legends: Some(["a", "b", "c", "d", "e"].map(str::to_owned).to_vec()),
        highlight_atom_lists: Some(vec![vec![0], vec![6], vec![1, 2], vec![]]),
        draw_options: Some(atom_indices()),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mols_to_grid_image_use_svg_with_params(&ms, &params).unwrap()),
        "473bd83698759e9613f19f200889447fab63b833"
    );
    assert_eq!(png_size(&mols_to_grid_image(&ms).unwrap()), (600, 400));
    let no_rows = MolsToGridImageParams {
        mols_per_row: 0,
        ..Default::default()
    };
    assert!(mols_to_grid_image_use_svg_with_params(&ms, &no_rows).is_err());
}

#[test]
fn draw_mol_with_matches_matches_python() {
    let m = aspirin();
    let q = mol_from_smarts("C(=O)O").unwrap();
    let matches = m.get_substruct_matches(&q).unwrap();
    assert_eq!(matches, [[3, 4, 5], [8, 9, 7]]);
    assert_eq!(
        sha1(&draw_mol_with_matches_as_svg(&m, &matches).unwrap()),
        "8e90eb551fce647f35ebed7bd00bd49f39a70bf4"
    );
    let params = DrawMolWithMatchesParams {
        qry: Some(q),
        colors: Some(vec![
            DrawColour::new(1.0, 0.0, 0.0),
            DrawColour::new(0.0, 0.0, 1.0),
        ]),
        ..Default::default()
    };
    assert_eq!(
        sha1(&draw_mol_with_matches_as_svg_with_params(&m, &matches, &params).unwrap()),
        "eb6bac7a96db96112c1c3837312afc669552d322"
    );
    let no_colors = DrawMolWithMatchesParams {
        colors: Some(Vec::new()),
        ..Default::default()
    };
    assert!(draw_mol_with_matches_as_svg_with_params(&m, &matches, &no_colors).is_err());
    assert!(draw_mol_with_matches_as_svg_with_params(&m, &[], &no_colors).is_ok());
}

#[test]
fn mol_to_image_is_a_png_of_its_size() {
    let params = MolToImageParams {
        size: (250, 200),
        ..Default::default()
    };
    assert_eq!(
        png_size(&mol_to_image_with_params(&aspirin(), &params).unwrap()),
        (250, 200)
    );
}

fn phenylcyclopropane() -> Mol {
    mol_from_smiles("c1ccccc1CC1CC1").unwrap()
}

fn morgan_bit_info(m: &Mol) -> MorganBitInfo {
    let mut ao = AdditionalOutput::default();
    ao.allocate_bit_info_map();
    let params = GetFingerprintParams {
        additional_output: Some(ao),
        ..Default::default()
    };
    let fpg = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 2,
        ..Default::default()
    })
    .unwrap();
    let (_, outputs) = fpg.get_fingerprint_with_outputs(m, &params).unwrap();
    outputs.get_bit_info_map().unwrap().clone()
}

fn rdkit_bit_info(m: &Mol) -> RDKitBitInfo {
    let params = RDKFingerprintParams {
        max_path: 5,
        ..Default::default()
    };
    rdk_fingerprint_with_outputs(m, &params).unwrap().1.bit_info
}

// UnitTestDraw.testDrawMorgan, with a Morgan generator's bit information.
#[test]
fn draw_morgan_matches_python() {
    let m = phenylcyclopropane();
    let bi = morgan_bit_info(&m);
    assert_eq!(bi[&872], [(6, 2)]);
    let svg1 = draw_morgan_bit(&m, 872, &bi).unwrap();
    assert_eq!(svg1, draw_morgan_env(&m, 6, 2).unwrap());
    assert_eq!(sha1(&svg1), "4c91713afb5df53a558a983fab1336f2d739318b");
    for colour in ["#CCCCCC", "#E5E533", "#9999E5"] {
        assert!(svg1.contains(&format!("style='fill:{colour};")));
    }
    let no_center = DrawMorganBitParams {
        env: DrawMorganEnvParams {
            center_color: None,
            ..Default::default()
        },
        ..Default::default()
    };
    let svg1 = draw_morgan_bit_with_params(&m, 872, &bi, &no_center).unwrap();
    assert_eq!(
        svg1,
        draw_morgan_env_with_params(&m, 6, 2, &no_center.env).unwrap()
    );
    assert_eq!(sha1(&svg1), "7b55f4c2546f49ec5670b7d87e6c19abf5b2cd93");
    assert!(!svg1.contains("style='fill:#9999E5;"));
    assert!(matches!(
        draw_morgan_bit(&m, 32, &bi),
        Err(Error::Key { .. })
    ));
    let second = DrawMorganBitParams {
        which_example: 1,
        ..Default::default()
    };
    assert!(matches!(
        draw_morgan_bit_with_params(&m, 872, &bi, &second),
        Err(Error::Index { .. })
    ));
    assert_eq!(
        sha1(&draw_morgan_env(&m, 0, 0).unwrap()),
        "369e0a9ccff3c67d97e7e757a2fe0d5b6eb53a78"
    );
    // An environment of no bonds is drawn from its SMILES.
    let salt = mol_from_smiles("[Na+].[Cl-]").unwrap();
    assert_eq!(
        sha1(&draw_morgan_env(&salt, 0, 1).unwrap()),
        "618e0a7d73e3475298d10d61189b3305e984055e"
    );
    let envs_params = DrawMorganEnvsParams {
        mols_per_row: 2,
        legends: Some(["a", "b", "c"].map(str::to_owned).to_vec()),
        ..Default::default()
    };
    assert_eq!(
        sha1(
            &draw_morgan_envs_with_params(&[(&m, 0, 0), (&m, 6, 2), (&m, 7, 1)], &envs_params)
                .unwrap()
        ),
        "ff74b0bcc2d1240207e1b9d836bee170369c66cf"
    );
    assert_eq!(
        sha1(&draw_morgan_bits(&[(&m, 872, &bi, 0), (&m, 872, &bi, 0)]).unwrap()),
        "89208c3abbcb8214bf9c90f834fb9e2babb04e88"
    );
    // Github #3796: the PNG carries no molecule metadata.
    let png = draw_morgan_bit_as_png(&m, 872, &bi).unwrap();
    assert_eq!(png_size(&png), (150, 150));
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { mol_from_png_string(&png) }.is_err());
    assert_eq!(
        png_size(&draw_morgan_bits_as_png(&[(&m, 872, &bi, 0)]).unwrap()),
        (450, 150)
    );
}

// UnitTestDraw.testDrawRDKit.
#[test]
fn draw_rdkit_matches_python() {
    let m = phenylcyclopropane();
    let bi = rdkit_bit_info(&m);
    assert_eq!(bi[&1553], [vec![0, 1, 9, 5, 4], vec![2, 3, 4, 9, 5]]);
    let svg1 = draw_rdkit_bit(&m, 1553, &bi).unwrap();
    assert_eq!(svg1, draw_rdkit_env(&m, &bi[&1553][0]).unwrap());
    assert_eq!(sha1(&svg1), "e7d08f5214ef9855a5bd7296651860d020ca6543");
    assert!(svg1.contains("style='fill:#E5E533;"));
    assert!(!svg1.contains("style='fill:#CCCCCC;"));
    assert!(!svg1.contains("style='fill:#9999E5;"));
    assert!(matches!(
        draw_rdkit_bit(&m, 32, &bi),
        Err(Error::Key { .. })
    ));
    let non_aromatic = DrawRDKitEnvParams {
        non_aromatic_color: Some(DrawColour::new(0.2, 0.8, 0.2)),
        ..Default::default()
    };
    assert_eq!(
        sha1(&draw_rdkit_env_with_params(&m, &[5, 6, 7], &non_aromatic).unwrap()),
        "392bab8cfa6710e8ba0006c2cb80aecac187df44"
    );
    let legends = DrawRDKitEnvsParams {
        legends: Some(["x", "y"].map(str::to_owned).to_vec()),
        ..Default::default()
    };
    assert_eq!(
        sha1(
            &draw_rdkit_envs_with_params(&[(&m, &bi[&1553][0]), (&m, &[5, 6])], &legends).unwrap()
        ),
        "c221de23347abe4f2ca8d3980bc1478d4efb69f2"
    );
    assert_eq!(
        sha1(&draw_rdkit_bits(&[(&m, 1553, &bi, 0)]).unwrap()),
        "d71a63e31d9f82d07d240c6bf5586d67f8546c33"
    );
    let png = draw_rdkit_bit_as_png(&m, 1553, &bi).unwrap();
    assert_eq!(png_size(&png), (150, 150));
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { mol_from_png_string(&png) }.is_err());
}

#[test]
fn debug_draw_matches_python() {
    let m = aspirin();
    assert_eq!(
        sha1(&debug_draw(&m).unwrap()),
        "7e8cbe59c15a6de6d001225dc4f15652ce75e746"
    );
    let params = DebugDrawParams {
        size: (300, 250),
        use_bw: false,
        include_h_labels: false,
        add_atom_indices: false,
        add_bond_indices: true,
    };
    assert_eq!(
        sha1(&debug_draw_with_params(&m, &params).unwrap()),
        "0337fde4bd0dcb70b4531d11a2de6b8cac48ce67"
    );
    // The caller's molecule gets no notes.
    assert!(!m.get_atom_with_idx(0).unwrap().has_prop("atomNote"));
    assert_eq!(png_size(&debug_draw_as_png(&m).unwrap()), (350, 350));
    let mut drawer = MolDraw2DSVG::new(350, 350).unwrap();
    debug_draw_with_drawer(&m, &mut drawer, &DebugDrawParams::default()).unwrap();
    drawer.finish_drawing();
    assert_eq!(
        sha1(&drawer.get_drawing_text().unwrap()),
        "7e8cbe59c15a6de6d001225dc4f15652ce75e746"
    );
}

// UnitTestDraw.testDrawReaction.
#[test]
fn reaction_to_image_matches_python() {
    let rxn = reaction_from_smarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]").unwrap();
    assert_eq!(
        sha1(&reaction_to_image_use_svg(&rxn).unwrap()),
        "2e9a0591896d3f927c8761191f17f73f42718091"
    );
    let params = ReactionToImageParams {
        sub_img_size: (150, 120),
        highlight_by_reactant: true,
        ..Default::default()
    };
    assert_eq!(
        sha1(&reaction_to_image_use_svg_with_params(&rxn, &params).unwrap()),
        "43325d1a274254a3055f4bd485e201855f0e7a90"
    );
    assert_eq!(png_size(&reaction_to_image(&rxn).unwrap()), (800, 200));
    let big = reaction_from_smarts(
        "[c;H1:3]1:[c:4]:[c:5]:[c;H1:6]:[c:7]2:[nH:8]:[c:9]:[c;H1:1]:[c:2]:1:2.O=[C:10]1[#6;H2:11][#6;H2:12][N:13][#6;H2:14][#6;H2:15]1>>[#6;H2:12]3[#6;H1:11]=[C:10]([c:1]1:[c:9]:[n:8]:[c:7]2:[c:6]:[c:5]:[c:4]:[c:3]:[c:2]:1:2)[#6;H2:15][#6;H2:14][N:13]3",
    )
    .unwrap();
    assert!(reaction_to_image_use_svg(&big).unwrap().contains("<svg"));
}

/// UnitTestDraw's grid: rows of oligomers, with gaps.
fn oligomer_matrix() -> Vec<Vec<Option<Mol>>> {
    [vec![1], vec![0, 2], vec![3, 0, 4]]
        .iter()
        .map(|row| {
            row.iter()
                .map(|&count| {
                    (count > 0).then(|| mol_from_smiles(&"NC(C)C(=O)".repeat(count)).unwrap())
                })
                .collect()
        })
        .collect()
}

fn borrowed(matrix: &[Vec<Option<Mol>>]) -> Vec<Vec<Option<&Mol>>> {
    matrix
        .iter()
        .map(|row| row.iter().map(Option::as_ref).collect())
        .collect()
}

// UnitTestDraw.testMolsMatrixToGridImage and testMolsMatrixToLinear.
#[test]
fn mols_matrix_to_grid_image_matches_python() {
    let owned = oligomer_matrix();
    let mm = borrowed(&owned);
    let repeats = [vec![1], vec![0, 2], vec![3, 0, 4]];
    let legends: Vec<Vec<String>> = repeats
        .iter()
        .map(|row| row.iter().map(|c| format!("{c} unit(s)")).collect())
        .collect();
    let lists = |offset: u32| -> Vec<Vec<Vec<u32>>> {
        repeats
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&c| (0..c).map(|n| n * 5 + offset).collect())
                    .collect()
            })
            .collect()
    };
    assert_eq!(
        sha1(&mols_matrix_to_grid_image_use_svg(&mm).unwrap()),
        "cefbc52fbf428de88f5c9829155b44cdcb2da104"
    );
    let all = MolsMatrixToGridImageParams {
        sub_img_size: (150, 150),
        legends_matrix: Some(legends.clone()),
        highlight_atom_lists_matrix: Some(lists(0)),
        highlight_bond_lists_matrix: Some(lists(1)),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mols_matrix_to_grid_image_use_svg_with_params(&mm, &all).unwrap()),
        "56ba9ff273e9d53c57f92167b18643d233b797a8"
    );
    let opts = MolsMatrixToGridImageParams {
        legends_matrix: Some(legends.clone()),
        draw_options: Some(atom_indices()),
        ..Default::default()
    };
    assert_eq!(
        sha1(&mols_matrix_to_grid_image_use_svg_with_params(&mm, &opts).unwrap()),
        "3418c9c8c778602e3b9202d02570912635ff2c01"
    );
    assert_eq!(
        png_size(&mols_matrix_to_grid_image(&mm).unwrap()),
        (600, 600)
    );
    assert_eq!(
        png_size(&mols_matrix_to_grid_image_with_params(&mm, &all).unwrap()),
        (450, 450)
    );
    // Matrices of another shape are refused, with Python's messages.
    let short = MolsMatrixToGridImageParams {
        legends_matrix: Some(legends[..1].to_vec()),
        ..Default::default()
    };
    let Err(Error::Value { message }) = mols_matrix_to_grid_image_use_svg_with_params(&mm, &short)
    else {
        panic!("a short legends matrix draws");
    };
    assert_eq!(
        message,
        "If legendsMatrix is provided it must be the same length (have the same number of \
         sub-iterables) as molsMatrix, 3; its length is 1."
    );
    let mut short_row = lists(0);
    short_row[0].pop();
    let short_row = MolsMatrixToGridImageParams {
        highlight_atom_lists_matrix: Some(short_row),
        ..Default::default()
    };
    assert!(mols_matrix_to_grid_image_use_svg_with_params(&mm, &short_row).is_err());
    assert!(mols_matrix_to_grid_image_use_svg(&[]).is_err());
    assert!(mols_matrix_to_grid_image_use_svg(&[vec![]]).is_err());
    // A grid of gaps draws nothing.
    assert!(
        mols_matrix_to_grid_image_use_svg(&[vec![None, None]])
            .unwrap()
            .contains("<svg")
    );
}

// UnitTestDraw.testCairoFile.
#[test]
fn mol_to_file_matches_python() {
    let dir = std::env::temp_dir().join(format!("rdkit-cxx-mol-to-file-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let m = aspirin();
    let svg = dir.join("m.svg");
    mol_to_file(&m, &svg).unwrap();
    assert_eq!(
        sha1(&std::fs::read_to_string(&svg).unwrap()),
        "91e470f6aa690aeb25fe9ae4b9ad737d50e9cd22"
    );
    let params = MolToFileParams {
        size: (200, 180),
        highlight_atoms: vec![1, 2],
        legend: "x".to_owned(),
        ..Default::default()
    };
    mol_to_file_with_params(&m, &svg, &params).unwrap();
    assert_eq!(
        sha1(&std::fs::read_to_string(&svg).unwrap()),
        "6349cd5e23c7515e9f4c04cfac3b093b5efc22f4"
    );
    let png = dir.join("m.png");
    mol_to_file(&m, &png).unwrap();
    assert_eq!(png_size(&std::fs::read(&png).unwrap()), (300, 300));
    // The image type overrides the extension.
    let typed = MolToFileParams {
        image_type: Some("png".to_owned()),
        ..Default::default()
    };
    let other = dir.join("m.img");
    mol_to_file_with_params(&m, &other, &typed).unwrap();
    assert_eq!(png_size(&std::fs::read(&other).unwrap()), (300, 300));
    // Python's extensions are case-sensitive.
    for bad in ["m.PNG", "m.jpg", "m", ".png"] {
        assert!(matches!(
            mol_to_file(&m, dir.join(bad)),
            Err(Error::Value { .. })
        ));
    }
    assert!(matches!(mol_to_file(&m, ""), Err(Error::Value { .. })));
    assert!(matches!(
        mol_to_file(&m, dir.join("missing").join("m.svg")),
        Err(Error::Os { .. })
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn set_comic_mode_sets_the_comic_font() {
    let mut opts = MolDrawOptions::default();
    set_comic_mode(&mut opts).unwrap();
    assert!(opts.comic_mode);
    assert!(
        opts.font_file.ends_with("Data/Fonts/ComicNeue-Regular.ttf"),
        "{}",
        opts.font_file
    );
    // RDKit draws with the font.
    let svg = mol_to_svg_with_params(
        &aspirin(),
        &MolToSVGParams {
            draw_options: Some(opts),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(svg.contains("<svg"));
}

#[test]
fn draw_params_defaults_match_python() {
    let env = DrawMorganEnvParams::default();
    assert_eq!(env.mol_size, (150, 150));
    assert_eq!(env.base_rad, 0.3);
    assert_eq!(env.aromatic_color, Some(DrawColour::new(0.9, 0.9, 0.2)));
    assert_eq!(env.ring_color, Some(DrawColour::new(0.8, 0.8, 0.8)));
    assert_eq!(env.center_color, Some(DrawColour::new(0.6, 0.6, 0.9)));
    assert_eq!(env.extra_color, DrawColour::new(0.9, 0.9, 0.9));
    assert_eq!(DrawMorganEnvsParams::default().mols_per_row, 3);
    assert_eq!(DrawRDKitEnvParams::default().non_aromatic_color, None);
    assert_eq!(DrawRDKitBitParams::default().which_example, 0);
    assert_eq!(DebugDrawParams::default().size, (350, 350));
    assert_eq!(ReactionToImageParams::default().sub_img_size, (200, 200));
    assert_eq!(
        MolsMatrixToGridImageParams::default().sub_img_size,
        (200, 200)
    );
    assert_eq!(MolToFileParams::default().image_type, None);
}
