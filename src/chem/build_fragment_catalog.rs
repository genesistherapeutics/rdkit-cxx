// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.BuildFragmentCatalog`, a Python script for
//! CASE-type analysis: fragment catalogs built from molecules, the
//! molecules scored by the fragments they have, and the fragments ranked by
//! how much they tell about the molecules' activities.
//!
//! The script's command-line interface is not ported, nor its progress
//! messages (each function's `reportFreq`), which Python writes to stdout.

use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::chem::Mol;
use crate::chem::fragment_catalog::{py_float, py_int};
use crate::chem::rdfragcatalog::{FragCatGenerator, FragCatParams, FragCatalog, FragFPGenerator};
use crate::data_structs::ExplicitBitVect;
use crate::ml::info_theory::{InfoBitRanker, InfoType};
use crate::{Error, Result};

/// The functional groups `build_catalog` removes by default: RDKit's
/// `Data/FunctionalGroups.txt`.
const DEFAULT_FUNCTIONAL_GROUPS: &str =
    include_str!("../../wrapper/third_party/rdkit/Data/FunctionalGroups.txt");

/// Where the molecules' activities come from. Python's `actName` and `acts`
/// arguments.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum Activities {
    /// The property of each molecule named by the first molecule's last
    /// property. Python's default.
    #[default]
    LastProp,
    /// A property of each molecule, as an integer. Python's `actName`.
    Prop(String),
    /// The activities, by molecule. Python's `acts`; as in Python, an empty
    /// list reads the first molecule's last property.
    Values(Vec<i64>),
}

/// The activities of molecules, read from `mols` (`None` for a molecule
/// that didn't parse) as Python reads them.
struct ActivityReader<'a> {
    name: Option<String>,
    values: &'a [i64],
}

impl<'a> ActivityReader<'a> {
    fn new(activities: &'a Activities, first: Option<&Mol>) -> Result<ActivityReader<'a>> {
        let (name, values) = match activities {
            Activities::Values(v) if !v.is_empty() => (None, v.as_slice()),
            Activities::Prop(name) if !name.is_empty() => (Some(name.clone()), &[][..]),
            _ => {
                let first = first.ok_or_else(|| Error::Value {
                    message: "no first molecule to name the activity".to_owned(),
                })?;
                let name = first.get_prop_names()?.pop().ok_or_else(|| Error::Value {
                    message: "the first molecule has no properties".to_owned(),
                })?;
                (Some(name), &[][..])
            }
        };
        Ok(ActivityReader { name, values })
    }

    /// The activity of the `i`-th molecule.
    fn get(&self, i: usize, mol: &Mol) -> Result<i64> {
        match &self.name {
            Some(name) => py_int(&mol.get_prop(name)?),
            None => self.values.get(i).copied().ok_or(Error::Index {
                index: to_u32(i),
                len: to_u32(self.values.len()),
            }),
        }
    }
}

fn to_u32(n: usize) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

/// The index `i` of a numpy axis of length `len`, negative from the end.
fn np_index(i: i64, len: usize) -> Result<usize> {
    let len_i = i64::try_from(len).unwrap_or(i64::MAX);
    let j = if i < 0 { i + len_i } else { i };
    usize::try_from(j)
        .ok()
        .filter(|&j| j < len)
        .ok_or(Error::Index {
            index: u32::try_from(i).unwrap_or(u32::MAX),
            len: to_u32(len),
        })
}

/// Parameters for `build_catalog_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BuildCatalogParams {
    /// The number of molecules after which to stop; zero or less for all.
    pub max_pts: i64,
    /// The file of functional groups to remove from fragments; `None` for
    /// RDKit's `Data/FunctionalGroups.txt`, a copy of which is embedded.
    pub group_file_name: Option<PathBuf>,
    /// The fewest bonds of a fragment with a bit.
    pub min_path: u32,
    /// The most bonds of a fragment with a bit.
    pub max_path: u32,
}

impl Default for BuildCatalogParams {
    fn default() -> BuildCatalogParams {
        BuildCatalogParams {
            max_pts: -1,
            group_file_name: None,
            min_path: 2,
            max_path: 6,
        }
    }
}

/// `FragCatParams` from the embedded functional groups, which RDKit reads
/// only from a file: written to a new file in the temporary directory,
/// removed once read.
fn default_frag_cat_params(min_path: u32, max_path: u32) -> Result<FragCatParams> {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let os_error = |e: std::io::Error| Error::Os {
        message: e.to_string(),
    };
    loop {
        let path = std::env::temp_dir().join(format!(
            "rdkit-cxx-FunctionalGroups-{}-{}.txt",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let mut file = match std::fs::File::create_new(&path) {
            Ok(file) => file,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(os_error(e)),
        };
        let written = file
            .write_all(DEFAULT_FUNCTIONAL_GROUPS.as_bytes())
            .and_then(|()| file.sync_all());
        drop(file);
        let params = written
            .map_err(os_error)
            .and_then(|()| FragCatParams::new(min_path, max_path, &path));
        let _ = std::fs::remove_file(&path);
        return params;
    }
}

/// A fragment catalog of molecules' fragments of 2 to 6 bonds. Mirrors
/// `BuildFragmentCatalog.BuildCatalog(suppl)`.
pub fn build_catalog(mols: &[Mol]) -> Result<FragCatalog> {
    build_catalog_with_params(mols, &BuildCatalogParams::default())
}

/// `build_catalog` with parameters. Mirrors
/// `BuildFragmentCatalog.BuildCatalog(suppl, maxPts, groupFileName,
/// minPath, maxPath)`. `Error::Os` when the default functional groups can't
/// be written to the temporary directory, and `FragCatParams::new`'s
/// errors.
pub fn build_catalog_with_params(mols: &[Mol], params: &BuildCatalogParams) -> Result<FragCatalog> {
    let fp_params = match &params.group_file_name {
        Some(file) => FragCatParams::new(params.min_path, params.max_path, file)?,
        None => default_frag_cat_params(params.min_path, params.max_path)?,
    };
    let mut catalog = FragCatalog::new(&fp_params);
    let fgen = FragCatGenerator;
    let n_pts = usize::try_from(params.max_pts)
        .ok()
        .filter(|&n| n > 0)
        .unwrap_or(mols.len());
    for mol in mols.iter().take(n_pts) {
        fgen.add_frags_from_mol(mol, &mut catalog)?;
    }
    Ok(catalog)
}

/// Per bit, whether a molecule lacks (0) or has (1) it, per activity, the
/// number of molecules: numpy's `nBits x 2 x nActs` array of `int32`, which
/// wraps.
pub type ScoreTable = Vec<[Vec<i32>; 2]>;

/// Parameters for `score_molecules_with_params`, with Python's defaults
/// (A2). Python's `maxPts` is not one: `ScoreMolecules` ignores it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoreMoleculesParams {
    /// Where the activities come from.
    pub activities: Activities,
    /// The number of activities.
    pub n_acts: usize,
}

impl Default for ScoreMoleculesParams {
    fn default() -> ScoreMoleculesParams {
        ScoreMoleculesParams {
            activities: Activities::default(),
            n_acts: 2,
        }
    }
}

/// Counts a molecule's bits into the table, as Python does: its bit `id`
/// in row `id - 1`, so bit 0 in the last row.
fn score_bits(
    table: &mut ScoreTable,
    bits: impl IntoIterator<Item = u32>,
    act: usize,
) -> Result<()> {
    for row in table.iter_mut() {
        row[0][act] = row[0][act].wrapping_add(1);
    }
    for id in bits {
        let row = np_index(i64::from(id) - 1, table.len())?;
        let row = &mut table[row];
        row[0][act] = row[0][act].wrapping_sub(1);
        row[1][act] = row[1][act].wrapping_add(1);
    }
    Ok(())
}

fn new_table(n_bits: u32, n_acts: usize) -> ScoreTable {
    vec![[vec![0; n_acts], vec![0; n_acts]]; n_bits as usize]
}

/// The molecules, `None` for those that didn't parse, scored by a catalog's
/// fragments: the table of counts and each molecule's bits, empty for
/// `None`. Mirrors `BuildFragmentCatalog.ScoreMolecules(suppl, catalog)`.
pub fn score_molecules(
    mols: &[Option<Mol>],
    catalog: &FragCatalog,
) -> Result<(ScoreTable, Vec<Vec<u32>>)> {
    score_molecules_with_params(mols, catalog, &ScoreMoleculesParams::default())
}

/// `score_molecules` with parameters. Mirrors
/// `BuildFragmentCatalog.ScoreMolecules(suppl, catalog, maxPts, actName,
/// acts, nActs)`, which counts bit `id` in row `id - 1`, as this does.
/// `Error::Index` for an activity outside `-n_acts..n_acts` (negative ones
/// count from the end, as numpy's indices do) or fewer activities than
/// molecules; `Error::Key` or `Error::Value` for a property missing or not
/// an integer.
pub fn score_molecules_with_params(
    mols: &[Option<Mol>],
    catalog: &FragCatalog,
    params: &ScoreMoleculesParams,
) -> Result<(ScoreTable, Vec<Vec<u32>>)> {
    let mut table = new_table(catalog.get_fp_length(), params.n_acts);
    let reader = ActivityReader::new(&params.activities, mols.first().and_then(Option::as_ref))?;
    let fpgen = FragFPGenerator;
    let mut obls = Vec::with_capacity(mols.len());
    for (i, mol) in mols.iter().enumerate() {
        let Some(mol) = mol else {
            obls.push(Vec::new());
            continue;
        };
        let act = np_index(reader.get(i, mol)?, params.n_acts)?;
        let bits = fpgen.get_fp_for_mol(mol, catalog)?.get_on_bits();
        score_bits(&mut table, bits.iter().copied(), act)?;
        obls.push(bits);
    }
    Ok((table, obls))
}

/// Parameters for `score_from_lists_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoreFromListsParams {
    /// The number of molecules to score; zero or less for all the lists.
    pub max_pts: i64,
    /// Where the activities come from.
    pub activities: Activities,
    /// The number of activities.
    pub n_acts: usize,
}

impl Default for ScoreFromListsParams {
    fn default() -> ScoreFromListsParams {
        ScoreFromListsParams {
            max_pts: -1,
            activities: Activities::default(),
            n_acts: 2,
        }
    }
}

/// `score_molecules`'s table from the molecules' bits, computed before;
/// the molecules give only their activities. Mirrors
/// `BuildFragmentCatalog.ScoreFromLists(bitLists, suppl, catalog)`.
pub fn score_from_lists(
    bit_lists: &[Vec<u32>],
    mols: &[Mol],
    catalog: &FragCatalog,
) -> Result<ScoreTable> {
    score_from_lists_with_params(bit_lists, mols, catalog, &ScoreFromListsParams::default())
}

/// `score_from_lists` with parameters. Mirrors
/// `BuildFragmentCatalog.ScoreFromLists(bitLists, suppl, catalog, maxPts,
/// actName, acts, nActs)`. `Error::Index` for fewer molecules or lists than
/// it scores and for a bit past the catalog's, and as
/// `score_molecules_with_params` for activities.
pub fn score_from_lists_with_params(
    bit_lists: &[Vec<u32>],
    mols: &[Mol],
    catalog: &FragCatalog,
    params: &ScoreFromListsParams,
) -> Result<ScoreTable> {
    let n_pts = usize::try_from(params.max_pts)
        .ok()
        .filter(|&n| n > 0)
        .unwrap_or(bit_lists.len());
    let mut table = new_table(catalog.get_fp_length(), params.n_acts);
    let reader = ActivityReader::new(&params.activities, mols.first())?;
    for i in 0..n_pts {
        let index_error = |len: usize| Error::Index {
            index: to_u32(i),
            len: to_u32(len),
        };
        let mol = mols.get(i).ok_or_else(|| index_error(mols.len()))?;
        let act = np_index(reader.get(i, mol)?, params.n_acts)?;
        let bits = bit_lists
            .get(i)
            .ok_or_else(|| index_error(bit_lists.len()))?;
        // Python counts each row once, as a set of `id - 1`.
        let mut ids: Vec<u32> = bits.clone();
        ids.sort_unstable();
        ids.dedup();
        score_bits(&mut table, ids, act)?;
    }
    Ok(table)
}

/// Parameters for `calc_gains_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CalcGainsParams {
    /// The number of bits to rank; negative for all.
    pub top_n: i64,
    /// Where the activities come from.
    pub activities: Activities,
    /// The number of activities.
    pub n_acts: u32,
    /// The activities whose bits rank first, by biased entropy; `None` or
    /// empty to rank by entropy.
    pub bias_list: Option<Vec<u32>>,
    /// Returns the molecules' fingerprints too.
    pub collect_fps: bool,
}

impl Default for CalcGainsParams {
    fn default() -> CalcGainsParams {
        CalcGainsParams {
            top_n: -1,
            activities: Activities::default(),
            n_acts: 2,
            bias_list: None,
            collect_fps: false,
        }
    }
}

/// Parameters for `calc_gains_from_fps_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CalcGainsFromFpsParams {
    /// The number of bits to rank; negative for all.
    pub top_n: i64,
    /// Where the activities come from.
    pub activities: Activities,
    /// The number of activities.
    pub n_acts: u32,
    /// The activities whose bits rank first, by biased entropy; `None` or
    /// empty to rank by entropy.
    pub bias_list: Option<Vec<u32>>,
}

impl Default for CalcGainsFromFpsParams {
    fn default() -> CalcGainsFromFpsParams {
        CalcGainsFromFpsParams {
            top_n: -1,
            activities: Activities::default(),
            n_acts: 2,
            bias_list: None,
        }
    }
}

fn new_ranker(n_bits: u32, n_acts: u32, bias_list: Option<&[u32]>) -> Result<InfoBitRanker> {
    Ok(match bias_list {
        Some(bias) if !bias.is_empty() => {
            let mut ranker =
                InfoBitRanker::new_with_info_type(n_bits, n_acts, InfoType::Biasentropy);
            ranker.set_bias_list(bias)?;
            ranker
        }
        _ => InfoBitRanker::new_with_info_type(n_bits, n_acts, InfoType::Entropy),
    })
}

fn top_n(top_n: i64, n_bits: u32) -> u32 {
    u32::try_from(top_n).unwrap_or(n_bits)
}

fn label(act: i64) -> Result<u32> {
    u32::try_from(act).map_err(|_| Error::Value {
        message: format!("invalid activity {act}"),
    })
}

/// The bits of a catalog ranked by their information gains over the
/// molecules' activities: per bit, its id, gain, and per activity the
/// number of molecules with it, as `InfoBitRanker::get_top_n` gives them.
/// Mirrors `BuildFragmentCatalog.CalcGains(suppl, catalog)`.
pub fn calc_gains(
    mols: &[Mol],
    catalog: &FragCatalog,
) -> Result<(Vec<Vec<f64>>, Vec<ExplicitBitVect>)> {
    calc_gains_with_params(mols, catalog, &CalcGainsParams::default())
}

/// `calc_gains` with parameters, returning the molecules' fingerprints too
/// for `collect_fps` (otherwise none). Mirrors
/// `BuildFragmentCatalog.CalcGains(suppl, catalog, topN, actName, acts,
/// nActs, reportFreq, biasList, collectFps)`. `Error::Value` for a negative
/// activity, `Error::Key` or `Error::Value` for a property missing or not
/// an integer, and the ranker's errors.
pub fn calc_gains_with_params(
    mols: &[Mol],
    catalog: &FragCatalog,
    params: &CalcGainsParams,
) -> Result<(Vec<Vec<f64>>, Vec<ExplicitBitVect>)> {
    let n_bits = catalog.get_fp_length();
    let reader = ActivityReader::new(&params.activities, mols.first())?;
    let mut ranker = new_ranker(n_bits, params.n_acts, params.bias_list.as_deref())?;
    let fpgen = FragFPGenerator;
    let mut fps = Vec::new();
    for (i, mol) in mols.iter().enumerate() {
        let act = label(reader.get(i, mol)?)?;
        let fp = fpgen.get_fp_for_mol(mol, catalog)?;
        ranker.accumulate_votes(&fp, act)?;
        if params.collect_fps {
            fps.push(fp);
        }
    }
    Ok((ranker.get_top_n(top_n(params.top_n, n_bits))?, fps))
}

/// `calc_gains` from the molecules' fingerprints, computed before; the
/// molecules give only their activities. Mirrors
/// `BuildFragmentCatalog.CalcGainsFromFps(suppl, fps)`.
pub fn calc_gains_from_fps(mols: &[Mol], fps: &[ExplicitBitVect]) -> Result<Vec<Vec<f64>>> {
    calc_gains_from_fps_with_params(mols, fps, &CalcGainsFromFpsParams::default())
}

/// `calc_gains_from_fps` with parameters. Mirrors
/// `BuildFragmentCatalog.CalcGainsFromFps(suppl, fps, topN, actName, acts,
/// nActs, reportFreq, biasList)`. `Error::Index` for no fingerprints or
/// fewer than molecules, and as `calc_gains_with_params`.
pub fn calc_gains_from_fps_with_params(
    mols: &[Mol],
    fps: &[ExplicitBitVect],
    params: &CalcGainsFromFpsParams,
) -> Result<Vec<Vec<f64>>> {
    let index_error = |i: usize| Error::Index {
        index: to_u32(i),
        len: to_u32(fps.len()),
    };
    let n_bits = fps.first().ok_or_else(|| index_error(0))?.get_num_bits();
    let reader = ActivityReader::new(&params.activities, mols.first())?;
    let mut ranker = new_ranker(n_bits, params.n_acts, params.bias_list.as_deref())?;
    for (i, mol) in mols.iter().enumerate() {
        let act = label(reader.get(i, mol)?)?;
        let fp = fps.get(i).ok_or_else(|| index_error(i))?;
        ranker.accumulate_votes(fp, act)?;
    }
    ranker.get_top_n(top_n(params.top_n, n_bits))
}

/// `'%.6f' % x` of a Python float.
fn percent_f6(x: f64) -> String {
    if x.is_nan() {
        "nan".to_owned()
    } else if x.is_infinite() {
        if x > 0.0 { "inf" } else { "-inf" }.to_owned()
    } else {
        format!("{x:.6}")
    }
}

/// `'%d' % x` of a Python float: truncated toward zero.
fn percent_d(x: f64) -> Result<String> {
    if !x.is_finite() {
        return Err(Error::Value {
            message: format!("cannot convert float {x} to integer"),
        });
    }
    let t = x.trunc();
    Ok(if t == 0.0 {
        "0".to_owned()
    } else {
        format!("{t:.0}")
    })
}

/// Writes gains, as `calc_gains` gives them, as CSV: a header, then per bit
/// its id, its description in `cat` (none for `None`), its gain to six
/// decimals, and its counts. Mirrors
/// `BuildFragmentCatalog.OutputGainsData(outF, gains, cat)`.
pub fn output_gains_data(
    out_f: &mut impl Write,
    gains: &[Vec<f64>],
    cat: Option<&FragCatalog>,
) -> Result<()> {
    output_gains_data_with_n_acts(out_f, gains, cat, 2)
}

/// `output_gains_data` with `n_acts` activities named in the header.
/// Mirrors `BuildFragmentCatalog.OutputGainsData(outF, gains, cat, nActs)`.
/// `Error::Os` when writing fails, `Error::Index` for a gain with fewer
/// than two values, and the catalog's errors for a bit past it.
pub fn output_gains_data_with_n_acts(
    out_f: &mut impl Write,
    gains: &[Vec<f64>],
    cat: Option<&FragCatalog>,
    n_acts: u32,
) -> Result<()> {
    let os_error = |e: std::io::Error| Error::Os {
        message: e.to_string(),
    };
    let act_headers: Vec<String> = (0..n_acts).map(|x| format!("Act-{x}")).collect();
    let header = if cat.is_some() {
        "id,Description,Gain"
    } else {
        "id,Gain"
    };
    writeln!(out_f, "{header},{}", act_headers.join(",")).map_err(os_error)?;
    for entry in gains {
        if entry.len() < 2 {
            return Err(Error::Index {
                index: to_u32(entry.len()),
                len: to_u32(entry.len()),
            });
        }
        // Python's `int(entry[0])`, then the catalog's unsigned index.
        let id = entry[0].trunc();
        let mut out_l = vec![percent_d(id)?];
        if let Some(cat) = cat {
            let idx = if (0.0..=f64::from(u32::MAX)).contains(&id) {
                id as u32
            } else {
                return Err(Error::Value {
                    message: format!("invalid bit id {id}"),
                });
            };
            out_l.push(cat.get_bit_description(idx)?);
        }
        out_l.push(percent_f6(entry[1]));
        for &x in &entry[2..] {
            out_l.push(percent_d(x)?);
        }
        writeln!(out_f, "{}", out_l.join(",")).map_err(os_error)?;
    }
    Ok(())
}

/// Parameters for `process_gains_data_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProcessGainsDataParams {
    /// What separates a line's fields.
    pub delim: String,
    /// The field with the bit's id.
    pub id_col: usize,
    /// The field with the bit's gain.
    pub gain_col: usize,
}

impl Default for ProcessGainsDataParams {
    fn default() -> ProcessGainsDataParams {
        ProcessGainsDataParams {
            delim: ",".to_owned(),
            id_col: 0,
            gain_col: 1,
        }
    }
}

/// The ids and gains of bits read from CSV after a header line, as
/// `output_gains_data` without a catalog writes them. Mirrors
/// `BuildFragmentCatalog.ProcessGainsData(inF)`.
pub fn process_gains_data(in_f: impl BufRead) -> Result<Vec<(String, f64)>> {
    process_gains_data_with_params(in_f, &ProcessGainsDataParams::default())
}

/// `process_gains_data` with parameters. Mirrors
/// `BuildFragmentCatalog.ProcessGainsData(inF, delim, idCol, gainCol)`.
/// `Error::Os` when reading fails, and `Error::Value` for a line without
/// the fields or a gain that isn't a number.
pub fn process_gains_data_with_params(
    mut in_f: impl BufRead,
    params: &ProcessGainsDataParams,
) -> Result<Vec<(String, f64)>> {
    let mut text = String::new();
    in_f.read_to_string(&mut text).map_err(|e| Error::Os {
        message: e.to_string(),
    })?;
    let missing = || Error::Value {
        message: "list index out of range".to_owned(),
    };
    let mut res = Vec::new();
    for line in text.split_inclusive('\n').skip(1) {
        let split: Vec<&str> = line.trim().split(params.delim.as_str()).collect();
        let id = split.get(params.id_col).ok_or_else(missing)?;
        let gain = py_float(split.get(params.gain_col).ok_or_else(missing)?)?;
        res.push((id.to_string(), gain));
    }
    Ok(res)
}
