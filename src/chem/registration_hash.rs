// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.RegistrationHash`, a Python module: a hash of a
//! molecule's chemistry, so that molecules that are chemically the same hash
//! the same.
//!
//! The hash improves on SMILES because it ignores features that are not
//! chemically meaningful (atom map numbers), and canonicalizes enhanced
//! stereo groups (`C[C@H](O)CC |&1:1|` and `C[C@@H](O)CC |&1:1|` hash the
//! same) and substance group data (polymer data, for example). A hash scheme
//! picks the layers that go into it; one considers tautomers equivalent.
//!
//! The port gives the hash, the layers, and the S-group layer's JSON exactly
//! as Python does. Python's `logger` and its private helpers are not public
//! here.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};

use crate::chem::rd_mol_hash::{HashFunction, MolHashParams, mol_hash, mol_hash_with_params};
use crate::chem::{
    CXSmilesFields, Mol, MolToCXSmilesParams, PropValue, RemoveHsParameters, StereoGroupType,
    SubstanceGroupRef, canonicalize_enhanced_stereo, get_mol_substance_groups,
    mol_to_cx_smiles_with_params, mol_to_smiles, remove_hs_with_params, remove_stereochemistry,
};
use crate::{Error, Result};

/// The atom property that holds an atom map number. Mirrors
/// `RegistrationHash.ATOM_PROP_MAP_NUMBER`.
pub fn atom_prop_map_number() -> &'static str {
    "molAtomMapNumber"
}

/// The CXSMILES extensions the hash includes by default. Mirrors
/// `RegistrationHash.DEFAULT_CXFLAG`.
pub fn default_cxflag() -> CXSmilesFields {
    CXSmilesFields::CX_ATOM_LABELS
        | CXSmilesFields::CX_ENHANCEDSTEREO
        | CXSmilesFields::CX_BOND_ATROPISOMER
}

/// The pattern of an enhanced stereo group in CXSMILES, uncompiled: the
/// crate has no regular expressions. The module does not use it. Mirrors
/// `RegistrationHash.ENHANCED_STEREO_GROUP_REGEX`.
pub fn enhanced_stereo_group_regex() -> &'static str {
    r"((?:a|[&o]\d+):\d+(?:,\d+)*)"
}

/// A weight for each kind of enhanced stereo group. The module does not use
/// them. Mirrors `RegistrationHash.ENHANCED_STEREO_GROUP_WEIGHTS`.
pub fn enhanced_stereo_group_weights() -> HashMap<StereoGroupType, u32> {
    HashMap::from([
        (StereoGroupType::StereoAnd, 1000),
        (StereoGroupType::StereoOr, 2000),
        (StereoGroupType::StereoAbsolute, 3000),
    ])
}

/// The tautomer hash of a molecule without atoms. Mirrors
/// `RegistrationHash.EMPTY_MOL_TAUTOMER_HASH`.
pub fn empty_mol_tautomer_hash() -> &'static str {
    "_0_0"
}

/// A layer of information about a molecule. Mirrors
/// `RegistrationHash.HashLayer`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[non_exhaustive]
pub enum HashLayer {
    /// RDKit's canonical CXSMILES, with the extensions the layers were made
    /// with. Mirrors `CANONICAL_SMILES`.
    CanonicalSmiles,
    /// Any other information to include. Mirrors `ESCAPE`.
    Escape,
    /// The molecular formula. Mirrors `FORMULA`.
    Formula,
    /// RDKit's canonical SMILES with all stereo removed. Mirrors
    /// `NO_STEREO_SMILES`.
    NoStereoSmiles,
    /// The tautomer hash with all stereo removed. Mirrors
    /// `NO_STEREO_TAUTOMER_HASH`.
    NoStereoTautomerHash,
    /// The canonicalized data of the substance groups, as JSON. Mirrors
    /// `SGROUP_DATA`.
    SgroupData,
    /// A SMILES-like form of a generic tautomer. Mirrors `TAUTOMER_HASH`.
    TautomerHash,
}

/// Which layers a hash uses. Mirrors `RegistrationHash.HashScheme`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum HashScheme {
    /// Every layer, the strictest scheme. Mirrors `ALL_LAYERS`.
    AllLayers,
    /// Leaves out the layers that see stereo. Mirrors
    /// `STEREO_INSENSITIVE_LAYERS`.
    StereoInsensitiveLayers,
    /// Leaves out the layers that see tautomers. Mirrors
    /// `TAUTOMER_INSENSITIVE_LAYERS`.
    TautomerInsensitiveLayers,
}

impl HashScheme {
    /// The scheme's layers, in the order they are hashed. Mirrors
    /// `HashScheme.value`.
    pub fn value(self) -> &'static [HashLayer] {
        use HashLayer::*;
        match self {
            HashScheme::AllLayers => &[
                CanonicalSmiles,
                Escape,
                Formula,
                NoStereoSmiles,
                NoStereoTautomerHash,
                SgroupData,
                TautomerHash,
            ],
            HashScheme::StereoInsensitiveLayers => &[
                Escape,
                Formula,
                NoStereoSmiles,
                NoStereoTautomerHash,
                SgroupData,
            ],
            HashScheme::TautomerInsensitiveLayers => &[
                Escape,
                Formula,
                NoStereoTautomerHash,
                SgroupData,
                TautomerHash,
            ],
        }
    }
}

/// Mirrors `RegistrationHash.EnhancedStereoUpdateMode`, which the module
/// does not use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum EnhancedStereoUpdateMode {
    /// Mirrors `ADD_WEIGHTS`.
    AddWeights,
    /// Mirrors `REMOVE_WEIGHTS`.
    RemoveWeights,
}

/// The SHA-1 of the `AllLayers` scheme's layers, as 40 hex digits. Mirrors
/// `RegistrationHash.GetMolHash`.
pub fn get_mol_hash(all_layers: &BTreeMap<HashLayer, String>) -> String {
    get_mol_hash_with_hash_scheme(all_layers, HashScheme::AllLayers)
}

/// Of the scheme's layers, skipping any `all_layers` lacks. Mirrors
/// `RegistrationHash.GetMolHash(all_layers, hash_scheme)`.
pub fn get_mol_hash_with_hash_scheme(
    all_layers: &BTreeMap<HashLayer, String>,
    hash_scheme: HashScheme,
) -> String {
    let mut h = sha1_smol::Sha1::new();
    for layer in hash_scheme.value() {
        if let Some(text) = all_layers.get(layer) {
            h.update(text.as_bytes());
        }
    }
    h.digest().to_string()
}

/// Parameters for `get_mol_layers_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetMolLayersParams {
    /// The names of the data S-group fields to include; `None` (or empty)
    /// for `Atrop` alone.
    pub data_field_names: Option<Vec<String>>,
    /// Any other information to include, as the `Escape` layer.
    pub escape: Option<String>,
    /// The CXSMILES extensions to include.
    pub cxflag: CXSmilesFields,
    /// Uses version 2 of the tautomer hash.
    pub enable_tautomer_hash_v2: bool,
}

impl Default for GetMolLayersParams {
    fn default() -> GetMolLayersParams {
        GetMolLayersParams {
            data_field_names: None,
            escape: None,
            cxflag: default_cxflag(),
            enable_tautomer_hash_v2: false,
        }
    }
}

/// The layers of information that identify a molecule. Mirrors
/// `RegistrationHash.GetMolLayers`.
pub fn get_mol_layers(original_molecule: &Mol) -> Result<BTreeMap<HashLayer, String>> {
    get_mol_layers_with_params(original_molecule, &GetMolLayersParams::default())
}

/// Mirrors `RegistrationHash.GetMolLayers(original_molecule,
/// data_field_names, escape, cxflag, enable_tautomer_hash_v2)`. A data
/// S-group with more than one data field is `Error::Value`, as are an SRU's
/// odd number of `XBCORR` bonds and a property of a type Python's code could
/// not use.
pub fn get_mol_layers_with_params(
    original_molecule: &Mol,
    params: &GetMolLayersParams,
) -> Result<BTreeMap<HashLayer, String>> {
    // Work on a copy with all non-stereogenic hydrogens removed.
    let mut mol = remove_unnecessary_hs(original_molecule, true)?;
    strip_atom_map_labels(&mut mol)?;

    // Canonicalize the stereo groups only when writing enhanced stereo.
    if params.cxflag & CXSmilesFields::CX_ENHANCEDSTEREO != CXSmilesFields::CX_NONE {
        canonicalize_enhanced_stereo(&mut mol)?;
    }

    let formula = mol_hash(&mol, HashFunction::MolFormula)?;

    let cxsmiles = mol_to_cx_smiles_with_params(
        &mol,
        &MolToCXSmilesParams {
            flags: params.cxflag,
            ..Default::default()
        },
    )?;

    let tautomer_hash = get_stereo_tautomer_hash_with_params(
        &mol,
        &GetStereoTautomerHashParams {
            cxflag: params.cxflag,
            enable_tautomer_hash_v2: params.enable_tautomer_hash_v2,
        },
    )?;

    let sgroup_data = canonicalize_sgroups(&mol, params.data_field_names.as_deref())?;

    let (no_stereo_tautomer_hash, no_stereo_smiles) =
        get_no_stereo_layers_with_enable_tautomer_hash_v2(&mol, params.enable_tautomer_hash_v2)?;

    Ok(BTreeMap::from([
        (HashLayer::CanonicalSmiles, cxsmiles),
        (HashLayer::Escape, params.escape.clone().unwrap_or_default()),
        (HashLayer::Formula, formula),
        (HashLayer::NoStereoSmiles, no_stereo_smiles),
        (HashLayer::NoStereoTautomerHash, no_stereo_tautomer_hash),
        (HashLayer::SgroupData, sgroup_data),
        (HashLayer::TautomerHash, tautomer_hash),
    ]))
}

fn strip_atom_map_labels(mol: &mut Mol) -> Result<()> {
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)?
            .clear_prop(atom_prop_map_number());
    }
    Ok(())
}

/// Removes the hydrogens the hash does not need, keeping hydrogen isotopes.
fn remove_unnecessary_hs(rdk_mol: &Mol, preserve_stereogenic_hs: bool) -> Result<Mol> {
    let remove_hs_params = RemoveHsParameters {
        update_explicit_count: true,
        remove_defining_bond_stereo: !preserve_stereogenic_hs,
        sanitize: false,
        ..Default::default()
    };
    let mut edited_mol = remove_hs_with_params(rdk_mol, &remove_hs_params)?;
    edited_mol.update_property_cache_with_strict(false)?;
    Ok(edited_mol)
}

/// Parameters for `get_stereo_tautomer_hash_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetStereoTautomerHashParams {
    /// The CXSMILES extensions to include.
    pub cxflag: CXSmilesFields,
    /// Uses version 2 of the tautomer hash.
    pub enable_tautomer_hash_v2: bool,
}

impl Default for GetStereoTautomerHashParams {
    fn default() -> GetStereoTautomerHashParams {
        GetStereoTautomerHashParams {
            cxflag: default_cxflag(),
            enable_tautomer_hash_v2: false,
        }
    }
}

/// The tautomer hash, with stereo and the CXSMILES extensions. Mirrors
/// `RegistrationHash.GetStereoTautomerHash`.
pub fn get_stereo_tautomer_hash(molecule: &Mol) -> Result<String> {
    get_stereo_tautomer_hash_with_params(molecule, &GetStereoTautomerHashParams::default())
}

/// Mirrors `RegistrationHash.GetStereoTautomerHash(molecule, cxflag,
/// enable_tautomer_hash_v2)`.
pub fn get_stereo_tautomer_hash_with_params(
    molecule: &Mol,
    params: &GetStereoTautomerHashParams,
) -> Result<String> {
    if molecule.get_num_atoms() == 0 {
        return Ok(empty_mol_tautomer_hash().to_owned());
    }

    // SHARED-7909: works around https://github.com/rdkit/rdkit/issues/4234.
    let no_h_mol = remove_unnecessary_hs(molecule, false)?;

    // Always writes CXSMILES, for the enhanced stereo.
    let cx_flags_to_skip = CXSmilesFields::CX_ALL.bits() ^ params.cxflag.bits();
    mol_hash_with_params(
        &no_h_mol,
        tautomer_hash_function(params.enable_tautomer_hash_v2),
        &MolHashParams {
            use_cx_smiles: true,
            cx_flags_to_skip,
        },
    )
}

/// The tautomer hash and canonical SMILES with all stereo removed. Mirrors
/// `RegistrationHash.GetNoStereoLayers`.
pub fn get_no_stereo_layers(mol: &Mol) -> Result<(String, String)> {
    get_no_stereo_layers_with_enable_tautomer_hash_v2(mol, false)
}

/// Mirrors `RegistrationHash.GetNoStereoLayers(mol,
/// enable_tautomer_hash_v2)`.
pub fn get_no_stereo_layers_with_enable_tautomer_hash_v2(
    mol: &Mol,
    enable_tautomer_hash_v2: bool,
) -> Result<(String, String)> {
    // SHARED-7909: strips all hydrogens, stereogenic ones too.
    let mut no_stereo_mol = remove_unnecessary_hs(mol, false)?;
    remove_stereochemistry(&mut no_stereo_mol)?;
    let no_stereo_tautomer_hash = mol_hash(
        &no_stereo_mol,
        tautomer_hash_function(enable_tautomer_hash_v2),
    )?;
    let no_stereo_smiles = mol_hash(&no_stereo_mol, HashFunction::CanonicalSmiles)?;
    Ok((no_stereo_tautomer_hash, no_stereo_smiles))
}

fn tautomer_hash_function(enable_tautomer_hash_v2: bool) -> HashFunction {
    if enable_tautomer_hash_v2 {
        HashFunction::HetAtomTautomerv2
    } else {
        HashFunction::HetAtomTautomer
    }
}

/// A bond as a pair of atom ranks or indices.
type BondPair = (i64, i64);

/// Each atom's position in the canonical SMILES, and the bonds as pairs of
/// those positions, lower first.
fn get_canonical_atom_ranks_and_bonds(mol: &Mol) -> Result<(Vec<i64>, Vec<BondPair>)> {
    if mol.get_num_atoms() == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    mol_to_smiles(mol)?;
    let ordertxt = mol.get_prop("_smilesAtomOutputOrder")?;
    let inner = ordertxt
        .strip_prefix('[')
        .and_then(|text| text.strip_suffix(']'))
        .ok_or_else(|| value_error(format!("unexpected atom output order {ordertxt:?}")))?;
    let smi_order = inner
        .split(',')
        .filter(|x| !x.is_empty())
        .map(|x| x.parse::<usize>())
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|_| value_error(format!("unexpected atom output order {ordertxt:?}")))?;
    let mut at_ranks = vec![0; smi_order.len()];
    for (i, &idx) in smi_order.iter().enumerate() {
        *at_ranks
            .get_mut(idx)
            .ok_or_else(|| value_error(format!("unexpected atom output order {ordertxt:?}")))? =
            i as i64;
    }
    let mut bnd_order = Vec::new();
    for bnd in mol.get_bonds() {
        let mut bo = lookup(&at_ranks, bnd.get_begin_atom_idx())?;
        let mut eo = lookup(&at_ranks, bnd.get_end_atom_idx())?;
        if bo > eo {
            std::mem::swap(&mut bo, &mut eo);
        }
        bnd_order.push((bo, eo));
    }
    Ok((at_ranks, bnd_order))
}

fn lookup<T: Copy>(items: &[T], idx: u32) -> Result<T> {
    items.get(idx as usize).copied().ok_or(Error::Index {
        index: idx,
        len: items.len() as u32,
    })
}

fn value_error(message: impl Into<String>) -> Error {
    Error::Value {
        message: message.into(),
    }
}

/// A canonicalized S-group: Python's dict, its keys in insertion order.
type SGroupData = Vec<(&'static str, Value)>;

fn sorted_ranks(at_ranks: &[i64], atoms: &[u32]) -> Result<Value> {
    let mut ats = atoms
        .iter()
        .map(|&x| lookup(at_ranks, x))
        .collect::<Result<Vec<_>>>()?;
    ats.sort();
    Ok(Value::List(ats.into_iter().map(Value::Int).collect()))
}

fn sorted_bonds(bnd_order: &[BondPair], bonds: &[u32]) -> Result<Value> {
    let mut bnds = bonds
        .iter()
        .map(|&x| lookup(bnd_order, x))
        .collect::<Result<Vec<_>>>()?;
    bnds.sort();
    Ok(Value::List(bnds.into_iter().map(Value::pair).collect()))
}

fn canonicalize_data_sgroup(
    sg: &SubstanceGroupRef<'_>,
    at_ranks: &[i64],
    bnd_order: &[BondPair],
    field_names: &[String],
) -> Result<Option<SGroupData>> {
    if sg.get_prop("TYPE")? != "DAT" || !sg.has_prop("FIELDNAME") {
        return Ok(None);
    }
    let field_name = sg.get_prop("FIELDNAME")?;
    if !field_names.contains(&field_name) {
        return Ok(None);
    }
    let data = sg.get_string_vect_prop("DATAFIELDS")?;
    if data.len() > 1 {
        return Err(value_error(
            "cannot canonicalize data groups with multiple data fields",
        ));
    }
    let data = data
        .into_iter()
        .next()
        .ok_or(Error::Index { index: 0, len: 0 })?;
    Ok(Some(vec![
        ("fieldName", Value::Str(field_name)),
        ("atom", sorted_ranks(at_ranks, &sg.get_atoms())?),
        ("bonds", sorted_bonds(bnd_order, &sg.get_bonds())?),
        ("value", Value::Str(data)),
    ]))
}

/// The bond's atoms, the lower-ranked first (the lower index on a tie).
fn get_canonical_bond_rep(mol: &Mol, bond_idx: u32, atom_ranks: &[i64]) -> Result<BondPair> {
    let bond = mol.get_bond_with_idx(bond_idx)?;
    let mut aid1 = bond.get_begin_atom_idx();
    let mut aid2 = bond.get_end_atom_idx();
    let (rank1, rank2) = (lookup(atom_ranks, aid1)?, lookup(atom_ranks, aid2)?);
    if rank1 > rank2 || (rank1 == rank2 && aid1 > aid2) {
        std::mem::swap(&mut aid1, &mut aid2);
    }
    Ok((aid1.into(), aid2.into()))
}

/// The bond indices in an `XBHEAD` or `XBCORR` property.
fn bond_indices(value: &PropValue, key: &str) -> Result<Vec<u32>> {
    match value {
        PropValue::UnsignedVec(values) => Ok(values.clone()),
        PropValue::IntVec(values) => values
            .iter()
            .map(|&x| u32::try_from(x).map_err(|_| value_error(format!("negative bond in {key}"))))
            .collect(),
        _ => Err(value_error(format!("{key} is not a list of bonds"))),
    }
}

fn canonicalize_sru_sgroup(
    mol: &Mol,
    sg: &SubstanceGroupRef<'_>,
    at_ranks: &[i64],
    bnd_order: &[BondPair],
) -> Result<Option<SGroupData>> {
    if sg.get_prop("TYPE")? != "SRU" {
        return Ok(None);
    }
    let props = sg.get_props_as_dict()?;
    let prop = |key: &str, default: Value| props.get(key).map_or(default, Value::from);
    let mut res = vec![
        ("type", Value::Str("SRU".to_owned())),
        ("atoms", sorted_ranks(at_ranks, &sg.get_atoms())?),
        ("bonds", sorted_bonds(bnd_order, &sg.get_bonds())?),
        ("index", prop("index", Value::Int(0))),
        ("connect", prop("CONNECT", Value::Str(String::new()))),
        ("label", prop("LABEL", Value::Str(String::new()))),
    ];
    if let Some(parent) = props.get("PARENT") {
        res.push(("PARENT", Value::from(parent)));
    }
    if let Some(xbhead) = props.get("XBHEAD") {
        let mut xbhbonds = bond_indices(xbhead, "XBHEAD")?
            .into_iter()
            .map(|x| get_canonical_bond_rep(mol, x, at_ranks))
            .collect::<Result<Vec<_>>>()?;
        xbhbonds.sort();
        res.push((
            "XBHEAD",
            Value::List(xbhbonds.into_iter().map(Value::pair).collect()),
        ));
    }
    if let Some(xbcorr) = props.get("XBCORR") {
        let xbcorrbonds = bond_indices(xbcorr, "XBCORR")?
            .into_iter()
            .map(|x| get_canonical_bond_rep(mol, x, at_ranks))
            .collect::<Result<Vec<_>>>()?;
        if xbcorrbonds.len() % 2 != 0 {
            return Err(value_error("XBCORR should have 2N bonds"));
        }
        // These are pairs, so sort them as such.
        let mut pairs: Vec<_> = xbcorrbonds
            .chunks(2)
            .map(|pair| (pair[0].min(pair[1]), pair[0].max(pair[1])))
            .collect();
        pairs.sort();
        res.push((
            "XBCORR",
            Value::List(
                pairs
                    .into_iter()
                    .map(|(b1, b2)| Value::List(vec![Value::pair(b1), Value::pair(b2)]))
                    .collect(),
            ),
        ));
    }
    Ok(Some(res))
}

fn canonicalize_cop_sgroup(
    sg: &SubstanceGroupRef<'_>,
    at_ranks: &[i64],
) -> Result<Option<SGroupData>> {
    if sg.get_prop("TYPE")? != "COP" {
        return Ok(None);
    }
    let props = sg.get_props_as_dict()?;
    // Leaves out the label on purpose: Marvin adds one, BIOVIA Draw does not.
    Ok(Some(vec![
        ("type", Value::Str("COP".to_owned())),
        ("atoms", sorted_ranks(at_ranks, &sg.get_atoms())?),
        (
            "index",
            props.get("index").map_or(Value::Int(0), Value::from),
        ),
    ]))
}

/// The canonicalized DAT, SRU, and COP S-groups, as Python's `json.dumps`
/// writes them.
fn canonicalize_sgroups(mol: &Mol, data_field_names: Option<&[String]>) -> Result<String> {
    let atrop = ["Atrop".to_owned()];
    let field_names = match data_field_names {
        Some(names) if !names.is_empty() => names,
        _ => &atrop,
    };
    let (at_ranks, bnd_order) = get_canonical_atom_ranks_and_bonds(mol)?;
    let mut res = Vec::new();
    for sg in get_mol_substance_groups(mol) {
        let lres = match sg.get_prop("TYPE")?.as_str() {
            "DAT" => canonicalize_data_sgroup(&sg, &at_ranks, &bnd_order, field_names)?,
            "SRU" => canonicalize_sru_sgroup(mol, &sg, &at_ranks, &bnd_order)?,
            "COP" => canonicalize_cop_sgroup(&sg, &at_ranks)?,
            _ => None,
        };
        res.extend(lres);
    }
    if res.len() > 1 {
        sort_sgroups(&mut res)?;
        // Update the "index" and "PARENT" properties.
        let mut idxmap: Vec<(Value, i64)> = Vec::new();
        for (i, itm) in res.iter_mut().enumerate() {
            if let Some(index) = field_mut(itm, "index") {
                let new_index = i as i64 + 1;
                match idxmap.iter_mut().find(|(key, _)| key.py_eq(index)) {
                    Some(entry) => entry.1 = new_index,
                    None => idxmap.push((index.clone(), new_index)),
                }
                *index = Value::Int(new_index);
            }
        }
        for itm in &mut res {
            if let Some(parent) = field_mut(itm, "PARENT") {
                let new_parent = idxmap
                    .iter()
                    .find(|(key, _)| key.py_eq(parent))
                    .ok_or_else(|| Error::Key {
                        key: parent.to_json(),
                    })?
                    .1;
                *parent = Value::Int(new_parent);
            }
        }
    }
    let mut json = String::from("[");
    for (i, itm) in res.iter().enumerate() {
        if i > 0 {
            json.push_str(", ");
        }
        json.push('{');
        for (j, (key, value)) in itm.iter().enumerate() {
            if j > 0 {
                json.push_str(", ");
            }
            write_json_string(key, &mut json);
            json.push_str(": ");
            value.write_json(&mut json);
        }
        json.push('}');
    }
    json.push(']');
    Ok(json)
}

fn field_mut<'a>(itm: &'a mut SGroupData, key: &str) -> Option<&'a mut Value> {
    itm.iter_mut()
        .find(|(name, _)| *name == key)
        .map(|(_, value)| value)
}

/// Sorts as Python's `sorted(tuple(x.items()) for x in res)`: each S-group
/// by its `(key, value)` pairs, in order. Values Python cannot order
/// (`TypeError`) are `Error::Value`.
fn sort_sgroups(res: &mut [SGroupData]) -> Result<()> {
    let as_tuple = |itm: &SGroupData| {
        Value::List(
            itm.iter()
                .map(|(key, value)| Value::List(vec![Value::Str((*key).to_owned()), value.clone()]))
                .collect(),
        )
    };
    let mut error = None;
    res.sort_by(|a, b| {
        let (a, b) = (as_tuple(a), as_tuple(b));
        match (a.py_lt(&b), b.py_lt(&a)) {
            (Ok(true), _) => Ordering::Less,
            (Ok(false), Ok(true)) => Ordering::Greater,
            (Ok(false), Ok(false)) => Ordering::Equal,
            (Err(err), _) | (_, Err(err)) => {
                error.get_or_insert(err);
                Ordering::Equal
            }
        }
    });
    error.map_or(Ok(()), Err)
}

/// A Python value in a canonicalized S-group.
#[derive(Clone, Debug)]
enum Value {
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
    List(Vec<Value>),
    None,
}

impl From<&PropValue> for Value {
    fn from(value: &PropValue) -> Value {
        match value {
            PropValue::Int(x) => Value::Int((*x).into()),
            PropValue::Unsigned(x) => Value::Int((*x).into()),
            PropValue::Double(x) => Value::Float(*x),
            PropValue::Bool(x) => Value::Bool(*x),
            PropValue::String(x) => Value::Str(x.clone()),
            PropValue::IntVec(xs) => {
                Value::List(xs.iter().map(|&x| Value::Int(x.into())).collect())
            }
            PropValue::UnsignedVec(xs) => {
                Value::List(xs.iter().map(|&x| Value::Int(x.into())).collect())
            }
            PropValue::DoubleVec(xs) => Value::List(xs.iter().map(|&x| Value::Float(x)).collect()),
            PropValue::StringVec(xs) => {
                Value::List(xs.iter().map(|x| Value::Str(x.clone())).collect())
            }
            PropValue::None => Value::None,
        }
    }
}

impl Value {
    fn pair((first, second): BondPair) -> Value {
        Value::List(vec![Value::Int(first), Value::Int(second)])
    }

    /// The value as a number, if Python treats it as one (`bool` is an
    /// `int`).
    fn number(&self) -> Option<f64> {
        match self {
            Value::Int(x) => Some(*x as f64),
            Value::Float(x) => Some(*x),
            Value::Bool(x) => Some(f64::from(u8::from(*x))),
            _ => None,
        }
    }

    /// Python's `==`.
    fn py_eq(&self, other: &Value) -> bool {
        match (self, other) {
            (Value::Int(a), Value::Int(b)) => a == b,
            (Value::Str(a), Value::Str(b)) => a == b,
            (Value::List(a), Value::List(b)) => {
                a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.py_eq(y))
            }
            (Value::None, Value::None) => true,
            _ => match (self.number(), other.number()) {
                (Some(a), Some(b)) => a == b,
                _ => false,
            },
        }
    }

    /// Python's `<`, which raises `TypeError` between types it cannot
    /// order.
    fn py_lt(&self, other: &Value) -> Result<bool> {
        match (self, other) {
            (Value::Int(a), Value::Int(b)) => Ok(a < b),
            (Value::Str(a), Value::Str(b)) => Ok(a < b),
            (Value::List(a), Value::List(b)) => match a.iter().zip(b).find(|(x, y)| !x.py_eq(y)) {
                Some((x, y)) => x.py_lt(y),
                None => Ok(a.len() < b.len()),
            },
            _ => match (self.number(), other.number()) {
                (Some(a), Some(b)) => Ok(a < b),
                _ => Err(value_error(format!(
                    "cannot order {} and {}",
                    self.to_json(),
                    other.to_json()
                ))),
            },
        }
    }

    fn to_json(&self) -> String {
        let mut json = String::new();
        self.write_json(&mut json);
        json
    }

    /// Writes the value as Python's `json.dumps` does, with its default
    /// separators and `ensure_ascii`.
    fn write_json(&self, json: &mut String) {
        match self {
            Value::Int(x) => json.push_str(&x.to_string()),
            Value::Float(x) => json.push_str(&python_float_json(*x)),
            Value::Bool(x) => json.push_str(if *x { "true" } else { "false" }),
            Value::Str(x) => write_json_string(x, json),
            Value::List(xs) => {
                json.push('[');
                for (i, x) in xs.iter().enumerate() {
                    if i > 0 {
                        json.push_str(", ");
                    }
                    x.write_json(json);
                }
                json.push(']');
            }
            Value::None => json.push_str("null"),
        }
    }
}

/// Writes a string as `json.dumps` does with `ensure_ascii`: everything
/// outside printable ASCII escaped, astral characters as surrogate pairs.
fn write_json_string(text: &str, json: &mut String) {
    json.push('"');
    for c in text.chars() {
        match c {
            '"' => json.push_str("\\\""),
            '\\' => json.push_str("\\\\"),
            '\n' => json.push_str("\\n"),
            '\r' => json.push_str("\\r"),
            '\t' => json.push_str("\\t"),
            '\u{8}' => json.push_str("\\b"),
            '\u{c}' => json.push_str("\\f"),
            ' '..='~' => json.push(c),
            _ => {
                let mut units = [0; 2];
                for unit in c.encode_utf16(&mut units) {
                    json.push_str(&format!("\\u{unit:04x}"));
                }
            }
        }
    }
    json.push('"');
}

/// Writes a float as `json.dumps` does: Python's `repr`, which writes the
/// shortest digits that round-trip, in exponent form outside
/// `1e-4 <= |x| < 1e16`.
fn python_float_json(x: f64) -> String {
    if x.is_nan() {
        return "NaN".to_owned();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Infinity" } else { "-Infinity" }.to_owned();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0" } else { "0.0" }.to_owned();
    }
    // Rust's `{:e}` writes the same shortest digits: "-1.2345e-7".
    let scientific = format!("{x:e}");
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("`{:e}` writes an exponent");
    let exponent: i32 = exponent.parse().expect("`{:e}` writes an integer exponent");
    let (sign, mantissa) = match mantissa.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", mantissa),
    };
    let digits = mantissa.replace('.', "");
    if (-4..16).contains(&exponent) {
        if exponent >= 0 {
            let point = exponent as usize + 1;
            if digits.len() > point {
                format!("{sign}{}.{}", &digits[..point], &digits[point..])
            } else {
                format!("{sign}{digits}{}.0", "0".repeat(point - digits.len()))
            }
        } else {
            format!("{sign}0.{}{digits}", "0".repeat((-exponent - 1) as usize))
        }
    } else {
        let exponent_sign = if exponent < 0 { '-' } else { '+' };
        format!("{sign}{mantissa}e{exponent_sign}{:02}", exponent.abs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_floats() {
        for (x, repr) in [
            (1.5, "1.5"),
            (100.0, "100.0"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (-1.25e-7, "-1.25e-07"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (1.2345678901234568e17, "1.2345678901234568e+17"),
            (123456.789, "123456.789"),
            (-0.0, "-0.0"),
        ] {
            assert_eq!(python_float_json(x), repr);
        }
    }

    #[test]
    fn json_strings() {
        let mut json = String::new();
        write_json_string("a\"b\\c\n\u{7f}é😀", &mut json);
        assert_eq!(json, r#""a\"b\\c\n\u007f\u00e9\ud83d\ude00""#);
    }

    #[test]
    fn python_ordering() {
        let list = |xs: Vec<Value>| Value::List(xs);
        assert!(
            list(vec![Value::Int(1)])
                .py_lt(&list(vec![Value::Int(1), Value::Int(0)]))
                .unwrap()
        );
        assert!(Value::Int(1).py_lt(&Value::Float(1.5)).unwrap());
        assert!(Value::Int(1).py_eq(&Value::Bool(true)));
        assert!(Value::Int(1).py_lt(&Value::Str("1".to_owned())).is_err());
        // Unequal types compare unequal before `<` is tried.
        assert!(!Value::Int(1).py_eq(&Value::Str("1".to_owned())));
    }
}
