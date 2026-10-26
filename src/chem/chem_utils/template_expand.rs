// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.ChemUtils.TemplateExpand`, a Python script: a
//! virtual library from a template, whose attachment points are dummies
//! labeled 1, 2, ..., and sidechains for each, written as an SD file. Its
//! command-line entry point is not ported.
//!
//! Python counts the molecules it writes in the module's `nDumped`; here
//! each expansion counts its own. Problems Python prints to standard error,
//! with their tracebacks, go to standard error here without them.

use std::io::Write;

use crate::chem::all_chem::{
    ConstrainedEmbedParams, compute_2d_coords, constrained_embed_with_params,
};
use crate::chem::chem_utils::align_depict::align_depict;
use crate::chem::{
    DeleteSubstructsParams, Mol, ReplaceSubstructsParams, SmilesWriteParams, add_hs,
    delete_substructs, delete_substructs_with_params, kekulize, mol_from_smarts, mol_from_smiles,
    mol_to_mol_block, mol_to_smiles, mol_to_smiles_with_params, replace_substructs_with_params,
    sanitize_mol,
};
use crate::rd_base::log_message;
use crate::{Error, Result};

/// A sidechain: its 1-based position in its supplier, and the molecule,
/// whose first atom attaches to the template.
/// A Rust addition, naming Python's sidechain tuples.
pub type Sidechain = (usize, Mol);

fn dummy() -> Mol {
    mol_from_smiles("[*]").expect("a dummy parses")
}

/// Copies the molecule's properties (not its private or computed ones) to
/// `to`.
fn copy_props(from: &Mol, to: &mut Mol) -> Result<()> {
    for prop_n in from.get_prop_names()? {
        to.set_prop(&prop_n, &from.get_prop(&prop_n)?);
    }
    Ok(())
}

/// The molecule without its dummies, written from the atom on each dummy,
/// so that atom comes first: one for the first dummy, or one for each with
/// `use_all`. Each keeps the molecule's properties. Mirrors
/// `TemplateExpand.MoveDummyNeighborsToBeginning(mol)`.
pub fn move_dummy_neighbors_to_beginning(mol: &Mol) -> Result<Vec<Mol>> {
    move_dummy_neighbors_to_beginning_with_use_all(mol, false)
}

/// Mirrors `TemplateExpand.MoveDummyNeighborsToBeginning(mol, useAll)`.
pub fn move_dummy_neighbors_to_beginning_with_use_all(
    mol: &Mol,
    use_all: bool,
) -> Result<Vec<Mol>> {
    let dummy_patt = dummy();
    let mut res = Vec::new();
    for match_ in mol.get_substruct_matches(&dummy_patt)? {
        let smi = mol_to_smiles_with_params(
            mol,
            &SmilesWriteParams {
                rooted_at_atom: match_[0] as i32,
                ..Default::default()
            },
        )?;
        let mut entry = delete_substructs(&mol_from_smiles(&smi)?, &dummy_patt)?;
        copy_props(mol, &mut entry)?;
        res.push(entry);
        if !use_all {
            break;
        }
    }
    Ok(res)
}

/// Parameters for `construct_sidechains_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ConstructSidechainsParams {
    /// SMARTS for the attachment point; `None` takes each molecule as it is,
    /// attached at its first atom.
    pub sma: Option<String>,
    /// Removes the atoms the pattern matches, attaching at their neighbor;
    /// otherwise attaches at the first matched atom.
    pub replace: bool,
    /// A sidechain for each match of the pattern, rather than the first.
    pub use_all: bool,
}

impl Default for ConstructSidechainsParams {
    fn default() -> ConstructSidechainsParams {
        ConstructSidechainsParams {
            sma: None,
            replace: true,
            use_all: false,
        }
    }
}

/// The sidechains of the supplier's molecules, each taken as it is,
/// attached at its first atom; `None` items, records that did not parse,
/// are skipped. Mirrors `TemplateExpand.ConstructSidechains(suppl)`.
pub fn construct_sidechains(
    suppl: impl IntoIterator<Item = Option<Mol>>,
) -> Result<Option<Vec<Sidechain>>> {
    construct_sidechains_with_params(suppl, &ConstructSidechainsParams::default())
}

/// As [`construct_sidechains`], at the atoms `params.sma` matches. SMARTS
/// RDKit cannot parse gives `None`, with an error in RDKit's log, as in
/// Python; a molecule the pattern does not match is logged as a warning.
/// Where removing the matched atoms leaves no dummy, Python keeps the
/// molecule without its position in its list; here it keeps the position.
/// Mirrors `TemplateExpand.ConstructSidechains(suppl, sma, replace,
/// useAll)`.
pub fn construct_sidechains_with_params(
    suppl: impl IntoIterator<Item = Option<Mol>>,
    params: &ConstructSidechainsParams,
) -> Result<Option<Vec<Sidechain>>> {
    let patt = match &params.sma {
        Some(sma) if !sma.is_empty() => match mol_from_smarts(sma) {
            Ok(patt) => Some(patt),
            Err(_) => {
                log_message(
                    "rdApp.error",
                    &format!("ERROR: could not construct pattern from smarts: {sma}\n"),
                );
                return Ok(None);
            }
        },
        _ => None,
    };
    let replacement = dummy();
    let mut res = Vec::new();
    for (idx, mol) in suppl.into_iter().enumerate() {
        let Some(mol) = mol else {
            continue;
        };
        let Some(patt) = &patt else {
            res.push((idx + 1, mol));
            continue;
        };
        if !mol.has_substruct_match(patt)? {
            log_message(
                "rdApp.warning",
                &format!(
                    "WARNING: The substructure pattern did not match sidechain {}. This may result in errors.\n",
                    idx + 1
                ),
            );
        }
        if params.replace {
            let mut tmp = replace_substructs_with_params(
                &mol,
                patt,
                &replacement,
                &ReplaceSubstructsParams::default(),
            )?;
            if !params.use_all {
                tmp.truncate(1);
            }
            for entry in tmp {
                match move_dummy_neighbors_to_beginning(&entry)?
                    .into_iter()
                    .next()
                {
                    Some(mut entry) => {
                        copy_props(&mol, &mut entry)?;
                        res.push((idx + 1, entry));
                    }
                    None => res.push((idx + 1, entry)),
                }
            }
        } else {
            for match_ in mol.get_substruct_matches(patt)? {
                let smi = mol_to_smiles_with_params(
                    &mol,
                    &SmilesWriteParams {
                        rooted_at_atom: match_[0] as i32,
                        ..Default::default()
                    },
                )?;
                let mut entry = mol_from_smiles(&smi)?;
                copy_props(&mol, &mut entry)?;
                res.push((idx + 1, entry));
            }
        }
    }
    Ok(Some(res))
}

/// Parameters for `explode_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ExplodeParams {
    /// Names each product `TemplateEnum: Mol_<n>`; otherwise the template's
    /// name and its building blocks'.
    pub auto_names: bool,
    /// Embeds each product in 3D, constrained to the template's
    /// coordinates, rather than depicting it aligned to the template.
    pub do_3d: bool,
    /// With `do_3d`, tethers the template's atoms while cleaning up.
    pub use_tethers: bool,
}

impl Default for ExplodeParams {
    fn default() -> ExplodeParams {
        ExplodeParams {
            auto_names: true,
            do_3d: false,
            use_tethers: false,
        }
    }
}

struct Exploder<'a> {
    sidechains: &'a [Vec<Sidechain>],
    core: Mol,
    template_name: String,
    params: &'a ExplodeParams,
    n_dumped: usize,
}

impl Exploder<'_> {
    fn explode(
        &mut self,
        mol: &Mol,
        depth: usize,
        chain_indices: &[(usize, usize)],
        out: &mut dyn FnMut(Mol) -> Result<()>,
    ) -> Result<()> {
        let patt = mol_from_smiles(&format!("[{}*]", depth + 1))?;
        for (i, (chain_idx, chain)) in self.sidechains[depth].iter().enumerate() {
            let mut tchain = chain_indices.to_vec();
            tchain.push((i, *chain_idx));
            let rs = replace_substructs_with_params(
                mol,
                &patt,
                chain,
                &ReplaceSubstructsParams {
                    replace_all: true,
                    ..Default::default()
                },
            )?;
            let Some(r) = rs.into_iter().next() else {
                continue;
            };
            if depth + 1 < self.sidechains.len() {
                self.explode(&r, depth + 1, &tchain, out)?;
                continue;
            }
            if let Some(r) = self.finish(r, &tchain)? {
                out(r)?;
            }
        }
        Ok(())
    }

    fn finish(&mut self, mut r: Mol, tchain: &[(usize, usize)]) -> Result<Option<Mol>> {
        if let Err(e) = sanitize_mol(&mut r) {
            eprintln!("{e}");
            return Ok(None);
        }
        if !self.params.do_3d {
            if r.has_substruct_match(&self.core)? {
                if let Err(e) = align_depict(&mut r, &self.core) {
                    eprintln!("{e}");
                    eprintln!("{}", mol_to_smiles(&r)?);
                }
            } else {
                eprintln!(">>> no match");
                compute_2d_coords(&mut r)?;
            }
        } else {
            r = add_hs(&r)?;
            constrained_embed_with_params(
                &mut r,
                &self.core,
                &ConstrainedEmbedParams {
                    use_tethers: self.params.use_tethers,
                    ..Default::default()
                },
            )?;
        }
        kekulize(&mut r)?;
        let bb_name = |bb_i: usize, bb: &(usize, usize)| -> Result<String> {
            let bb_mol = &self.sidechains[bb_i][bb.0].1;
            Ok(if bb_mol.has_prop("_Name") {
                bb_mol.get_prop("_Name")?
            } else {
                bb.1.to_string()
            })
        };
        let t_name = if self.params.auto_names {
            format!("TemplateEnum: Mol_{}", self.n_dumped + 1)
        } else {
            let mut t_name = self.template_name.clone();
            for (bb_i, bb) in tchain.iter().enumerate() {
                t_name.push('_');
                t_name.push_str(&bb_name(bb_i, bb)?);
            }
            t_name
        };
        r.set_prop("_Name", &t_name);
        r.set_prop("seq_num", &(self.n_dumped + 1).to_string());
        let indices: Vec<String> = tchain.iter().map(|x| x.1.to_string()).collect();
        r.set_prop("reagent_indices", &indices.join("_"));
        for (bb_i, bb) in tchain.iter().enumerate() {
            let bb_mol = &self.sidechains[bb_i][bb.0].1;
            r.set_prop(&format!("building_block_{}", bb_i + 1), &bb_name(bb_i, bb)?);
            r.set_int_prop(&format!("_idx_building_block_{}", bb_i + 1), bb.1 as i32);
            for prop_n in bb_mol.get_prop_names()? {
                r.set_prop(
                    &format!("building_block_{}_{prop_n}", bb_i + 1),
                    &bb_mol.get_prop(&prop_n)?,
                );
            }
        }
        self.n_dumped += 1;
        if self.n_dumped.is_multiple_of(100) {
            log_message(
                "rdApp.info",
                &format!("INFO: Done {} molecules\n", self.n_dumped),
            );
        }
        Ok(Some(r))
    }
}

/// Writes the library of the template and the sidechains, one list for
/// each attachment point, to `out_f` as SD records, with properties naming
/// each product's building blocks. A product RDKit cannot sanitize is
/// skipped, as in Python, and a write that fails is `Error::Os`. Mirrors
/// `TemplateExpand.Explode(template, sidechains, outF)`.
pub fn explode(
    template: &Mol,
    sidechains: &[Vec<Sidechain>],
    out_f: &mut impl Write,
) -> Result<()> {
    explode_with_params(template, sidechains, out_f, &ExplodeParams::default())
}

/// As [`explode`]. Mirrors `TemplateExpand.Explode(template, sidechains,
/// outF, autoNames, do3D, useTethers)`.
pub fn explode_with_params(
    template: &Mol,
    sidechains: &[Vec<Sidechain>],
    out_f: &mut impl Write,
    params: &ExplodeParams,
) -> Result<()> {
    let core =
        delete_substructs_with_params(template, &dummy(), &DeleteSubstructsParams::default())?;
    let template_name = if template.has_prop("_Name") {
        template.get_prop("_Name")?
    } else {
        "template".to_owned()
    };
    if sidechains.is_empty() {
        // Python indexes the first list.
        return Err(Error::index(0, 0));
    }
    let mut exploder = Exploder {
        sidechains,
        core,
        template_name,
        params,
        n_dumped: 0,
    };
    let os = |e: std::io::Error| Error::Os {
        message: e.to_string(),
    };
    let mut write = |mol: Mol| -> Result<()> {
        out_f
            .write_all(mol_to_mol_block(&mol)?.as_bytes())
            .map_err(os)?;
        for p_n in mol.get_prop_names()? {
            write!(out_f, ">  <{p_n}>\n{}\n\n", mol.get_prop(&p_n)?).map_err(os)?;
        }
        writeln!(out_f, "$$$$").map_err(os)
    };
    exploder.explode(template, 0, &[], &mut write)
}
