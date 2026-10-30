// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The other formats: TPL, Mol2, Marvin, XYZ, PDB, CML, RDKit SVG, ChemDraw,
//! and biopolymer sequences.

use std::path::Path;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::{expect_ok, to_string};
use crate::ffi::molfiles::ffi;
use crate::util::path_bytes;

/// Defines a parameter struct with Python's defaults (CONVENTIONS.md A2).
macro_rules! defaults {
    (
        $(#[$meta:meta])*
        $name:ident { $($(#[$fmeta:meta])* $field:ident: $ty:ty = $default:expr),* $(,)? }
    ) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq)]
        pub struct $name {
            $($(#[$fmeta])* pub $field: $ty,)*
        }

        impl Default for $name {
            fn default() -> $name {
                $name { $($field: $default,)* }
            }
        }
    };
}

defaults! {
    /// Parameters for `mol_from_tpl_block_with_params`, with Python's defaults.
    MolFromTPLBlockParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Skips the first conformer, as for TPLs written by CombiCode.
        skip_first_conf: bool = false,
    }
}

defaults! {
    /// Parameters for `mol_from_tpl_file_with_params`, with Python's defaults.
    MolFromTPLFileParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Skips the first conformer, as for TPLs written by CombiCode.
        skip_first_conf: bool = false,
    }
}

/// Parses a TPL block. Mirrors `rdmolfiles.MolFromTPLBlock`, which returns
/// `None` where this returns an error.
pub fn mol_from_tpl_block(tpl_block: &str) -> Result<Mol> {
    mol_from_tpl_block_with_params(tpl_block, &MolFromTPLBlockParams::default())
}

/// Mirrors `rdmolfiles.MolFromTPLBlock(tplBlock, sanitize, skipFirstConf)`.
pub fn mol_from_tpl_block_with_params(
    tpl_block: &str,
    params: &MolFromTPLBlockParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_tpl_block(
        tpl_block,
        params.sanitize,
        params.skip_first_conf,
    )?))
}

/// Reads a TPL file. Mirrors `rdmolfiles.MolFromTPLFile`.
pub fn mol_from_tpl_file(file_name: impl AsRef<Path>) -> Result<Mol> {
    mol_from_tpl_file_with_params(file_name, &MolFromTPLFileParams::default())
}

/// Mirrors `rdmolfiles.MolFromTPLFile(fileName, sanitize, skipFirstConf)`.
pub fn mol_from_tpl_file_with_params(
    file_name: impl AsRef<Path>,
    params: &MolFromTPLFileParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_tpl_file(
        path_bytes(file_name.as_ref()),
        params.sanitize,
        params.skip_first_conf,
    )?))
}

defaults! {
    /// Parameters for `mol_to_tpl_block_with_params`, with Python's defaults.
    MolToTPLBlockParams {
        /// The atom property that holds partial charges.
        partial_charge_prop: String = "_GasteigerCharge".into(),
        /// Writes the first conformer twice, for CombiCode.
        write_first_conf_twice: bool = false,
    }
}

defaults! {
    /// Parameters for `mol_to_tpl_file_with_params`, with Python's defaults.
    MolToTPLFileParams {
        /// The atom property that holds partial charges.
        partial_charge_prop: String = "_GasteigerCharge".into(),
        /// Writes the first conformer twice, for CombiCode.
        write_first_conf_twice: bool = false,
    }
}

/// A TPL block; empty for a molecule without conformers. Mirrors
/// `rdmolfiles.MolToTPLBlock`.
pub fn mol_to_tpl_block(mol: &Mol) -> Result<String> {
    mol_to_tpl_block_with_params(mol, &MolToTPLBlockParams::default())
}

/// Mirrors `rdmolfiles.MolToTPLBlock(mol, partialChargeProp,
/// writeFirstConfTwice)`.
pub fn mol_to_tpl_block_with_params(mol: &Mol, params: &MolToTPLBlockParams) -> Result<String> {
    to_string(ffi::mol_to_tpl_block(
        mol.as_ffi(),
        &params.partial_charge_prop,
        params.write_first_conf_twice,
    )?)
}

/// Writes a TPL file. Mirrors `rdmolfiles.MolToTPLFile`.
pub fn mol_to_tpl_file(mol: &Mol, file_name: impl AsRef<Path>) -> Result<()> {
    mol_to_tpl_file_with_params(mol, file_name, &MolToTPLFileParams::default())
}

/// Mirrors `rdmolfiles.MolToTPLFile(mol, fileName, partialChargeProp,
/// writeFirstConfTwice)`.
pub fn mol_to_tpl_file_with_params(
    mol: &Mol,
    file_name: impl AsRef<Path>,
    params: &MolToTPLFileParams,
) -> Result<()> {
    Ok(ffi::mol_to_tpl_file(
        mol.as_ffi(),
        path_bytes(file_name.as_ref()),
        &params.partial_charge_prop,
        params.write_first_conf_twice,
    )?)
}

defaults! {
    /// Parameters for `mol_from_mol2_block_with_params`, with Python's
    /// defaults.
    MolFromMol2BlockParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
        /// Standardizes substructures Mol2 files write oddly.
        cleanup_substructures: bool = true,
    }
}

defaults! {
    /// Parameters for `mol_from_mol2_file_with_params`, with Python's
    /// defaults.
    MolFromMol2FileParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
        /// Standardizes substructures Mol2 files write oddly.
        cleanup_substructures: bool = true,
    }
}

/// Parses a Tripos Mol2 block with Corina atom types. Mirrors
/// `rdmolfiles.MolFromMol2Block`, which returns `None` where this returns an
/// error.
pub fn mol_from_mol2_block(mol2_block: &str) -> Result<Mol> {
    mol_from_mol2_block_with_params(mol2_block, &MolFromMol2BlockParams::default())
}

/// Mirrors `rdmolfiles.MolFromMol2Block(mol2Block, sanitize, removeHs,
/// cleanupSubstructures)`.
pub fn mol_from_mol2_block_with_params(
    mol2_block: &str,
    params: &MolFromMol2BlockParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_mol2_block(
        mol2_block,
        params.sanitize,
        params.remove_hs,
        params.cleanup_substructures,
    )?))
}

/// Reads a Tripos Mol2 file. Mirrors `rdmolfiles.MolFromMol2File`.
pub fn mol_from_mol2_file(mol2_file_name: impl AsRef<Path>) -> Result<Mol> {
    mol_from_mol2_file_with_params(mol2_file_name, &MolFromMol2FileParams::default())
}

/// Mirrors `rdmolfiles.MolFromMol2File(mol2FileName, sanitize, removeHs,
/// cleanupSubstructures)`.
pub fn mol_from_mol2_file_with_params(
    mol2_file_name: impl AsRef<Path>,
    params: &MolFromMol2FileParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_mol2_file(
        path_bytes(mol2_file_name.as_ref()),
        params.sanitize,
        params.remove_hs,
        params.cleanup_substructures,
    )?))
}

defaults! {
    /// Parameters for `mol_from_mrv_block_with_params`, with Python's
    /// defaults.
    MolFromMrvBlockParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
    }
}

defaults! {
    /// Parameters for `mol_from_mrv_file_with_params`, with Python's
    /// defaults.
    MolFromMrvFileParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
    }
}

/// Parses a Marvin document. A bond with fewer than two atom references is
/// `Error::Parse`. Mirrors `rdmolfiles.MolFromMrvBlock`, which returns `None`
/// where this returns an error.
pub fn mol_from_mrv_block(mrv_block: &str) -> Result<Mol> {
    mol_from_mrv_block_with_params(mrv_block, &MolFromMrvBlockParams::default())
}

/// Mirrors `rdmolfiles.MolFromMrvBlock(mrvBlock, sanitize, removeHs)`.
pub fn mol_from_mrv_block_with_params(
    mrv_block: &str,
    params: &MolFromMrvBlockParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_mrv_block(
        mrv_block,
        params.sanitize,
        params.remove_hs,
    )?))
}

/// Reads a Marvin file, which it reads once and checks as
/// `mol_from_mrv_block` does. Mirrors `rdmolfiles.MolFromMrvFile`.
pub fn mol_from_mrv_file(mol_file_name: impl AsRef<Path>) -> Result<Mol> {
    mol_from_mrv_file_with_params(mol_file_name, &MolFromMrvFileParams::default())
}

/// Mirrors `rdmolfiles.MolFromMrvFile(molFileName, sanitize, removeHs)`.
pub fn mol_from_mrv_file_with_params(
    mol_file_name: impl AsRef<Path>,
    params: &MolFromMrvFileParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_mrv_file(
        path_bytes(mol_file_name.as_ref()),
        params.sanitize,
        params.remove_hs,
    )?))
}

defaults! {
    /// Parameters for `mol_to_mrv_block_with_params`, with Python's defaults.
    MolToMrvBlockParams {
        /// Writes stereochemistry.
        include_stereo: bool = true,
        /// The conformer to write, or -1 for the default.
        conf_id: i32 = -1,
        /// Kekulizes the molecule before writing.
        kekulize: bool = true,
        /// Indents the XML.
        pretty_print: bool = false,
    }
}

defaults! {
    /// Parameters for `mol_to_mrv_file_with_params`, with Python's defaults.
    MolToMrvFileParams {
        /// Writes stereochemistry.
        include_stereo: bool = true,
        /// The conformer to write, or -1 for the default.
        conf_id: i32 = -1,
        /// Kekulizes the molecule before writing.
        kekulize: bool = true,
        /// Indents the XML.
        pretty_print: bool = false,
    }
}

/// A Marvin document for a molecule. Mirrors `rdmolfiles.MolToMrvBlock`.
pub fn mol_to_mrv_block(mol: &Mol) -> Result<String> {
    mol_to_mrv_block_with_params(mol, &MolToMrvBlockParams::default())
}

/// Mirrors `rdmolfiles.MolToMrvBlock(mol, includeStereo, confId, kekulize,
/// prettyPrint)`.
pub fn mol_to_mrv_block_with_params(mol: &Mol, params: &MolToMrvBlockParams) -> Result<String> {
    to_string(ffi::mol_to_mrv_block(
        mol.as_ffi(),
        params.include_stereo,
        params.conf_id,
        params.kekulize,
        params.pretty_print,
    )?)
}

/// Writes a Marvin file. Mirrors `rdmolfiles.MolToMrvFile`.
pub fn mol_to_mrv_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_mrv_file_with_params(mol, filename, &MolToMrvFileParams::default())
}

/// Mirrors `rdmolfiles.MolToMrvFile(mol, filename, includeStereo, confId,
/// kekulize, prettyPrint)`.
pub fn mol_to_mrv_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolToMrvFileParams,
) -> Result<()> {
    Ok(ffi::mol_to_mrv_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        params.include_stereo,
        params.conf_id,
        params.kekulize,
        params.pretty_print,
    )?)
}

/// Parses an XYZ block into atoms with coordinates and no bonds. Mirrors
/// `rdmolfiles.MolFromXYZBlock`, which returns `None` where this returns an
/// error.
pub fn mol_from_xyz_block(xyz_block: &str) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_xyz_block(xyz_block)?))
}

/// Reads an XYZ file. Mirrors `rdmolfiles.MolFromXYZFile`.
pub fn mol_from_xyz_file(xyz_file_name: impl AsRef<Path>) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_xyz_file(path_bytes(
        xyz_file_name.as_ref(),
    ))?))
}

defaults! {
    /// Parameters for `mol_to_xyz_block_with_params`, with Python's defaults.
    MolToXYZBlockParams {
        /// The conformer to write, or -1 for the default.
        conf_id: i32 = -1,
        /// Digits after the decimal point, at most 1000 (`Error::Value`).
        precision: u32 = 6,
    }
}

defaults! {
    /// Parameters for `mol_to_xyz_file_with_params`, with Python's defaults.
    MolToXYZFileParams {
        /// The conformer to write, or -1 for the default.
        conf_id: i32 = -1,
        /// Digits after the decimal point, at most 1000 (`Error::Value`).
        precision: u32 = 6,
    }
}

/// An XYZ block; empty for a molecule without conformers. Mirrors
/// `rdmolfiles.MolToXYZBlock`.
pub fn mol_to_xyz_block(mol: &Mol) -> Result<String> {
    mol_to_xyz_block_with_params(mol, &MolToXYZBlockParams::default())
}

/// Mirrors `rdmolfiles.MolToXYZBlock(mol, confId, precision)`.
pub fn mol_to_xyz_block_with_params(mol: &Mol, params: &MolToXYZBlockParams) -> Result<String> {
    to_string(ffi::mol_to_xyz_block(
        mol.as_ffi(),
        params.conf_id,
        params.precision,
    )?)
}

/// Writes an XYZ file. Mirrors `rdmolfiles.MolToXYZFile`.
pub fn mol_to_xyz_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_xyz_file_with_params(mol, filename, &MolToXYZFileParams::default())
}

/// Mirrors `rdmolfiles.MolToXYZFile(mol, filename, confId, precision)`.
pub fn mol_to_xyz_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolToXYZFileParams,
) -> Result<()> {
    Ok(ffi::mol_to_xyz_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        params.conf_id,
        params.precision,
    )?)
}

defaults! {
    /// Parameters for `mol_from_pdb_block_with_params`, with Python's
    /// defaults.
    MolFromPDBBlockParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
        /// RDKit's PDB flavor bits.
        flavor: u32 = 0,
        /// Bonds atoms by distance. Coordinates must then be finite and at
        /// most 1e7 in magnitude (`Error::Value`): RDKit hashes them in int.
        proximity_bonding: bool = true,
    }
}

defaults! {
    /// Parameters for `mol_from_pdb_file_with_params`, with Python's
    /// defaults.
    MolFromPDBFileParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
        /// RDKit's PDB flavor bits.
        flavor: u32 = 0,
        /// Bonds atoms by distance. Coordinates must then be finite and at
        /// most 1e7 in magnitude (`Error::Value`): RDKit hashes them in int.
        proximity_bonding: bool = true,
    }
}

/// Parses a PDB block, up to its `END` record. Mirrors
/// `rdmolfiles.MolFromPDBBlock`, which returns `None` where this returns an
/// error.
pub fn mol_from_pdb_block(mol_block: &str) -> Result<Mol> {
    mol_from_pdb_block_with_params(mol_block, &MolFromPDBBlockParams::default())
}

/// Mirrors `rdmolfiles.MolFromPDBBlock(molBlock, sanitize, removeHs, flavor,
/// proximityBonding)`.
pub fn mol_from_pdb_block_with_params(
    mol_block: &str,
    params: &MolFromPDBBlockParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_pdb_block(
        mol_block,
        params.sanitize,
        params.remove_hs,
        params.flavor,
        params.proximity_bonding,
    )?))
}

/// Reads a PDB file, which it reads once and checks as `mol_from_pdb_block`
/// does. Mirrors `rdmolfiles.MolFromPDBFile`.
pub fn mol_from_pdb_file(pdb_file_name: impl AsRef<Path>) -> Result<Mol> {
    mol_from_pdb_file_with_params(pdb_file_name, &MolFromPDBFileParams::default())
}

/// Mirrors `rdmolfiles.MolFromPDBFile(pdbFileName, sanitize, removeHs,
/// flavor, proximityBonding)`.
pub fn mol_from_pdb_file_with_params(
    pdb_file_name: impl AsRef<Path>,
    params: &MolFromPDBFileParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_pdb_file(
        path_bytes(pdb_file_name.as_ref()),
        params.sanitize,
        params.remove_hs,
        params.flavor,
        params.proximity_bonding,
    )?))
}

defaults! {
    /// Parameters for `mol_to_pdb_block_with_params`, with Python's defaults.
    MolToPDBBlockParams {
        /// The conformer to write, or -1 for all of them as models.
        conf_id: i32 = -1,
        /// RDKit's PDB flavor bits (MODEL lines, CONECT records, and so on).
        flavor: u32 = 0,
    }
}

defaults! {
    /// Parameters for `mol_to_pdb_file_with_params`, with Python's defaults.
    MolToPDBFileParams {
        /// The conformer to write, or -1 for all of them as models.
        conf_id: i32 = -1,
        /// RDKit's PDB flavor bits (MODEL lines, CONECT records, and so on).
        flavor: u32 = 0,
    }
}

/// A PDB block. Mirrors `rdmolfiles.MolToPDBBlock`.
pub fn mol_to_pdb_block(mol: &Mol) -> Result<String> {
    mol_to_pdb_block_with_params(mol, &MolToPDBBlockParams::default())
}

/// Mirrors `rdmolfiles.MolToPDBBlock(mol, confId, flavor)`.
pub fn mol_to_pdb_block_with_params(mol: &Mol, params: &MolToPDBBlockParams) -> Result<String> {
    to_string(ffi::mol_to_pdb_block(
        mol.as_ffi(),
        params.conf_id,
        params.flavor,
    )?)
}

/// Writes a PDB file. Mirrors `rdmolfiles.MolToPDBFile`.
pub fn mol_to_pdb_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_pdb_file_with_params(mol, filename, &MolToPDBFileParams::default())
}

/// Mirrors `rdmolfiles.MolToPDBFile(mol, filename, confId, flavor)`.
pub fn mol_to_pdb_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolToPDBFileParams,
) -> Result<()> {
    Ok(ffi::mol_to_pdb_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        params.conf_id,
        params.flavor,
    )?)
}

defaults! {
    /// Parameters for `mol_to_cml_block_with_params`, with Python's defaults.
    MolToCMLBlockParams {
        /// The conformer to write, or -1 for the default.
        conf_id: i32 = -1,
        /// Kekulizes the molecule before writing.
        kekulize: bool = true,
    }
}

defaults! {
    /// Parameters for `mol_to_cml_file_with_params`, with Python's defaults.
    MolToCMLFileParams {
        /// The conformer to write, or -1 for the default.
        conf_id: i32 = -1,
        /// Kekulizes the molecule before writing.
        kekulize: bool = true,
    }
}

/// A Chemical Markup Language document for a molecule. Formal charges that
/// sum past the `i32` range are `Error::Value`, since RDKit sums them in
/// `int`. Mirrors `rdmolfiles.MolToCMLBlock`.
pub fn mol_to_cml_block(mol: &Mol) -> Result<String> {
    mol_to_cml_block_with_params(mol, &MolToCMLBlockParams::default())
}

/// Mirrors `rdmolfiles.MolToCMLBlock(mol, confId, kekulize)`.
pub fn mol_to_cml_block_with_params(mol: &Mol, params: &MolToCMLBlockParams) -> Result<String> {
    to_string(ffi::mol_to_cml_block(
        mol.as_ffi(),
        params.conf_id,
        params.kekulize,
    )?)
}

/// Writes a CML file. Mirrors `rdmolfiles.MolToCMLFile`.
pub fn mol_to_cml_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_cml_file_with_params(mol, filename, &MolToCMLFileParams::default())
}

/// Mirrors `rdmolfiles.MolToCMLFile(mol, filename, confId, kekulize)`.
pub fn mol_to_cml_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolToCMLFileParams,
) -> Result<()> {
    Ok(ffi::mol_to_cml_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        params.conf_id,
        params.kekulize,
    )?)
}

defaults! {
    /// Parameters for `mol_from_rdkit_svg_with_params`, with Python's
    /// defaults.
    MolFromRDKitSVGParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// Removes hydrogens (only with `sanitize`).
        remove_hs: bool = true,
    }
}

/// Reads the molecule an RDKit drawing records in its SVG metadata. Mirrors
/// `rdmolfiles.MolFromRDKitSVG`, which returns `None` where this returns an
/// error.
///
/// # Safety
///
/// RDKit builds the atoms and bonds with `AtomFromSmiles` and
/// `BondFromSmiles`, which write the SMILES parser's debug level without a
/// lock, so no other thread may parse SMILES during the call
/// (CONVENTIONS.md U7).
pub unsafe fn mol_from_rdkit_svg(svg: &str) -> Result<Mol> {
    // SAFETY: the caller upholds the contract above.
    unsafe { mol_from_rdkit_svg_with_params(svg, &MolFromRDKitSVGParams::default()) }
}

/// Mirrors `rdmolfiles.MolFromRDKitSVG(svg, sanitize, removeHs)`.
///
/// # Safety
///
/// As `mol_from_rdkit_svg`.
pub unsafe fn mol_from_rdkit_svg_with_params(
    svg: &str,
    params: &MolFromRDKitSVGParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_rdkit_svg(
        svg,
        params.sanitize,
        params.remove_hs,
    )?))
}

defaults! {
    /// Parameters for `mol_from_sequence_with_params`, with Python's
    /// defaults.
    MolFromSequenceParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// 0 for L-protein, 1 for D-protein, 2 to 5 for RNA and 6 to 9 for DNA
        /// (no caps, 5' cap, 3' cap, both caps).
        flavor: i32 = 0,
    }
}

defaults! {
    /// Parameters for `mol_from_fasta_with_params`, with Python's defaults.
    MolFromFASTAParams {
        /// Sanitizes the molecule.
        sanitize: bool = true,
        /// As `MolFromSequenceParams::flavor`.
        flavor: i32 = 0,
    }
}

/// Builds a peptide or nucleic acid from its one-letter sequence. Mirrors
/// `rdmolfiles.MolFromSequence`, which returns `None` where this returns an
/// error.
pub fn mol_from_sequence(text: &str) -> Result<Mol> {
    mol_from_sequence_with_params(text, &MolFromSequenceParams::default())
}

/// Mirrors `rdmolfiles.MolFromSequence(text, sanitize, flavor)`.
pub fn mol_from_sequence_with_params(text: &str, params: &MolFromSequenceParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_sequence(
        text,
        params.sanitize,
        params.flavor,
    )?))
}

/// Builds a biopolymer from FASTA, named for its title line. Mirrors
/// `rdmolfiles.MolFromFASTA`, which returns `None` where this returns an
/// error.
pub fn mol_from_fasta(text: &str) -> Result<Mol> {
    mol_from_fasta_with_params(text, &MolFromFASTAParams::default())
}

/// Mirrors `rdmolfiles.MolFromFASTA(text, sanitize, flavor)`.
pub fn mol_from_fasta_with_params(text: &str, params: &MolFromFASTAParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_fasta(
        text,
        params.sanitize,
        params.flavor,
    )?))
}

/// Builds a biopolymer from HELM. Mirrors `rdmolfiles.MolFromHELM`, which
/// returns `None` where this returns an error.
pub fn mol_from_helm(text: &str) -> Result<Mol> {
    mol_from_helm_with_sanitize(text, true)
}

/// Mirrors `rdmolfiles.MolFromHELM(text, sanitize)`.
pub fn mol_from_helm_with_sanitize(text: &str, sanitize: bool) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_helm(text, sanitize)?))
}

/// The one-letter sequence of a molecule's residues, read from their PDB
/// residue info. Mirrors `rdmolfiles.MolToSequence`.
pub fn mol_to_sequence(mol: &Mol) -> Result<String> {
    to_string(ffi::mol_to_sequence(mol.as_ffi())?)
}

/// FASTA for a molecule's residues. Mirrors `rdmolfiles.MolToFASTA`.
pub fn mol_to_fasta(mol: &Mol) -> Result<String> {
    to_string(ffi::mol_to_fasta(mol.as_ffi())?)
}

/// HELM for a molecule's residues; empty unless every atom has PDB residue
/// info. Mirrors `rdmolfiles.MolToHELM`.
pub fn mol_to_helm(mol: &Mol) -> Result<String> {
    to_string(ffi::mol_to_helm(mol.as_ffi())?)
}

rdkit_enum! {
    /// A ChemDraw format. Mirrors `rdmolfiles.CDXMLFormat`.
    CDXMLFormat: ffi::CDXMLFormat {
        Cdxml = CDXML,
        Cdx = CDX,
        Auto = Auto,
    }
}

/// Whether RDKit was built with ChemDraw's library, which reads and writes
/// binary CDX. Mirrors `rdmolfiles.HasChemDrawCDXSupport`.
pub fn has_chem_draw_cdx_support() -> bool {
    expect_ok(ffi::has_chem_draw_cdx_support())
}

/// A CDXML document for a molecule. RDKit writes coordinates as 32-bit fixed
/// point, scaling 2D ones to ChemDraw's bond length, so coordinates that do
/// not fit, or a 2D molecule without bonds of nonzero length, are
/// `Error::Value`. Mirrors `rdmolfiles.MolToCDXMLBlock`.
pub fn mol_to_cdxml_block(mol: &Mol) -> Result<String> {
    let bytes = mol_to_cdxml_block_format(mol, CDXMLFormat::Cdxml)?;
    String::from_utf8(bytes).map_err(|err| crate::Error::Value {
        message: format!("invalid UTF-8 in a string from RDKit: {err}"),
    })
}

/// A ChemDraw document in `format`, as bytes: binary for `Cdx`, text
/// otherwise. Mirrors `rdmolfiles.MolToCDXMLBlock(mol, format)`, which returns
/// `bytes` for `CDX` and `str` otherwise.
pub fn mol_to_cdxml_block_format(mol: &Mol, format: CDXMLFormat) -> Result<Vec<u8>> {
    Ok(ffi::mol_to_cdxml_block(mol.as_ffi(), format.to_ffi())?)
}
