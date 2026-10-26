// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolProcessing`: fingerprints for the molecules in a
//! file.
//!
//! RDKit's `getFingerprintsForMolsInFile` opens the file itself, so the
//! wrapper does what it does with RDKit's own suppliers and generators, over
//! a stream it checks as RDKit reads it (CONVENTIONS.md U12). RDKit's parsing
//! threads fingerprint the molecules they read, as in Python, except with a
//! generator that may construct RDKit's Ctrl-C handler (feature invariants,
//! or the CIP labeler under the new stereo perception): then RDKit's threads
//! only parse, and the calling thread, which holds the crate's lock,
//! fingerprints the molecules (U4).

use std::path::Path;

use crate::Result;
use crate::chem::rd_fingerprint_generator::FingerprintGenerator64;
use crate::data_structs::ExplicitBitVect;
use crate::error::expect_ok;
use crate::ffi::fingerprint_generator::ffi as fp_ffi;
use crate::ffi::mol_processing::ffi;
use crate::util::path_bytes;

/// How `get_fingerprints_for_mols_in_file_with_params` reads a file. Mirrors
/// `rdMolProcessing.SupplierOptions`; `Default` is a default-constructed
/// RDKit `SupplierOptions` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct SupplierOptions {
    /// The threads that parse and fingerprint the molecules of SD and SMILES
    /// files: 0 for one per core, and fewer counting back from it; one thread
    /// reads in order on the calling thread. More than 1024 is
    /// `Error::Value`: RDKit starts every one before it reads.
    pub num_threads: i32,
    /// Sanitizes each molecule.
    pub sanitize: bool,
    /// Removes hydrogens (SD files, with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed SD records rather than reading what it can.
    pub strict_parsing: bool,
    /// The separator of a SMILES file's columns.
    pub delimiter: String,
    /// The column of a SMILES file that holds the SMILES.
    pub smiles_column: i32,
    /// The column of a SMILES file that holds the name, or -1 for none; below
    /// -1 is `Error::Value`, since RDKit would index before the columns.
    pub name_column: i32,
    /// Whether a SMILES file's first line names its columns.
    pub title_line: bool,
    /// The TDT field that holds the name.
    pub name_record: String,
    /// The id a TDT file's 2D coordinates get, or -1 to skip them.
    pub conf_id_2d: i32,
    /// The id a TDT file's 3D coordinates get, or -1 to skip them.
    pub conf_id_3d: i32,
}

impl Default for SupplierOptions {
    fn default() -> SupplierOptions {
        let data = expect_ok(ffi::supplier_options_default());
        SupplierOptions {
            num_threads: data.num_threads,
            sanitize: data.sanitize,
            remove_hs: data.remove_hs,
            strict_parsing: data.strict_parsing,
            delimiter: data.delimiter,
            smiles_column: data.smiles_column,
            name_column: data.name_column,
            title_line: data.title_line,
            name_record: data.name_record,
            conf_id_2d: data.conf_id_2d,
            conf_id_3d: data.conf_id_3d,
        }
    }
}

impl SupplierOptions {
    fn to_ffi(&self) -> ffi::SupplierOptionsData {
        ffi::SupplierOptionsData {
            num_threads: self.num_threads,
            sanitize: self.sanitize,
            remove_hs: self.remove_hs,
            strict_parsing: self.strict_parsing,
            delimiter: self.delimiter.clone(),
            smiles_column: self.smiles_column,
            name_column: self.name_column,
            title_line: self.title_line,
            name_record: self.name_record.clone(),
            conf_id_2d: self.conf_id_2d,
            conf_id_3d: self.conf_id_3d,
        }
    }
}

/// Parameters for `get_fingerprints_for_mols_in_file_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, Default)]
pub struct GetFingerprintsForMolsInFileParams<'a> {
    /// The generator; `None` for RDKit's default, Morgan fingerprints of
    /// radius 3 (from a 32-bit generator, as Python's first overload uses).
    pub generator: Option<&'a FingerprintGenerator64>,
    /// How the file is read.
    pub options: SupplierOptions,
}

/// Morgan fingerprints (radius 3, 2048 bits) of the molecules in a file. See
/// `get_fingerprints_for_mols_in_file_with_params`. Mirrors
/// `rdMolProcessing.GetFingerprintsForMolsInFile`.
///
/// # Examples
///
/// ```no_run
/// use rdkit_cxx::chem::rd_mol_processing::get_fingerprints_for_mols_in_file;
///
/// for fp in get_fingerprints_for_mols_in_file("molecules.smi")? {
///     println!("{:?}", fp.map(|fp| fp.get_num_on_bits()));
/// }
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub fn get_fingerprints_for_mols_in_file(
    filename: impl AsRef<Path>,
) -> Result<Vec<Option<ExplicitBitVect>>> {
    get_fingerprints_for_mols_in_file_with_params(
        filename,
        &GetFingerprintsForMolsInFileParams::default(),
    )
}

/// The fingerprints of the molecules in a file, in the file's order. The
/// extension picks the format, as for RDKit's `GeneralMolSupplier`: `.sdf`,
/// `.smi`, `.csv`, `.txt`, `.tsv`, or `.tdt`, each optionally gzipped
/// (`.gz`, or `.sdfgz`). A file that cannot be opened or is empty, or an
/// extension RDKit does not know, is `Error::Os`; MAE files (`.mae`,
/// `.maegz`), which RDKit reads with the `MaeMolSupplier` this crate does not
/// bind, are `Error::Value`, as is input past the limits U12 sets.
///
/// A record RDKit cannot read has no fingerprint (`None`). As in Python, how
/// the result ends depends on the threads: one thread (when
/// `options.num_threads` comes to one, and for TDT files) reads in order and
/// keeps every record but a final unreadable one, and an error fingerprinting
/// a molecule is the call's error; several threads keep records only up to
/// the last one that has a fingerprint, and a molecule whose fingerprint
/// fails has none. RDKit's one-thread SMILES and TDT suppliers seek, which a
/// gzip stream cannot, so they read nothing from a gzipped file, as in
/// Python. Mirrors
/// `rdMolProcessing.GetFingerprintsForMolsInFile(filename, generator,
/// options)`; Python's overloads for 32- and 64-bit generators are one
/// function here, since every generator Python makes is 64-bit.
pub fn get_fingerprints_for_mols_in_file_with_params(
    filename: impl AsRef<Path>,
    params: &GetFingerprintsForMolsInFileParams<'_>,
) -> Result<Vec<Option<ExplicitBitVect>>> {
    let path = path_bytes(filename.as_ref());
    let options = params.options.to_ffi();
    let mut list = match params.generator {
        Some(generator) => ffi::get_fingerprints_for_mols_in_file_with_generator(
            path,
            generator.as_ffi(),
            &options,
        )?,
        None => ffi::get_fingerprints_for_mols_in_file(path, &options)?,
    };
    let size = expect_ok(fp_ffi::fingerprint_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    Ok((0..size)
        .map(|i| {
            let fp = expect_ok(fp_ffi::fingerprint_list_take_bits(list.pin_mut(), i));
            (!fp.is_null()).then(|| ExplicitBitVect::from_ffi(fp))
        })
        .collect())
}
