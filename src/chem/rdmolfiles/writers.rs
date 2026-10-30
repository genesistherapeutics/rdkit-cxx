// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecule writers: SD, SMILES, TDT, PDB, and Maestro files, one molecule at
//! a time.
//!
//! A writer writes to a file, to standard output (the file name `"-"`, except
//! for `MaeWriter`), or to any `Write` value in place of Python's file
//! objects, which it may borrow:
//!
//! ```no_run
//! use rdkit_cxx::chem::*;
//!
//! let mol = mol_from_smiles("CCO")?;
//! let mut out = Vec::new();
//! let mut writer = SDWriter::from_file_obj(&mut out)?;
//! writer.write(&mol)?;
//! writer.close()?;
//! drop(writer);
//! let text = String::from_utf8(out).unwrap();
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```
//!
//! The borrow outlives the writer, which writes to it when dropped:
//!
//! ```compile_fail
//! use rdkit_cxx::chem::*;
//!
//! let _writer;
//! {
//!     let mut out = Vec::new();
//!     _writer = SDWriter::from_file_obj(&mut out).unwrap();
//! }
//! ```
//!
//! Output that cannot be written, a failed `Write` included, is `Error::Os`
//! from the call that finds it and from every call after. Writes are
//! buffered: `flush` and `close` write what is left, and dropping a writer
//! closes it, as Python's `__exit__` does, without reporting errors.

use std::fmt;
use std::io::Write;
use std::marker::PhantomData;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::writers::ffi;
use crate::util::path_bytes;

/// A Rust `Write` value that C++ writes through, in place of a Python file
/// object.
pub(crate) struct RustWriter(Box<dyn Write + Send>);

impl RustWriter {
    /// Boxes `writer` for C++, which may keep it no longer than `'w`.
    ///
    /// The box's lifetime is erased. Each writer below owns the C++ object
    /// that owns the box and never lends it, and carries `'w` with an
    /// explicit `Drop`, so the borrow checker keeps `'w` alive until the
    /// writer, and with it the box, is dropped; a writer that is leaked
    /// leaks the box with it, unused.
    fn boxed<'w, W: Write + Send + 'w>(writer: W) -> Box<RustWriter> {
        let writer: Box<dyn Write + Send + 'w> = Box::new(writer);
        // SAFETY: as above. The two types differ only in the trait object's
        // lifetime bound, so they have one layout.
        let writer: Box<dyn Write + Send + 'static> = unsafe { std::mem::transmute(writer) };
        Box::new(RustWriter(writer))
    }

    /// Writes all of `buf`, retrying interrupted writes; a panic is an
    /// error, so it does not unwind into C++.
    pub(crate) fn write_all(&mut self, buf: &[u8]) -> std::result::Result<(), String> {
        match catch_unwind(AssertUnwindSafe(|| self.0.write_all(buf))) {
            Ok(result) => result.map_err(|err| err.to_string()),
            Err(_) => Err("the writer panicked".into()),
        }
    }

    pub(crate) fn flush(&mut self) -> std::result::Result<(), String> {
        match catch_unwind(AssertUnwindSafe(|| self.0.flush())) {
            Ok(result) => result.map_err(|err| err.to_string()),
            Err(_) => Err("the writer panicked".into()),
        }
    }
}

fn non_null<T: cxx::memory::UniquePtrTarget>(ptr: UniquePtr<T>) -> UniquePtr<T> {
    assert!(!ptr.is_null(), "RDKit returned a null writer");
    ptr
}

fn strings<S: AsRef<str>>(props: &[S]) -> Vec<String> {
    props.iter().map(|prop| prop.as_ref().to_owned()).collect()
}

// SAFETY (for each `unsafe impl Send` below): RDKit's writers keep no
// thread-local state, so one may be used from any thread, and a Rust writer
// is `Send`. Standard output, which several writers may share, is
// synchronized. The Rust API takes `&mut self` to write, so one thread uses a
// writer at a time.

/// The methods every writer has, stamped on each.
macro_rules! writer_methods {
    ($python:literal, $write:ident, $flush:ident, $close:ident, $num_mols:ident) => {
        /// Writes a molecule's default conformer. Mirrors
        #[doc = concat!("`", $python, ".write`.")]
        pub fn write(&mut self, mol: &Mol) -> Result<()> {
            self.write_with_conf_id(mol, -1)
        }

        /// Writes a molecule's conformer `conf_id` (-1 for the default). Mirrors
        #[doc = concat!("`", $python, ".write(mol, confId)`.")]
        pub fn write_with_conf_id(&mut self, mol: &Mol, conf_id: i32) -> Result<()> {
            Ok(ffi::$write(self.ptr.pin_mut(), mol.as_ffi(), conf_id)?)
        }

        /// Writes what is buffered. After `close`, this is `Error::Runtime`,
        /// as in Python. Mirrors
        #[doc = concat!("`", $python, ".flush`.")]
        pub fn flush(&mut self) -> Result<()> {
            Ok(ffi::$flush(self.ptr.pin_mut())?)
        }

        /// Writes what is buffered and closes the output; a file the writer
        /// opened is closed, and a `Write` value is kept until the writer is
        /// dropped. Writing afterwards is `Error::Runtime`, as in Python, and
        /// closing again does nothing. Mirrors
        #[doc = concat!("`", $python, ".close`.")]
        pub fn close(&mut self) -> Result<()> {
            Ok(ffi::$close(self.ptr.pin_mut())?)
        }

        /// The number of molecules written. Mirrors
        #[doc = concat!("`", $python, ".NumMols`.")]
        pub fn num_mols(&self) -> u32 {
            expect_ok(ffi::$num_mols(&self.ptr))
        }
    };
}

/// Parameters for `SDWriter::get_text_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SDWriterGetTextParams {
    /// The conformer to write, or -1 for the default.
    pub conf_id: i32,
    /// Kekulizes the molecule.
    pub kekulize: bool,
    /// Writes a V3000 mol block.
    pub force_v3000: bool,
    /// The molecule's id, which numbers its properties as id + 1, or -1 for
    /// none. `i32::MAX` is `Error::Value`: RDKit adds 1 in `int`.
    pub molid: i32,
}

impl Default for SDWriterGetTextParams {
    fn default() -> SDWriterGetTextParams {
        SDWriterGetTextParams {
            conf_id: -1,
            kekulize: true,
            force_v3000: false,
            molid: -1,
        }
    }
}

/// Writes SD files. Mirrors `rdmolfiles.SDWriter`.
///
/// Each record is the molecule's mol block (with 2D coordinates computed for
/// one without a conformer) and its properties: those `set_props` names, or
/// else every property that is not computed. RDKit numbers the properties of
/// a writer's molecules in `int`, so a writer that has written `i32::MAX`
/// molecules refuses more with `Error::Value`.
pub struct SDWriter<'w> {
    /// Never null.
    ptr: UniquePtr<ffi::SDWriter>,
    output: PhantomData<&'w mut ()>,
}

unsafe impl Send for SDWriter<'_> {}

/// Closes the writer, which may write to its borrowed output, so the
/// borrow must outlive it (an explicit `Drop` makes the borrow checker
/// require that).
impl Drop for SDWriter<'_> {
    fn drop(&mut self) {}
}

impl SDWriter<'static> {
    /// A writer to a file, or standard output for `"-"`; a file that cannot
    /// be opened is `Error::Os`. Mirrors `SDWriter(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<SDWriter<'static>> {
        Ok(SDWriter {
            ptr: non_null(ffi::sd_writer_new(path_bytes(file_name.as_ref()))?),
            output: PhantomData,
        })
    }
}

impl<'w> SDWriter<'w> {
    /// A writer to `file_obj`, which it keeps until it is dropped. Mirrors
    /// `SDWriter(fileObj)`.
    pub fn from_file_obj<W: Write + Send + 'w>(file_obj: W) -> Result<SDWriter<'w>> {
        Ok(SDWriter {
            ptr: non_null(ffi::sd_writer_from_writer(RustWriter::boxed(file_obj))?),
            output: PhantomData,
        })
    }

    writer_methods!(
        "SDWriter",
        sd_writer_write,
        sd_writer_flush,
        sd_writer_close,
        sd_writer_num_mols
    );

    /// Writes only these properties of each molecule, rather than every one
    /// that is not computed; RDKit warns when molecules have been written.
    /// Mirrors `SDWriter.SetProps`.
    pub fn set_props<S: AsRef<str>>(&mut self, props: &[S]) {
        expect_ok(ffi::sd_writer_set_props(
            self.ptr.pin_mut(),
            &strings(props),
        ))
    }

    /// Whether mol blocks are written as V3000. Mirrors
    /// `SDWriter.GetForceV3000`.
    pub fn get_force_v3000(&self) -> bool {
        expect_ok(ffi::sd_writer_get_force_v3000(&self.ptr))
    }

    /// Mirrors `SDWriter.SetForceV3000`.
    pub fn set_force_v3000(&mut self, val: bool) {
        expect_ok(ffi::sd_writer_set_force_v3000(self.ptr.pin_mut(), val))
    }

    /// Whether molecules are kekulized. Mirrors `SDWriter.GetKekulize`.
    pub fn get_kekulize(&self) -> bool {
        expect_ok(ffi::sd_writer_get_kekulize(&self.ptr))
    }

    /// Mirrors `SDWriter.SetKekulize`.
    pub fn set_kekulize(&mut self, val: bool) {
        expect_ok(ffi::sd_writer_set_kekulize(self.ptr.pin_mut(), val))
    }
}

impl SDWriter<'_> {
    /// The SD record a writer would write for a molecule, with every property
    /// that is not computed. Mirrors `SDWriter.GetText`.
    pub fn get_text(mol: &Mol) -> Result<String> {
        SDWriter::get_text_with_params(mol, &SDWriterGetTextParams::default())
    }

    /// Mirrors `SDWriter.GetText(mol, confId, kekulize, force_v3000, molid)`.
    pub fn get_text_with_params(mol: &Mol, params: &SDWriterGetTextParams) -> Result<String> {
        to_string(ffi::sd_writer_get_text(
            mol.as_ffi(),
            params.conf_id,
            params.kekulize,
            params.force_v3000,
            params.molid,
        )?)
    }
}

impl fmt::Debug for SDWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SDWriter")
            .field("num_mols", &self.num_mols())
            .finish_non_exhaustive()
    }
}

/// Parameters for `SmilesWriter::new_with_params` and
/// `SmilesWriter::from_file_obj_with_params`, which take the same optional
/// arguments, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SmilesWriterParams {
    /// Separates the columns.
    pub delimiter: String,
    /// The name column's header; empty for no names.
    pub name_header: String,
    /// Writes a header line.
    pub include_header: bool,
    /// Writes isomeric SMILES.
    pub isomeric_smiles: bool,
    /// Writes Kekulé SMILES.
    pub kekule_smiles: bool,
}

impl Default for SmilesWriterParams {
    fn default() -> SmilesWriterParams {
        SmilesWriterParams {
            delimiter: " ".into(),
            name_header: "Name".into(),
            include_header: true,
            isomeric_smiles: true,
            kekule_smiles: false,
        }
    }
}

impl SmilesWriterParams {
    fn to_ffi(&self) -> ffi::SmilesWriterParamsData {
        ffi::SmilesWriterParamsData {
            delimiter: self.delimiter.clone(),
            name_header: self.name_header.clone(),
            include_header: self.include_header,
            isomeric_smiles: self.isomeric_smiles,
            kekule_smiles: self.kekule_smiles,
        }
    }
}

/// Writes SMILES files: a SMILES, a name, and the properties `set_props`
/// names on each line. Mirrors `rdmolfiles.SmilesWriter`.
pub struct SmilesWriter<'w> {
    /// Never null.
    ptr: UniquePtr<ffi::SmilesWriter>,
    output: PhantomData<&'w mut ()>,
}

unsafe impl Send for SmilesWriter<'_> {}

/// Closes the writer, which may write to its borrowed output, so the
/// borrow must outlive it (an explicit `Drop` makes the borrow checker
/// require that).
impl Drop for SmilesWriter<'_> {
    fn drop(&mut self) {}
}

impl SmilesWriter<'static> {
    /// A writer to a file, or standard output for `"-"`; a file that cannot
    /// be opened is `Error::Os`. Mirrors `SmilesWriter(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<SmilesWriter<'static>> {
        SmilesWriter::new_with_params(file_name, &SmilesWriterParams::default())
    }

    /// Mirrors `SmilesWriter(fileName, delimiter, nameHeader, includeHeader,
    /// isomericSmiles, kekuleSmiles)`.
    pub fn new_with_params<P: AsRef<Path>>(
        file_name: P,
        params: &SmilesWriterParams,
    ) -> Result<SmilesWriter<'static>> {
        Ok(SmilesWriter {
            ptr: non_null(ffi::smiles_writer_new(
                path_bytes(file_name.as_ref()),
                &params.to_ffi(),
            )?),
            output: PhantomData,
        })
    }
}

impl<'w> SmilesWriter<'w> {
    /// A writer to `file_obj`, which it keeps until it is dropped. Mirrors
    /// `SmilesWriter(fileObj)`.
    pub fn from_file_obj<W: Write + Send + 'w>(file_obj: W) -> Result<SmilesWriter<'w>> {
        SmilesWriter::from_file_obj_with_params(file_obj, &SmilesWriterParams::default())
    }

    /// Mirrors `SmilesWriter(fileObj, delimiter, nameHeader, includeHeader,
    /// isomericSmiles, kekuleSmiles)`.
    pub fn from_file_obj_with_params<W: Write + Send + 'w>(
        file_obj: W,
        params: &SmilesWriterParams,
    ) -> Result<SmilesWriter<'w>> {
        Ok(SmilesWriter {
            ptr: non_null(ffi::smiles_writer_from_writer(
                RustWriter::boxed(file_obj),
                &params.to_ffi(),
            )?),
            output: PhantomData,
        })
    }

    writer_methods!(
        "SmilesWriter",
        smiles_writer_write,
        smiles_writer_flush,
        smiles_writer_close,
        smiles_writer_num_mols
    );

    /// The properties to write after each name. Once a molecule has been
    /// written, RDKit logs an error and keeps the old ones. Mirrors
    /// `SmilesWriter.SetProps`.
    pub fn set_props<S: AsRef<str>>(&mut self, props: &[S]) {
        expect_ok(ffi::smiles_writer_set_props(
            self.ptr.pin_mut(),
            &strings(props),
        ))
    }
}

impl fmt::Debug for SmilesWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SmilesWriter")
            .field("num_mols", &self.num_mols())
            .finish_non_exhaustive()
    }
}

/// Writes Daylight TDT files. Mirrors `rdmolfiles.TDTWriter`.
///
/// `close` writes the last record's closing `|` line.
pub struct TDTWriter<'w> {
    /// Never null.
    ptr: UniquePtr<ffi::TDTWriter>,
    output: PhantomData<&'w mut ()>,
}

unsafe impl Send for TDTWriter<'_> {}

/// Closes the writer, which may write to its borrowed output, so the
/// borrow must outlive it (an explicit `Drop` makes the borrow checker
/// require that).
impl Drop for TDTWriter<'_> {
    fn drop(&mut self) {}
}

impl TDTWriter<'static> {
    /// A writer to a file, or standard output for `"-"`; a file that cannot
    /// be opened is `Error::Os`. Mirrors `TDTWriter(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<TDTWriter<'static>> {
        Ok(TDTWriter {
            ptr: non_null(ffi::tdt_writer_new(path_bytes(file_name.as_ref()))?),
            output: PhantomData,
        })
    }
}

impl<'w> TDTWriter<'w> {
    /// A writer to `file_obj`, which it keeps until it is dropped. Mirrors
    /// `TDTWriter(fileObj)`.
    pub fn from_file_obj<W: Write + Send + 'w>(file_obj: W) -> Result<TDTWriter<'w>> {
        Ok(TDTWriter {
            ptr: non_null(ffi::tdt_writer_from_writer(RustWriter::boxed(file_obj))?),
            output: PhantomData,
        })
    }

    writer_methods!(
        "TDTWriter",
        tdt_writer_write,
        tdt_writer_flush,
        tdt_writer_close,
        tdt_writer_num_mols
    );

    /// Writes only these properties of each molecule, rather than every one
    /// that is not computed. Mirrors `TDTWriter.SetProps`.
    pub fn set_props<S: AsRef<str>>(&mut self, props: &[S]) {
        expect_ok(ffi::tdt_writer_set_props(
            self.ptr.pin_mut(),
            &strings(props),
        ))
    }

    /// The significant digits of each coordinate. Mirrors
    /// `TDTWriter.GetNumDigits`.
    pub fn get_num_digits(&self) -> u32 {
        expect_ok(ffi::tdt_writer_get_num_digits(&self.ptr))
    }

    /// Above 1000 is `Error::Value`: C++ formats a number in a stack buffer
    /// as long as its digits. Mirrors `TDTWriter.SetNumDigits`.
    pub fn set_num_digits(&mut self, num_digits: u32) -> Result<()> {
        Ok(ffi::tdt_writer_set_num_digits(
            self.ptr.pin_mut(),
            num_digits,
        )?)
    }

    /// Whether coordinates are written in 2D. Mirrors `TDTWriter.GetWrite2D`.
    pub fn get_write_2d(&self) -> bool {
        expect_ok(ffi::tdt_writer_get_write_2d(&self.ptr))
    }

    /// Writes 2D coordinates. Mirrors `TDTWriter.SetWrite2D`.
    pub fn set_write_2d(&mut self) {
        self.set_write_2d_with_state(true)
    }

    /// Mirrors `TDTWriter.SetWrite2D(state)`.
    pub fn set_write_2d_with_state(&mut self, state: bool) {
        expect_ok(ffi::tdt_writer_set_write_2d(self.ptr.pin_mut(), state))
    }

    /// Whether each molecule's name is written. Mirrors
    /// `TDTWriter.GetWriteNames`.
    pub fn get_write_names(&self) -> bool {
        expect_ok(ffi::tdt_writer_get_write_names(&self.ptr))
    }

    /// Writes names. Mirrors `TDTWriter.SetWriteNames`.
    pub fn set_write_names(&mut self) {
        self.set_write_names_with_state(true)
    }

    /// Mirrors `TDTWriter.SetWriteNames(state)`.
    pub fn set_write_names_with_state(&mut self, state: bool) {
        expect_ok(ffi::tdt_writer_set_write_names(self.ptr.pin_mut(), state))
    }
}

impl fmt::Debug for TDTWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TDTWriter")
            .field("num_mols", &self.num_mols())
            .finish_non_exhaustive()
    }
}

/// Writes PDB files. Mirrors `rdmolfiles.PDBWriter`.
///
/// `flavor` is `MolToPDBBlock`'s; flavor 1 writes each molecule as a model.
pub struct PDBWriter<'w> {
    /// Never null.
    ptr: UniquePtr<ffi::PDBWriter>,
    output: PhantomData<&'w mut ()>,
}

unsafe impl Send for PDBWriter<'_> {}

/// Closes the writer, which may write to its borrowed output, so the
/// borrow must outlive it (an explicit `Drop` makes the borrow checker
/// require that).
impl Drop for PDBWriter<'_> {
    fn drop(&mut self) {}
}

impl PDBWriter<'static> {
    /// A writer to a file, or standard output for `"-"`; a file that cannot
    /// be opened is `Error::Os`. Mirrors `PDBWriter(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<PDBWriter<'static>> {
        PDBWriter::new_with_flavor(file_name, 0)
    }

    /// Mirrors `PDBWriter(fileName, flavor)`.
    pub fn new_with_flavor<P: AsRef<Path>>(
        file_name: P,
        flavor: u32,
    ) -> Result<PDBWriter<'static>> {
        Ok(PDBWriter {
            ptr: non_null(ffi::pdb_writer_new(path_bytes(file_name.as_ref()), flavor)?),
            output: PhantomData,
        })
    }
}

impl<'w> PDBWriter<'w> {
    /// A writer to `file_obj`, which it keeps until it is dropped. Mirrors
    /// `PDBWriter(fileObj)`.
    pub fn from_file_obj<W: Write + Send + 'w>(file_obj: W) -> Result<PDBWriter<'w>> {
        PDBWriter::from_file_obj_with_flavor(file_obj, 0)
    }

    /// Mirrors `PDBWriter(fileObj, flavor)`.
    pub fn from_file_obj_with_flavor<W: Write + Send + 'w>(
        file_obj: W,
        flavor: u32,
    ) -> Result<PDBWriter<'w>> {
        Ok(PDBWriter {
            ptr: non_null(ffi::pdb_writer_from_writer(
                RustWriter::boxed(file_obj),
                flavor,
            )?),
            output: PhantomData,
        })
    }

    writer_methods!(
        "PDBWriter",
        pdb_writer_write,
        pdb_writer_flush,
        pdb_writer_close,
        pdb_writer_num_mols
    );
}

impl fmt::Debug for PDBWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PDBWriter")
            .field("num_mols", &self.num_mols())
            .finish_non_exhaustive()
    }
}

/// Parameters for `MaeWriter::get_text_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MaeWriterGetTextParams {
    /// The conformer to write, or -1 for the default.
    pub conf_id: i32,
    /// The properties to write; empty for every one that is not computed.
    pub props_list: Vec<String>,
}

impl Default for MaeWriterGetTextParams {
    fn default() -> MaeWriterGetTextParams {
        MaeWriterGetTextParams {
            conf_id: -1,
            props_list: Vec::new(),
        }
    }
}

/// Writes Maestro files. Mirrors `rdmolfiles.MaeWriter`.
///
/// RDKit's Maestro writer casts every monomer info to a PDB residue info and
/// matches each property name it writes with a regex that recurses once per
/// character, so a molecule with another kind of monomer info, or a property
/// name longer than 1000 characters, is `Error::Value`. A molecule without
/// atoms is `Error::Value`, as in Python.
pub struct MaeWriter<'w> {
    /// Never null.
    ptr: UniquePtr<ffi::MaeWriter>,
    output: PhantomData<&'w mut ()>,
}

unsafe impl Send for MaeWriter<'_> {}

/// Closes the writer, which may write to its borrowed output, so the
/// borrow must outlive it (an explicit `Drop` makes the borrow checker
/// require that).
impl Drop for MaeWriter<'_> {
    fn drop(&mut self) {}
}

impl MaeWriter<'static> {
    /// A writer to a file (`"-"` is a file name); a file that cannot be
    /// opened is `Error::Os`. Mirrors `MaeWriter(filename)`.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<MaeWriter<'static>> {
        Ok(MaeWriter {
            ptr: non_null(ffi::mae_writer_new(path_bytes(filename.as_ref()))?),
            output: PhantomData,
        })
    }
}

impl<'w> MaeWriter<'w> {
    /// A writer to `fileobj`, which it keeps until it is dropped. Mirrors
    /// `MaeWriter(fileobj)`.
    pub fn from_fileobj<W: Write + Send + 'w>(fileobj: W) -> Result<MaeWriter<'w>> {
        Ok(MaeWriter {
            ptr: non_null(ffi::mae_writer_from_writer(RustWriter::boxed(fileobj))?),
            output: PhantomData,
        })
    }

    writer_methods!(
        "MaeWriter",
        mae_writer_write,
        mae_writer_flush,
        mae_writer_close,
        mae_writer_num_mols
    );

    /// Writes only these properties of each molecule, atom, and bond, rather
    /// than every one that is not computed. Mirrors `MaeWriter.SetProps`.
    pub fn set_props<S: AsRef<str>>(&mut self, props_list: &[S]) {
        expect_ok(ffi::mae_writer_set_props(
            self.ptr.pin_mut(),
            &strings(props_list),
        ))
    }
}

impl MaeWriter<'_> {
    /// The Maestro block a writer would write for a molecule, without the
    /// file's opening block. Mirrors `MaeWriter.GetText`.
    pub fn get_text(mol: &Mol) -> Result<String> {
        MaeWriter::get_text_with_params(mol, &MaeWriterGetTextParams::default())
    }

    /// Mirrors `MaeWriter.GetText(mol, confId, props_list)`.
    pub fn get_text_with_params(mol: &Mol, params: &MaeWriterGetTextParams) -> Result<String> {
        to_string(ffi::mae_writer_get_text(
            mol.as_ffi(),
            params.conf_id,
            &params.props_list,
        )?)
    }
}

impl fmt::Debug for MaeWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MaeWriter")
            .field("num_mols", &self.num_mols())
            .finish_non_exhaustive()
    }
}
