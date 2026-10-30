// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecule suppliers: SD, SMILES, and TDT files, one record at a time.
//!
//! A supplier reads its input through a checked stream (CONVENTIONS.md U12):
//! input past RDKit's limits, or a file that changes while it is read, is
//! `Error::Value` from the call that reaches it and from every call after. A
//! record RDKit cannot read is an `Err` where Python yields `None`; RDKit logs
//! why. Iterators are named for the molecules they yield, and end after an
//! exception's `Err`. A record index past the end is `Error::Index`; negative
//! indices are not bound.

use std::fmt;
use std::io::Read;
use std::marker::PhantomData;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::suppliers::ffi;
use crate::util::path_bytes;
use crate::{Error, Result};

/// A Rust `Read` value that C++ reads through `rust_reader_read`, in place
/// of a Python file object.
pub(crate) struct RustReader(Box<dyn Read + Send>);

impl RustReader {
    /// Boxes `reader` for C++, which may keep it no longer than `'r`.
    ///
    /// The box's lifetime is erased, as `RustWriter::boxed` erases a
    /// writer's: the supplier that owns the C++ object that owns the box
    /// never lends it, and carries `'r` with an explicit `Drop`, so the
    /// borrow checker keeps `'r` alive until the supplier, and with it the
    /// box, is dropped.
    fn boxed<'r, R: Read + Send + 'r>(reader: R) -> Box<RustReader> {
        let reader: Box<dyn Read + Send + 'r> = Box::new(reader);
        // SAFETY: as above. The two types differ only in the trait object's
        // lifetime bound, so they have one layout.
        let reader: Box<dyn Read + Send + 'static> = unsafe { std::mem::transmute(reader) };
        Box::new(RustReader(reader))
    }

    /// Reads into `buf`, retrying interrupted reads; a panic is an error, so
    /// it does not unwind into C++. `Read` is a safe trait, so a count larger
    /// than `buf`, which C++ would read past the buffer by, is an error too.
    pub(crate) fn read(&mut self, buf: &mut [u8]) -> std::result::Result<usize, String> {
        loop {
            match catch_unwind(AssertUnwindSafe(|| self.0.read(buf))) {
                Ok(Ok(n)) if n > buf.len() => {
                    return Err(format!(
                        "the reader said it read {n} bytes into a buffer of {}",
                        buf.len()
                    ));
                }
                Ok(Ok(n)) => return Ok(n),
                Ok(Err(err)) if err.kind() == std::io::ErrorKind::Interrupted => continue,
                Ok(Err(err)) => return Err(err.to_string()),
                Err(_) => return Err("the reader panicked".into()),
            }
        }
    }
}

fn record(mol: UniquePtr<crate::ffi::mol::ffi::RWMol>) -> Result<Mol> {
    if mol.is_null() {
        Err(Error::Parse {
            message: "RDKit could not read the record; its log says why".into(),
        })
    } else {
        Ok(Mol::from_ffi(mol))
    }
}

/// One step of an iterator over a supplier. An exception from RDKit ends the
/// iteration after its `Err`, as it ends Python's (which reads some others as
/// `None` and goes on, though a supplier may then never advance).
fn step(
    done: &mut bool,
    mol: std::result::Result<UniquePtr<crate::ffi::mol::ffi::RWMol>, cxx::Exception>,
) -> Result<Mol> {
    match mol {
        Ok(mol) => record(mol),
        Err(err) => {
            *done = true;
            Err(err.into())
        }
    }
}

fn non_null<T: cxx::memory::UniquePtrTarget>(ptr: UniquePtr<T>) -> UniquePtr<T> {
    assert!(!ptr.is_null(), "RDKit returned a null supplier");
    ptr
}

// SAFETY (for each `unsafe impl Send` below): RDKit's suppliers keep no
// thread-local state, so one may be used from any thread; the Rust API takes
// `&mut self` to read, so one thread uses it at a time. The multithreaded
// suppliers synchronize their own threads, and a Rust reader is `Send`.

/// Parameters for `SDMolSupplier::new_with_params` and
/// `SDMolSupplier::set_data_with_params`, which take the same optional
/// arguments, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SDMolSupplierParams {
    /// Sanitizes each molecule.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed records rather than reading what it can.
    pub strict_parsing: bool,
}

impl Default for SDMolSupplierParams {
    fn default() -> SDMolSupplierParams {
        SDMolSupplierParams {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
        }
    }
}

/// The molecules of an SD file, read as they are asked for, with random
/// access. Mirrors `rdmolfiles.SDMolSupplier`. Each record's SD data become
/// the molecule's properties.
///
/// `Default` is a supplier without records, to which `set_data` gives some.
pub struct SDMolSupplier {
    /// Never null.
    ptr: UniquePtr<ffi::SDMolSupplier>,
}

unsafe impl Send for SDMolSupplier {}

impl Default for SDMolSupplier {
    /// Mirrors `SDMolSupplier()`.
    fn default() -> SDMolSupplier {
        SDMolSupplier {
            ptr: non_null(expect_ok(ffi::sd_mol_supplier_new())),
        }
    }
}

impl SDMolSupplier {
    fn raw(&self) -> &ffi::SDMolSupplier {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::SDMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier for an SD file; one that cannot be opened, or is empty, is
    /// `Error::Os`. Mirrors `SDMolSupplier(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<SDMolSupplier> {
        SDMolSupplier::new_with_params(file_name, &SDMolSupplierParams::default())
    }

    /// Mirrors `SDMolSupplier(fileName, sanitize, removeHs, strictParsing)`.
    pub fn new_with_params<P: AsRef<Path>>(
        file_name: P,
        params: &SDMolSupplierParams,
    ) -> Result<SDMolSupplier> {
        Ok(SDMolSupplier {
            ptr: non_null(ffi::sd_mol_supplier_from_file(
                path_bytes(file_name.as_ref()),
                params.sanitize,
                params.remove_hs,
                params.strict_parsing,
            )?),
        })
    }

    /// Reads `data` from its start instead, with the default parameters,
    /// which also resets `get_process_property_lists`. Mirrors
    /// `SDMolSupplier.SetData`.
    pub fn set_data(&mut self, data: &str) -> Result<()> {
        self.set_data_with_params(data, &SDMolSupplierParams::default())
    }

    /// Mirrors `SDMolSupplier.SetData(data, sanitize, removeHs,
    /// strictParsing)`.
    pub fn set_data_with_params(&mut self, data: &str, params: &SDMolSupplierParams) -> Result<()> {
        self.set_data_bytes_with_params(data.as_bytes(), params)
    }

    /// As `set_data`, for bytes. Mirrors `SDMolSupplier.SetData(data:
    /// bytes)`.
    pub fn set_data_bytes(&mut self, data: &[u8]) -> Result<()> {
        self.set_data_bytes_with_params(data, &SDMolSupplierParams::default())
    }

    /// Mirrors `SDMolSupplier.SetData(data: bytes, sanitize, removeHs,
    /// strictParsing)`.
    pub fn set_data_bytes_with_params(
        &mut self,
        data: &[u8],
        params: &SDMolSupplierParams,
    ) -> Result<()> {
        Ok(ffi::sd_mol_supplier_set_data(
            self.raw_mut(),
            data,
            params.sanitize,
            params.remove_hs,
            params.strict_parsing,
        )?)
    }

    /// The number of records, which reads the whole input once. Mirrors
    /// `len(SDMolSupplier)`.
    pub fn len(&mut self) -> Result<u32> {
        Ok(ffi::sd_mol_supplier_length(self.raw_mut())?)
    }

    /// Whether there are no records.
    pub fn is_empty(&mut self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// Whether the last record has been read, which it also is once reading
    /// the input failed; iterating then yields the error. Mirrors
    /// `SDMolSupplier.atEnd`.
    pub fn at_end(&mut self) -> bool {
        self.try_at_end().unwrap_or(true)
    }

    /// As `at_end`, with the error reading the input failed with.
    fn try_at_end(&mut self) -> Result<bool> {
        Ok(ffi::sd_mol_supplier_at_end(self.raw_mut())?)
    }

    /// Returns to the first record. Mirrors `SDMolSupplier.reset`.
    pub fn reset(&mut self) {
        expect_ok(ffi::sd_mol_supplier_reset(self.raw_mut()))
    }

    /// Record `idx`; reading continues after it. Mirrors
    /// `SDMolSupplier[idx]`.
    pub fn get(&mut self, idx: u32) -> Result<Mol> {
        record(ffi::sd_mol_supplier_get(self.raw_mut(), idx)?)
    }

    /// The text of record `index`. Mirrors `SDMolSupplier.GetItemText`.
    pub fn get_item_text(&mut self, index: u32) -> Result<String> {
        to_string(ffi::sd_mol_supplier_get_item_text(self.raw_mut(), index)?)
    }

    /// Whether `atom.prop.*` and `bond.prop.*` data become atom and bond
    /// properties. Mirrors `SDMolSupplier.GetProcessPropertyLists`.
    pub fn get_process_property_lists(&self) -> bool {
        expect_ok(ffi::sd_mol_supplier_get_process_property_lists(self.raw()))
    }

    /// Mirrors `SDMolSupplier.SetProcessPropertyLists`.
    pub fn set_process_property_lists(&mut self, val: bool) {
        expect_ok(ffi::sd_mol_supplier_set_process_property_lists(
            self.raw_mut(),
            val,
        ))
    }

    /// The molecules, from the first. Mirrors `iter(SDMolSupplier)`, which
    /// resets the supplier.
    pub fn mols(&mut self) -> SDMols<'_> {
        self.reset();
        self.mols_from_here()
    }

    /// The molecules from where the supplier stands, as Python's `next`
    /// continues.
    pub fn mols_from_here(&mut self) -> SDMols<'_> {
        SDMols {
            suppl: self,
            done: false,
        }
    }
}

impl fmt::Debug for SDMolSupplier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SDMolSupplier").finish_non_exhaustive()
    }
}

/// The molecules of an `SDMolSupplier`, from where it stands.
/// A Rust addition: Python iterates the supplier itself.
#[derive(Debug)]
pub struct SDMols<'s> {
    suppl: &'s mut SDMolSupplier,
    done: bool,
}

/// The next molecule. Mirrors `next(SDMolSupplier)`.
impl Iterator for SDMols<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match self.suppl.try_at_end() {
            Ok(false) => {}
            Ok(true) => return None,
            Err(err) => {
                self.done = true;
                return Some(Err(err));
            }
        }
        Some(step(
            &mut self.done,
            ffi::sd_mol_supplier_next(self.suppl.raw_mut()),
        ))
    }
}

/// Parameters for `ForwardSDMolSupplier::new_with_params` and
/// `ForwardSDMolSupplier::from_fileobj_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ForwardSDMolSupplierParams {
    /// Sanitizes each molecule.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed records rather than reading what it can.
    pub strict_parsing: bool,
}

impl Default for ForwardSDMolSupplierParams {
    fn default() -> ForwardSDMolSupplierParams {
        ForwardSDMolSupplierParams {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
        }
    }
}

/// The molecules of an SD stream, read once, in order. Mirrors
/// `rdmolfiles.ForwardSDMolSupplier`.
///
/// A supplier reading a Rust reader borrows it for `'r`; one reading a file
/// is `ForwardSDMolSupplier<'static>`.
pub struct ForwardSDMolSupplier<'r> {
    /// Never null.
    ptr: UniquePtr<ffi::ForwardSDMolSupplier>,
    input: PhantomData<&'r mut ()>,
}

unsafe impl Send for ForwardSDMolSupplier<'_> {}

/// Drops the borrowed input with the supplier, so the borrow must outlive
/// it (an explicit `Drop` makes the borrow checker require that).
impl Drop for ForwardSDMolSupplier<'_> {
    fn drop(&mut self) {}
}

impl ForwardSDMolSupplier<'static> {
    /// A supplier for an SD file; one that cannot be opened is `Error::Os`.
    /// Mirrors `ForwardSDMolSupplier(filename)`.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<ForwardSDMolSupplier<'static>> {
        ForwardSDMolSupplier::new_with_params(filename, &ForwardSDMolSupplierParams::default())
    }

    /// Mirrors `ForwardSDMolSupplier(filename, sanitize, removeHs,
    /// strictParsing)`.
    pub fn new_with_params<P: AsRef<Path>>(
        filename: P,
        params: &ForwardSDMolSupplierParams,
    ) -> Result<ForwardSDMolSupplier<'static>> {
        Ok(ForwardSDMolSupplier {
            ptr: non_null(ffi::forward_sd_mol_supplier_from_file(
                path_bytes(filename.as_ref()),
                params.sanitize,
                params.remove_hs,
                params.strict_parsing,
            )?),
            input: PhantomData,
        })
    }
}

impl<'r> ForwardSDMolSupplier<'r> {
    fn raw(&self) -> &ffi::ForwardSDMolSupplier {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::ForwardSDMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier reading `fileobj`, such as a decompressing reader, which
    /// it keeps until it is dropped. A read error ends the input, and the
    /// call that met it is `Error::Os`. Mirrors
    /// `ForwardSDMolSupplier(fileobj)`.
    pub fn from_fileobj<R: Read + Send + 'r>(fileobj: R) -> Result<ForwardSDMolSupplier<'r>> {
        ForwardSDMolSupplier::from_fileobj_with_params(
            fileobj,
            &ForwardSDMolSupplierParams::default(),
        )
    }

    /// Mirrors `ForwardSDMolSupplier(fileobj, sanitize, removeHs,
    /// strictParsing)`.
    pub fn from_fileobj_with_params<R: Read + Send + 'r>(
        fileobj: R,
        params: &ForwardSDMolSupplierParams,
    ) -> Result<ForwardSDMolSupplier<'r>> {
        Ok(ForwardSDMolSupplier {
            ptr: non_null(ffi::forward_sd_mol_supplier_from_reader(
                RustReader::boxed(fileobj),
                params.sanitize,
                params.remove_hs,
                params.strict_parsing,
            )?),
            input: PhantomData,
        })
    }

    /// Whether the input is exhausted, which it also is once reading it
    /// failed; iterating then yields the error. Mirrors
    /// `ForwardSDMolSupplier.atEnd`.
    pub fn at_end(&mut self) -> bool {
        self.try_at_end().unwrap_or(true)
    }

    /// As `at_end`, with the error reading the input failed with.
    fn try_at_end(&mut self) -> Result<bool> {
        Ok(ffi::forward_sd_mol_supplier_at_end(self.raw_mut())?)
    }

    /// Whether the last read met the end of the input instead of a record.
    /// Mirrors `ForwardSDMolSupplier.GetEOFHitOnRead`.
    pub fn get_eof_hit_on_read(&self) -> bool {
        expect_ok(ffi::forward_sd_mol_supplier_get_eof_hit_on_read(self.raw()))
    }

    /// Mirrors `ForwardSDMolSupplier.GetProcessPropertyLists`.
    pub fn get_process_property_lists(&self) -> bool {
        expect_ok(ffi::forward_sd_mol_supplier_get_process_property_lists(
            self.raw(),
        ))
    }

    /// Mirrors `ForwardSDMolSupplier.SetProcessPropertyLists`.
    pub fn set_process_property_lists(&mut self, val: bool) {
        expect_ok(ffi::forward_sd_mol_supplier_set_process_property_lists(
            self.raw_mut(),
            val,
        ))
    }

    /// The molecules not yet read. Mirrors `iter(ForwardSDMolSupplier)`.
    pub fn mols(&mut self) -> ForwardSDMols<'_> {
        ForwardSDMols {
            suppl: self.raw_mut(),
            done: false,
        }
    }
}

impl fmt::Debug for ForwardSDMolSupplier<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ForwardSDMolSupplier")
            .finish_non_exhaustive()
    }
}

/// The molecules of a `ForwardSDMolSupplier` not yet read.
/// A Rust addition: Python iterates the supplier itself.
pub struct ForwardSDMols<'s> {
    suppl: Pin<&'s mut ffi::ForwardSDMolSupplier>,
    done: bool,
}

impl fmt::Debug for ForwardSDMols<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ForwardSDMols")
            .field("done", &self.done)
            .finish_non_exhaustive()
    }
}

/// The next molecule. Mirrors `next(ForwardSDMolSupplier)`, which also stops
/// when a read meets the end of the input.
impl Iterator for ForwardSDMols<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match ffi::forward_sd_mol_supplier_at_end(self.suppl.as_mut()) {
            Ok(false) => {}
            Ok(true) => return None,
            Err(err) => {
                self.done = true;
                return Some(Err(err.into()));
            }
        }
        let mol = ffi::forward_sd_mol_supplier_next(self.suppl.as_mut());
        if mol.is_ok()
            && matches!(
                ffi::forward_sd_mol_supplier_at_end(self.suppl.as_mut()),
                Ok(true)
            )
            && expect_ok(ffi::forward_sd_mol_supplier_get_eof_hit_on_read(
                &self.suppl,
            ))
        {
            return None;
        }
        Some(step(&mut self.done, mol))
    }
}

/// Parameters for `SmilesMolSupplier::new_with_params`, `set_data_with_params`,
/// and `smiles_mol_supplier_from_text_with_params`, which take the same
/// optional arguments, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SmilesMolSupplierParams {
    /// The characters that separate a line's fields.
    pub delimiter: String,
    /// The field that holds the SMILES.
    pub smiles_column: i32,
    /// The field that holds the name, or -1 to name each molecule by its line
    /// number; below -1 is `Error::Value`, since RDKit would read before a
    /// line's fields.
    pub name_column: i32,
    /// Reads the first line as the fields' names.
    pub title_line: bool,
    /// Sanitizes each molecule.
    pub sanitize: bool,
}

impl Default for SmilesMolSupplierParams {
    fn default() -> SmilesMolSupplierParams {
        SmilesMolSupplierParams {
            delimiter: " ".into(),
            smiles_column: 0,
            name_column: 1,
            title_line: true,
            sanitize: true,
        }
    }
}

impl SmilesMolSupplierParams {
    fn to_ffi(&self) -> ffi::SmilesMolSupplierParamsData {
        ffi::SmilesMolSupplierParamsData {
            delimiter: self.delimiter.clone(),
            smiles_column: self.smiles_column,
            name_column: self.name_column,
            title_line: self.title_line,
            sanitize: self.sanitize,
        }
    }
}

/// The molecules of a SMILES table, one per line, with random access.
/// Mirrors `rdmolfiles.SmilesMolSupplier`; the other fields become
/// properties, named by the title line.
///
/// `Default` is a supplier without records, to which `set_data` gives some.
pub struct SmilesMolSupplier {
    /// Never null.
    ptr: UniquePtr<ffi::SmilesMolSupplier>,
}

unsafe impl Send for SmilesMolSupplier {}

impl Default for SmilesMolSupplier {
    /// Mirrors `SmilesMolSupplier()`.
    fn default() -> SmilesMolSupplier {
        SmilesMolSupplier {
            ptr: non_null(expect_ok(ffi::smiles_mol_supplier_new())),
        }
    }
}

impl SmilesMolSupplier {
    fn raw_mut(&mut self) -> Pin<&mut ffi::SmilesMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier for a SMILES file; one that cannot be opened, or is empty,
    /// is `Error::Os`. Mirrors `SmilesMolSupplier(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<SmilesMolSupplier> {
        SmilesMolSupplier::new_with_params(file_name, &SmilesMolSupplierParams::default())
    }

    /// Mirrors `SmilesMolSupplier(fileName, delimiter, smilesColumn,
    /// nameColumn, titleLine, sanitize)`.
    pub fn new_with_params<P: AsRef<Path>>(
        file_name: P,
        params: &SmilesMolSupplierParams,
    ) -> Result<SmilesMolSupplier> {
        Ok(SmilesMolSupplier {
            ptr: non_null(ffi::smiles_mol_supplier_from_file(
                path_bytes(file_name.as_ref()),
                &params.to_ffi(),
            )?),
        })
    }

    /// Reads `data` from its start instead. Mirrors
    /// `SmilesMolSupplier.SetData`.
    pub fn set_data(&mut self, data: &str) -> Result<()> {
        self.set_data_with_params(data, &SmilesMolSupplierParams::default())
    }

    /// Mirrors `SmilesMolSupplier.SetData(data, delimiter, smilesColumn,
    /// nameColumn, titleLine, sanitize)`.
    pub fn set_data_with_params(
        &mut self,
        data: &str,
        params: &SmilesMolSupplierParams,
    ) -> Result<()> {
        Ok(ffi::smiles_mol_supplier_set_data(
            self.raw_mut(),
            data,
            &params.to_ffi(),
        )?)
    }

    /// The number of records, which reads the whole input once. Mirrors
    /// `len(SmilesMolSupplier)`.
    pub fn len(&mut self) -> Result<u32> {
        Ok(ffi::smiles_mol_supplier_length(self.raw_mut())?)
    }

    /// Whether there are no records.
    pub fn is_empty(&mut self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// Whether the last record has been read, which it also is once reading
    /// the input failed; iterating then yields the error.
    pub fn at_end(&mut self) -> bool {
        self.try_at_end().unwrap_or(true)
    }

    /// As `at_end`, with the error reading the input failed with.
    fn try_at_end(&mut self) -> Result<bool> {
        Ok(ffi::smiles_mol_supplier_at_end(self.raw_mut())?)
    }

    /// Returns to the first record. Mirrors `SmilesMolSupplier.reset`.
    pub fn reset(&mut self) {
        expect_ok(ffi::smiles_mol_supplier_reset(self.raw_mut()))
    }

    /// Record `idx`; reading continues after it. Mirrors
    /// `SmilesMolSupplier[idx]`.
    pub fn get(&mut self, idx: u32) -> Result<Mol> {
        record(ffi::smiles_mol_supplier_get(self.raw_mut(), idx)?)
    }

    /// The line of record `index`. Mirrors `SmilesMolSupplier.GetItemText`.
    pub fn get_item_text(&mut self, index: u32) -> Result<String> {
        to_string(ffi::smiles_mol_supplier_get_item_text(
            self.raw_mut(),
            index,
        )?)
    }

    /// The molecules, from the first. Mirrors `iter(SmilesMolSupplier)`,
    /// which resets the supplier.
    pub fn mols(&mut self) -> SmilesMols<'_> {
        self.reset();
        self.mols_from_here()
    }

    /// The molecules from where the supplier stands, as Python's `next`
    /// continues.
    pub fn mols_from_here(&mut self) -> SmilesMols<'_> {
        SmilesMols {
            suppl: self,
            done: false,
        }
    }
}

impl fmt::Debug for SmilesMolSupplier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SmilesMolSupplier").finish_non_exhaustive()
    }
}

/// A supplier for SMILES text. Mirrors `rdmolfiles.SmilesMolSupplierFromText`.
pub fn smiles_mol_supplier_from_text(text: &str) -> Result<SmilesMolSupplier> {
    smiles_mol_supplier_from_text_with_params(text, &SmilesMolSupplierParams::default())
}

/// Mirrors `rdmolfiles.SmilesMolSupplierFromText(text, delimiter,
/// smilesColumn, nameColumn, titleLine, sanitize)`.
pub fn smiles_mol_supplier_from_text_with_params(
    text: &str,
    params: &SmilesMolSupplierParams,
) -> Result<SmilesMolSupplier> {
    let mut suppl = SmilesMolSupplier::default();
    suppl.set_data_with_params(text, params)?;
    Ok(suppl)
}

/// The molecules of a `SmilesMolSupplier`, from where it stands.
/// A Rust addition: Python iterates the supplier itself.
#[derive(Debug)]
pub struct SmilesMols<'s> {
    suppl: &'s mut SmilesMolSupplier,
    done: bool,
}

/// The next molecule. Mirrors `next(SmilesMolSupplier)`.
impl Iterator for SmilesMols<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match self.suppl.try_at_end() {
            Ok(false) => {}
            Ok(true) => return None,
            Err(err) => {
                self.done = true;
                return Some(Err(err));
            }
        }
        Some(step(
            &mut self.done,
            ffi::smiles_mol_supplier_next(self.suppl.raw_mut()),
        ))
    }
}

/// Parameters for `TDTMolSupplier::new_with_params` and
/// `TDTMolSupplier::set_data_with_params`, which take the same optional
/// arguments, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct TDTMolSupplierParams {
    /// The record whose value names each molecule; empty for none.
    pub name_record: String,
    /// The conformer id to give 2D coordinates, or -1 to skip them.
    pub conf_id_2d: i32,
    /// The conformer id to give 3D coordinates, or -1 to skip them.
    pub conf_id_3d: i32,
    /// Sanitizes each molecule.
    pub sanitize: bool,
}

impl Default for TDTMolSupplierParams {
    fn default() -> TDTMolSupplierParams {
        TDTMolSupplierParams {
            name_record: String::new(),
            conf_id_2d: -1,
            conf_id_3d: -1,
            sanitize: true,
        }
    }
}

impl TDTMolSupplierParams {
    fn to_ffi(&self) -> ffi::TDTMolSupplierParamsData {
        ffi::TDTMolSupplierParamsData {
            name_record: self.name_record.clone(),
            conf_id_2d: self.conf_id_2d,
            conf_id_3d: self.conf_id_3d,
            sanitize: self.sanitize,
        }
    }
}

/// The molecules of a Daylight TDT file, with random access. Mirrors
/// `rdmolfiles.TDTMolSupplier`.
///
/// The supplier reads its input into memory and finds its records, the lines
/// that start with `$SMI<`, itself; RDKit reads each record. RDKit's own
/// supplier trusts the record positions it collects, reading past them for
/// some inputs, so records are read as it reads them once it has counted
/// them: Python's supplier may differ on input whose records lack their `|`
/// lines. `Default` is a supplier without records, to which `set_data` gives
/// some.
pub struct TDTMolSupplier {
    /// Never null.
    ptr: UniquePtr<ffi::TDTMolSupplier>,
}

unsafe impl Send for TDTMolSupplier {}

impl Default for TDTMolSupplier {
    /// Mirrors `TDTMolSupplier()`.
    fn default() -> TDTMolSupplier {
        TDTMolSupplier {
            ptr: non_null(expect_ok(ffi::tdt_mol_supplier_new())),
        }
    }
}

impl TDTMolSupplier {
    fn raw_mut(&mut self) -> Pin<&mut ffi::TDTMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier for a TDT file; one that cannot be opened, or is empty, is
    /// `Error::Os`. Mirrors `TDTMolSupplier(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<TDTMolSupplier> {
        TDTMolSupplier::new_with_params(file_name, &TDTMolSupplierParams::default())
    }

    /// Mirrors `TDTMolSupplier(fileName, nameRecord, confId2D, confId3D,
    /// sanitize)`.
    pub fn new_with_params<P: AsRef<Path>>(
        file_name: P,
        params: &TDTMolSupplierParams,
    ) -> Result<TDTMolSupplier> {
        Ok(TDTMolSupplier {
            ptr: non_null(ffi::tdt_mol_supplier_from_file(
                path_bytes(file_name.as_ref()),
                &params.to_ffi(),
            )?),
        })
    }

    /// Reads `data` from its start instead. Mirrors `TDTMolSupplier.SetData`.
    pub fn set_data(&mut self, data: &str) -> Result<()> {
        self.set_data_with_params(data, &TDTMolSupplierParams::default())
    }

    /// Mirrors `TDTMolSupplier.SetData(data, nameRecord, confId2D, confId3D,
    /// sanitize)`.
    pub fn set_data_with_params(
        &mut self,
        data: &str,
        params: &TDTMolSupplierParams,
    ) -> Result<()> {
        Ok(ffi::tdt_mol_supplier_set_data(
            self.raw_mut(),
            data,
            &params.to_ffi(),
        )?)
    }

    /// The number of records. Mirrors `len(TDTMolSupplier)`.
    pub fn len(&mut self) -> Result<u32> {
        Ok(ffi::tdt_mol_supplier_length(self.raw_mut())?)
    }

    /// Whether there are no records.
    pub fn is_empty(&mut self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// Whether the last record has been read.
    pub fn at_end(&mut self) -> bool {
        expect_ok(ffi::tdt_mol_supplier_at_end(self.raw_mut()))
    }

    /// As `at_end`: a TDT supplier reads text it holds, so it cannot fail.
    fn try_at_end(&mut self) -> Result<bool> {
        Ok(self.at_end())
    }

    /// Returns to the first record. Mirrors `TDTMolSupplier.reset`.
    pub fn reset(&mut self) {
        expect_ok(ffi::tdt_mol_supplier_reset(self.raw_mut()))
    }

    /// Record `idx`; reading continues after it. Mirrors
    /// `TDTMolSupplier[idx]`.
    pub fn get(&mut self, idx: u32) -> Result<Mol> {
        record(ffi::tdt_mol_supplier_get(self.raw_mut(), idx)?)
    }

    /// The text of record `index`. Mirrors `TDTMolSupplier.GetItemText`.
    pub fn get_item_text(&mut self, index: u32) -> Result<String> {
        to_string(ffi::tdt_mol_supplier_get_item_text(self.raw_mut(), index)?)
    }

    /// The molecules, from the first. Mirrors `iter(TDTMolSupplier)`, which
    /// resets the supplier.
    pub fn mols(&mut self) -> TDTMols<'_> {
        self.reset();
        self.mols_from_here()
    }

    /// The molecules from where the supplier stands, as Python's `next`
    /// continues.
    pub fn mols_from_here(&mut self) -> TDTMols<'_> {
        TDTMols {
            suppl: self,
            done: false,
        }
    }
}

impl fmt::Debug for TDTMolSupplier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TDTMolSupplier").finish_non_exhaustive()
    }
}

/// The molecules of a `TDTMolSupplier`, from where it stands.
/// A Rust addition: Python iterates the supplier itself.
#[derive(Debug)]
pub struct TDTMols<'s> {
    suppl: &'s mut TDTMolSupplier,
    done: bool,
}

/// The next molecule. Mirrors `next(TDTMolSupplier)`.
impl Iterator for TDTMols<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match self.suppl.try_at_end() {
            Ok(false) => {}
            Ok(true) => return None,
            Err(err) => {
                self.done = true;
                return Some(Err(err));
            }
        }
        Some(step(
            &mut self.done,
            ffi::tdt_mol_supplier_next(self.suppl.raw_mut()),
        ))
    }
}

fn mt_params(
    num_writer_threads: u32,
    size_input_queue: usize,
    size_output_queue: usize,
) -> ffi::MultithreadedParamsData {
    ffi::MultithreadedParamsData {
        num_writer_threads,
        size_input_queue,
        size_output_queue,
    }
}

/// The next record of a multithreaded supplier: `None` at the end.
fn mt_record(
    done: &mut bool,
    mol: std::result::Result<UniquePtr<crate::ffi::mol::ffi::RWMol>, cxx::Exception>,
    id: u32,
) -> Option<Result<Mol>> {
    match mol {
        Ok(mol) if mol.is_null() && id == 0 => None,
        mol => Some(step(done, mol)),
    }
}

/// Parameters for `MultithreadedSDMolSupplier::new_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MultithreadedSDMolSupplierParams {
    /// Sanitizes each molecule.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed records rather than reading what it can.
    pub strict_parsing: bool,
    /// The parsing threads: 0 for one per core. More than 1024 is
    /// `Error::Value`: RDKit starts every one at the first read. If the
    /// system cannot start them all, that read is `Error::Runtime`, and so is
    /// every read after it.
    pub num_writer_threads: u32,
    /// The records read ahead; 0, which never lets RDKit's threads finish, and
    /// sizes above 2^20 are `Error::Value`: RDKit 2025.09 allocates each queue
    /// whole, and dereferences null if that allocation fails.
    pub size_input_queue: usize,
    /// The molecules parsed ahead, limited as `size_input_queue`.
    pub size_output_queue: usize,
}

impl Default for MultithreadedSDMolSupplierParams {
    fn default() -> MultithreadedSDMolSupplierParams {
        MultithreadedSDMolSupplierParams {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
            num_writer_threads: 1,
            size_input_queue: 5,
            size_output_queue: 5,
        }
    }
}

/// The molecules of an SD file, parsed on threads of its own, in the order
/// they finish. Mirrors `rdmolfiles.MultithreadedSDMolSupplier`.
///
/// The threads start at the first molecule asked for and stop when the
/// supplier is dropped. `Default` is a supplier without records.
pub struct MultithreadedSDMolSupplier {
    /// Never null.
    ptr: UniquePtr<ffi::MultithreadedSDMolSupplier>,
}

unsafe impl Send for MultithreadedSDMolSupplier {}

impl Default for MultithreadedSDMolSupplier {
    /// Mirrors `MultithreadedSDMolSupplier()`.
    fn default() -> MultithreadedSDMolSupplier {
        MultithreadedSDMolSupplier {
            ptr: non_null(expect_ok(ffi::mt_sd_mol_supplier_new())),
        }
    }
}

impl MultithreadedSDMolSupplier {
    fn raw(&self) -> &ffi::MultithreadedSDMolSupplier {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::MultithreadedSDMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier for an SD file; one that cannot be opened, or is empty, is
    /// `Error::Os`. Mirrors `MultithreadedSDMolSupplier(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<MultithreadedSDMolSupplier> {
        MultithreadedSDMolSupplier::new_with_params(
            file_name,
            &MultithreadedSDMolSupplierParams::default(),
        )
    }

    /// Mirrors `MultithreadedSDMolSupplier(fileName, sanitize, removeHs,
    /// strictParsing, numWriterThreads, sizeInputQueue, sizeOutputQueue)`.
    pub fn new_with_params<P: AsRef<Path>>(
        file_name: P,
        params: &MultithreadedSDMolSupplierParams,
    ) -> Result<MultithreadedSDMolSupplier> {
        Ok(MultithreadedSDMolSupplier {
            ptr: non_null(ffi::mt_sd_mol_supplier_from_file(
                path_bytes(file_name.as_ref()),
                params.sanitize,
                params.remove_hs,
                params.strict_parsing,
                &mt_params(
                    params.num_writer_threads,
                    params.size_input_queue,
                    params.size_output_queue,
                ),
            )?),
        })
    }

    /// Whether every record has been returned, which they also have once
    /// reading the input failed; iterating then yields the error. Mirrors
    /// `MultithreadedSDMolSupplier.atEnd`.
    pub fn at_end(&mut self) -> bool {
        self.try_at_end().unwrap_or(true)
    }

    /// As `at_end`, with the error reading the input failed with.
    fn try_at_end(&mut self) -> Result<bool> {
        Ok(ffi::mt_sd_mol_supplier_at_end(self.raw_mut())?)
    }

    /// The 1-based record number of the last molecule returned, or 0.
    /// Mirrors `MultithreadedSDMolSupplier.GetLastRecordId`.
    pub fn get_last_record_id(&self) -> u32 {
        expect_ok(ffi::mt_sd_mol_supplier_get_last_record_id(self.raw()))
    }

    /// The text of the last record returned. Mirrors
    /// `MultithreadedSDMolSupplier.GetLastItemText`.
    pub fn get_last_item_text(&self) -> Result<String> {
        to_string(expect_ok(ffi::mt_sd_mol_supplier_get_last_item_text(
            self.raw(),
        )))
    }

    /// Mirrors `MultithreadedSDMolSupplier.GetProcessPropertyLists`.
    pub fn get_process_property_lists(&self) -> bool {
        expect_ok(ffi::mt_sd_mol_supplier_get_process_property_lists(
            self.raw(),
        ))
    }

    /// Once reading has started this is `Error::Value`: RDKit's parsing
    /// threads read the setting without a lock. Mirrors
    /// `MultithreadedSDMolSupplier.SetProcessPropertyLists`.
    pub fn set_process_property_lists(&mut self, val: bool) -> Result<()> {
        Ok(ffi::mt_sd_mol_supplier_set_process_property_lists(
            self.raw_mut(),
            val,
        )?)
    }

    /// The molecules not yet returned. Mirrors
    /// `iter(MultithreadedSDMolSupplier)`.
    pub fn mols(&mut self) -> MultithreadedSDMols<'_> {
        MultithreadedSDMols {
            suppl: self,
            done: false,
        }
    }
}

impl fmt::Debug for MultithreadedSDMolSupplier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MultithreadedSDMolSupplier")
            .field("last_record_id", &self.get_last_record_id())
            .finish_non_exhaustive()
    }
}

/// The molecules of a `MultithreadedSDMolSupplier` not yet returned.
/// A Rust addition: Python iterates the supplier itself.
#[derive(Debug)]
pub struct MultithreadedSDMols<'s> {
    suppl: &'s mut MultithreadedSDMolSupplier,
    done: bool,
}

/// The next molecule to finish. Mirrors `next(MultithreadedSDMolSupplier)`.
impl Iterator for MultithreadedSDMols<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match self.suppl.try_at_end() {
            Ok(false) => {}
            Ok(true) => return None,
            Err(err) => {
                self.done = true;
                return Some(Err(err));
            }
        }
        let mut id = 0;
        let mol = ffi::mt_sd_mol_supplier_next(self.suppl.raw_mut(), &mut id);
        mt_record(&mut self.done, mol, id)
    }
}

/// Parameters for `MultithreadedSmilesMolSupplier::new_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MultithreadedSmilesMolSupplierParams {
    /// The characters that separate a line's fields.
    pub delimiter: String,
    /// The field that holds the SMILES.
    pub smiles_column: i32,
    /// As `SmilesMolSupplierParams::name_column`.
    pub name_column: i32,
    /// Reads the first line as the fields' names.
    pub title_line: bool,
    /// Sanitizes each molecule.
    pub sanitize: bool,
    /// As `MultithreadedSDMolSupplierParams::num_writer_threads`.
    pub num_writer_threads: u32,
    /// As `MultithreadedSDMolSupplierParams::size_input_queue`.
    pub size_input_queue: usize,
    /// As `MultithreadedSDMolSupplierParams::size_output_queue`.
    pub size_output_queue: usize,
}

impl Default for MultithreadedSmilesMolSupplierParams {
    fn default() -> MultithreadedSmilesMolSupplierParams {
        MultithreadedSmilesMolSupplierParams {
            delimiter: " \t".into(),
            smiles_column: 0,
            name_column: 1,
            title_line: true,
            sanitize: true,
            num_writer_threads: 1,
            size_input_queue: 5,
            size_output_queue: 5,
        }
    }
}

/// The molecules of a SMILES file, parsed on threads of its own, in the order
/// they finish. Mirrors `rdmolfiles.MultithreadedSmilesMolSupplier`.
///
/// RDKit reads the end of the file as a last, empty record, which is an
/// `Err`, as Python yields `None` for it. The threads start at the first
/// molecule asked for and stop when the supplier is dropped. `Default` is a
/// supplier without records.
pub struct MultithreadedSmilesMolSupplier {
    /// Never null.
    ptr: UniquePtr<ffi::MultithreadedSmilesMolSupplier>,
}

unsafe impl Send for MultithreadedSmilesMolSupplier {}

impl Default for MultithreadedSmilesMolSupplier {
    /// Mirrors `MultithreadedSmilesMolSupplier()`.
    fn default() -> MultithreadedSmilesMolSupplier {
        MultithreadedSmilesMolSupplier {
            ptr: non_null(expect_ok(ffi::mt_smiles_mol_supplier_new())),
        }
    }
}

impl MultithreadedSmilesMolSupplier {
    fn raw(&self) -> &ffi::MultithreadedSmilesMolSupplier {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::MultithreadedSmilesMolSupplier> {
        self.ptr.pin_mut()
    }

    /// A supplier for a SMILES file; one that cannot be opened, or is empty,
    /// is `Error::Os`. Mirrors `MultithreadedSmilesMolSupplier(fileName)`.
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<MultithreadedSmilesMolSupplier> {
        MultithreadedSmilesMolSupplier::new_with_params(
            file_name,
            &MultithreadedSmilesMolSupplierParams::default(),
        )
    }

    /// Mirrors `MultithreadedSmilesMolSupplier(fileName, delimiter,
    /// smilesColumn, nameColumn, titleLine, sanitize, numWriterThreads,
    /// sizeInputQueue, sizeOutputQueue)`.
    pub fn new_with_params<P: AsRef<Path>>(
        file_name: P,
        params: &MultithreadedSmilesMolSupplierParams,
    ) -> Result<MultithreadedSmilesMolSupplier> {
        let smiles = ffi::SmilesMolSupplierParamsData {
            delimiter: params.delimiter.clone(),
            smiles_column: params.smiles_column,
            name_column: params.name_column,
            title_line: params.title_line,
            sanitize: params.sanitize,
        };
        Ok(MultithreadedSmilesMolSupplier {
            ptr: non_null(ffi::mt_smiles_mol_supplier_from_file(
                path_bytes(file_name.as_ref()),
                &smiles,
                &mt_params(
                    params.num_writer_threads,
                    params.size_input_queue,
                    params.size_output_queue,
                ),
            )?),
        })
    }

    /// Whether every record has been returned, which they also have once
    /// reading the input failed; iterating then yields the error. Mirrors
    /// `MultithreadedSmilesMolSupplier.atEnd`.
    pub fn at_end(&mut self) -> bool {
        self.try_at_end().unwrap_or(true)
    }

    /// As `at_end`, with the error reading the input failed with.
    fn try_at_end(&mut self) -> Result<bool> {
        Ok(ffi::mt_smiles_mol_supplier_at_end(self.raw_mut())?)
    }

    /// The 1-based record number of the last molecule returned, or 0.
    /// Mirrors `MultithreadedSmilesMolSupplier.GetLastRecordId`.
    pub fn get_last_record_id(&self) -> u32 {
        expect_ok(ffi::mt_smiles_mol_supplier_get_last_record_id(self.raw()))
    }

    /// The line of the last record returned. Mirrors
    /// `MultithreadedSmilesMolSupplier.GetLastItemText`.
    pub fn get_last_item_text(&self) -> Result<String> {
        to_string(expect_ok(ffi::mt_smiles_mol_supplier_get_last_item_text(
            self.raw(),
        )))
    }

    /// The molecules not yet returned. Mirrors
    /// `iter(MultithreadedSmilesMolSupplier)`.
    pub fn mols(&mut self) -> MultithreadedSmilesMols<'_> {
        MultithreadedSmilesMols {
            suppl: self,
            done: false,
        }
    }
}

impl fmt::Debug for MultithreadedSmilesMolSupplier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MultithreadedSmilesMolSupplier")
            .field("last_record_id", &self.get_last_record_id())
            .finish_non_exhaustive()
    }
}

/// The molecules of a `MultithreadedSmilesMolSupplier` not yet returned.
/// A Rust addition: Python iterates the supplier itself.
#[derive(Debug)]
pub struct MultithreadedSmilesMols<'s> {
    suppl: &'s mut MultithreadedSmilesMolSupplier,
    done: bool,
}

/// The next molecule to finish. Mirrors
/// `next(MultithreadedSmilesMolSupplier)`.
impl Iterator for MultithreadedSmilesMols<'_> {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if self.done {
            return None;
        }
        match self.suppl.try_at_end() {
            Ok(false) => {}
            Ok(true) => return None,
            Err(err) => {
                self.done = true;
                return Some(Err(err));
            }
        }
        let mut id = 0;
        let mol = ffi::mt_smiles_mol_supplier_next(self.suppl.raw_mut(), &mut id);
        mt_record(&mut self.done, mol, id)
    }
}
