// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::chem::rdmolfiles::RustWriter;

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `SmilesWriter`'s constructor arguments.
    struct SmilesWriterParamsData {
        delimiter: String,
        name_header: String,
        include_header: bool,
        isomeric_smiles: bool,
        kekule_smiles: bool,
    }

    extern "Rust" {
        type RustWriter;
        fn rust_writer_write(writer: &mut RustWriter, buf: &[u8]) -> Result<()>;
        fn rust_writer_flush(writer: &mut RustWriter) -> Result<()>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/writers.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type SDWriter;
        type SmilesWriter;
        type TDTWriter;
        type PDBWriter;
        type MaeWriter;

        fn sd_writer_new(path: &[u8]) -> Result<UniquePtr<SDWriter>>;
        fn sd_writer_from_writer(writer: Box<RustWriter>) -> Result<UniquePtr<SDWriter>>;
        fn sd_writer_write(w: Pin<&mut SDWriter>, mol: &RWMol, conf_id: i32) -> Result<()>;
        fn sd_writer_flush(w: Pin<&mut SDWriter>) -> Result<()>;
        fn sd_writer_close(w: Pin<&mut SDWriter>) -> Result<()>;
        fn sd_writer_num_mols(w: &SDWriter) -> Result<u32>;
        fn sd_writer_set_props(w: Pin<&mut SDWriter>, props: &Vec<String>) -> Result<()>;
        fn sd_writer_get_force_v3000(w: &SDWriter) -> Result<bool>;
        fn sd_writer_set_force_v3000(w: Pin<&mut SDWriter>, val: bool) -> Result<()>;
        fn sd_writer_get_kekulize(w: &SDWriter) -> Result<bool>;
        fn sd_writer_set_kekulize(w: Pin<&mut SDWriter>, val: bool) -> Result<()>;
        fn sd_writer_get_text(
            mol: &RWMol,
            conf_id: i32,
            kekulize: bool,
            force_v3000: bool,
            molid: i32,
        ) -> Result<UniquePtr<CxxString>>;

        fn smiles_writer_new(
            path: &[u8],
            params: &SmilesWriterParamsData,
        ) -> Result<UniquePtr<SmilesWriter>>;
        fn smiles_writer_from_writer(
            writer: Box<RustWriter>,
            params: &SmilesWriterParamsData,
        ) -> Result<UniquePtr<SmilesWriter>>;
        fn smiles_writer_write(w: Pin<&mut SmilesWriter>, mol: &RWMol, conf_id: i32) -> Result<()>;
        fn smiles_writer_flush(w: Pin<&mut SmilesWriter>) -> Result<()>;
        fn smiles_writer_close(w: Pin<&mut SmilesWriter>) -> Result<()>;
        fn smiles_writer_num_mols(w: &SmilesWriter) -> Result<u32>;
        fn smiles_writer_set_props(w: Pin<&mut SmilesWriter>, props: &Vec<String>) -> Result<()>;

        fn tdt_writer_new(path: &[u8]) -> Result<UniquePtr<TDTWriter>>;
        fn tdt_writer_from_writer(writer: Box<RustWriter>) -> Result<UniquePtr<TDTWriter>>;
        fn tdt_writer_write(w: Pin<&mut TDTWriter>, mol: &RWMol, conf_id: i32) -> Result<()>;
        fn tdt_writer_flush(w: Pin<&mut TDTWriter>) -> Result<()>;
        fn tdt_writer_close(w: Pin<&mut TDTWriter>) -> Result<()>;
        fn tdt_writer_num_mols(w: &TDTWriter) -> Result<u32>;
        fn tdt_writer_set_props(w: Pin<&mut TDTWriter>, props: &Vec<String>) -> Result<()>;
        fn tdt_writer_get_num_digits(w: &TDTWriter) -> Result<u32>;
        fn tdt_writer_set_num_digits(w: Pin<&mut TDTWriter>, num_digits: u32) -> Result<()>;
        fn tdt_writer_get_write_2d(w: &TDTWriter) -> Result<bool>;
        fn tdt_writer_set_write_2d(w: Pin<&mut TDTWriter>, state: bool) -> Result<()>;
        fn tdt_writer_get_write_names(w: &TDTWriter) -> Result<bool>;
        fn tdt_writer_set_write_names(w: Pin<&mut TDTWriter>, state: bool) -> Result<()>;

        fn pdb_writer_new(path: &[u8], flavor: u32) -> Result<UniquePtr<PDBWriter>>;
        fn pdb_writer_from_writer(
            writer: Box<RustWriter>,
            flavor: u32,
        ) -> Result<UniquePtr<PDBWriter>>;
        fn pdb_writer_write(w: Pin<&mut PDBWriter>, mol: &RWMol, conf_id: i32) -> Result<()>;
        fn pdb_writer_flush(w: Pin<&mut PDBWriter>) -> Result<()>;
        fn pdb_writer_close(w: Pin<&mut PDBWriter>) -> Result<()>;
        fn pdb_writer_num_mols(w: &PDBWriter) -> Result<u32>;

        fn mae_writer_new(path: &[u8]) -> Result<UniquePtr<MaeWriter>>;
        fn mae_writer_from_writer(writer: Box<RustWriter>) -> Result<UniquePtr<MaeWriter>>;
        fn mae_writer_write(w: Pin<&mut MaeWriter>, mol: &RWMol, conf_id: i32) -> Result<()>;
        fn mae_writer_flush(w: Pin<&mut MaeWriter>) -> Result<()>;
        fn mae_writer_close(w: Pin<&mut MaeWriter>) -> Result<()>;
        fn mae_writer_num_mols(w: &MaeWriter) -> Result<u32>;
        fn mae_writer_set_props(w: Pin<&mut MaeWriter>, props: &Vec<String>) -> Result<()>;
        fn mae_writer_get_text(
            mol: &RWMol,
            conf_id: i32,
            props: &Vec<String>,
        ) -> Result<UniquePtr<CxxString>>;
    }
}

fn rust_writer_write(writer: &mut RustWriter, buf: &[u8]) -> Result<(), String> {
    writer.write_all(buf)
}

fn rust_writer_flush(writer: &mut RustWriter) -> Result<(), String> {
    writer.flush()
}
