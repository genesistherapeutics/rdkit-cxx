// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::chem::rdmolfiles::RustReader;

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `SmilesMolSupplierParams`.
    struct SmilesMolSupplierParamsData {
        delimiter: String,
        smiles_column: i32,
        name_column: i32,
        title_line: bool,
        sanitize: bool,
    }

    /// `TDTMolSupplierParams`.
    struct TDTMolSupplierParamsData {
        name_record: String,
        conf_id_2d: i32,
        conf_id_3d: i32,
        sanitize: bool,
    }

    /// `MultithreadedMolSupplier::Parameters`.
    struct MultithreadedParamsData {
        num_writer_threads: u32,
        size_input_queue: usize,
        size_output_queue: usize,
    }

    extern "Rust" {
        type RustReader;
        fn rust_reader_read(reader: &mut RustReader, buf: &mut [u8]) -> Result<usize>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/suppliers.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type SDMolSupplier;
        type ForwardSDMolSupplier;
        type SmilesMolSupplier;
        type TDTMolSupplier;
        type MultithreadedSDMolSupplier;
        type MultithreadedSmilesMolSupplier;

        fn sd_mol_supplier_new() -> Result<UniquePtr<SDMolSupplier>>;
        fn sd_mol_supplier_from_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<SDMolSupplier>>;
        fn sd_mol_supplier_set_data(
            s: Pin<&mut SDMolSupplier>,
            data: &[u8],
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<()>;
        fn sd_mol_supplier_next(s: Pin<&mut SDMolSupplier>) -> Result<UniquePtr<RWMol>>;
        fn sd_mol_supplier_get(s: Pin<&mut SDMolSupplier>, idx: u32) -> Result<UniquePtr<RWMol>>;
        fn sd_mol_supplier_get_item_text(
            s: Pin<&mut SDMolSupplier>,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn sd_mol_supplier_length(s: Pin<&mut SDMolSupplier>) -> Result<u32>;
        fn sd_mol_supplier_at_end(s: Pin<&mut SDMolSupplier>) -> Result<bool>;
        fn sd_mol_supplier_reset(s: Pin<&mut SDMolSupplier>) -> Result<()>;
        fn sd_mol_supplier_get_process_property_lists(s: &SDMolSupplier) -> Result<bool>;
        fn sd_mol_supplier_set_process_property_lists(
            s: Pin<&mut SDMolSupplier>,
            val: bool,
        ) -> Result<()>;

        fn forward_sd_mol_supplier_from_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<ForwardSDMolSupplier>>;
        fn forward_sd_mol_supplier_from_reader(
            reader: Box<RustReader>,
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<ForwardSDMolSupplier>>;
        fn forward_sd_mol_supplier_next(
            s: Pin<&mut ForwardSDMolSupplier>,
        ) -> Result<UniquePtr<RWMol>>;
        fn forward_sd_mol_supplier_at_end(s: Pin<&mut ForwardSDMolSupplier>) -> Result<bool>;
        fn forward_sd_mol_supplier_get_eof_hit_on_read(s: &ForwardSDMolSupplier) -> Result<bool>;
        fn forward_sd_mol_supplier_get_process_property_lists(
            s: &ForwardSDMolSupplier,
        ) -> Result<bool>;
        fn forward_sd_mol_supplier_set_process_property_lists(
            s: Pin<&mut ForwardSDMolSupplier>,
            val: bool,
        ) -> Result<()>;

        fn smiles_mol_supplier_new() -> Result<UniquePtr<SmilesMolSupplier>>;
        fn smiles_mol_supplier_from_file(
            path: &[u8],
            params: &SmilesMolSupplierParamsData,
        ) -> Result<UniquePtr<SmilesMolSupplier>>;
        fn smiles_mol_supplier_set_data(
            s: Pin<&mut SmilesMolSupplier>,
            data: &str,
            params: &SmilesMolSupplierParamsData,
        ) -> Result<()>;
        fn smiles_mol_supplier_next(s: Pin<&mut SmilesMolSupplier>) -> Result<UniquePtr<RWMol>>;
        fn smiles_mol_supplier_get(
            s: Pin<&mut SmilesMolSupplier>,
            idx: u32,
        ) -> Result<UniquePtr<RWMol>>;
        fn smiles_mol_supplier_get_item_text(
            s: Pin<&mut SmilesMolSupplier>,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn smiles_mol_supplier_length(s: Pin<&mut SmilesMolSupplier>) -> Result<u32>;
        fn smiles_mol_supplier_at_end(s: Pin<&mut SmilesMolSupplier>) -> Result<bool>;
        fn smiles_mol_supplier_reset(s: Pin<&mut SmilesMolSupplier>) -> Result<()>;

        fn tdt_mol_supplier_new() -> Result<UniquePtr<TDTMolSupplier>>;
        fn tdt_mol_supplier_from_file(
            path: &[u8],
            params: &TDTMolSupplierParamsData,
        ) -> Result<UniquePtr<TDTMolSupplier>>;
        fn tdt_mol_supplier_set_data(
            s: Pin<&mut TDTMolSupplier>,
            data: &str,
            params: &TDTMolSupplierParamsData,
        ) -> Result<()>;
        fn tdt_mol_supplier_next(s: Pin<&mut TDTMolSupplier>) -> Result<UniquePtr<RWMol>>;
        fn tdt_mol_supplier_get(s: Pin<&mut TDTMolSupplier>, idx: u32) -> Result<UniquePtr<RWMol>>;
        fn tdt_mol_supplier_get_item_text(
            s: Pin<&mut TDTMolSupplier>,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn tdt_mol_supplier_length(s: Pin<&mut TDTMolSupplier>) -> Result<u32>;
        fn tdt_mol_supplier_at_end(s: Pin<&mut TDTMolSupplier>) -> Result<bool>;
        fn tdt_mol_supplier_reset(s: Pin<&mut TDTMolSupplier>) -> Result<()>;

        fn mt_sd_mol_supplier_new() -> Result<UniquePtr<MultithreadedSDMolSupplier>>;
        fn mt_sd_mol_supplier_from_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
            params: &MultithreadedParamsData,
        ) -> Result<UniquePtr<MultithreadedSDMolSupplier>>;
        fn mt_sd_mol_supplier_next(
            s: Pin<&mut MultithreadedSDMolSupplier>,
            record: &mut u32,
        ) -> Result<UniquePtr<RWMol>>;
        fn mt_sd_mol_supplier_at_end(s: Pin<&mut MultithreadedSDMolSupplier>) -> Result<bool>;
        fn mt_sd_mol_supplier_get_last_record_id(s: &MultithreadedSDMolSupplier) -> Result<u32>;
        fn mt_sd_mol_supplier_get_last_item_text(
            s: &MultithreadedSDMolSupplier,
        ) -> Result<UniquePtr<CxxString>>;
        fn mt_sd_mol_supplier_get_process_property_lists(
            s: &MultithreadedSDMolSupplier,
        ) -> Result<bool>;
        fn mt_sd_mol_supplier_set_process_property_lists(
            s: Pin<&mut MultithreadedSDMolSupplier>,
            val: bool,
        ) -> Result<()>;

        fn mt_smiles_mol_supplier_new() -> Result<UniquePtr<MultithreadedSmilesMolSupplier>>;
        fn mt_smiles_mol_supplier_from_file(
            path: &[u8],
            smiles: &SmilesMolSupplierParamsData,
            params: &MultithreadedParamsData,
        ) -> Result<UniquePtr<MultithreadedSmilesMolSupplier>>;
        fn mt_smiles_mol_supplier_next(
            s: Pin<&mut MultithreadedSmilesMolSupplier>,
            record: &mut u32,
        ) -> Result<UniquePtr<RWMol>>;
        fn mt_smiles_mol_supplier_at_end(
            s: Pin<&mut MultithreadedSmilesMolSupplier>,
        ) -> Result<bool>;
        fn mt_smiles_mol_supplier_get_last_record_id(
            s: &MultithreadedSmilesMolSupplier,
        ) -> Result<u32>;
        fn mt_smiles_mol_supplier_get_last_item_text(
            s: &MultithreadedSmilesMolSupplier,
        ) -> Result<UniquePtr<CxxString>>;
    }
}

fn rust_reader_read(reader: &mut RustReader, buf: &mut [u8]) -> Result<usize, String> {
    reader.read(buf)
}
