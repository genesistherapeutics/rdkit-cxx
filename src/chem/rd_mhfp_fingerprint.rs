// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMHFPFingerprint`: MinHash fingerprints (MHFP)
//! and their folded form (SECFP).

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::data_structs::ExplicitBitVect;
use crate::error::expect_ok;
use crate::ffi::mhfp::ffi;

/// The shingling arguments of every `MHFPEncoder` method but the SECFP
/// ones, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EncodeParams {
    /// The largest radius of the circular substructures.
    pub radius: u8,
    /// Adds the SSSR rings to the shingling.
    pub rings: bool,
    /// Writes isomeric SMILES.
    pub isomeric: bool,
    /// Kekulizes the molecule first; `Error::Kekulize` when it cannot be.
    pub kekulize: bool,
    /// The smallest radius; 0 adds each atom.
    pub min_radius: u8,
}

impl Default for EncodeParams {
    fn default() -> EncodeParams {
        EncodeParams {
            radius: 3,
            rings: true,
            isomeric: false,
            kekulize: true,
            min_radius: 1,
        }
    }
}

impl EncodeParams {
    fn to_ffi(&self) -> ffi::ShinglingOptions {
        ffi::ShinglingOptions {
            radius: self.radius,
            rings: self.rings,
            isomeric: self.isomeric,
            kekulize: self.kekulize,
            min_radius: self.min_radius,
        }
    }
}

/// The arguments of the SECFP methods, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EncodeSECFPParams {
    /// The shingling arguments.
    pub shingling: EncodeParams,
    /// The fingerprint's size in bits; `Error::Value` for 0 or more than
    /// `i32::MAX` (RDKit folds modulo it, U3).
    pub length: usize,
}

impl Default for EncodeSECFPParams {
    fn default() -> EncodeSECFPParams {
        EncodeSECFPParams {
            shingling: EncodeParams::default(),
            length: 2048,
        }
    }
}

/// Parameters for `MHFPEncoder::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MHFPEncoderParams {
    /// The number of hash permutations: the MinHash's length.
    pub n_permutations: u32,
    /// The permutations' random seed.
    pub seed: u32,
}

impl Default for MHFPEncoderParams {
    fn default() -> MHFPEncoderParams {
        MHFPEncoderParams {
            n_permutations: 2048,
            seed: 42,
        }
    }
}

/// Encodes molecules as MinHashes of their circular substructures' SMILES.
/// Mirrors `rdMHFPFingerprint.MHFPEncoder`. Encoders with different
/// permutation counts or seeds give incomparable MinHashes.
///
/// The SMILES forms return `Error::Parse` for SMILES that do not parse.
pub struct MHFPEncoder {
    /// Never null.
    ptr: UniquePtr<ffi::MHFPEncoder>,
}

impl fmt::Debug for MHFPEncoder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MHFPEncoder").finish_non_exhaustive()
    }
}

impl Default for MHFPEncoder {
    fn default() -> MHFPEncoder {
        MHFPEncoder::new()
    }
}

impl MHFPEncoder {
    /// An encoder with 2048 permutations and seed 42. Mirrors
    /// `MHFPEncoder()`.
    pub fn new() -> MHFPEncoder {
        MHFPEncoder::new_with_params(&MHFPEncoderParams::default())
    }

    /// Mirrors `MHFPEncoder(n_permutations, seed)`.
    pub fn new_with_params(params: &MHFPEncoderParams) -> MHFPEncoder {
        MHFPEncoder {
            ptr: expect_ok(ffi::mhfp_encoder_new(params.n_permutations, params.seed)),
        }
    }

    /// The MinHash of arbitrary strings. Mirrors
    /// `MHFPEncoder.FromStringArray`.
    pub fn from_string_array<S: AsRef<str>>(&self, vec: &[S]) -> Vec<u32> {
        let strs: Vec<&str> = vec.iter().map(AsRef::as_ref).collect();
        expect_ok(ffi::mhfp_from_string_array(&self.ptr, &strs))
    }

    /// The MinHash of arbitrary integers. Mirrors `MHFPEncoder.FromArray`.
    pub fn from_array(&self, vec: &[u32]) -> Vec<u32> {
        expect_ok(ffi::mhfp_from_array(&self.ptr, vec))
    }

    /// The SMILES of the molecule's circular substructures (and rings).
    /// Mirrors `MHFPEncoder.CreateShinglingFromSmiles`.
    pub fn create_shingling_from_smiles(&self, smiles: &str) -> Result<Vec<String>> {
        self.create_shingling_from_smiles_with_params(smiles, &EncodeParams::default())
    }

    /// Mirrors `MHFPEncoder.CreateShinglingFromSmiles(smiles, radius, rings,
    /// isomeric, kekulize, min_radius)`.
    pub fn create_shingling_from_smiles_with_params(
        &self,
        smiles: &str,
        params: &EncodeParams,
    ) -> Result<Vec<String>> {
        Ok(ffi::mhfp_shingling_from_smiles(
            &self.ptr,
            smiles,
            &params.to_ffi(),
        )?)
    }

    /// Mirrors `MHFPEncoder.CreateShinglingFromMol`.
    pub fn create_shingling_from_mol(&self, mol: &Mol) -> Result<Vec<String>> {
        self.create_shingling_from_mol_with_params(mol, &EncodeParams::default())
    }

    /// Mirrors `MHFPEncoder.CreateShinglingFromMol(mol, radius, rings,
    /// isomeric, kekulize, min_radius)`.
    pub fn create_shingling_from_mol_with_params(
        &self,
        mol: &Mol,
        params: &EncodeParams,
    ) -> Result<Vec<String>> {
        Ok(ffi::mhfp_shingling_from_mol(
            &self.ptr,
            mol.as_ffi(),
            &params.to_ffi(),
        )?)
    }

    /// The MinHash of a molecule's shingling. Mirrors
    /// `MHFPEncoder.EncodeSmiles`.
    pub fn encode_smiles(&self, smiles: &str) -> Result<Vec<u32>> {
        self.encode_smiles_with_params(smiles, &EncodeParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeSmiles(smiles, radius, rings, isomeric,
    /// kekulize, min_radius)`.
    pub fn encode_smiles_with_params(
        &self,
        smiles: &str,
        params: &EncodeParams,
    ) -> Result<Vec<u32>> {
        Ok(ffi::mhfp_encode_smiles(
            &self.ptr,
            smiles,
            &params.to_ffi(),
        )?)
    }

    /// Mirrors `MHFPEncoder.EncodeMol`.
    pub fn encode_mol(&self, mol: &Mol) -> Result<Vec<u32>> {
        self.encode_mol_with_params(mol, &EncodeParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeMol(mol, radius, rings, isomeric,
    /// kekulize, min_radius)`.
    pub fn encode_mol_with_params(&self, mol: &Mol, params: &EncodeParams) -> Result<Vec<u32>> {
        Ok(ffi::mhfp_encode_mol(
            &self.ptr,
            mol.as_ffi(),
            &params.to_ffi(),
        )?)
    }

    /// Each SMILES's MinHash. Mirrors `MHFPEncoder.EncodeSmilesBulk`.
    pub fn encode_smiles_bulk<S: AsRef<str>>(&self, smiles: &[S]) -> Result<Vec<Vec<u32>>> {
        self.encode_smiles_bulk_with_params(smiles, &EncodeParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeSmilesBulk(smiles, radius, rings,
    /// isomeric, kekulize, min_radius)`.
    pub fn encode_smiles_bulk_with_params<S: AsRef<str>>(
        &self,
        smiles: &[S],
        params: &EncodeParams,
    ) -> Result<Vec<Vec<u32>>> {
        smiles
            .iter()
            .map(|s| self.encode_smiles_with_params(s.as_ref(), params))
            .collect()
    }

    /// Each molecule's MinHash. Mirrors `MHFPEncoder.EncodeMolsBulk`.
    pub fn encode_mols_bulk(&self, mols: &[Mol]) -> Result<Vec<Vec<u32>>> {
        self.encode_mols_bulk_with_params(mols, &EncodeParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeMolsBulk(mols, radius, rings, isomeric,
    /// kekulize, min_radius)`.
    pub fn encode_mols_bulk_with_params(
        &self,
        mols: &[Mol],
        params: &EncodeParams,
    ) -> Result<Vec<Vec<u32>>> {
        mols.iter()
            .map(|mol| self.encode_mol_with_params(mol, params))
            .collect()
    }

    /// The shingling's hashes folded into a bit vector (SECFP). Mirrors
    /// `MHFPEncoder.EncodeSECFPSmiles`.
    pub fn encode_secfp_smiles(&self, smiles: &str) -> Result<ExplicitBitVect> {
        self.encode_secfp_smiles_with_params(smiles, &EncodeSECFPParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeSECFPSmiles(smiles, radius, rings,
    /// isomeric, kekulize, min_radius, length)`.
    pub fn encode_secfp_smiles_with_params(
        &self,
        smiles: &str,
        params: &EncodeSECFPParams,
    ) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(ffi::mhfp_encode_secfp_smiles(
            &self.ptr,
            smiles,
            &params.shingling.to_ffi(),
            params.length,
        )?))
    }

    /// Mirrors `MHFPEncoder.EncodeSECFPMol`.
    pub fn encode_secfp_mol(&self, mol: &Mol) -> Result<ExplicitBitVect> {
        self.encode_secfp_mol_with_params(mol, &EncodeSECFPParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeSECFPMol(mol, radius, rings, isomeric,
    /// kekulize, min_radius, length)`.
    pub fn encode_secfp_mol_with_params(
        &self,
        mol: &Mol,
        params: &EncodeSECFPParams,
    ) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(ffi::mhfp_encode_secfp_mol(
            &self.ptr,
            mol.as_ffi(),
            &params.shingling.to_ffi(),
            params.length,
        )?))
    }

    /// Each SMILES's SECFP. Mirrors `MHFPEncoder.EncodeSECFPSmilesBulk`.
    pub fn encode_secfp_smiles_bulk<S: AsRef<str>>(
        &self,
        smiles: &[S],
    ) -> Result<Vec<ExplicitBitVect>> {
        self.encode_secfp_smiles_bulk_with_params(smiles, &EncodeSECFPParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeSECFPSmilesBulk(smiles, radius, rings,
    /// isomeric, kekulize, min_radius, length)`.
    pub fn encode_secfp_smiles_bulk_with_params<S: AsRef<str>>(
        &self,
        smiles: &[S],
        params: &EncodeSECFPParams,
    ) -> Result<Vec<ExplicitBitVect>> {
        smiles
            .iter()
            .map(|s| self.encode_secfp_smiles_with_params(s.as_ref(), params))
            .collect()
    }

    /// Each molecule's SECFP. Mirrors `MHFPEncoder.EncodeSECFPMolsBulk`.
    pub fn encode_secfp_mols_bulk(&self, mols: &[Mol]) -> Result<Vec<ExplicitBitVect>> {
        self.encode_secfp_mols_bulk_with_params(mols, &EncodeSECFPParams::default())
    }

    /// Mirrors `MHFPEncoder.EncodeSECFPMolsBulk(mols, radius, rings,
    /// isomeric, kekulize, min_radius, length)`.
    pub fn encode_secfp_mols_bulk_with_params(
        &self,
        mols: &[Mol],
        params: &EncodeSECFPParams,
    ) -> Result<Vec<ExplicitBitVect>> {
        mols.iter()
            .map(|mol| self.encode_secfp_mol_with_params(mol, params))
            .collect()
    }

    /// The fraction of positions where two MinHashes differ. Mirrors
    /// `MHFPEncoder.Distance`. `Error::Value` for MinHashes of different
    /// lengths (RDKit reads the second at every index of the first); NaN
    /// for two empty ones.
    pub fn distance(a: &[u32], b: &[u32]) -> Result<f64> {
        Ok(ffi::mhfp_distance(a, b)?)
    }
}
