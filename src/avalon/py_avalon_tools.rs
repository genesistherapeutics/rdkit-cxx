// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Avalon.pyAvalonTools`: the Avalon toolkit's
//! fingerprints, canonical SMILES, 2D coordinates, and structure checker.
//!
//! Only the module's constants are bound. Its functions run the Avalon
//! toolkit (`External/AvalonTools`, about 34,000 lines of C with unlocked
//! global state), which parses SMILES and mol blocks without validating
//! them: malformed SMILES such as `ClC(Cl)=C1=CC=CC..1` corrupt its heap,
//! and malformed mol blocks are read out of bounds (crashes in Python).
//! Molecules reach the same code as the SMILES or mol blocks RDKit writes
//! for them, and the fingerprints and the layout write past an atom's
//! table of 20 neighbors, and go on, for an atom of more bonds, which a
//! sanitized molecule can have. The structure checker reads a file it
//! could not open through a null `FILE*`. This binding has not audited the
//! toolkit for the input it trusts (CONVENTIONS.md U1); the manifest
//! (`inventory/manifest.toml`) gives each function's reason.

/// The fingerprint bits for substructure screening, as a fingerprint's
/// `bitFlags`. Mirrors `pyAvalonTools.avalonSSSBits`.
pub fn avalon_sss_bits() -> i32 {
    0x007FFF
}

/// The fingerprint bits for similarity, as a fingerprint's `bitFlags`.
/// Mirrors `pyAvalonTools.avalonSimilarityBits`.
pub fn avalon_similarity_bits() -> i32 {
    0xF07FFF
}

/// Defines a flag type (CONVENTIONS.md F1) whose constants keep Python's
/// spelling, with the values of the toolkit's `struchk.h`.
macro_rules! struchk_flags {
    ($(#[$meta:meta])* $type:ident { $($name:ident = $value:expr),* $(,)? }) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $type(u32);

        #[allow(non_upper_case_globals)]
        impl $type {
            $(
                #[doc = concat!("Mirrors `", stringify!($type), ".", stringify!($name), "`.")]
                pub const $name: $type = $type($value);
            )*

            const NAMES: &[(&str, $type)] = &[$((stringify!($name), $type::$name)),*];
        }

        crate::flags::flag_ops!($type);
    };
}

struchk_flags!(
    /// The corrections and problems the structure checker reports: bits,
    /// combined. Mirrors `pyAvalonTools.StruChkFlag` (F1); the constants
    /// keep Python's spelling.
    StruChkFlag {
        bad_molecule = 1,
        alias_conversion_failed = 2,
        transformed = 4,
        fragments_found = 8,
        either_warning = 16,
        stereo_error = 32,
        dubious_stereo_removed = 64,
        atom_clash = 128,
        atom_check_failed = 256,
        size_check_failed = 512,
        recharged = 1024,
        stereo_forced_bad = 2048,
        stereo_transformed = 4096,
        template_transformed = 8192,
    }
);

struchk_flags!(
    /// Masks of the structure checker's flags: those of a bad structure
    /// (`bad_set`) and of a transformed one (`transformed_set`, which also
    /// has struchk's 16384, for a tautomer transformed). Mirrors
    /// `pyAvalonTools.StruChkResult` (F1); the constants keep Python's
    /// spelling.
    StruChkResult {
        success = 0,
        bad_set = 2979,
        transformed_set = 29788,
    }
);
