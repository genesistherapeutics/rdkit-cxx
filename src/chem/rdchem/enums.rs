// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::ffi::{atom, bond, stereo};

/// Defines a public enum mirroring a Python enum, with conversions to and from
/// the cxx-checked C++ enum (CONVENTIONS.md N4).
macro_rules! rdkit_enum {
    (
        $(#[$meta:meta])*
        $name:ident: $ffi:ty { $($variant:ident = $cxx:ident),* $(,)? }
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            $(#[doc = concat!("Mirrors `", stringify!($cxx), "`.")] $variant,)*
        }

        impl $name {
            #[allow(dead_code)] // not every enum is set from Rust yet
            pub(crate) fn from_ffi(value: $ffi) -> $name {
                $(if value == <$ffi>::$cxx {
                    return $name::$variant;
                })*
                panic!(concat!("unknown ", stringify!($name), " value {}"), value.repr)
            }

            #[allow(dead_code)] // not every enum is set from Rust yet
            pub(crate) fn to_ffi(self) -> $ffi {
                match self {
                    $($name::$variant => <$ffi>::$cxx,)*
                }
            }
        }
    };
}

pub(crate) use rdkit_enum;

rdkit_enum! {
    /// A bond's type. Mirrors `rdchem.BondType`.
    BondType: bond::ffi::BondType {
        Unspecified = UNSPECIFIED,
        Single = SINGLE,
        Double = DOUBLE,
        Triple = TRIPLE,
        Quadruple = QUADRUPLE,
        Quintuple = QUINTUPLE,
        Hextuple = HEXTUPLE,
        Oneandahalf = ONEANDAHALF,
        Twoandahalf = TWOANDAHALF,
        Threeandahalf = THREEANDAHALF,
        Fourandahalf = FOURANDAHALF,
        Fiveandahalf = FIVEANDAHALF,
        Aromatic = AROMATIC,
        Ionic = IONIC,
        Hydrogen = HYDROGEN,
        Threecenter = THREECENTER,
        Dativeone = DATIVEONE,
        Dative = DATIVE,
        Dativel = DATIVEL,
        Dativer = DATIVER,
        Other = OTHER,
        Zero = ZERO,
    }
}

rdkit_enum! {
    /// An atom's hybridization. Mirrors `rdchem.HybridizationType`.
    HybridizationType: atom::ffi::HybridizationType {
        Unspecified = UNSPECIFIED,
        S = S,
        Sp = SP,
        Sp2 = SP2,
        Sp3 = SP3,
        Sp2d = SP2D,
        Sp3d = SP3D,
        Sp3d2 = SP3D2,
        Other = OTHER,
    }
}

rdkit_enum! {
    /// An atom's chirality. Mirrors `rdchem.ChiralType`.
    ChiralType: atom::ffi::ChiralType {
        ChiUnspecified = CHI_UNSPECIFIED,
        ChiTetrahedralCw = CHI_TETRAHEDRAL_CW,
        ChiTetrahedralCcw = CHI_TETRAHEDRAL_CCW,
        ChiOther = CHI_OTHER,
        ChiTetrahedral = CHI_TETRAHEDRAL,
        ChiAllene = CHI_ALLENE,
        ChiSquareplanar = CHI_SQUAREPLANAR,
        ChiTrigonalbipyramidal = CHI_TRIGONALBIPYRAMIDAL,
        ChiOctahedral = CHI_OCTAHEDRAL,
    }
}

rdkit_enum! {
    /// Which valence `get_valence` returns. Mirrors `rdchem.ValenceType`.
    ValenceType: atom::ffi::ValenceType {
        Implicit = IMPLICIT,
        Explicit = EXPLICIT,
    }
}

rdkit_enum! {
    /// A bond's drawing direction, for stereochemistry. Mirrors
    /// `rdchem.BondDir`.
    BondDir: bond::ffi::BondDir {
        None = NONE,
        Beginwedge = BEGINWEDGE,
        Begindash = BEGINDASH,
        Enddownright = ENDDOWNRIGHT,
        Endupright = ENDUPRIGHT,
        Eitherdouble = EITHERDOUBLE,
        Unknown = UNKNOWN,
    }
}

rdkit_enum! {
    /// A bond's stereochemistry. Mirrors `rdchem.BondStereo`.
    BondStereo: bond::ffi::BondStereo {
        Stereonone = STEREONONE,
        Stereoany = STEREOANY,
        Stereoz = STEREOZ,
        Stereoe = STEREOE,
        Stereocis = STEREOCIS,
        Stereotrans = STEREOTRANS,
        Stereoatropcw = STEREOATROPCW,
        Stereoatropccw = STEREOATROPCCW,
    }
}

rdkit_enum! {
    /// How a stereo group's atoms relate. Mirrors `rdchem.StereoGroupType`.
    StereoGroupType: stereo::ffi::StereoGroupType {
        StereoAbsolute = STEREO_ABSOLUTE,
        StereoOr = STEREO_OR,
        StereoAnd = STEREO_AND,
    }
}

rdkit_enum! {
    /// The kind of a potential stereo element. Mirrors `rdchem.StereoType`.
    StereoType: stereo::ffi::StereoType {
        Unspecified = Unspecified,
        AtomTetrahedral = Atom_Tetrahedral,
        AtomSquarePlanar = Atom_SquarePlanar,
        AtomTrigonalBipyramidal = Atom_TrigonalBipyramidal,
        AtomOctahedral = Atom_Octahedral,
        BondDouble = Bond_Double,
        BondCumuleneEven = Bond_Cumulene_Even,
        BondAtropisomer = Bond_Atropisomer,
    }
}

rdkit_enum! {
    /// Whether a stereo element is specified. Mirrors `rdchem.StereoSpecified`.
    StereoSpecified: stereo::ffi::StereoSpecified {
        Unspecified = Unspecified,
        Specified = Specified,
        Unknown = Unknown,
    }
}

rdkit_enum! {
    /// A stereo element's configuration. Mirrors `rdchem.StereoDescriptor`;
    /// Python lacks the atropisomer values. RDKit's `None` is Python's
    /// `NoValue`.
    StereoDescriptor: stereo::ffi::StereoDescriptor {
        NoValue = None,
        TetCw = Tet_CW,
        TetCcw = Tet_CCW,
        BondCis = Bond_Cis,
        BondTrans = Bond_Trans,
        BondAtropCw = Bond_AtropCW,
        BondAtropCcw = Bond_AtropCCW,
    }
}
