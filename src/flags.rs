// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! What every flag newtype shares (CONVENTIONS.md F1).

/// Implements `bits`, `from_bits_retain`, the bit operators, and a `Debug`
/// that names the value from `$type::NAMES`, for a `u32` newtype.
macro_rules! flag_ops {
    ($type:ident) => {
        impl $type {
            /// The raw bits.
            pub const fn bits(self) -> u32 {
                self.0
            }

            /// Flags from raw bits, keeping bits no constant names.
            pub const fn from_bits_retain(bits: u32) -> $type {
                $type(bits)
            }
        }

        impl std::fmt::Debug for $type {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match $type::NAMES.iter().find(|(_, flags)| flags == self) {
                    Some((name, _)) => write!(f, concat!(stringify!($type), "::{}"), name),
                    None => write!(f, concat!(stringify!($type), "({:#x})"), self.0),
                }
            }
        }

        impl std::ops::BitOr for $type {
            type Output = $type;

            fn bitor(self, rhs: $type) -> $type {
                $type(self.0 | rhs.0)
            }
        }

        impl std::ops::BitOrAssign for $type {
            fn bitor_assign(&mut self, rhs: $type) {
                self.0 |= rhs.0;
            }
        }

        impl std::ops::BitAnd for $type {
            type Output = $type;

            fn bitand(self, rhs: $type) -> $type {
                $type(self.0 & rhs.0)
            }
        }

        impl std::ops::BitAndAssign for $type {
            fn bitand_assign(&mut self, rhs: $type) {
                self.0 &= rhs.0;
            }
        }

        impl std::ops::Not for $type {
            type Output = $type;

            fn not(self) -> $type {
                $type(!self.0)
            }
        }
    };
}

pub(crate) use flag_ops;
