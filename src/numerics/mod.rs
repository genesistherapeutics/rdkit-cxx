// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Numerics`, which exposes `rdAlignment` as `Alignment`.

pub mod rd_alignment;

pub use rd_alignment as alignment;
