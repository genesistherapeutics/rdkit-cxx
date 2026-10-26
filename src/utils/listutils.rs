// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.utils.listutils`, a pure-Python module: utilities for
//! lists.

use std::fmt::Debug;

/// A compact expression for a list, with runs of equal items multiplied
/// out: `[0, 1, 1, 1, 1, 0]` gives `[0]+[1]*4+[0]`, and an empty list
/// `[]`. Mirrors `listutils.CompactListRepr`.
///
/// Items are written with `Debug`, which for integers and most floats is
/// Python's `repr`; strings come out in double quotes, where Python uses
/// single ones.
pub fn compact_list_repr<T: PartialEq + Debug>(lst: &[T]) -> String {
    let Some(first) = lst.first() else {
        return "[]".to_owned();
    };
    let label = |item: &T, count: usize| {
        if count > 1 {
            format!("[{item:?}]*{count}")
        } else {
            format!("[{item:?}]")
        }
    };
    let mut components = Vec::new();
    let mut last = first;
    let mut count = 1;
    for item in &lst[1..] {
        if item != last {
            components.push(label(last, count));
            count = 1;
            last = item;
        } else {
            count += 1;
        }
    }
    components.push(label(last, count));
    components.join("+")
}
