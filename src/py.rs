// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! What the ports of RDKit's Python-only functionality need from Python and
//! NumPy to compute what Python computes, to the last bit and character
//! (CONVENTIONS.md R2): its parsing and printing of numbers and strings,
//! and its floating-point reductions.

use crate::{Error, Result};

// Parsing, as Python's `str` methods and `int()` and `float()` parse.

/// Whether Python's `str.isspace` counts the character as whitespace.
pub(crate) fn is_py_space(c: char) -> bool {
    c.is_whitespace() || ('\x1c'..='\x1f').contains(&c)
}

/// Python's `str.strip()`.
pub(crate) fn py_strip(s: &str) -> &str {
    s.trim_matches(is_py_space)
}

/// Python's `str.split()`: the runs of non-whitespace.
pub(crate) fn py_split(s: &str) -> Vec<&str> {
    s.split(is_py_space).filter(|w| !w.is_empty()).collect()
}

/// Python's `int(s)` for a string: an optional sign and decimal digits,
/// single underscores allowed between digits, whitespace around.
pub(crate) fn py_int(s: &str) -> Option<i64> {
    let t = py_strip(s);
    let (sign, digits) = match t.strip_prefix('-') {
        Some(rest) => (-1, rest),
        None => (1, t.strip_prefix('+').unwrap_or(t)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
        || !digits.chars().all(|c| c.is_ascii_digit() || c == '_')
    {
        return None;
    }
    let value: i64 = digits.replace('_', "").parse().ok()?;
    Some(sign * value)
}

/// Python's `float(s)` for a string: Rust's parse, after stripping
/// whitespace, with Python's spellings of infinity and NaN, and single
/// underscores between digits.
pub(crate) fn py_float(s: &str) -> Option<f64> {
    let t = py_strip(s);
    let unsigned = t.trim_start_matches(['+', '-']);
    if t.len() - unsigned.len() > 1 {
        return None;
    }
    let lower = unsigned.to_ascii_lowercase();
    if ["inf", "infinity", "nan"].contains(&lower.as_str()) {
        let negative = t.starts_with('-');
        let v = if lower == "nan" {
            f64::NAN
        } else {
            f64::INFINITY
        };
        return Some(if negative { -v } else { v });
    }
    if !unsigned.starts_with(|c: char| c.is_ascii_digit() || c == '.') {
        return None;
    }
    let chars: Vec<char> = t.chars().collect();
    for (i, &c) in chars.iter().enumerate() {
        if c == '_' {
            let digit = |j: Option<&char>| j.is_some_and(char::is_ascii_digit);
            if !digit(i.checked_sub(1).and_then(|j| chars.get(j))) || !digit(chars.get(i + 1)) {
                return None;
            }
        }
    }
    t.replace('_', "").parse().ok()
}

/// Python's `int(s)`, refusing what Python refuses with its message.
pub(crate) fn parse_int(s: &str) -> Result<i64> {
    py_int(s).ok_or_else(|| Error::Value {
        message: format!("invalid literal for int() with base 10: '{s}'"),
    })
}

/// Python's `float(s)`, refusing what Python refuses with its message.
pub(crate) fn parse_float(s: &str) -> Result<f64> {
    py_float(s).ok_or_else(|| Error::Value {
        message: format!("could not convert string to float: '{s}'"),
    })
}

// Printing, as Python's `repr` and `%` formatting print.

/// `repr(x)` (and `str(x)`) of a Python `float`: the shortest digits that
/// round-trip, in fixed notation for exponents from -4 to 15 and in
/// scientific notation, with at least two exponent digits, otherwise.
pub(crate) fn float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_owned();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_owned();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0" } else { "0.0" }.to_owned();
    }
    // Rust's `{:e}` writes the shortest digits that round-trip, as Python's
    // repr does.
    let sci = format!("{x:e}");
    let (mantissa, exp) = sci.split_once('e').expect("{:e} has an exponent");
    let exp: i32 = exp.parse().expect("{:e} has an integer exponent");
    let (sign, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => ("-", m),
        None => ("", mantissa),
    };
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let decpt = exp + 1;
    let body = if -4 < decpt && decpt <= 16 {
        if decpt <= 0 {
            format!("0.{}{digits}", "0".repeat(decpt.unsigned_abs() as usize))
        } else if decpt as usize >= digits.len() {
            format!("{digits}{}.0", "0".repeat(decpt as usize - digits.len()))
        } else {
            let (int, frac) = digits.split_at(decpt as usize);
            format!("{int}.{frac}")
        }
    } else {
        let (first, rest) = digits.split_at(1);
        let mantissa = if rest.is_empty() {
            first.to_owned()
        } else {
            format!("{first}.{rest}")
        };
        format!(
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.unsigned_abs()
        )
    };
    format!("{sign}{body}")
}

/// `repr(s)` of a Python `str`: single quotes unless the string holds one
/// and no double quote, with backslashes, the quote, and characters that
/// are not printable escaped.
pub(crate) fn str_repr(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            // Python's str.isprintable, near enough: ASCII and Latin-1
            // controls are escaped.
            c if (c as u32) < 0x20 || (0x7f..0xa1).contains(&(c as u32)) || c == '\u{ad}' => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// The lines of `s` as Python's `str.splitlines()` splits them, without
/// their line breaks.
pub(crate) fn splitlines(s: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = s.char_indices().peekable();
    while let Some((i, c)) = chars.next() {
        let is_break = matches!(
            c,
            '\n' | '\r'
                | '\x0b'
                | '\x0c'
                | '\x1c'
                | '\x1d'
                | '\x1e'
                | '\u{85}'
                | '\u{2028}'
                | '\u{2029}'
        );
        if is_break {
            lines.push(&s[start..i]);
            let mut end = i + c.len_utf8();
            if c == '\r' && chars.peek().is_some_and(|&(_, n)| n == '\n') {
                chars.next();
                end += 1;
            }
            start = end;
        }
    }
    if start < s.len() {
        lines.push(&s[start..]);
    }
    lines
}

/// Python's `%f`: six decimals, with Python's `nan` and `inf`.
pub(crate) fn format_f(x: f64) -> String {
    if x.is_nan() {
        "nan".to_owned()
    } else {
        format!("{x:.6}")
    }
}

// Reductions and built-ins.

/// Python's `sum` of floats (3.12 and later): Neumaier's compensated sum,
/// left to right from 0.
pub(crate) fn py_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut c) = (0.0f64, 0.0f64);
    for x in values {
        let t = sum + x;
        if sum.abs() >= x.abs() {
            c += (sum - t) + x;
        } else {
            c += (x - t) + sum;
        }
        sum = t;
    }
    if c != 0.0 && c.is_finite() {
        sum += c;
    }
    sum
}

/// NumPy's pairwise sum of a contiguous `float64` array (`np.add.reduce`).
pub(crate) fn np_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - n % 8 {
            for (j, rj) in r.iter_mut().enumerate() {
                *rj += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for &x in &a[i..] {
            res += x;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        np_sum(&a[..n2]) + np_sum(&a[n2..])
    }
}

/// NumPy's `np.std` of a list of floats: the population standard deviation.
pub(crate) fn np_std(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    let mean = (0.0 + np_sum(a)) / n;
    let squares: Vec<f64> = a.iter().map(|&x| (x - mean) * (x - mean)).collect();
    ((0.0 + np_sum(&squares)) / n).sqrt()
}

/// Python's `math.exp`, which raises OverflowError where the result
/// overflows: `Error::Value` then.
pub(crate) fn py_exp(x: f64) -> Result<f64> {
    let y = x.exp();
    if y.is_infinite() && x.is_finite() {
        return Err(Error::Value {
            message: "math range error".to_owned(),
        });
    }
    Ok(y)
}

/// Python's float division, which raises ZeroDivisionError for a zero
/// divisor: `Error::Value` then.
pub(crate) fn py_div(a: f64, b: f64) -> Result<f64> {
    if b == 0.0 {
        return Err(Error::Value {
            message: "float division by zero".to_owned(),
        });
    }
    Ok(a / b)
}

/// Python's `bisect.bisect_right`: the index after the last element of
/// `bins` that is not greater than `x`, found by Python's binary search, so
/// unsorted bins and NaN give Python's answers.
pub(crate) fn bisect_right(bins: &[f64], x: f64) -> usize {
    bisect_right_by(bins, |b| x < *b)
}

/// `bisect_right` for any elements, with `is_below(e)` for Python's
/// `x < e`.
pub(crate) fn bisect_right_by<T>(a: &[T], is_below: impl Fn(&T) -> bool) -> usize {
    let (mut lo, mut hi) = (0, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if is_below(&a[mid]) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Python's `max(values)`: the first value no later one exceeds, as Python
/// compares them (so a leading NaN wins). `Error::Value` for no values, as
/// Python raises.
pub(crate) fn py_max(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut values = values.into_iter();
    let first = values.next().ok_or_else(empty)?;
    Ok(values.fold(first, |best, x| if x > best { x } else { best }))
}

/// Python's `min(values)`, as `py_max`.
pub(crate) fn py_min(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut values = values.into_iter();
    let first = values.next().ok_or_else(empty)?;
    Ok(values.fold(first, |best, x| if x < best { x } else { best }))
}

fn empty() -> Error {
    Error::Value {
        message: "max() or min() of an empty sequence".to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::{float_repr, parse_float, parse_int, py_float, py_int, splitlines, str_repr};

    #[test]
    fn floats_parse_as_python_parses_them() {
        for (s, value) in [
            ("1.5", 1.5),
            (" 2\n", 2.0),
            ("+.5", 0.5),
            ("1.", 1.0),
            ("-1e3", -1000.0),
            ("1_0", 10.0),
            ("1_0.2_5", 10.25),
            ("1e1_0", 1e10),
            ("Infinity", f64::INFINITY),
            ("-inf", f64::NEG_INFINITY),
            ("\u{2003}3\u{1f}", 3.0),
        ] {
            assert_eq!(py_float(s), Some(value), "{s:?}");
            assert_eq!(parse_float(s).unwrap(), value, "{s:?}");
        }
        assert!(py_float("nan").unwrap().is_nan());
        assert!(py_float("-NaN").unwrap().is_nan());
        for s in [
            "", "_1", "1_", "1__0", "1_.5", "1._5", "1e_5", "+-1", "0x10", "e5", "in f", "1 0",
        ] {
            assert_eq!(py_float(s), None, "{s:?}");
            assert_eq!(
                parse_float(s).unwrap_err().to_string(),
                format!("could not convert string to float: '{s}'")
            );
        }
    }

    #[test]
    fn ints_parse_as_python_parses_them() {
        for (s, value) in [("12", 12), (" -3 ", -3), ("+4", 4), ("1_000", 1000)] {
            assert_eq!(py_int(s), Some(value), "{s:?}");
            assert_eq!(parse_int(s).unwrap(), value, "{s:?}");
        }
        for s in ["", "1.0", "_1", "1_", "1__0", "--1", "0x10"] {
            assert_eq!(py_int(s), None, "{s:?}");
            assert!(parse_int(s).is_err(), "{s:?}");
        }
    }

    #[test]
    fn str_reprs_match_python() {
        for (s, repr) in [
            ("abc", "'abc'"),
            ("a'b", "\"a'b\""),
            ("a'b\"", "'a\\'b\"'"),
            ("a\tb\n", "'a\\tb\\n'"),
            ("back\\slash", "'back\\\\slash'"),
            ("\x01é", "'\\x01é'"),
        ] {
            assert_eq!(str_repr(s), repr);
        }
    }

    #[test]
    fn lines_split_as_python_splits_them() {
        assert_eq!(splitlines("a\nb\r\nc\rd\x0ce\n"), ["a", "b", "c", "d", "e"]);
        assert_eq!(splitlines("\na\n\n"), ["", "a", ""]);
        assert_eq!(splitlines(""), Vec::<&str>::new());
    }

    #[test]
    fn float_reprs_match_python() {
        for (x, repr) in [
            (1.0, "1.0"),
            (0.1, "0.1"),
            (-2.5, "-2.5"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (1e15, "1000000000000000.0"),
            (123456789012345678.0, "1.2345678901234568e+17"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (1.2345e-7, "1.2345e-07"),
            (0.30000000000000004, "0.30000000000000004"),
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (f64::INFINITY, "inf"),
            (f64::NAN, "nan"),
            (1e100, "1e+100"),
            (12.0, "12.0"),
        ] {
            assert_eq!(float_repr(x), repr);
        }
    }
}
