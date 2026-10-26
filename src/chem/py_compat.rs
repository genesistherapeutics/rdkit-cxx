// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! What the ports of RDKit's pure-Python modules need from Python itself to
//! compute what Python computes (CONVENTIONS.md R2).

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

#[cfg(test)]
mod tests {
    use super::{float_repr, splitlines, str_repr};

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
