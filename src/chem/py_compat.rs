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

#[cfg(test)]
mod tests {
    use super::float_repr;

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
