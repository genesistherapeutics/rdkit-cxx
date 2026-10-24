//! Mirrors `rdkit.rdBase`.

use crate::error::expect_ok;
use crate::ffi::rd_base::ffi;

/// The RDKit version. Mirrors `rdBase.rdkitVersion`.
pub fn rdkit_version() -> String {
    let version = expect_ok(ffi::rdkit_version());
    version
        .to_str()
        .expect("RDKit's version is UTF-8")
        .to_owned()
}

/// The Boost version RDKit was built with. Mirrors `rdBase.boostVersion`.
pub fn boost_version() -> String {
    let version = expect_ok(ffi::boost_version());
    version
        .to_str()
        .expect("Boost's version is UTF-8")
        .to_owned()
}

/// Disables the log named by `spec`, such as `"rdApp.error"`. Mirrors
/// `rdBase.DisableLog`.
pub fn disable_log(spec: &str) {
    expect_ok(ffi::disable_log(spec));
}

/// Enables the log named by `spec`, such as `"rdApp.error"`. Mirrors
/// `rdBase.EnableLog`.
pub fn enable_log(spec: &str) {
    expect_ok(ffi::enable_log(spec));
}
