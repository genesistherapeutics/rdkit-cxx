use rdkit_cxx::rd_base;

#[test]
fn rdkit_version_is_the_supported_series() {
    assert!(rd_base::rdkit_version().starts_with("2026.09."));
}

#[test]
fn boost_version_is_empty_or_boost_lib_version() {
    // RDKit fills this from CMake's Boost_LIB_VERSION ("1_90"), which a
    // CONFIG-mode find_package(Boost) leaves unset, as in conda-forge's build.
    let version = rd_base::boost_version();
    assert!(version.is_empty() || version.contains('_'), "{version}");
}

#[test]
fn log_can_be_disabled_and_enabled() {
    rd_base::disable_log("rdApp.error");
    rd_base::enable_log("rdApp.error");
}
